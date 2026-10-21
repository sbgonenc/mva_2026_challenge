=head1 LICENSE

Copyright [2026] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME

 SpliceAI2

=head1 SYNOPSIS

 mv SpliceAI2.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin SpliceAI2,snv=/path/to/snv/chr*.tsv.gz,indel=/path/to/indel/chr*.tsv.gz
 ./vep -i variations.vcf --plugin SpliceAI2,file=/path/to/variants.tsv.spliceai2.gz,cols=all

=head1 DESCRIPTION

 An Ensembl VEP plugin that retrieves pre-calculated SpliceAI2 predictions.
 SpliceAI2, developed by Illumina, Inc., is a deep learning model that predicts
 variant effects on splicing through quantitative modelling of splice sites,
 splice junctions and transcripts. It is the successor of SpliceAI.

 SpliceAI2 predictions are summarised into a single spliceai2_summary_score
 ranging from 0 to 1 (the maximum of the donor/acceptor gain/loss delta scores).
 The author-suggested cutoffs (with their SpliceAI equivalents) are:
   * 0.1  (SpliceAI 0.2, high recall)
   * 0.25 (SpliceAI 0.5, balance of precision and recall)
   * 0.5  (SpliceAI 0.8, high precision)

 This plugin is available for GRCh38.

 More information can be found at:
 https://github.com/Illumina/SpliceAI2

 Please cite the SpliceAI2 publication alongside Ensembl VEP if you use this resource:
 Jaganathan K, Chen J, Liu X, Zhang Y, et al. A unified framework for quantitative
 splicing and transcript prediction. bioRxiv (2026). https://doi.org/10.64898/2026.10.08.757673

 Options are passed to the plugin as key=value pairs:

   snv      : Tabix-indexed SpliceAI2 file(s) with SNV scores
   indel    : Tabix-indexed SpliceAI2 file(s) with indel scores
   file     : Any other tabix-indexed SpliceAI2 file(s), e.g. the output of the
              spliceai2 command line tool on your own variants
              At least one of snv, indel or file is mandatory. Each value can be
              a single file, a directory (all *.gz files in it are loaded) or a
              glob pattern (e.g. /path/to/snv/chr*.tsv.gz); the last two load the
              per-chromosome files distributed by Illumina.

   cols     : (optional) Colon-separated list of columns to return from the
              SpliceAI2 files, or 'all' to return all score columns
              (default: "spliceai2_summary_score" plus the delta scores and
              distances of the strongest donor/acceptor gain/loss, which
              correspond to SpliceAI's DS_* and DP_* fields)

   match_to : (optional) How SpliceAI2 predictions are matched to transcripts.
              SpliceAI2 scores a variant separately for each gene strand.
                "strand" (default): only return predictions computed for the
                                    strand of the overlapping transcript
                "any"             : match on genomic position and alleles only

 Output:
  Each returned column is prefixed with "SpliceAI2_", e.g. "SpliceAI2_spliceai2_summary_score"
  and "SpliceAI2_donor_gain_delta_score_0". Column names follow the
  {effect_type}_{field}_{k} scheme of SpliceAI2, where k ranks the strongest changes per
  effect type (0 being the strongest). The pre-calculated files keep the 3 strongest
  changes (k = 0-2); the spliceai2 tool itself reports 10 (k = 0-9). If more than one prediction matches
  a transcript, the one with the highest spliceai2_summary_score is returned.

  For JSON output the selected columns are returned (without prefix) as a hash under
  the key "SpliceAI2".

 The following steps are necessary to run this plugin:

 1. Download the pre-calculated scores (GRCh38, GENCODE v48 protein-coding genes;
    ~350 GB) from Hugging Face. Access is restricted to non-commercial use and
    requires accepting the SpliceAI2 terms of use:
    https://huggingface.co/datasets/illumina-ai/SpliceAI2-data

    hf download illumina-ai/SpliceAI2-data --repo-type dataset \
      --include "precomputed_scores_v2.0/*" --local-dir SpliceAI2_data

 2. Index the files. They are already sorted and bgzip-compressed, with the column
    names on the first line, which is skipped by tabix (-S 1) and read by the plugin:

    for f in SpliceAI2_data/precomputed_scores_v2.0/GRCh38/{snv,indel}/chr*.tsv.gz; do
      tabix -s 1 -b 2 -e 2 -S 1 $f
    done

    Files produced by running the spliceai2 tool yourself (columns chrom, pos, ref,
    alt, strand followed by the predictions) must first be sorted by chromosome and
    position, compressed with bgzip and then indexed in the same way.

 3. Run the plugin:

    ./vep -i variations.vcf --plugin SpliceAI2,snv=/path/to/snv/chr*.tsv.gz,indel=/path/to/indel/chr*.tsv.gz
    ./vep -i variations.vcf --plugin SpliceAI2,snv=/path/to/snv,indel=/path/to/indel,cols=all

 You must have the Bio::DB::HTS module or the tabix utility must be installed
 in your path to use this plugin.

=cut

package SpliceAI2;

use strict;
use warnings;

use File::Basename;
use IO::Uncompress::Gunzip qw($GunzipError);

use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my $PREFIX = 'SpliceAI2_';

# Columns identifying the variant; these are not returned as annotations
my @ID_COLS = qw(chrom pos ref alt strand);

my @DEFAULT_COLS = qw(
  spliceai2_summary_score
  donor_gain_delta_score_0 donor_loss_delta_score_0
  acceptor_gain_delta_score_0 acceptor_loss_delta_score_0
  donor_gain_dist_0 donor_loss_dist_0
  acceptor_gain_dist_0 acceptor_loss_dist_0
);

my %EFFECT_DESCRIPTION = (
  donor_gain    => 'new or strengthened donor',
  donor_loss    => 'lost or weakened donor',
  acceptor_gain => 'new or strengthened acceptor',
  acceptor_loss => 'lost or weakened acceptor',
  jxn_gain      => 'new or strengthened junction',
  jxn_loss      => 'lost or weakened junction',
);

my %FIELD_DESCRIPTION = (
  delta_score   => 'absolute change in splicing probability',
  ref_score     => 'splicing probability for the reference allele',
  alt_score     => 'splicing probability for the alternate allele',
  dist          => 'distance (bp) between variant and splice site',
  donor_dist    => 'distance (bp) between variant and donor site',
  acceptor_dist => 'distance (bp) between variant and acceptor site',
);

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  my $param_hash = $self->params_to_hash();

  my @files;
  for my $key (qw(snv indel file)) {
    next unless defined $param_hash->{$key};
    my $value = $param_hash->{$key};
    my @matched = -d $value ? glob("$value/*.gz") : grep { -e $_ } glob($value);
    die "\n  ERROR: No SpliceAI2 file found for $key=$param_hash->{$key}\n" unless @matched;
    push @files, @matched;
  }
  die "\n  ERROR: No SpliceAI2 file specified\nTry using 'snv=/path/to/snv/chr*.tsv.gz,indel=/path/to/indel/chr*.tsv.gz'\n"
    unless @files;
  $self->add_file($_) for @files;

  $self->_parse_colnames($param_hash->{cols} || join(':', @DEFAULT_COLS));

  my $match_to = $param_hash->{match_to} || 'strand';
  die "\n  ERROR: match_to must be 'strand' or 'any'\n"
    unless $match_to eq 'strand' || $match_to eq 'any';
  if ($match_to eq 'strand' && !defined $self->{col_index}->{strand}) {
    warn "WARNING: SpliceAI2 file has no 'strand' column; matching on position and alleles only\n";
    $match_to = 'any';
  }
  $self->{match_to} = $match_to;

  # Chromosome naming of the SpliceAI2 files (Illumina uses 'chr1', VEP uses '1')
  $self->{chr_prefix} = $self->_files_use_chr_prefix();

  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub _read_header {
  my ($self, $file) = @_;

  # Header prefixed with '#'
  open(my $fh, "tabix -H $file |") or die "ERROR: cannot open tabix file $file: $!\n";
  my $last;
  $last = $_ while <$fh>;
  close $fh;

  # Header on the first line, skipped by tabix (Illumina files indexed with -S 1)
  unless (defined $last) {
    my $gz = IO::Uncompress::Gunzip->new($file) or die "ERROR: cannot open $file: $GunzipError\n";
    $last = $gz->getline;
    $gz->close;
  }

  die "\n  ERROR: No column header found in $file\n" .
      "  The first line of the file must contain the column names (see plugin documentation)\n"
    unless defined $last && $last =~ /^#?chrom\t/;

  $last =~ s/(^#|\s+$)//g;
  return [ split /\t/, $last ];
}

sub _parse_colnames {
  my ($self, $cols) = @_;

  my @files = @{ $self->files };
  my $colnames = $self->_read_header($files[0]);

  # All files (SNV, indel, per chromosome) must share the same layout
  my $expected = join("\t", @$colnames);
  for my $file (@files[1..$#files]) {
    die "\n  ERROR: Column header of $file differs from $files[0]\n"
      unless join("\t", @{ $self->_read_header($file) }) eq $expected;
  }

  my %col_index;
  @col_index{@$colnames} = (0..$#$colnames);

  my @missing = grep { !defined $col_index{$_} } qw(chrom pos ref alt);
  die "\n  ERROR: Required columns not found in file header: ", join(", ", @missing), "\n" if @missing;

  my %is_id_col = map { $_ => 1 } @ID_COLS;
  my @score_cols = grep { !$is_id_col{$_} } @$colnames;

  my @selected;
  if ($cols eq 'all') {
    @selected = @score_cols;
  } else {
    @selected = split /:/, $cols;
    my @invalid = grep { !defined $col_index{$_} } @selected;
    die "\n  ERROR: The following columns were not found in file header: ", join(", ", @invalid), "\n"
      if @invalid;
  }

  $self->{colnames}  = $colnames;
  $self->{col_index} = \%col_index;
  $self->{cols}      = \@selected;
}

sub _files_use_chr_prefix {
  my $self = shift;

  open(my $fh, "tabix -l " . $self->files->[0] . " |") or return 0;
  my $first = <$fh>;
  close $fh;

  return (defined $first && $first =~ /^chr/) ? 1 : 0;
}

sub _describe_col {
  my $col = shift;

  return "SpliceAI2 summary score (0-1): max of donor/acceptor gain/loss delta scores. " .
         "Guide cutoffs: >=0.1 high recall; >=0.25 balanced; >=0.5 high precision"
    if $col eq 'spliceai2_summary_score';

  if ($col =~ /^(donor_gain|donor_loss|acceptor_gain|acceptor_loss|jxn_gain|jxn_loss)_(.+)_(\d+)$/
      && $FIELD_DESCRIPTION{$2}) {
    my $rank = $3 == 0 ? 'strongest change' : 'change ranked ' . ($3 + 1);
    return "SpliceAI2 $FIELD_DESCRIPTION{$2} for $EFFECT_DESCRIPTION{$1} ($rank)";
  }

  return "SpliceAI2 column $col";
}

sub get_header_info {
  my $self = shift;

  my @files = @{ $self->files };
  my $source = basename $files[0];
  $source .= " and " . (@files - 1) . " other file(s)" if @files > 1;

  my %header = map { $PREFIX . $_ => _describe_col($_) . " (from $source)" } @{ $self->{cols} };
  return \%header;
}

sub run {
  my ($self, $tva) = @_;

  my $vf = $tva->variation_feature;
  return {} unless defined $vf->{chr};

  (my $chr = $vf->{chr}) =~ s/^chr//i;
  $chr = 'chr' . $chr if $self->{chr_prefix};

  my ($start, $end) = ($vf->{start}, $vf->{end});
  ($start, $end) = ($end, $start) if $start > $end;

  # SpliceAI2 indels are VCF-style, anchored on the preceding base
  my @data = @{ $self->get_data($chr, $start - 1, $end) };
  return {} unless @data;

  my $tr_strand = $tva->transcript->strand;
  my $vf_ref = $vf->ref_allele_string;
  my $vf_alt = $tva->variation_feature_seq;

  my $best;
  for my $candidate (@data) {
    next if $self->{match_to} eq 'strand' && $candidate->{strand} != $tr_strand;

    # SpliceAI2 alleles are always given on the forward strand
    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf_ref,
        alts   => [ $vf_alt ],
        pos    => $vf->{start},
        strand => $vf->strand,
      },
      {
        ref    => $candidate->{ref},
        alts   => [ $candidate->{alt} ],
        pos    => $candidate->{start},
        strand => 1,
      }
    );
    next unless @$matches;

    $best = $candidate if !defined $best || _summary($candidate) > _summary($best);
  }

  return {} unless defined $best;

  my %result = map { $_ => $best->{result}->{$_} } @{ $self->{cols} };

  return { SpliceAI2 => \%result }
    if $self->{config}->{output_format} eq 'json' || $self->{config}->{rest};

  return { map { $PREFIX . $_ => $result{$_} } keys %result };
}

sub _summary {
  my $score = $_[0]->{result}->{spliceai2_summary_score};
  return (defined $score && $score ne '') ? $score : -1;
}

sub parse_data {
  my ($self, $line) = @_;

  chomp $line;
  my @values = split /\t/, $line, -1;
  my $idx = $self->{col_index};

  my ($chr, $pos, $ref, $alt) = @values[ @$idx{qw(chrom pos ref alt)} ];
  my $strand = defined $idx->{strand} ? $values[ $idx->{strand} ] : '+';
  $strand = ($strand eq '-' || $strand eq '-1') ? -1 : 1;

  my $end = $pos + length($ref) - 1;

  # VCF-like adjustment of anchored indels for comparison with VEP
  if (length($ref) != length($alt) && substr($ref, 0, 1) eq substr($alt, 0, 1)) {
    $pos++;
    $ref = substr($ref, 1) || '-';
    $alt = substr($alt, 1) || '-';
  }

  my %result;
  @result{ @{ $self->{colnames} } } = @values;

  return {
    chr    => $chr,
    start  => $pos,
    end    => $end,
    ref    => $ref,
    alt    => $alt,
    strand => $strand,
    result => \%result,
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
