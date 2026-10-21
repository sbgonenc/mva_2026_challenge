#!/bin/bash

# Default values
reference=""
gff_annotation=""
known_variants_vcf=""
n_cores=1
tmp_path="/tmp/annotate_vcf/"

# Help function
usage() {
    echo "Usage: $0 -i <input.vcf> -r <reference.fasta> -g <annotation.gff> -k <known_variants.vcf> -o <output.vcf> [-c <cores>]"
    echo "Options:"
    echo "  -i    Input VCF file (required)"
    echo "  -r    Reference genome FASTA file (required)"
    echo "  -g    GFF annotation file (required)"
    echo "  -k    Known variants VCF file (required)"
    echo "  -o    Output VCF file (required)"
    echo "  -c    Number of CPU cores to use (default: 1)"
    echo "  -t    temp folder path default ${tmp_path}"
    echo "  -h    Display this help message"
    exit 1
}

# Parse command line arguments
while getopts "i:r:g:k:o:c:t:h" opt; do
    case $opt in
        i) in_vcf="$OPTARG";;
        r) reference="$OPTARG";;
        g) gff_annotation="$OPTARG";;
        k) known_variants_vcf="$OPTARG";;
        o) out_vcf="$OPTARG";;
        c) n_cores="$OPTARG";;
        t) tmp_path="$OPTARG";;
        h) usage;;
        \?) usage;;
    esac
done

# Check if required arguments are provided
if [ -z "$in_vcf" ] || [ -z "$reference" ] || [ -z "$gff_annotation" ] || [ -z "$known_variants_vcf" ] || [ -z "$out_vcf" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Check if input files exist
for file in "$in_vcf" "$reference" "$gff_annotation" "$known_variants_vcf"; do
    if [ ! -f "$file" ]; then
        echo "Error: File $file does not exist"
        exit 1
    fi
done

mkdir -p ${tmp_path}

### Check VCF
if [[ $in_vcf == *.gz ]]; then
    # If input is gzipped, use zcat to read it
    zgrep '#' $in_vcf > ${tmp_path}head.vcf
    zgrep -v '^#' $in_vcf  > ${tmp_path}tmp.vcf
else
    # If input is not gzipped, use regular grep
    grep '#' $in_vcf > ${tmp_path}head.vcf
    grep -v '^#' $in_vcf > ${tmp_path}tmp.vcf
fi

# Check and convert chromosome names if they start with 'chr'
sed 's/ID=chr/ID=/' ${tmp_path}head.vcf > ${tmp_path}no_chrom_tmp.vcf
sed 's/^chr//' ${tmp_path}tmp.vcf | sort -V -k1,1 -k2,2 >> ${tmp_path}no_chrom_tmp.vcf


bgzip ${tmp_path}no_chrom_tmp.vcf
bcftools index --tbi ${tmp_path}no_chrom_tmp.vcf

bcftools norm \
            --atomize \
            --check-ref ws \
            --fasta-ref $reference \
            -m - \
            --gff-annot $gff_annotation \
            -O v ${tmp_path}no_chrom_tmp.vcf.gz \
            -o ${tmp_path}norm_tmp.vcf \
            --threads ${n_cores} \
    && \
bcftools csq --fasta-ref $reference \
                 --gff-annot $gff_annotation \
                 --phase s \
                 --unify-chr-names 1 \
                 --threads ${n_cores} \
                 --regions-overlap 2 \
                 -O z ${tmp_path}norm_tmp.vcf \
                 -o ${tmp_path}csq_tmp.vcf.gz \
                 --write-index=tbi \
    && \
bcftools annotate --annotations $known_variants_vcf \
    		      --columns INFO \
    		      -O z ${tmp_path}csq_tmp.vcf.gz \
                  --write-index=tbi \
    		      -o ${out_vcf}



# Clean up temporary files
rm -fr ${tmp_path}