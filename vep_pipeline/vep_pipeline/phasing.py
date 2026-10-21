"""Read-backed (WhatsHap) phase sidecar: `{phasing.results_dir}/{sample}.readback.tsv`, header `key<TAB>gt<TAB>ps`.

It holds only phased records (GT with `|`, PS set), keyed by the normalised VEP key (`jobs._key`); `prioritise` loads it
as the `whatshap` phase source.

The `phase` command is one SLURM job (`phase.sbatch`): `--prepare` finds the candidate-gene regions and extracts the het
calls in them, WhatsHap phases those against the BAM, `--finish` writes the sidecar. All work files live in
`phasing.results_dir`.
"""
from __future__ import annotations
import shlex
from itertools import combinations
from pathlib import Path

import yaml

from .candidates import build_alleles, collect
from .config import ConfigError
from .exomiser import load_remm
from .genedb import g2p_path, load_g2p
from .jobs import PKG_ROOT, _key, _out_dir, _prefix, _sbatch_header, _tool, final_vcf_path
from .rank_settings import exomiser_settings, phasing_settings, prioritise_settings

HEADER = "key\tgt\tps"


def readback_path(cfg: dict, ph: dict) -> Path:
    return Path(ph["results_dir"]) / f"{cfg['sample']}.readback.tsv"


def load_readback(path) -> dict:
    """{key: (gt, ps)} of the phased records in a read-back sidecar; unphased rows are skipped. A row without exactly
    three fields (a partial or foreign file) is a ConfigError."""
    out = {}
    with open(path) as fh:
        if fh.readline().rstrip("\n") != HEADER:
            raise ConfigError(f"{path}: expected the header {HEADER!r} (key, gt, ps tab-separated)")
        for n, line in enumerate(fh, 2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ConfigError(f"{path} line {n}: expected 3 tab-separated fields (key, gt, ps), got "
                                  f"{len(fields)}; re-run `vep-pipeline phase --submit`")
            key, gt, ps = fields
            if "|" in gt and ps not in ("", "."):
                out[key] = (gt, ps)
    return out


def _write_atomic(path: Path, text: str) -> None:
    """Write a temporary file next to `path`, then move it into place: `path` is never partial."""
    tmp = path.with_name(f"{path.name}.tmp")
    try:
        tmp.write_text(text)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


# ----------------------------------------------------------------------------- the `phase` job

def phase_regions(by_gene: dict, pad: int) -> list:
    """Merged [(chrom, start, end)] (1-based, inclusive): per gene with >= 2 het alleles, from the first allele start
    - pad to the last allele end + pad. Overlapping or touching regions of one contig are merged."""
    spans = []
    for alleles in by_gene.values():
        het = [a for a in alleles if a.zyg == "het"]
        if len(het) < 2:
            continue
        chrom = het[0].first.call.chrom
        start = min(a.first.call.pos for a in het)
        end = max(a.pieces[-1].call.pos + len(a.pieces[-1].call.ref) - 1 for a in het)
        spans.append((chrom, max(1, start - pad), end + pad))
    merged = []
    for chrom, start, end in sorted(spans):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2] + 1:
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((chrom, start, end))
    return merged


def write_regions_bed(regions, path):
    """BED (0-based start, end exclusive) of 1-based inclusive regions."""
    Path(path).write_text("".join(f"{c}\t{s - 1}\t{e}\n" for c, s, e in regions))


def bam_contigs(header_sam) -> list:
    """Contig names of the @SQ lines of `samtools view -H` output."""
    names = []
    for line in Path(header_sam).read_text().splitlines():
        if line.startswith("@SQ"):
            names += [f[3:] for f in line.split("\t") if f.startswith("SN:")]
    return names


def contig_map(bam_names, vcf_names) -> dict:
    """{vcf contig: bam contig}: the same name, `chr` added, `chr` stripped, or MT <-> chrM; else ConfigError."""
    have, out = set(bam_names), {}
    for name in vcf_names:
        stripped = name[3:] if name.startswith("chr") else name
        mito = {"MT": "chrM", "chrM": "MT"}.get(name)
        found = next((c for c in (name, f"chr{name}", stripped, mito) if c in have), None)
        if found is None:
            raise ConfigError(f"contig {name!r} of the candidate regions is not in the BAM (contigs: "
                              f"{', '.join(sorted(have)[:5])}, ...)")
        out[name] = found
    return out


def _work(ph) -> Path:
    return Path(ph["results_dir"])


def _plural(n, word) -> str:
    return f"{n} {word}" + ("" if n == 1 else "s")


def render_phase(cfg, ph) -> Path:
    """Freeze the config to {results_dir}/phase.config.yaml and write phase.sbatch, the job that runs the whole command."""
    wd = _work(ph)
    wd.mkdir(parents=True, exist_ok=True)
    Path(cfg["slurm"]["log_dir"]).mkdir(parents=True, exist_ok=True)
    cfg_path = wd / "phase.config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    q, tools, env = shlex.quote, shlex.quote(cfg["tools"]["env"]), shlex.quote(ph["env"])
    job_cfg = {**cfg, "slurm": {**cfg["slurm"], "mem": ph["mem"], "time": ph["time"]}}
    phase = f"micromamba run -n {tools} python -m vep_pipeline.cli phase"
    whatshap = (f"micromamba run -n {env} whatshap phase --reference {q(ph['reference'])} --ignore-read-groups "
                + ("--indels " if ph["indels"] else "")
                + f"-o {q(str(wd / 'phased.vcf.gz'))} {q(str(wd / 'subset.vcf.gz'))} {q(ph['bam'])}")
    script = wd / "phase.sbatch"
    stale = (readback_path(cfg, ph), wd / "phased.vcf.gz", wd / "phased.vcf.gz.tbi")   # an earlier run's outputs
    script.write_text(
        _sbatch_header(job_cfg, "phase", "phase_%j.out", 4)
        + f"export PYTHONPATH={q(str(PKG_ROOT))}${{PYTHONPATH:+:$PYTHONPATH}}\n"
        + f"rm -f -- {' '.join(q(str(p)) for p in stale)}\n"
        + f"micromamba run -n {env} samtools view -H {q(ph['bam'])} > {q(str(wd / 'bam_header.sam'))}\n"
        + f"{phase} --prepare -c {q(str(cfg_path))}\n"
        + f"if [[ -s {q(str(wd / 'regions.bed'))} ]]; then\n"
        + f"  {whatshap}\n"
        + f"  micromamba run -n {tools} bcftools index -t {q(str(wd / 'phased.vcf.gz'))}\n"
        + f"  {phase} --finish -c {q(str(cfg_path))}\n"
        + "fi\n")
    return script


def _write_map(path, pairs):
    Path(path).write_text("".join(f"{a}\t{b}\n" for a, b in sorted(pairs)))


def prepare(cfg, ph) -> None:
    """Regions of the candidate genes, region_alleles.tsv, the contig maps and subset.vcf.gz (het PASS calls inside the
    regions, contigs renamed to the BAM's). With no region: an empty sidecar and an empty regions.bed. The sidecar and
    phased.vcf.gz of an earlier run are removed first, so a job that fails later leaves no sidecar to read."""
    wd = _work(ph)
    for old in (readback_path(cfg, ph), wd / "phased.vcf.gz", wd / "phased.vcf.gz.tbi"):
        old.unlink(missing_ok=True)
    header = wd / "bam_header.sam"
    if not header.is_file():
        raise ConfigError(f"{header} not found: the job's `samtools view -H` step writes it")
    settings, ex = prioritise_settings(cfg), exomiser_settings(cfg)
    vcf = final_vcf_path(cfg)
    if not vcf.is_file():
        raise ConfigError(f"{vcf} not found: run `vep-pipeline run` first")
    genes = load_g2p(g2p_path(cfg, settings), settings["g2p"]["confidence"])
    remm = load_remm(cfg, ex) if ex["enabled"] else {}
    cands = collect(vcf, cfg["sample"], settings, genes, remm, [])
    by_gene = build_alleles(cands, genes, settings, ex["sex"])
    regions = phase_regions(by_gene, ph["pad"])
    write_regions_bed(regions, wd / "regions.bed")
    (wd / "region_alleles.tsv").write_text("".join(
        f"{g}\t{a.key}\n" for g, alleles in by_gene.items() if sum(a.zyg == "het" for a in alleles) > 1
        for a in alleles if a.zyg == "het"))
    if not regions:
        _write_atomic(readback_path(cfg, ph), HEADER + "\n")
        print("no gene has two het alleles: nothing to phase")
        return
    to_bam = contig_map(bam_contigs(header), sorted({r[0] for r in regions}))
    _write_map(wd / "contigs_to_bam.txt", to_bam.items())
    _write_map(wd / "contigs_from_bam.txt", [(b, v) for v, b in to_bam.items()])
    raw, subset = wd / "subset.raw.bcf", wd / "subset.vcf.gz"
    _tool(["bcftools", "view", "-R", wd / "regions.bed", "-s", cfg["sample"], "-f", "PASS,.", "-i", 'GT="het"',
           "-Ob", "-o", raw, _out_dir(cfg) / f"{_prefix(cfg)}.input.vcf.gz"])
    _tool(["bcftools", "annotate", "--rename-chrs", wd / "contigs_to_bam.txt", "-Oz", "-o", subset, raw])
    raw.unlink()
    _tool(["bcftools", "index", "-f", "-t", subset])
    print(f"{_plural(len(regions), 'region')} for phasing; het calls extracted to {subset}")


def finish(cfg, ph) -> Path:
    """Write the read-back sidecar from phased.vcf.gz (contigs mapped back to the VCF's names) and print the counts."""
    wd = _work(ph)
    phased = wd / "phased.vcf.gz"
    if not phased.is_file():
        raise ConfigError(f"{phased} not found: the WhatsHap step of the job writes it")
    back = dict(line.split("\t") for line in (wd / "contigs_from_bam.txt").read_text().splitlines())
    rows = []
    for line in _tool(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%PS]\n", phased]).splitlines():
        chrom, pos, ref, alt, gt, ps = line.split("\t")
        if "|" in gt and ps != ".":
            chrom = back.get(chrom, chrom)
            rows.append((chrom, _key(chrom, pos, ref, alt), gt, ps))
    side = readback_path(cfg, ph)
    _write_atomic(side, HEADER + "\n" + "".join(f"{key}\t{gt}\t{ps}\n" for _c, key, gt, ps in rows))
    block = {key: (chrom, ps) for chrom, key, _gt, ps in rows}
    keys_by_gene = {}
    for line in (wd / "region_alleles.tsv").read_text().splitlines():
        gene, key = line.split("\t")
        keys_by_gene.setdefault(gene, []).append(key)
    pairs = sum(1 for keys in keys_by_gene.values() for a, b in combinations(keys, 2)
                if a in block and b in block and block[a] == block[b])
    print(f"{_plural(len(rows), 'phased record')} in {_plural(len(set(block.values())), 'block')}; "
          f"{_plural(pairs, 'candidate pair')} with both alleles in one block; sidecar {side}")
    return side
