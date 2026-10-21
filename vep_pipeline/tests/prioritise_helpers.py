"""Shared builders for the prioritise tests (not a test module: no test_ prefix)."""
import subprocess
from pathlib import Path

from vep_pipeline.evidence import Evidence


def make_ev(**kw) -> Evidence:
    base = dict(symbol="G", consequence="missense_variant", hgvs="", popmax=0.0, clinvar=None, stars=0, revel=None,
                am=None, splice=None, cadd=None, nmd_escape=False, loeuf=None, phaplo=None, noncoding_signal=False)
    base.update(kw)
    return Evidence(**base)


FORMAT_TYPES = {"GT": ("1", "String"), "AD": ("R", "Integer"), "DP": ("1", "Integer"), "GQ": ("1", "Integer"),
                "PS": ("1", "Integer"), "PGT": ("1", "String"), "PID": ("1", "String")}
DEFAULT_SAMPLE = {"GT": "0/1", "AD": "10,10", "DP": "20", "GQ": "50"}


def write_vcf(path: Path, records: list, csq_fields: list, fmt=("GT", "AD", "DP", "GQ"), sample="S1") -> Path:
    """bgzipped single-sample VCF. record keys: pos, csq (list of {field: value}); optional chrom, ref, alt,
    filter, orig (EXO_ORIG) and FORMAT values by tag name."""
    head = ["##fileformat=VCFv4.2", "##contig=<ID=1,length=1000000>", "##contig=<ID=X,length=156040895>",
            '##FILTER=<ID=PASS,Description="All filters passed">', '##FILTER=<ID=LowQual,Description="low">',
            '##INFO=<ID=EXO_ORIG,Number=1,Type=String,Description="original record">',
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
            f'Format: {"|".join(csq_fields)}">']
    head += [f'##FORMAT=<ID={t},Number={FORMAT_TYPES[t][0]},Type={FORMAT_TYPES[t][1]},Description="{t}">' for t in fmt]
    head.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]))
    lines = []
    for r in records:
        csq = ",".join("|".join(e.get(f, "") for f in csq_fields) for e in r["csq"])
        info = f"CSQ={csq}" + (f";EXO_ORIG={r['orig']}" if r.get("orig") else "")
        values = ":".join(r.get(t, DEFAULT_SAMPLE.get(t, ".")) for t in fmt)
        lines.append("\t".join([r.get("chrom", "1"), str(r["pos"]), ".", r.get("ref", "A"), r.get("alt", "G"), ".",
                                r.get("filter", "PASS"), info, ":".join(fmt), values]))
    plain = Path(str(path).removesuffix(".gz"))
    plain.write_text("\n".join(head + lines) + "\n")
    subprocess.run(["bgzip", "-f", str(plain)], check=True)
    return Path(f"{plain}.gz")


from vep_pipeline.candidates import Call, Candidate


def make_call(pos, gt="0/1", chrom="1", ref="A", alt="G", orig=".", ad="10,10", dp="20", gq="50", ps=".",
              pgt=".", pid=".", filt="PASS") -> Call:
    return Call(chrom, pos, ref, alt, filt, orig, gt, ad, dp, gq, ps, pgt, pid)


def make_cand(pos, symbol="G", call=None, **ev_kw) -> Candidate:
    return Candidate(call or make_call(pos), make_ev(symbol=symbol, **ev_kw))


GENES_TSV_HEADER = ["#RANK", "ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "MOI", "P-VALUE", "EXOMISER_GENE_COMBINED_SCORE",
                    "EXOMISER_GENE_PHENO_SCORE"]


def write_genes_tsv(path: Path, rows: list) -> Path:
    """rows: (rank, symbol, moi, combined, pheno)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(GENES_TSV_HEADER)]
    lines += ["\t".join(map(str, (r, f"{g}_{m}", g, 1, m, 0.0, c, p))) for r, g, m, c, p in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


VARIANTS_TSV_HEADER = ["#RANK", "ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "MOI", "CONTIG", "START", "END", "REF", "ALT",
                       "ALL_PATH"]


def write_variants_tsv(path: Path, rows: list) -> Path:
    """Exomiser TSV_VARIANT stand-in. rows: (contig, start, ref, alt, all_path, moi='AD')."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(VARIANTS_TSV_HEADER)]
    for i, (contig, start, ref, alt, all_path, *moi) in enumerate(rows, 1):
        lines.append("\t".join(map(str, (i, f"{contig}-{start}-{ref}-{alt}_{moi[0] if moi else 'AD'}", "G", 1,
                                         moi[0] if moi else "AD", contig, start, start + len(ref) - 1, ref, alt,
                                         all_path))))
    path.write_text("\n".join(lines) + "\n")
    return path
