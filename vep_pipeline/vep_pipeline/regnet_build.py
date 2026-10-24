"""Compile the scDoRI eGRN/GRN output and the S-LDSC M1 results into the small TSVs `regnet` reads.

The eGRN edge lists are huge (about 112 M lines in total), so they are streamed once and only edges on a peak that is
linked to a gene are kept. Gene names may contain underscores, so a gene is parsed from the links file name by
stripping the fixed suffix, never by splitting on `_`. Stdlib only.
"""
from __future__ import annotations
import csv
import datetime
import json
import re
import statistics
from pathlib import Path

from .config import ConfigError

LINK_SUFFIX = "_enhancers_peak_links_.tsv"
EGRN_SUFFIX = "_eGRN_edge_list.tsv"
GRN_SUFFIX = "_GRN_edge_list.tsv"
PEAK = re.compile(r"^(chr[0-9XYM]+):(\d+)-(\d+)$")
PREFIX = {"M1": "TF_", "M1bin": "TFbin_"}
SLDSC_COLS = ("Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value")
SCDORI_DIRS = ("gene_peak_interactions", "eGRN_act_peak_gene", "GRN_act_TF_gene")


def bh(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values (monotone, capped at 1), in input order."""
    n = len(pvals)
    order = sorted(range(n), key=pvals.__getitem__)
    out = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running = min(running, pvals[i] * n / rank)
        out[i] = running
    return out


def _rows(path):
    """Yield the header, then the data rows, of a tab-separated file."""
    with open(path, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        yield next(r)
        yield from r


def _need_dir(path: Path, hint: str = "") -> Path:
    if not path.is_dir():
        raise ConfigError(f"{path}: directory missing{hint}")
    return path


def read_links(scdori: Path) -> dict:
    """{peak: {gene: max link score over cell types}} from gene_peak_interactions."""
    out = {}
    for p in sorted((scdori / "gene_peak_interactions").glob(f"*{LINK_SUFFIX}")):
        gene = p.name[: -len(LINK_SUFFIX)]
        rows = _rows(p)
        next(rows)
        for row in rows:
            if not PEAK.match(row[0]):
                raise ConfigError(f"{p}: malformed peak name {row[0]!r} (expected chrN:start-end)")
            vals = [float(x) for x in row[1:] if x not in ("", "nan")]
            out.setdefault(row[0], {})[gene] = max(vals) if vals else 0.0
    return out


def read_peak_tfs(scdori: Path, peaks: set) -> dict:
    """{peak: {tf: [celltypes, max_weight]}} over all eGRN edge lists, linked peaks only."""
    out = {}
    for p in sorted((scdori / "eGRN_act_peak_gene").glob(f"*{EGRN_SUFFIX}")):
        ct = p.name[: -len(EGRN_SUFFIX)]
        rows = _rows(p)
        next(rows)
        for tf, peak, w in rows:
            if peak in peaks:
                rec = out.setdefault(peak, {}).setdefault(tf, [set(), 0.0])
                rec[0].add(ct)
                rec[1] = max(rec[1], float(w))
    return out


def read_regulons(scdori: Path) -> dict:
    """{tf: {gene: set of cell types}} over all TF->gene GRN edge lists."""
    out = {}
    for p in sorted((scdori / "GRN_act_TF_gene").glob(f"*{GRN_SUFFIX}")):
        ct = p.name[: -len(GRN_SUFFIX)]
        rows = _rows(p)
        next(rows)
        for tf, gene, _w in rows:
            out.setdefault(tf, {}).setdefault(gene, set()).add(ct)
    return out


def gene_universe(scdori: Path) -> list:
    """Sorted union of the first column (gene_symbol) of the *_GRN.tsv matrices."""
    genes = set()
    for p in sorted((scdori / "GRN_act_TF_gene").glob("*_GRN.tsv")):
        rows = _rows(p)
        next(rows)
        genes.update(row[0] for row in rows)
    return sorted(genes)


def read_sldsc(dirs: list, model: str) -> list:
    """[(trait, source, tf, coef, se, p, fdr)] for traits listed in full/passing.txt, BH within each trait.

    A line of passing.txt with a comma is an rg pair and is ignored. A trait seen in an earlier dir is not read from a
    later one. `source` is the directory name."""
    prefix = PREFIX[model]
    out, seen = [], set()
    for d in dirs:
        passing = [ln.strip() for ln in (d / "full" / "passing.txt").read_text().splitlines()]
        for trait in sorted({t for t in passing if t and "," not in t} - seen):
            res = d / "h2_cts" / model / f"{trait}.cell_type_results.txt"
            if not res.is_file():
                continue
            with open(res, newline="") as fh:
                rd = csv.DictReader(fh, delimiter="\t")
                missing = [c for c in SLDSC_COLS if c not in (rd.fieldnames or ())]
                if missing:
                    raise ConfigError(f"{res}: columns {missing} missing")
                rows = [r for r in rd if r["Name"].startswith(prefix)]
            seen.add(trait)
            fdr = bh([float(r["Coefficient_P_value"]) for r in rows])
            out += [(trait, d.name, r["Name"][len(prefix):], float(r["Coefficient"]),
                     float(r["Coefficient_std_error"]), float(r["Coefficient_P_value"]), q)
                    for r, q in zip(rows, fdr)]
    return out


def read_labels(path: Path | None) -> dict:
    """{trait: label} from a `trait\\tlabel` TSV; empty when path is None."""
    if path is None:
        return {}
    if not path.is_file():
        raise ConfigError(f"trait label file missing: {path} (pass an existing trait<TAB>label TSV or omit it)")
    with open(path, newline="") as fh:
        return {r["trait"]: r["label"] for r in csv.DictReader(fh, delimiter="\t")}


def _write(path: Path, header: list, rows) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def build(scdori: Path, sldsc_dirs: list, model: str, out: Path, labels: Path | None) -> dict:
    """Write linked_peaks/peak_tfs/regulons/sldsc/traits TSVs, gene_universe.txt and build.json into `out`.

    Returns the counts: linked_peaks and genes_linked are unique peaks and genes; the rest are row counts."""
    if model not in PREFIX:
        raise ConfigError(f"unknown S-LDSC model {model!r} (expected one of {sorted(PREFIX)})")
    scdori = _need_dir(Path(scdori), " (set the scDoRI output directory)")
    for sub in SCDORI_DIRS:
        _need_dir(scdori / sub)
    sldsc_dirs = [Path(d) for d in sldsc_dirs]
    for d in sldsc_dirs:
        _need_dir(d, " (set the S-LDSC result directory)")
        _need_dir(d / "h2_cts" / model)
        if not (d / "full" / "passing.txt").is_file():
            raise ConfigError(f"{d / 'full' / 'passing.txt'}: file missing")
    label_map = read_labels(None if labels is None else Path(labels))
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)

    links = read_links(scdori)
    link_rows = []
    for peak, genes in links.items():
        chrom, start, end = PEAK.match(peak).groups()
        link_rows += [(peak, chrom, int(start), int(end), g, f"{v:.6g}") for g, v in genes.items()]
    link_rows.sort(key=lambda r: (r[1], r[2], r[3], r[4]))
    _write(out / "linked_peaks.tsv", ["peak", "chrom", "start", "end", "gene", "link_max"], link_rows)

    peak_tfs = read_peak_tfs(scdori, set(links))
    pt_rows = sorted((peak, tf, len(cts), f"{w:.6g}", ";".join(sorted(cts)))
                     for peak, tfs in peak_tfs.items() for tf, (cts, w) in tfs.items())
    _write(out / "peak_tfs.tsv", ["peak", "tf", "n_celltypes", "max_weight", "celltypes"], pt_rows)

    regulons = read_regulons(scdori)
    reg_rows = sorted((tf, g, len(cts)) for tf, genes in regulons.items() for g, cts in genes.items())
    _write(out / "regulons.tsv", ["tf", "gene", "n_celltypes"], reg_rows)

    universe = gene_universe(scdori)
    (out / "gene_universe.txt").write_text("".join(f"{g}\n" for g in universe))

    sl = read_sldsc(sldsc_dirs, model)
    _write(out / "sldsc.tsv", ["trait", "source", "tf", "coef", "se", "p", "fdr"],
           [(t, s, tf, f"{c:.6g}", f"{se:.6g}", f"{p:.6g}", f"{q:.6g}") for t, s, tf, c, se, p, q in sl])
    n_tf, source = {}, {}
    for t, s, *_ in sl:
        n_tf[t] = n_tf.get(t, 0) + 1
        source[t] = s
    _write(out / "traits.tsv", ["trait", "source", "label", "n_tf"],
           [(t, source[t], label_map.get(t, ""), n_tf[t]) for t in sorted(source)])

    counts = {"linked_peaks": len(links), "genes_linked": len({r[4] for r in link_rows}),
              "peak_tf_edges": len(pt_rows), "regulon_edges": len(reg_rows), "universe_genes": len(universe),
              "sldsc_rows": len(sl), "traits": len(source)}
    sizes = [len(g) for g in regulons.values()] or [0]
    meta = {"utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "scdori": str(scdori), "sldsc_dirs": [str(d) for d in sldsc_dirs], "model": model,
            "labels": None if labels is None else str(labels),
            "files": {"gene_peak_interactions": len(list((scdori / "gene_peak_interactions").glob(f"*{LINK_SUFFIX}"))),
                      "eGRN_edge_lists": len(list((scdori / "eGRN_act_peak_gene").glob(f"*{EGRN_SUFFIX}"))),
                      "GRN_edge_lists": len(list((scdori / "GRN_act_TF_gene").glob(f"*{GRN_SUFFIX}")))},
            "counts": counts,
            "regulon_size": {"min": min(sizes), "median": statistics.median(sizes), "max": max(sizes)}}
    (out / "build.json").write_text(json.dumps(meta, indent=2) + "\n")
    return counts
