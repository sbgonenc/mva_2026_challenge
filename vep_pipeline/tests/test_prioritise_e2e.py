import argparse
import csv

import pytest
from prioritise_helpers import write_genes_tsv, write_variants_tsv, write_vcf

from vep_pipeline import cli, prioritise as pr
from vep_pipeline.rank_settings import prioritise_settings

CSQ = ["Allele", "Consequence", "SYMBOL", "Feature", "MANE_SELECT", "HGVSp", "REVEL", "NMD", "ClinVar_CLNSIG",
       "ClinVar_CLNREVSTAT", "gnomADg_NFE_AF", "gnomADg_FIN_AF"]


def e(symbol, consequence="missense_variant", **kw):
    return {"SYMBOL": symbol, "Consequence": consequence, "Feature": f"ENST_{symbol}", "MANE_SELECT": "NM_1", **kw}


RECORDS = [
    {"pos": 1000, "ref": "C", "alt": "T", "csq": [e("GENEAR", "stop_gained")]},
    {"pos": 2000, "csq": [e("GENEAR", REVEL="0.95")]},
    {"pos": 5000, "orig": "1|5000|AC|GT|1", "csq": [e("GENEAR2", REVEL="0.95")]},
    {"pos": 5001, "ref": "C", "alt": "T", "orig": "1|5000|AC|GT|1", "csq": [e("GENEAR2", REVEL="0.95")]},
    {"pos": 8000, "alt": "C", "orig": "1|8000|A|C,G|1", "csq": [e("GENEAR3", REVEL="0.95")]},
    {"pos": 8000, "alt": "G", "orig": "1|8000|A|C,G|2", "csq": [e("GENEAR3", REVEL="0.95")]},
    {"pos": 9000, "csq": [e("TP53", REVEL="0.1", ClinVar_CLNSIG="Pathogenic",
                            ClinVar_CLNREVSTAT="criteria_provided&_multiple_submitters&_no_conflicts")]},
    {"pos": 10000, "csq": [e("GENEAD", REVEL="0.80")]},
    {"pos": 11000, "csq": [e("GENEAD2", REVEL="0.99", gnomADg_NFE_AF="0.2")]},
    {"pos": 12000, "filter": "LowQual", "csq": [e("GENEAD3", "stop_gained")]},
    {"pos": 13000, "csq": [e("GENEAD4", REVEL="0.99", gnomADg_FIN_AF="0.3")]},
]


@pytest.fixture
def e2e(minimal_cfg, tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    write_vcf(out / "S1.vep.vcf.gz", RECORDS, CSQ)
    g2p = tmp_path / "g2p.csv"
    g2p.write_text('"gene symbol","allelic requirement","confidence","molecular mechanism"\n'
                   + "".join(f"{g},biallelic_autosomal,definitive,loss of function\n"
                             for g in ("GENEAR", "GENEAR2", "GENEAR3"))
                   + "GENEAD,monoallelic_autosomal,definitive,loss of function\n"
                   + "TP53,monoallelic_autosomal,definitive,loss of function\n")
    res = out / "exomiser"
    pheno = [("GENEAR3", 0.95), ("GENEAR", 0.90), ("GENEAR2", 0.85), ("FILL1", 0.80), ("FILL2", 0.79),
             ("GENEAD", 0.70), ("TP53", 0.40)]
    write_genes_tsv(res / "phenotype" / "S1.genes.tsv", [(i, g, "ANY", p, p) for i, (g, p) in enumerate(pheno, 1)])
    for sub in ("exome", "genome"):
        write_genes_tsv(res / sub / "S1.genes.tsv", [(1, "GENEAR3", "AR", 0.9, 0.95)])
    write_variants_tsv(res / "genome" / "S1.variants.tsv", [("1", 1000, "C", "T", "CADD=0.9,REMM=0.5")])
    minimal_cfg["exomiser"] = {"hpo": ["HP:0001250"], "sex": "MALE"}
    minimal_cfg["prioritise"] = {"g2p": {"file": str(g2p)}, "pheno_rank": {"strong": 3, "moderate": 6}}
    return minimal_cfg


def test_end_to_end_csv(e2e):
    out, cand = pr.run(e2e, prioritise_settings(e2e), sensitivity=True)
    rows = list(csv.reader(open(out)))
    assert tuple(rows[0]) == pr.CSV_HEADER
    body = rows[1:]
    genes = [r[11].split("; ")[1] if r[10] == "primary" else "TP53" for r in body]
    assert genes == ["GENEAR3", "GENEAR", "GENEAD", "GENEAR2", "GENEAD4", "TP53"]
    assert body[0][1:9] == ["chr1", "8000", "A", "C", "chr1", "8000", "A", "G"]      # trans pair, original record
    assert body[3][1:9] == ["chr1", "5000", "AC", "GT", "", "", "", ""]              # atomised MNV = one allele
    assert sum(float(r[9]) for r in body if r[10] == "primary") == pytest.approx(0.95, abs=0.005)
    assert body[5][10] == "secondary" and body[5][9] == "0.300"
    assert body[5][11].startswith("Incidental finding unrelated to primary phenotype")
    assert "phase=unknown" in body[1][11]
    assert cand.is_file() and out.with_suffix(".sensitivity.tsv").is_file()
    assert out.with_suffix(".excluded.tsv").read_text().splitlines()[0] == "variant\tregion"
    scen = {l.split("\t")[0] for l in out.with_suffix(".sensitivity.tsv").read_text().splitlines()[1:]}
    assert scen == {name for name, _ in pr.SCENARIOS}


def test_run_is_deterministic(e2e):
    out, _ = pr.run(e2e, prioritise_settings(e2e))
    first = out.read_text()
    out, _ = pr.run(e2e, prioritise_settings(e2e))
    assert out.read_text() == first


def test_exomiser_off_marks_phenotype_none(e2e):
    e2e["exomiser"] = {"enabled": False}
    out, _ = pr.run(e2e, prioritise_settings(e2e))
    body = list(csv.reader(open(out)))[1:]
    assert all("phenotype=none" in r[11] for r in body)


def test_cli_slurm_renders_job(e2e, capsys):
    args = argparse.Namespace(top_n=5, output=None, sensitivity=True, slurm=True, submit=False)
    assert cli.prioritise_cmd(e2e, args) == 0
    script = capsys.readouterr().out.split()[2]
    text = open(script).read()
    assert "python -m vep_pipeline.cli prioritise -c" in text and "--top-n=5 --sensitivity" in text


def test_sensitivity_scenarios_bracket_p_trans():
    scen = dict(pr.SCENARIOS)
    assert "p_trans_1" not in scen
    assert scen["p_trans_0.15"] == {"p_trans_unphased": 0.15} and scen["p_trans_0.9"] == {"p_trans_unphased": 0.9}


def test_sensitivity_scenarios_isolate_panel_alternatives_and_read_backed_phase():
    scen = dict(pr.SCENARIOS)
    assert scen["no_gene_panel"] == {"gene_panel": {"name": "", "genes": [], "phenotype_floor": None,
                                                    "deep_second_hit": False}}
    assert scen["one_finding_per_gene"] == {"alternative_pairs": {"max_per_gene": 1, "min_phenotype": "strong"}}
    assert scen["no_read_backed_phase"] == {"phase_sources": ["multiallelic", "ps", "pid"]}


def _enable_phasing(cfg, tmp_path):
    for f in ("s.bam", "s.bam.bai", "ref.fa", "ref.fa.fai"):
        (tmp_path / f).write_text("")
    cfg["phasing"] = {"enabled": True, "bam": str(tmp_path / "s.bam"), "reference": str(tmp_path / "ref.fa"),
                      "env": "whatshap"}


def test_phasing_enabled_without_sidecar_is_an_error(e2e, tmp_path):
    _enable_phasing(e2e, tmp_path)
    with pytest.raises(pr.ConfigError, match="run `vep-pipeline phase --submit` first"):
        pr.run(e2e, prioritise_settings(e2e))


def test_readback_sidecar_phases_a_pair_and_shows_the_source(e2e, tmp_path):
    _enable_phasing(e2e, tmp_path)
    side = tmp_path / "out" / "phasing" / "S1.readback.tsv"
    side.parent.mkdir()
    side.write_text("key\tgt\tps\n1:1000:C>T\t0|1\t1000\n1:2000:A>G\t1|0\t1000\n")
    out, cand = pr.run(e2e, prioritise_settings(e2e), sensitivity=True)
    rows = {r[11].split("; ")[1]: r for r in list(csv.reader(open(out)))[1:] if r[10] == "primary"}
    assert "phase=trans (whatshap)" in rows["GENEAR"][11]
    cand_rows = [l.split("\t") for l in cand.read_text().splitlines()[1:]]
    ar = [r for r in cand_rows if r[0] == "GENEAR" and r[10]]
    assert ar and {r[8] for r in ar} == {"trans"} and "phase=trans (whatshap)" in ar[0][11]
    e2e["prioritise"]["phase_sources"] = ["ps"]
    out2, _ = pr.run(e2e, prioritise_settings(e2e))
    assert "phase=unknown" in [r for r in csv.reader(open(out2)) if "GENEAR;" in r[11]][0][11]


# ----------------------------------------------------------------------------- gene panel with alternative pairs

# AR panel gene: a stop (null), a 5'UTR deep-search candidate, a REVEL 0.80 missense and a splice-region candidate
# phased cis to the missense by PS. Not scored by Exomiser ('unknown'); the panel floor makes it 'strong'.
PANEL_RECORDS = [
    {"pos": 20000, "ref": "C", "alt": "T", "csq": [e("GENEMVA", "stop_gained", HGVSp="ENSP_M:p.Arg100Ter")]},
    {"pos": 20500, "csq": [e("GENEMVA", "5_prime_UTR_variant")]},
    {"pos": 21000, "GT": "0|1", "PS": "21000", "csq": [e("GENEMVA", REVEL="0.80", HGVSp="ENSP_M:p.Gly200Asp")]},
    {"pos": 21500, "GT": "0|1", "PS": "21000", "csq": [e("GENEMVA", "splice_region_variant&intron_variant")]},
]
CH_HEADER = ("gene\tvariant_1\tvariant_2\thgvs_1\thgvs_2\tlevel_1\tlevel_2\treasons_1\treasons_2\tphase\t"
             "phase_sources\tdistance\tcollapsed_into\tfit\ttier\tepcr_before_cap\treportable\tcsv_row")


@pytest.fixture
def e2e_panel(e2e, tmp_path):
    write_vcf(tmp_path / "out" / "S1.vep.vcf.gz", RECORDS + PANEL_RECORDS, CSQ, fmt=("GT", "AD", "DP", "GQ", "PS"))
    with open(tmp_path / "g2p.csv", "a") as fh:
        fh.write("GENEMVA,biallelic_autosomal,definitive,loss of function\n")
    e2e["prioritise"]["gene_panel"] = {"name": "MVA", "genes": ["GENEMVA"], "phenotype_floor": "strong",
                                       "deep_second_hit": True}
    e2e["prioritise"]["alternative_pairs"] = {"max_per_gene": 3, "min_phenotype": "strong"}
    return e2e


def _key(cols):
    chrom, pos, ref, alt = cols
    return f"{chrom.removeprefix('chr')}:{pos}:{ref}>{alt}" if pos else ""


def test_panel_gene_alternative_pairs_share_its_budget(e2e_panel):
    out, cand = pr.run(e2e_panel, prioritise_settings(e2e_panel))
    body = list(csv.reader(open(out)))[1:]
    genes = [r[11].split("; ")[1] if r[10] == "primary" else "TP53" for r in body]
    assert genes == ["GENEAR3", "GENEAR", "GENEAD", "GENEMVA", "GENEAR2", "GENEMVA", "GENEMVA", "GENEAD4", "TP53"]
    mva = [r for r in body if "; GENEMVA; " in r[11]]
    assert [r[1:9] for r in mva] == [["chr1", "20000", "C", "T", "chr1", "21000", "A", "G"],    # best pair first
                                     ["chr1", "20000", "C", "T", "chr1", "20500", "A", "G"],
                                     ["chr1", "20500", "A", "G", "chr1", "21000", "A", "G"]]
    assert [r[11].split("; ")[0] for r in mva] == ["T2", "T4", "T4"]
    assert "recessive/full" in mva[0][11] and all("recessive/panel_partner" in r[11] for r in mva[1:])
    assert all("panel=MVA (phenotype floor)" in r[11] and "phase=unknown" in r[11] for r in mva)
    assert float(mva[0][9]) / float(mva[1][9]) == pytest.approx(6, rel=0.02)            # 0.45 : 0.075 after the cap
    assert sum(float(r[9]) for r in body if r[10] == "primary") == pytest.approx(0.95, abs=0.005)
    shared = [float(r[3]) for r in csv.reader(open(cand), delimiter="\t") if r[0] == "GENEMVA"]
    assert shared == pytest.approx([0.225, 0.0375, 0.0375], abs=6e-4)
    assert sum(shared) == pytest.approx(0.60 * 0.5, abs=2e-3)       # B = T2 0.60, all x p_trans 0.5 (3-decimal output)


def test_compound_het_tsv(e2e_panel):
    out, _ = pr.run(e2e_panel, prioritise_settings(e2e_panel))
    body = list(csv.reader(open(out)))[1:]
    lines = out.with_suffix(".compound_het.tsv").read_text().splitlines()
    assert lines[0] == CH_HEADER
    rows = [dict(zip(CH_HEADER.split("\t"), l.split("\t"))) for l in lines[1:]]
    assert [r["gene"] for r in rows] == ["GENEAR", "GENEAR3"] + ["GENEMVA"] * 6
    mva = {(r["variant_1"].split(":")[1], r["variant_2"].split(":")[1]): r for r in rows if r["gene"] == "GENEMVA"}
    assert list(mva) == [("20000", "20500"), ("20000", "21000"), ("20000", "21500"), ("20500", "21000"),
                         ("20500", "21500"), ("21000", "21500")]
    assert {k: (r["reportable"], r["csv_row"]) for k, r in mva.items()} == {
        ("20000", "20500"): ("yes", "6"), ("20000", "21000"): ("yes", "4"), ("20000", "21500"): ("no (collapsed)", ""),
        ("20500", "21000"): ("yes", "7"), ("20500", "21500"): ("no (two second-hit-only alleles)", ""),
        ("21000", "21500"): ("no (cis)", "")}
    best = mva[("20000", "21000")]
    assert (best["hgvs_1"], best["hgvs_2"], best["level_1"], best["level_2"]) == (
        "ENSP_M:p.Arg100Ter", "ENSP_M:p.Gly200Asp", "strong", "moderate")
    assert (best["fit"], best["tier"], best["epcr_before_cap"], best["phase"], best["distance"]) == (
        "full", "T2", "0.225", "unknown", "1000")
    cis = mva[("21000", "21500")]
    assert (cis["phase"], cis["phase_sources"], cis["fit"], cis["epcr_before_cap"]) == ("cis", "ps", "", "")
    assert mva[("20000", "21500")]["collapsed_into"] == "1:20000:C>T|1:21000:A>G"
    assert mva[("20500", "21500")]["reasons_2"] == "panel candidate (splice)"
    for r in rows:                                                       # csv_row points at the same two variants
        if r["csv_row"]:
            row = body[int(r["csv_row"]) - 1]
            assert {_key(row[1:5]), _key(row[5:9])} == {r["variant_1"], r["variant_2"]}
    assert all(not r["csv_row"] for r in rows if r["reportable"] != "yes")


# ----------------------------------------------------------------------------- regulatory network (regnet)

REG_RECORD = {"pos": 30100, "csq": [{"Consequence": "intergenic_variant"}]}     # in a peak linked to GENEAR


def _regnet_build(tmp_path, gene="GENEAR"):
    from test_regnet_build import _w, make_sldsc
    from vep_pipeline.regnet_build import build
    sc = tmp_path / "scdori"
    for sub in ("gene_peak_interactions", "eGRN_act_peak_gene", "GRN_act_TF_gene"):
        (sc / sub).mkdir(parents=True)
    _w(sc / f"gene_peak_interactions/{gene}_enhancers_peak_links_.tsv", [["peak_name2", "Ct a"], ["chr1:30000-30500", 0.4]])
    _w(sc / "eGRN_act_peak_gene/Ct_a_eGRN_edge_list.tsv", [["source", "target", "weight"], ["TF1", "chr1:30000-30500", 0.1]])
    _w(sc / "GRN_act_TF_gene/Ct_a_GRN_edge_list.tsv", [["source", "target", "weight"], ["TF1", gene, 0.1]])
    _w(sc / "GRN_act_TF_gene/Ct_a_GRN.tsv", [["gene_symbol", "TF1"], [gene, 0.1], ["GENEAD", 0.0]])
    out = tmp_path / "regnet_build"
    build(sc, [make_sldsc(tmp_path, "srcA")], "M1", out, None)
    return out


def _outputs(out):
    return {p.name: p.read_bytes() for p in sorted(out.parent.iterdir()) if p.is_file() and p.name.startswith(out.stem)}


def test_regnet_writes_regulatory_tsv_and_summary(e2e, tmp_path):
    write_vcf(tmp_path / "out" / "S1.vep.vcf.gz", RECORDS + [REG_RECORD], CSQ)
    out, _ = pr.run(e2e, prioritise_settings(e2e))
    plain = _outputs(out)
    assert not out.with_suffix(".regulatory.tsv").exists()
    assert not out.with_suffix(".regulatory.summary.tsv").exists()
    for p in out.parent.iterdir():
        if p.is_file() and p.name.startswith(out.stem):
            p.unlink()
    e2e["regnet"] = {"enabled": True, "build_dir": str(_regnet_build(tmp_path)), "traits": ["T1"], "null_draws": 20}
    out, _ = pr.run(e2e, prioritise_settings(e2e))
    with open(out.with_suffix(".regulatory.tsv"), newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert [(r["chrom"], r["pos"], r["gene"], r["tier"], r["pheno"], r["sldsc_tfs"], r["v"], r["second_hit"])
            for r in rows] == [("chr1", "30100", "GENEAR", "R2", "strong", "TF1(T1:0.002)", "no", "no")]
    kv = dict(l.split("\t") for l in out.with_suffix(".regulatory.summary.tsv").read_text().splitlines())
    assert kv["hits"] == "1" and kv["rows_R2"] == "1" and kv["declared_traits"] == "T1:"
    assert kv["flagged_tfs_T1"] == "TF1" and kv["null_draws"] == "20"
    with_reg = _outputs(out)
    assert set(with_reg) - set(plain) == {out.with_suffix(".regulatory.tsv").name,
                                          out.with_suffix(".regulatory.summary.tsv").name}
    assert {k: v for k, v in with_reg.items() if k in plain} == plain     # the CSV and every TSV are unchanged


def test_regnet_disabled_block_changes_nothing(e2e, tmp_path):
    out, _ = pr.run(e2e, prioritise_settings(e2e))
    before = _outputs(out)
    e2e["regnet"] = {"enabled": False, "build_dir": str(tmp_path / "missing")}
    out, _ = pr.run(e2e, prioritise_settings(e2e))
    assert _outputs(out) == before and not any(".regulatory" in n for n in before)


# AR gene GENEREG (strong phenotype): a stop (strong null) and an intergenic CADD 26 variant in a peak linked to it.
# CADD 26 is >= 25.3, a non-coding CADD second hit on its own; the variant has no CSQ entry for GENEREG (link-only),
# so it must carry only the network reason and vanish under no_regulatory_network.
REG_CSQ = CSQ + ["CADD_PHRED"]
REG_SECOND = [{"pos": 30100, "csq": [{"Consequence": "intergenic_variant", "CADD_PHRED": "26"}]},
              {"pos": 31000, "ref": "C", "alt": "T", "csq": [e("GENEREG", "stop_gained", HGVSp="ENSP_R:p.Arg10Ter")]}]


@pytest.fixture
def e2e_reg(e2e, tmp_path):
    write_vcf(tmp_path / "out" / "S1.vep.vcf.gz", RECORDS + REG_SECOND, REG_CSQ)
    with open(tmp_path / "g2p.csv", "a") as fh:
        fh.write("GENEREG,biallelic_autosomal,definitive,loss of function\n")
    pheno = [("GENEAR3", 0.95), ("GENEREG", 0.92), ("GENEAR", 0.90), ("GENEAR2", 0.85), ("FILL1", 0.80),
             ("FILL2", 0.79), ("GENEAD", 0.70), ("TP53", 0.40)]
    write_genes_tsv(tmp_path / "out" / "exomiser" / "phenotype" / "S1.genes.tsv",
                    [(i, g, "ANY", p, p) for i, (g, p) in enumerate(pheno, 1)])
    e2e["prioritise"]["pheno_rank"] = {"strong": 4, "moderate": 7}
    return e2e


def _reg_on(cfg, tmp_path, **kw):
    cfg["regnet"] = {"enabled": True, "build_dir": str(_regnet_build(tmp_path, "GENEREG")), "traits": ["T1"],
                     "null_draws": 20, **kw}


def _tsv_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _sens(out):
    rows = [l.split("\t") for l in out.with_suffix(".sensitivity.tsv").read_text().splitlines()[1:]]
    by = {}
    for r in rows:
        by.setdefault(r[0], []).append(r[1:])
    return by


def test_regnet_second_hit_pairs_with_a_strong_null(e2e_reg, tmp_path):
    out, cand = pr.run(e2e_reg, prioritise_settings(e2e_reg), sensitivity=True)
    plain = _sens(out)
    assert [(r["tier"], r["variant_2"]) for r in _tsv_rows(cand) if r["gene"] == "GENEREG"] == [("T4", "")]
    _reg_on(e2e_reg, tmp_path)
    out, cand = pr.run(e2e_reg, prioritise_settings(e2e_reg), sensitivity=True)
    (reg,) = _tsv_rows(out.with_suffix(".regulatory.tsv"))
    assert (reg["gene"], reg["tier"], reg["second_hit"], reg["v"]) == ("GENEREG", "R1", "yes", "yes")
    (row,) = [r for r in _tsv_rows(cand) if r["gene"] == "GENEREG"]
    assert (row["tier"], row["epcr_before_cap"], row["fit"], row["level"], row["phase"]) == (
        "T3", "0.150", "full", "supporting", "unknown")                          # T3 0.30 x p_trans 0.5
    assert (row["variant_1"], row["variant_2"]) == ("1:30100:A>G", "1:31000:C>T")
    assert ("v1: chr1:30000-30500 het level=supporting [regulatory network: peak chr1:30000-30500 -> GENEREG; R1; "
            "S-LDSC TFs TF1 (second hit only)]") in row["notes"]
    sens = _sens(out)
    assert set(sens) == {name for name, _ in pr.SCENARIOS} | {"no_regulatory_network"}
    assert sens["no_regulatory_network"] == plain["base"]                        # the hook removed = the plain run
    assert sens["base"] != plain["base"]


def test_regnet_disabled_or_without_second_hit_keeps_every_output(e2e_reg, tmp_path):
    out, _ = pr.run(e2e_reg, prioritise_settings(e2e_reg), sensitivity=True)
    before = _outputs(out)
    assert "no_regulatory_network" not in _sens(out)
    _reg_on(e2e_reg, tmp_path, enabled=False)
    out, _ = pr.run(e2e_reg, prioritise_settings(e2e_reg), sensitivity=True)
    assert _outputs(out) == before                                               # byte-identical, sensitivity too
    e2e_reg["regnet"].update(enabled=True, second_hit=False)
    out, _ = pr.run(e2e_reg, prioritise_settings(e2e_reg), sensitivity=True)
    after = _outputs(out)
    assert {k: v for k, v in after.items() if ".regulatory" not in k} == before
    assert (out.with_suffix(".regulatory.tsv").name in after
            and "\tyes\t" in out.with_suffix(".regulatory.tsv").read_text())   # a second-hit row, not used
