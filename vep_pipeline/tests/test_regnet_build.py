import csv
import json

import pytest

from vep_pipeline.config import ConfigError
from vep_pipeline.regnet_build import bh, build


def _w(path, rows):
    path.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")


def make_scdori(root):
    (root / "gene_peak_interactions").mkdir(parents=True)
    (root / "eGRN_act_peak_gene").mkdir()
    (root / "GRN_act_TF_gene").mkdir()
    _w(root / "gene_peak_interactions/GENEA_enhancers_peak_links_.tsv",
       [["peak_name2", "Ct a", "Ct b"], ["chr1:100-600", 0.5, -0.2], ["chr1:1000-1500", -0.1, 1.2]])
    _w(root / "gene_peak_interactions/GENE_B_enhancers_peak_links_.tsv",
       [["peak_name2", "Ct a", "Ct b"], ["chr1:1000-1500", 0.3, 0.0]])
    _w(root / "eGRN_act_peak_gene/Ct_a_eGRN_edge_list.tsv",
       [["source", "target", "weight"], ["TF1", "chr1:100-600", 0.01], ["TF2", "chr2:5-505", 0.2]])
    _w(root / "eGRN_act_peak_gene/Ct_b_eGRN_edge_list.tsv",
       [["source", "target", "weight"], ["TF1", "chr1:100-600", 0.03]])
    _w(root / "GRN_act_TF_gene/Ct_a_GRN_edge_list.tsv", [["source", "target", "weight"], ["TF1", "GENEA", 0.1]])
    _w(root / "GRN_act_TF_gene/Ct_a_GRN.tsv", [["gene_symbol", "TF1"], ["GENEA", 0.1], ["GENE_B", 0], ["GENEC", 0]])


def make_sldsc(root, name):
    d = root / name
    (d / "full").mkdir(parents=True)
    (d / "h2_cts/M1").mkdir(parents=True)
    (d / "full/passing.txt").write_text("T1,T2\nT1\nT2\n")
    head = ["Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value"]
    _w(d / "h2_cts/M1/T1.cell_type_results.txt",
       [head, ["TF_TF1", 1e-9, 2e-9, "1e-3"], ["TF_TF2", 1e-10, 3e-9, "5e-1"]])
    _w(d / "h2_cts/M1/T2.cell_type_results.txt", [head, ["TF_TF1", 5e-10, 1e-9, "9e-1"]])
    _w(d / "h2_cts/M1/T3.cell_type_results.txt", [head, ["TF_TF1", 5e-10, 1e-9, "1e-5"]])
    (d / "h2_cts/M1/T1.runlog").write_text("log\n")
    return d


def tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@pytest.fixture
def built(tmp_path):
    scdori = tmp_path / "scdori"
    make_scdori(scdori)
    a = make_sldsc(tmp_path, "srcA")
    out = tmp_path / "out"
    counts = build(scdori, [a], "M1", out, None)
    return out, counts, tmp_path


def test_bh():
    assert bh([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.04, 0.04])
    assert bh([]) == []
    assert bh([0.9, 0.9]) == pytest.approx([0.9, 0.9])
    assert bh([0.5, 0.0001])[0] == pytest.approx(0.5)


def test_links_keep_gene_names_with_underscores(built):
    rows = tsv(built[0] / "linked_peaks.tsv")
    assert len(rows) == 3
    got = {(r["peak"], r["gene"]): r for r in rows}
    assert float(got[("chr1:1000-1500", "GENE_B")]["link_max"]) == pytest.approx(0.3)
    r = got[("chr1:1000-1500", "GENEA")]
    assert float(r["link_max"]) == pytest.approx(1.2)
    assert (r["chrom"], r["start"], r["end"]) == ("chr1", "1000", "1500")
    assert float(got[("chr1:100-600", "GENEA")]["link_max"]) == pytest.approx(0.5)


def test_peak_tfs_only_for_linked_peaks_and_counts_celltypes(built):
    rows = tsv(built[0] / "peak_tfs.tsv")
    assert [(r["peak"], r["tf"]) for r in rows] == [("chr1:100-600", "TF1")]
    assert rows[0]["n_celltypes"] == "2"
    assert float(rows[0]["max_weight"]) == pytest.approx(0.03)
    assert rows[0]["celltypes"] == "Ct_a;Ct_b"


def test_regulons_and_universe(built):
    out, counts, _ = built
    assert tsv(out / "regulons.tsv") == [{"tf": "TF1", "gene": "GENEA", "n_celltypes": "1"}]
    assert (out / "gene_universe.txt").read_text().split("\n") == ["GENEA", "GENEC", "GENE_B", ""]
    assert counts == {"linked_peaks": 2, "genes_linked": 2, "peak_tf_edges": 1, "regulon_edges": 1,
                      "universe_genes": 3, "sldsc_rows": 3, "traits": 2}
    meta = json.loads((out / "build.json").read_text())
    assert meta["counts"] == counts
    assert meta["regulon_size"] == {"min": 1, "median": 1, "max": 1}
    assert meta["model"] == "M1" and meta["utc"].endswith("Z")


def test_sldsc_passing_only_bh_within_trait(built):
    rows = tsv(built[0] / "sldsc.tsv")
    assert {r["trait"] for r in rows} == {"T1", "T2"}
    by = {(r["trait"], r["tf"]): r for r in rows}
    assert set(by) == {("T1", "TF1"), ("T1", "TF2"), ("T2", "TF1")}
    assert by[("T1", "TF1")]["source"] == "srcA"
    assert float(by[("T1", "TF1")]["p"]) == pytest.approx(1e-3)
    assert float(by[("T1", "TF1")]["fdr"]) == pytest.approx(2e-3)
    assert float(by[("T1", "TF2")]["fdr"]) == pytest.approx(0.5)
    assert float(by[("T2", "TF1")]["fdr"]) == pytest.approx(0.9)


def test_traits_first_source_wins_and_labels_join(tmp_path):
    scdori = tmp_path / "scdori"
    make_scdori(scdori)
    a = make_sldsc(tmp_path, "srcA")
    b = make_sldsc(tmp_path, "srcB")
    (b / "full/passing.txt").write_text("T1\nT2\nT4\n")
    _w(b / "h2_cts/M1/T4.cell_type_results.txt",
       [["Name", "Coefficient", "Coefficient_std_error", "Coefficient_P_value"], ["TF_TF9", 1, 1, "0.2"]])
    labels = tmp_path / "labels.tsv"
    _w(labels, [["trait", "label"], ["T1", "first trait"]])
    build(scdori, [a, b], "M1", tmp_path / "out", labels)
    traits = {r["trait"]: r for r in tsv(tmp_path / "out/traits.tsv")}
    assert set(traits) == {"T1", "T2", "T4"}
    assert traits["T1"]["source"] == "srcA" and traits["T4"]["source"] == "srcB"
    assert traits["T1"]["label"] == "first trait" and traits["T2"]["label"] == ""
    assert traits["T1"]["n_tf"] == "2"
    assert {r["source"] for r in tsv(tmp_path / "out/sldsc.tsv") if r["trait"] == "T1"} == {"srcA"}


def test_missing_dir_raises_configerror(tmp_path):
    scdori = tmp_path / "scdori"
    make_scdori(scdori)
    a = make_sldsc(tmp_path, "srcA")
    with pytest.raises(ConfigError, match="nowhere"):
        build(tmp_path / "nowhere", [a], "M1", tmp_path / "out", None)
    with pytest.raises(ConfigError, match="srcZ"):
        build(scdori, [a, tmp_path / "srcZ"], "M1", tmp_path / "out", None)
    with pytest.raises(ConfigError, match="eGRN_act_peak_gene"):
        (scdori / "eGRN_act_peak_gene" / "Ct_a_eGRN_edge_list.tsv").unlink()
        (scdori / "eGRN_act_peak_gene" / "Ct_b_eGRN_edge_list.tsv").unlink()
        (scdori / "eGRN_act_peak_gene").rmdir()
        build(scdori, [a], "M1", tmp_path / "out", None)


def test_malformed_peak_name_raises_configerror(tmp_path):
    scdori = tmp_path / "scdori"
    make_scdori(scdori)
    a = make_sldsc(tmp_path, "srcA")
    _w(scdori / "gene_peak_interactions/GENEA_enhancers_peak_links_.tsv",
       [["peak_name2", "Ct a"], ["chr1:100", 0.5]])
    with pytest.raises(ConfigError, match="chr1:100"):
        build(scdori, [a], "M1", tmp_path / "out", None)


def test_unknown_model_raises_configerror(tmp_path):
    scdori = tmp_path / "scdori"
    make_scdori(scdori)
    with pytest.raises(ConfigError, match="M2"):
        build(scdori, [make_sldsc(tmp_path, "srcA")], "M2", tmp_path / "out", None)
