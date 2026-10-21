from pathlib import Path

import pytest

from vep_pipeline.config import ConfigError
from vep_pipeline.genedb import GeneInfo, load_g2p, load_sf, moi_codes
from vep_pipeline.rank_settings import PRIORITISE_DEFAULTS

G2P = ('"g2p id","gene symbol","allelic requirement","confidence","molecular mechanism"\n'
       'G2P1,BUB1B,biallelic_autosomal,definitive,loss of function\n'
       'G2P2,GENEAD,monoallelic_autosomal,strong,gain of function\n'
       'G2P3,GENEAD,biallelic_autosomal,limited,undetermined\n')


def test_moi_codes():
    assert moi_codes("biallelic_autosomal") == {"AR"}
    assert moi_codes("monoallelic_autosomal") == {"AD"}
    assert moi_codes("monoallelic_X_hemizygous") == {"XR"}
    assert moi_codes("monoallelic_X_heterozygous") == {"XD"}
    assert moi_codes("mitochondrial") == {"MT"}
    assert moi_codes("biallelic_autosomal; monoallelic_autosomal") == {"AR", "AD"}
    assert moi_codes("monoallelic_X") == {"XD", "XR"}
    assert moi_codes("monoallelic_Y_hemizygous") == {"AD"}     # one copy suffices; treated like AD


def test_load_g2p_filters_confidence(tmp_path):
    p = tmp_path / "g2p.csv"
    p.write_text(G2P)
    genes = load_g2p(p, ["definitive", "strong", "moderate"])
    assert genes["BUB1B"] == GeneInfo(frozenset({"AR"}), frozenset({"loss of function"}))
    assert genes["GENEAD"].modes == frozenset({"AD"})          # the 'limited' AR row is dropped


def test_load_g2p_missing_column(tmp_path):
    p = tmp_path / "g2p.csv"
    p.write_text('"gene symbol","confidence"\nX,definitive\n')
    with pytest.raises(ConfigError, match="allelic requirement"):
        load_g2p(p, ["definitive"])


def test_real_g2p_header_is_accepted():
    head = Path(__file__).parent / "data" / "prioritise" / "g2p_head.csv"
    assert isinstance(load_g2p(head, ["definitive"]), dict)    # Task 1 fixture: the real header parses


def test_shipped_sf_table():
    sf = load_sf(PRIORITISE_DEFAULTS["secondary"]["genes"])
    assert len(sf) == 84
    assert sf["TP53"].rule == "all_plp" and sf["TTN"].rule == "truncating_only"
    assert sf["HFE"].rule == "hfe_c282y_hom" and sf["MUTYH"].rule == "biallelic"
    assert sf["ABCD1"].rule == "xl_biallelic" and sf["GLA"].rule == "xl_any"
    assert sf["LDLR"].inheritance == "SD"
