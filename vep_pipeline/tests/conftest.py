import pytest, yaml
from vep_pipeline.config import resolve_placeholders

@pytest.fixture
def minimal_cfg(tmp_path):
    """Config dict shaped like the shipped config.yaml (placeholders resolved); every path under tmp_path."""
    t = str(tmp_path)
    return resolve_placeholders({
        "assembly": "GRCh38",
        "cache_version": 116,
        "sample": "S1",
        "tools": {
            "env": "vep116",
            "vep_repo": f"{t}/ensembl-vep",
            "plugins_dir": f"{t}/VEP_plugins",
            "prepare_vcf": f"{t}/prepare-patient-vcf.sh",
            "prepare_env": "exomiser",
            "pget": f"{t}/pget.sh",
        },
        "input": {
            "vcf": f"{t}/in.vcf.gz",
            "normalize": False,
            "reference_fasta": f"{t}/ref.fa",
        },
        "sv": {"enabled": False},
        "output": {"dir": f"{t}/out", "prefix": "{sample}"},
        "data": {
            "root": f"{t}/data",
            "cache_dir": "{data.root}/cache",
            "plugin_data_dir": "{data.root}/plugins",
            "custom_data_dir": "{data.root}/custom",
        },
        "run": {"fork": 1, "buffer_size": 100, "extra_flags": []},
        "slurm": {"partition": "normal", "cpus": 1, "mem": "1G", "time": "01:00:00",
                  "log_dir": "{output.dir}/logs"},
        "smoke": {"region": "17:43044295-43125483"},
        "custom": {},
        "plugins": {},
    })

MINI_REGISTRY = {
    "Toy": {
        "module": "Toy",
        "category": "public",
        "default": False,
        "files": {"file": {"dest": "Toy/toy.tsv.gz",
                           "downloads": [{"url": "file:///dev/null", "as": "Toy/toy.tsv.gz"}],
                           "index": "tbi"}},
        "args": "Toy,file={file}",
        "output_fields": ["TOY"],
        "evidence_class": "unsupervised_predictor",
    },
    "Nodata": {
        "module": "Nodata",
        "category": "no_data",
        "args": "Nodata",
        "output_fields": ["NODATA"],
        "evidence_class": "consequence",
    },
}

@pytest.fixture
def mini_registry(tmp_path):
    """Two-entry registry (Toy, Nodata) as a dict of PluginSpec."""
    from vep_pipeline.registry import load_registry
    p = tmp_path / "plugins.yaml"
    p.write_text(yaml.safe_dump(MINI_REGISTRY, sort_keys=False))
    return load_registry(p)


# ----------------------------------------------------------------------------- preflight-check fixtures

import subprocess
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[2]
TESTDATA = SCRIPTS / "ensembl-vep" / "t" / "testdata"
TESTCACHE = TESTDATA / "cache"
TESTCACHE_FASTA = TESTCACHE / "homo_sapiens" / "84_GRCh38" / "test.fa"
PLUGINS_DIR = SCRIPTS / "VEP_plugins"


def _sh(*cmd):
    subprocess.run(cmd, check=True, capture_output=True, text=True)


@pytest.fixture
def vcf_factory(tmp_path):
    """vcf_factory(contigs={name: length}, records=[(chrom, pos, ref, alt[, gt])], samples=("S1",)) -> bgzipped,
    tabix-indexed VCF path. Records are written in contig-header order; gt defaults to 0/1 for every sample."""
    n = [0]
    def make(contigs, records=(), samples=("S1",)):
        n[0] += 1
        order = {c: i for i, c in enumerate(contigs)}
        lines = ["##fileformat=VCFv4.2"]
        lines += [f"##contig=<ID={c},length={length}>" for c, length in contigs.items()]
        lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]))
        for rec in sorted(records, key=lambda r: (order[r[0]], r[1])):
            chrom, pos, ref, alt = rec[:4]
            gt = rec[4] if len(rec) > 4 else "0/1"
            lines.append("\t".join([chrom, str(pos), ".", ref, alt, ".", ".", ".", "GT", *[gt] * len(samples)]))
        path = tmp_path / f"v{n[0]}.vcf"
        path.write_text("\n".join(lines) + "\n")
        _sh("bgzip", "-f", str(path))
        _sh("tabix", "-f", "-p", "vcf", f"{path}.gz")
        return Path(f"{path}.gz")
    return make


@pytest.fixture
def fasta_factory(tmp_path):
    """fasta_factory({name: length}, seq={name: sequence}) -> FASTA path with a .fai.

    A contig without `seq` is N (sparse and with a hand-written .fai when it is large, so a 249 Mb contig is cheap)."""
    n = [0]
    def make(lengths, seq=None):
        n[0] += 1
        seq = seq or {}
        path = tmp_path / f"f{n[0]}.fa"
        if all(c in seq for c in lengths) or sum(lengths.values()) <= 1_000_000:
            with open(path, "w") as fh:
                for c, length in lengths.items():
                    s = seq.get(c, "N" * length)
                    assert len(s) == length, (c, len(s), length)
                    fh.write(f">{c}\n" + "".join(s[i:i + 60] + "\n" for i in range(0, length, 60)))
            _sh("samtools", "faidx", str(path))
            return path
        fai, off = [], 0
        for c, length in lengths.items():
            off += len(c) + 2
            fai.append(f"{c}\t{length}\t{off}\t60\t61")
            off += length + (length + 59) // 60
        with open(path, "wb") as fh:
            fh.truncate(off)
        Path(f"{path}.fai").write_text("\n".join(fai) + "\n")
        return path
    return make


@pytest.fixture
def testcache_cfg(minimal_cfg, tmp_path):
    """minimal_cfg wired to the VEP test cache that ships with the repo (version 84, GRCh38, chr21 only)."""
    cfg = minimal_cfg
    cfg["cache_version"] = 84
    cfg["data"]["cache_dir"] = str(TESTCACHE)
    cfg["data"]["plugin_data_dir"] = str(tmp_path / "data" / "plugins")
    cfg["tools"]["plugins_dir"] = str(PLUGINS_DIR)
    cfg["input"]["reference_fasta"] = str(TESTCACHE_FASTA)
    cfg["output"]["dir"] = str(tmp_path / "out")
    return cfg
