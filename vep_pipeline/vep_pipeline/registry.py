"""Plugin registry: load plugins.yaml, resolve the plugins enabled in config.yaml, write fields.tsv."""
from __future__ import annotations
import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from .config import ConfigError, _normalise

CATEGORIES = ("no_data", "public", "existing", "manual", "user_file",
              "db_online", "grch37_only", "broken", "filter")
RUNNABLE = frozenset({"no_data", "public", "existing", "manual", "user_file"})
CONTROL_SOURCES = ("data", "vcf", "clinvar", "header")

# Why a category cannot be enabled, used when an entry has no `reason` of its own.
CATEGORY_REASONS = {
    "db_online": "needs a live Ensembl/web database; this pipeline runs VEP --offline",
    "grch37_only": "data/plugin is GRCh37 only; this pipeline is GRCh38",
    "broken": "plugin does not work in release/116",
    "filter": "filter plugin: drops consequences instead of annotating",
}

_PLACEHOLDER = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)\}")
_RENDER_KEYS = ("chunk", "workdir", "report_dir")

# Core VEP / custom-annotation fields: (field, evidence class). Trailing * = field family.
CORE_FIELDS = (
    [(f, "consequence") for f in ("Consequence", "IMPACT", "SYMBOL", "Gene", "Feature", "BIOTYPE",
                                  "HGVSc", "HGVSp", "HGVSg", "CANONICAL", "MANE_SELECT", "PICK")]
    + [(f, "known_variant") for f in ("Existing_variation", "CLIN_SIG", "SOMATIC", "PHENO",
                                      "PUBMED", "ClinVar_*")]
    + [(f, "population") for f in ("gnomADe_*", "gnomADg_*", "AF", "MAX_AF", "MAX_AF_POPS")]
    + [("SIFT", "clinically_trained_predictor"), ("PolyPhen", "clinically_trained_predictor"),
       ("ENCODE_cCRE", "regulatory_overlap"), ("GENE_PHENO", "gene_phenotype")]
)


@dataclass(frozen=True)
class Download:
    url: str
    as_: str
    size: int | None = None


@dataclass(frozen=True)
class FileSpec:
    key: str
    dest: str
    downloads: tuple = ()
    index: str | None = None
    index_url: str | None = None
    checksums_url: str | None = None
    prep: tuple = ()
    manual: str | None = None
    validate: dict = field(default_factory=dict)
    glob: str | None = None

    def members(self, path) -> list:
        """Data files behind `path`: the file itself, or for a `glob` entry the matching files in the directory."""
        path = Path(path)
        return sorted(path.glob(self.glob)) if self.glob else [path]


@dataclass(frozen=True)
class PluginSpec:
    name: str
    module: str
    category: str
    default: bool = False
    args: str | None = None
    files: tuple = ()
    params: dict = field(default_factory=dict)
    output_fields: tuple = ()
    evidence_class: str = ""
    field_classes: dict = field(default_factory=dict)
    eligible_consequences: tuple = ()
    conflicts_with: tuple = ()
    control: dict = field(default_factory=dict)
    reason: str | None = None


@dataclass(frozen=True)
class ResolvedPlugin:
    spec: PluginSpec
    paths: dict
    params: dict
    arg_template: str
    report_dir: Path | None = None

    def render(self, chunk: str, workdir: Path) -> str:
        values = {"chunk": str(chunk), "workdir": str(workdir)}
        if self.report_dir is not None:
            values["report_dir"] = str(self.report_dir)
        return _fill(self.arg_template, values)


# ----------------------------------------------------------------------------- loading

_ENTRY_KEYS = {"module", "category", "default", "args", "files", "params", "output_fields",
               "evidence_class", "field_classes", "eligible_consequences", "conflicts_with",
               "control", "reason"}
_FILE_KEYS = {"dest", "downloads", "index", "index_url", "checksums_url", "prep", "manual", "validate", "glob"}


def _fill(template: str, values: dict) -> str:
    """Replace {key} for known keys only; every other brace sequence is left untouched."""
    return _PLACEHOLDER.sub(lambda m: str(values[m.group(1)]) if m.group(1) in values else m.group(0),
                            template)


def fill_prep(cmd: str, path: Path, **extra) -> str:
    """Fill {path}, {dir} (and extras such as empty_indel) in a prep command; leave other braces (awk, sed)."""
    path = Path(path)
    return _fill(cmd, {"path": str(path), "dir": str(path.parent), **{k: v for k, v in extra.items() if v is not None}})


def _check_control(name, control):
    if not control:
        return
    src = control.get("source")
    if src not in CONTROL_SOURCES:
        raise ConfigError(f"{name}: control.source must be one of {CONTROL_SOURCES}, got {src!r}")
    if src == "data":
        cols = control.get("cols")
        if not (isinstance(cols, list) and len(cols) == 4 and all(isinstance(c, int) and c >= 1 for c in cols)):
            raise ConfigError(f"{name}: control.cols must be four 1-based column indices (CHROM, POS, REF, ALT)")
    if src == "clinvar" and not (isinstance(control.get("mc"), str) and control["mc"]):
        raise ConfigError(f"{name}: control.mc must be a non-empty string for source 'clinvar'")


def _parse_file(name, category, key, raw) -> FileSpec:
    bad = set(raw) - _FILE_KEYS
    if bad:
        raise ConfigError(f"{name}.files.{key}: unknown key(s) {sorted(bad)}")
    if "dest" not in raw:
        raise ConfigError(f"{name}.files.{key}: missing dest")
    dls = raw.get("downloads") or []
    downloads = []
    for d in dls:
        as_ = d.get("as", raw["dest"] if len(dls) == 1 else None)
        if not as_:
            raise ConfigError(f"{name}.files.{key}: download {d.get('url')} needs `as` (several downloads)")
        downloads.append(Download(d["url"], as_, d.get("size")))
    if raw.get("index") not in (None, "tbi", "csi"):
        raise ConfigError(f"{name}.files.{key}: index must be tbi or csi")
    if category in ("manual", "existing") and not raw.get("manual"):
        raise ConfigError(f"{name}.files.{key}: {category} entries need `manual` instructions")
    if raw.get("glob") and (dls or raw.get("prep") or raw.get("index_url") or not raw.get("manual")):
        raise ConfigError(f"{name}.files.{key}: a `glob` entry is a user-supplied directory: it needs `manual` "
                          f"instructions and cannot have downloads, prep or index_url")
    return FileSpec(key=key, dest=raw["dest"], downloads=tuple(downloads), index=raw.get("index"),
                    index_url=raw.get("index_url"), checksums_url=raw.get("checksums_url"),
                    prep=tuple(raw.get("prep") or ()), manual=raw.get("manual"),
                    validate=dict(raw.get("validate") or {}), glob=raw.get("glob"))


def _parse_entry(name, raw) -> PluginSpec:
    bad = set(raw) - _ENTRY_KEYS
    if bad:
        raise ConfigError(f"{name}: unknown registry key(s) {sorted(bad)}")
    for req in ("module", "category"):
        if req not in raw:
            raise ConfigError(f"{name}: missing {req}")
    if raw["category"] not in CATEGORIES:
        raise ConfigError(f"{name}: category must be one of {CATEGORIES}, got {raw['category']!r}")
    control = dict(raw.get("control") or {})
    _check_control(name, control)
    return PluginSpec(
        name=name, module=raw["module"], category=raw["category"],
        default=bool(raw.get("default", False)), args=raw.get("args"),
        files=tuple(_parse_file(name, raw["category"], k, v) for k, v in (raw.get("files") or {}).items()),
        params=dict(raw.get("params") or {}), output_fields=tuple(raw.get("output_fields") or ()),
        evidence_class=raw.get("evidence_class", ""), field_classes=dict(raw.get("field_classes") or {}),
        eligible_consequences=tuple(raw.get("eligible_consequences") or ()),
        conflicts_with=tuple(raw.get("conflicts_with") or ()), control=control,
        reason=raw.get("reason"))


def load_registry(path) -> dict:
    """plugins.yaml -> {name: PluginSpec}. Top-level keys starting with `_` are YAML anchors only."""
    raw = yaml.safe_load(Path(path).read_text()) or {}
    return {name: _parse_entry(name, entry) for name, entry in raw.items() if not str(name).startswith("_")}


def unsupported_reason(spec: PluginSpec) -> str:
    return spec.reason or CATEGORY_REASONS.get(spec.category, "")


# ----------------------------------------------------------------------------- resolving

def _blank(v):
    return v is None or v == ""


def _fmt_param(name, key, v):
    if isinstance(v, (list, tuple)):
        if not v:
            raise ConfigError(f"{name}: parameter '{key}' must not be empty")
        return ",".join(str(x) for x in v)
    return str(v)


def _resolve_one(name, spec, user, data_dir: Path, report_dir: Path) -> ResolvedPlugin:
    args = user.get("args")
    if spec.category not in RUNNABLE:
        if not (args and user.get("allow_unsupported")):
            raise ConfigError(
                f"{name}: category '{spec.category}' cannot be enabled ({unsupported_reason(spec)}). "
                f"To force it anyway set both `args: \"<plugin args>\"` and `allow_unsupported: true`.")
    if args:  # escape hatch: the user supplies the whole --plugin argument string
        return ResolvedPlugin(spec, {}, {}, args, report_dir)

    allowed = {"enabled", "args", "allow_unsupported"} | {f.key for f in spec.files} | set(spec.params)
    bad = set(user) - allowed
    if bad:
        raise ConfigError(f"{name}: unknown config key(s) {sorted(bad)}; valid: {sorted(allowed - {'enabled'})}")
    if spec.args is None:
        raise ConfigError(f"{name}: the registry has no recipe for this plugin"
                          f"{' (' + spec.reason + ')' if spec.reason else ''}; "
                          f"enable it with `args: \"{spec.module},...\"` (see {spec.module}.pm)")

    paths = {}
    for fs in spec.files:
        val = user.get(fs.key)
        if not _blank(val):
            paths[fs.key] = Path(str(val))
        elif fs.manual:
            raise ConfigError(f"{name}: data file '{fs.key}' is not configured "
                              f"(set plugins.{name}.{fs.key}). {' '.join(fs.manual.split())}")
        else:
            paths[fs.key] = data_dir / fs.dest
    params = {k: (user[k] if not _blank(user.get(k)) else v) for k, v in spec.params.items()}

    values = {**{k: str(p) for k, p in paths.items()},
              **{k: _fmt_param(name, k, v) for k, v in params.items()}}
    template = _fill(spec.args, values)
    left = [k for k in _PLACEHOLDER.findall(template) if k not in _RENDER_KEYS]
    if left:
        raise ConfigError(f"{name}: registry args reference undefined placeholder(s) {left}")
    return ResolvedPlugin(spec, paths, params, template, report_dir)


def report_dir(cfg: dict) -> Path:
    """`{output.dir}/reports` of the cfg being rendered: the one place {report_dir} (plugin args) and the chunk
    script's `mkdir -p` come from."""
    return Path(cfg["output"]["dir"]) / "reports"


def enabled_plugins(cfg: dict, registry: dict) -> list:
    """Plugins enabled in cfg['plugins'], resolved to paths/arguments, in config order."""
    section = _normalise(cfg.get("plugins"))
    wanted = [(n, e) for n, e in section.items() if e["enabled"]]
    for name, _ in wanted:
        if name not in registry:
            raise ConfigError(f"unknown plugin: {name} (see `list-plugins`)")
    names = [n for n, _ in wanted]
    for name in names:
        for other in registry[name].conflicts_with:
            if other in names:
                raise ConfigError(f"plugin conflict: {name} and {other} cannot be enabled together")
    data_dir = Path(cfg["data"]["plugin_data_dir"])
    return [_resolve_one(n, registry[n], e, data_dir, report_dir(cfg)) for n, e in wanted]


# ----------------------------------------------------------------------------- fields.tsv

def _plugin_fields(rp: ResolvedPlugin):
    spec = rp.spec
    names = rp.params["columns"] if isinstance(rp.params.get("columns"), list) else spec.output_fields
    for f in names:
        yield f, spec.field_classes.get(f, spec.evidence_class) or "unclassified"


def write_fields_tsv(plugins, out: Path) -> None:
    """field<TAB>evidence_class<TAB>source for the core VEP fields and every enabled plugin."""
    rows, seen = [], set()
    def add(f, cls, src):
        if f not in seen:
            seen.add(f)
            rows.append(f"{f}\t{cls}\t{src}")
    for f, cls in CORE_FIELDS:
        add(f, cls, "vep")
    for rp in plugins:
        for f, cls in _plugin_fields(rp):
            add(f, cls, rp.spec.name)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("field\tevidence_class\tsource\n" + "\n".join(rows) + "\n")
