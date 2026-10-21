from __future__ import annotations
import copy, re
from pathlib import Path
import yaml

_PH = re.compile(r"\{([A-Za-z_][A-Za-z0-9_.]*)\}")
REQUIRED = ["sample", "input.vcf", "input.reference_fasta", "data.root", "data.cache_dir",
            "data.plugin_data_dir", "data.custom_data_dir", "output.dir", "output.prefix",
            "tools.env", "tools.vep_repo", "tools.plugins_dir", "tools.prepare_vcf", "tools.pget"]

class ConfigError(ValueError):
    pass

def _get(cfg, dotted):
    node = cfg
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            raise ConfigError(f"unknown placeholder {{{dotted}}}")
        node = node[part]
    return node

def resolve_placeholders(cfg: dict, max_passes: int = 10) -> dict:
    cur = copy.deepcopy(cfg)
    def walk(node, root):
        if isinstance(node, dict): return {k: walk(v, root) for k, v in node.items()}
        if isinstance(node, list): return [walk(v, root) for v in node]
        if isinstance(node, str): return _PH.sub(lambda m: str(_get(root, m.group(1))), node)
        return node
    for _ in range(max_passes):
        nxt = walk(cur, cur)
        if nxt == cur:
            if any(_PH.search(s) for s in _strings(nxt)):
                raise ConfigError("placeholder cycle detected")
            return nxt
        cur = nxt
    raise ConfigError("placeholder cycle detected")

def _strings(node):
    if isinstance(node, dict):
        for v in node.values(): yield from _strings(v)
    elif isinstance(node, list):
        for v in node: yield from _strings(v)
    elif isinstance(node, str):
        yield node

def _normalise(section):
    out = {}
    for name, val in (section or {}).items():
        if isinstance(val, bool): out[name] = {"enabled": val}
        elif isinstance(val, dict): out[name] = {**val, "enabled": bool(val.get("enabled", True))}
        else: raise ConfigError(f"{name}: expected true/false or a mapping")
    return out

def load_config(path) -> dict:
    cfg = resolve_placeholders(yaml.safe_load(Path(path).read_text()) or {})
    for key in REQUIRED:
        try: val = _get(cfg, key)
        except ConfigError: val = None
        if val in (None, ""): raise ConfigError(f"missing required config key: {key}")
    if cfg.get("assembly") != "GRCh38":
        raise ConfigError("assembly must be GRCh38 (only GRCh38 data is configured)")
    for s in _strings({k: cfg[k] for k in ("input", "output", "data", "tools")}):
        if re.search(r"\s", s): raise ConfigError(f"whitespace not allowed in paths: {s!r}")
    cfg["plugins"] = _normalise(cfg.get("plugins"))
    cfg["custom"] = _normalise(cfg.get("custom"))
    cfg.setdefault("sv", {"enabled": False})
    return cfg
