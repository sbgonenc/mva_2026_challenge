# Local VEP plugin: SpliceAI2

`SpliceAI2.pm` annotates variants with Illumina's precomputed SpliceAI2 scores (v2.0, GENCODE v48 protein-coding genes; non-commercial terms).

It is not part of Ensembl's VEP_plugins `release/116`. To use it:

1. Copy it into the VEP_plugins checkout that `tools.plugins_dir` in `vep_pipeline/config.yaml` points to.
2. Apply `plugin_config.SpliceAI2.diff` to that checkout's `plugin_config.txt`: `git apply plugin_config.SpliceAI2.diff`.
3. Set `plugins.SpliceAI2.snv` and `.indel` in `vep_pipeline/config.yaml` to the directories of the per-chromosome score files. The files must be tabix-indexed with `-S 1`.
