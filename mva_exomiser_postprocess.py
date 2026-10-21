import marimo

__generated_with = "0.24.2"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import pandas as pd
    import os
    import numpy as np

    return (pd,)


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Clinvar Annotated Exomiser Run ##
    """)
    return


@app.cell
def _(pd):
    exomiser_in_path = "/fast/AG_Bunina/Berk/projects/mva_challenge_26/exomiser/clinvar_annotated_run/exome/SAMPLE_ID.variants.tsv"

    exome_vars = pd.read_csv(exomiser_in_path, sep="\t")
    exome_vars
    return (exome_vars,)


@app.cell
def _(exome_vars):
    fltrd_exome_vars = exome_vars[
        ((exome_vars["P-VALUE"] <= 0.05) | (exome_vars["CLINVAR_STAR_RATING"] >= 2))
        & (exome_vars["EXOMISER_VARIANT_SCORE"] > 0.1)
        & ((exome_vars["MAX_FREQ"] < 0.1) | exome_vars["MAX_FREQ"].isnull())
        & (~exome_vars["CLINVAR_PRIMARY_INTERPRETATION"].str.contains("^BENIGN"))
        & (~exome_vars["EXOMISER_ACMG_CLASSIFICATION"].str.contains("BENIGN"))
        & (exome_vars["#RANK"] < 100)
        & (exome_vars["FUNCTIONAL_CLASS"] != "synonymous_variant")
        ]#.drop_duplicates(subset=["CONTIG", "START"])
    fltrd_exome_vars
    return (fltrd_exome_vars,)


@app.cell
def _(pd):
    genomiser_in_path = "/fast/AG_Bunina/Berk/projects/mva_challenge_26/exomiser/clinvar_annotated_run/genome/SAMPLE_ID.variants.tsv"

    genome_vars = pd.read_csv(genomiser_in_path, sep="\t")
    genome_vars
    return (genome_vars,)


@app.cell
def _(exome_vars, genome_vars):
    fltrd_genome_vars = genome_vars[
        ((genome_vars["P-VALUE"] <= 0.05) | (genome_vars["CLINVAR_STAR_RATING"] >= 2))
        & (genome_vars["EXOMISER_VARIANT_SCORE"] > 0.1)
        & ((genome_vars["MAX_FREQ"] < 0.1) | genome_vars["MAX_FREQ"].isnull())
        & (~genome_vars["CLINVAR_PRIMARY_INTERPRETATION"].str.contains("^BENIGN"))
        & (genome_vars["#RANK"] < 100)
        & (genome_vars["FUNCTIONAL_CLASS"] != "synonymous_variant")
        & (~exome_vars["EXOMISER_ACMG_CLASSIFICATION"].str.contains("BENIGN"))
        ]#.drop_duplicates(subset=["CONTIG", "START"])
    fltrd_genome_vars
    return (fltrd_genome_vars,)


@app.cell
def _(fltrd_exome_vars, fltrd_genome_vars, pd):
    exome_genome_concat = pd.concat([fltrd_genome_vars, fltrd_exome_vars])#.drop_duplicates(subset=["CONTIG", "START"])
    output_df = exome_genome_concat[[
        "#RANK",
        "CONTIG",
        "START", 
        "REF", 
        "ALT",
        "GENE_SYMBOL",
        "EXOMISER_ACMG_DISEASE_NAME",
        "EXOMISER_VARIANT_SCORE", 
        "MAX_PATH", 
        "MAX_FREQ", 
        "EXOMISER_ACMG_EVIDENCE",
        "GENOTYPE",
        "FUNCTIONAL_CLASS"
    ]].copy()
    output_df
    return (output_df,)


@app.cell
def _(output_df):
    "proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes"

    output_df["proband_id"] = "SAMPLE_ID"

    output_df.rename({"CONTIG":"chrom_1", "START": "pos_1", "REF":"ref_1", "ALT":"alt_1", "EXOMISER_VARIANT_SCORE":"epcr"}, axis=1, inplace=True)
    return


@app.cell
def _(output_df):
    output_df["finding_type"] = output_df["GENE_SYMBOL"].apply(lambda x: "primary" if x=="BUB1B" else "secondary")
    return


@app.cell
def _(output_df):
    output_df
    return


@app.cell
def _(output_df):
    # 2. Define a lambda function to extract the first mode, handling all-NaN groups
    fill_mode = lambda x: x.fillna(x.mode()[0]) if not x.mode().empty else x

    # 3. Apply the transformation to columns outside the grouping keys
    cols = output_df.columns.difference(['chrom_1', 'pos_1'])
    output_df[cols] = output_df.groupby(['chrom_1', 'pos_1'])[cols].transform(fill_mode)
    return


@app.cell
def _(output_df):
    output_df
    return


@app.cell
def _(output_df):
    output_df["notes"] = output_df["EXOMISER_ACMG_DISEASE_NAME"] +";"+ output_df["EXOMISER_ACMG_EVIDENCE"] + ";" + output_df["GENOTYPE"]
    output_df
    return


@app.cell
def _(output_df):
    final_output_df = output_df.loc[:,~output_df.columns.str.isupper()].drop_duplicates(subset=["chrom_1", "pos_1", "ref_1", "alt_1"]).nlargest(10, "epcr")
    final_output_df
    return (final_output_df,)


@app.cell
def _(final_output_df):
    final_output_df[["chrom_2","pos_2","ref_2","alt_2"]] = None
    return


@app.cell
def _(final_output_df):
    final_output_df[["proband_id","chrom_1","pos_1","ref_1","alt_1","chrom_2","pos_2","ref_2","alt_2","epcr","finding_type","notes"]].to_csv("sbgonenc96-exomiser+genomiser_run.csv", index=False)
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
