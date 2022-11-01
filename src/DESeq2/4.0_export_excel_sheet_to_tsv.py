import pandas as pd
import fl_utils
import os

fl_utils.set_project_path()
df = pd.read_excel(snakemake.input["ball_xlsx"])
df = df.iloc[:-5]
df.to_csv(
    snakemake.output["ball_tsv"],
    sep="\t",
    index=False,
    decimal=",",
    na_rep="NA"
)
