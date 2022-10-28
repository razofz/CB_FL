import pandas as pd
import fl_utils
import os

fl_utils.set_project_path()
os.chdir("data/external/BALL-1988S-HTSeq")
df = pd.read_excel("B-ALL-subtyping.xlsx")
df = df.iloc[:-5]
df.to_csv("subtypes.tsv", sep="\t", index=False, decimal=",", na_rep="NA")
