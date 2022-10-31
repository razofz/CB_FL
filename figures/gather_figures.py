import re
import shutil
import pandas as pd
import fl_utils

df = pd.read_excel("figure_renaming.ods", header=1).dropna(how="all")
df.to_csv("figure_renaming_auto_converted_from_ods.csv")
del df


with open("figure_renaming_auto_converted_from_ods.csv", "r") as f:
    lines = [line[:-1] for line in f.readlines() if not line.startswith(",,")]
    # lines = lines[1:]
    lines = [re.sub(",{2,}", "", line) for line in lines]


fl_utils.set_project_path()

figures = {}

for line in lines:
    split = line.split(",")
    d = {"origin": split[2], "filenames": split[3:]}
    figures[split[1]] = d

NOTEBOOKS_PATH = "data/processed/notebooks/Figures"
DESEQ2_PATH = "data/processed/DEseq2/images"
GLOBAL_FIGURES_PATH = "figures"

for figure in list(figures.keys()):
    print(figure)
    fl_utils.print_files(figures[figure])
    if figures[figure]["origin"] not in ["py", "DEseq"]:
        print(f"Skipping: {figure} {figures[figure]=}")
    else:
        if figures[figure]["origin"] == "py":
            path = NOTEBOOKS_PATH
            ft = ".svg"
        elif figures[figure]["origin"] == "DEseq":
            path = DESEQ2_PATH
            ft = ".pdf"

        if isinstance(figures[figure]["filenames"], list):
            i = 1
            for fn in figures[figure]["filenames"]:
                try:
                    shutil.copy(
                        path + "/" + fn + ft,
                        GLOBAL_FIGURES_PATH
                        + "/"
                        + figure
                        + (
                            str(i)
                            if len(figures[figure]["filenames"]) > 1
                            else ""
                        )
                        + "_"
                        + fn
                        + ft,
                    )
                except:
                    print(f"file {fn} not existing")
                i += 1
        else:
            print(f'{type(figures[figure]["filenames"])=}')
