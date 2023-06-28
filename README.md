# Single-cell multiomics of human fetal hematopoiesis defines a developmental specific population and a fetal signature

**N. B.** Before running any scripts or notebooks, please edit the file
`project_path` to contain your local path to this folder. For example, if your
username is `bioinf` and this folder (directory) is called `Sommarin_et_al` and
resides in your home directory, the `project_path` file should contain:

    /home/bioinf/Sommarin_et_al

The scripts and notebooks read this file to set working directory, output
directories etc, so please do this first. You can also use
[direnv](https://direnv.net/), see the included `.envrc` file (which also reads
the `project_path` file). This will set the environment variable PROJECT_PATH
to be what you specified in the project_path file.

**N. B. II** The procedure to run this project is an amalgamation of
approaches, due to such things as different code authors with different
backgrounds, several revision rounds etc. To abstract, the pipeline for running
this project is as follows:

1. Run the Seurat part
2. Run the 1st notebooks part
3. Run the DESeq2 part
4. Run the 2nd notebooks part

One level down abstraction-wise:

1. Activate the conda environment specified in `envs/seurat.yaml` and run the
   `src/seurat-analysis/run_all.sh` script (working directory/where to execute this script from:
   `$PROJECT_PATH/src/seurat-analysis`).
2. Activate the conda environment specified in `envs/notebooks.yaml` and run
   the associated notebooks with this command:
   `jupyter nbconvert --inplace --execute --config config_for_running_notebooks_pre_DESeq2.py`
   (working directory `$PROJECT_PATH/notebooks`)
3. Install [snakemake](https://snakemake.readthedocs.io/en/stable/) (e.g. in a
   conda environment) and run the DESeq2 part with this command:
   `snakemake -j all --use-conda` (working directory `$PROJECT_PATH`).
4. Activate the conda environment specified in `envs/notebooks.yaml` and run
   the associated notebooks with this command:
   `jupyter nbconvert --inplace --execute --config config_for_running_notebooks_post_DESeq2.py`
   (working directory `$PROJECT_PATH/notebooks`)
5. Gather the produced figures used for the paper by running snakemake again,
   letting it find the scattered plots and copying them into the
   `figures/figures` folder. Use the same command as above:
   `snakemake -j all --use-conda` (working directory `$PROJECT_PATH`).

### Directory structure

Guide for folder (directory) names:

- `data`
    - contains the data generated in and needed for this study.
    - `external`: external data used
    - `processed`: the output of this project (data and figures)
    - `interim`: "middle-steps" that are used to generate processed data but is
      in itself not part of the final results
    - `raw`: the raw data generated in this study
- `envs`
    - contains conda environment specifications for the three different parts
      of this project. Also used by snakemake for the DESeq2 part of this
      study.
    - install [conda](https://docs.conda.io/en/latest/) (or even better,
      [mamba](https://github.com/mamba-org/mamba)), and then each environment
      can be generated with the command `conda env create -f
      <environment_file.yaml>` (substitute `conda` with `mamba` if you
      installed mamba).  (Don't forget to also activate each environment after
      creation before running any scripts or notebooks.)
        - (note I: if you already have an environment with the same name, you
          can e.g. change the name in the first line of corresponding .yaml
          file to avoid conflict errors.)
        - (note II: snakemake will create the environments it needs for itself
          from these files, but you will still need to create them for the
          seurat and notebooks parts of the study.)
- `figures`
    - directory where all the project's plots that are used as figures for the
      accompanying paper are gathered. This is done by snakemake. The naming
      convention for the figures files is: *\<figure index\>_\<original filename\>*.
      E.g. `Figure1b_UMAP_FL.svg`
- `notebooks`
    - jupyter notebooks (if you generate environments like specified above,
      jupyterlab is installed)
    - both .ipynb and .py representations are present, they are synced with
      [jupytext](https://jupytext.readthedocs.io/en/latest/). You can either
      run the .py files or the .ipynb files. For executing the .ipynb
      (notebook) files, including generating cell output, there are two files
      for that: `config_for_running_notebooks_pre_DESeq2.py` and
      `config_for_running_notebooks_post_DESeq2.py`. The files are not meant to
      be run themselves, see the first lines of comments in the files for a
      command to run (which uses the respective file as a config file). Also
      included is a command to trigger the syncing of .py and .ipynb files
      after changes made to any of the formats (syncing will nevertheless
      happen automatically if you use jupyter from the environment included in
      this project (`envs/notebooks.yaml`) to open an .ipynb file).
- `src`
    - scripts used in this study
    - `seurat-analysis`: scripts for the part using the R package `Seurat`
      (amongst others). There is a shell script for running all R scripts:
      `run_all.sh`. Activate the seurat environment and run this. You might
      have to allow it to run as well (e. g. `chmod +x run_all.sh`).
    - `DESeq2`: scripts for the part using the R package `DESeq2` (amongst
      others). Use snakemake to run this part.

This project is runnable in its current form on Unix-based systems (Linux or
macOS), i. e. not Windows. The paths are specified with a forward slash. If you
are using Windows and want to reproduce these results, one suggestion would be
to use WSL2 (Windows Subsystem for Linux).

### Docker container

Another approach is to use a [Docker](https://docs.docker.com/) container. We
have prepared one that should enable all parts of this analysis to run out of
the box. You can find it on dockerhub at
[https://hub.docker.com/r/razofz/fetal-liver](https://hub.docker.com/r/razofz/fetal-liver).
When you have docker installed you can pull (download) the container with this command:

    docker pull razofz/fetal-liver:0.4

Depending on how you configured Docker you might have to preface the command
with `sudo`. See e. g.
[here](https://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo)
for some pointers on how to configure Docker for use without sudo.  
When you have pulled the container you can start an instance of it with the
following command. If you use this exact command make sure you have changed the
`project_path` file to contain simply `/fl`. Note that this will reflect the
path on the container and not on your computer. You also need to have navigated
to the directory on your computer in which this project resides for this to
work (due to the `PWD` command in the command). In other words, if you have
downloaded this project as in the example above, you need to have navigated
to `/home/bioinf/Sommarin_et_al` to run the following command:

    docker run --interactive --tty --volume ${PWD}:/fl --rm --name Sommarin_et_al razofz/fetal-liver:0.4 /bin/bash

This will give you a bash session in the Docker container. The project files
are in the `/fl` dir, and you can follow the instructions above on activating
the conda environments (they have already been created in this container, see
the included `Dockerfile` for details).

/Rasmus

P.S.
Note that our adult bone marrow samples have the label yBM in many of the
scripts. This is due to that data being from a study where there was even older
bone marrow data, so in relation to that these were the younger bone marrow
data.
Also note that our CS22 data have label hpc at times.
D.S.

---

(Project structure is a slimmed-down version of the [Data Science
Cookiecutter](https://drivendata.github.io/cookiecutter-data-science/)) 
