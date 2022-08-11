Source code and input files related to the figures in the rMATS-turbo manuscript "rMATS-turbo: an efficient and flexible computational tool for alternative splicing analysis of large-scale RNA-seq data"

## Figure 3
  - The static HTML format of figure generating process is demonstrated in `plot_Fig_PC3E_GS689.html`
  - Source code for generating Figure 3 is demonstrated in `plot_Fig_PC3E_GS689.ipynb`
  - Input files for generating Figure 3 is located in folder `PC3E-GS689`
  - The output files are created in `./PC3E-GS689/plot_rmats2sashimi/output/Sashimi_plot/` and `./PC3E-GS689/plot_venn_volcano`
## Figure 4
  - The static HTML format of figure generating process is demonstrated in `plot_Fig_CCLE.html`
  - Source code for generating Figure 4 is demonstrated in `plot_Fig_CCLE.ipynb`
  - Input files for generating Figure 4 is located in folder `CCLE`
  - The output file is `./CCLE/plot_heatmap/plot_heatmap_sortedByEMT_usingCorrelation.pdf`

## Requirements
  - [Jupyter Notebook](https://jupyter.org/)
  - [rmats2sashimiplot](https://github.com/Xinglab/rmats2sashimiplot)
  - Python related dependencies:
    - rpy2
  - R packages
    - zeallot
    - [ComplexHeatmap](https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
    - circlize
    - ggplot2
    - ggrepel
  - GNU sed

## Install
* MAC users: install GNU sed
  + `brew install gnu-sed`
  + `PATH="/usr/local/opt/gnu-sed/libexec/gnubin:$PATH"`
* install most dependencies to a conda environment
  + `conda create --prefix ./conda_env`
  + `conda activate ./conda_env`
  + `conda install -c conda-forge -c bioconda --file ./conda_requirements.txt`
  + `CFLAGS="-I$(pwd)/conda_env/include" pip install rpy2==2.8.6`
  + `conda deactivate`
* install rmats2sashimiplot and samtools
  + `conda install -c conda-forge -c bioconda --file ./conda_requirements_rmats2sashimiplot.txt`

## Run
* Activate the main conda_env
  - `conda activate ./conda_env`
* The code in `plot_Fig_PC3E_GS689.ipynb` and `plot_Fig_CCLE.ipynb` can be run interactively in web browser via Jupyter Notebook.
  - Launching Jupyter Notebook
    + `jupyter notebook`
* It can also be run as an entire script
  - Convert notebooks to scripts
    + `jupyter nbconvert --to python ./plot_Fig_PC3E_GS689.ipynb`
    + `jupyter nbconvert --to python ./plot_Fig_CCLE.ipynb`
  - For the PC3E notebook:
    + `ipython ./plot_Fig_PC3E_GS689.py`
  - For the CCLE notebook:
    + `ipython ./plot_Fig_CCLE.py`
