# SkipGuide Analysis
The code used for the analysis and production of the results described in the paper **TBA**.

## Dependencies
The analysis was performed using Python 3.7.5 and Jupyter. The dependencies are listed in [environment.yml](environment.yml). We recommend using the conda package manager from [Anaconda Python](https://www.anaconda.com/distribution/) to create an environment for running the analysis:

`conda env create -f environment.yml`

Activate the environment by:

`conda activate skipguide_data_processing`

## Data Files
The provided Jupyter notebooks (see [Usage](#usage) section) can produce all the results starting from the [raw sequencing data](#raw-data-files). However, computations can take a very long time, on the order of hours or days. The notebooks are configured to skip certain long computations if pre-computed files are available. We recommend you instead download the [pre-computed files](#pre-computed-files) before running the notebooks.

### Raw Data Files
If you opt to not use the pre-computed files, the raw sequencing data needs to be available. Download them from **TBA**, and place them in the [`data/reads`](data/reads) directory before running the provided notebooks. Running all the notebooks may take on the order of hours or days.

### Pre-Computed Files
If you opt to use the pre-computed files, the raw sequencing data is not necessary. Download the pre-computed files from **TBA**, extract, and replace the [`cache`](cache) directory with the extracted `cache` directory. Running all the notebooks should then take less than half an hour.

## Usage
Note that the output images and tables ([`output` directory](output)) are already uploaded to this repository for your convenience. You can also open the provided notebooks and view the results. This section details how you can run the notebooks from scratch.

See [Data Files](#data-files) section to include the necessary data files.

If pre-computed files are not used, modify the `NUM_PROCESSES` variable in [config.py](src/config.py) to specify the number of cores for multiprocessing.

Start a Jupyter notebook server, e.g.:

`jupyter notebook --port=8888`

Run the provided notebooks under [`src`](src) in the following order:
1. [Sequence_Extraction.ipynb](src/Sequence_Extraction.ipynb)
2. [Barcode_Sequence_Lookup_Tables.ipynb](src/Barcode_Sequence_Lookup_Tables.ipynb)
3. [datA_Characterize_Sequences_Indels.ipynb](src/datA_Characterize_Sequences_Indels.ipynb)
4. [inDelphi_Evaluation.ipynb](src/inDelphi_Evaluation.ipynb)
5. [datB_Characterize_Skipping.ipynb](src/datB_Characterize_Skipping.ipynb)
6. [MMSplice_Predict_Skipping.ipynb](src/MMSplice_Predict_Skipping.ipynb)
7. [wMMSplice_SkipGuide_Evaluation.ipynb](src/wMMSplice_SkipGuide_Evaluation.ipynb)

Inspect the comments and markdown in the notebooks for more context.
