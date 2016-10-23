# Single-cell RNA-Seq TCC prep

This repository contains scripts needed to generate transcript compatibility count (TCC) matrices from single-cell RNA-Seq data. Included is error-correction of barcodes, collapsing of UMIs and pseudoalignment of reads to a transcriptome to obtain transcript compatibility counts. The scripts utilize [kallisto](http://pachterlab.github.io/kallisto) for pseudoalignment.

We currently support the 10X Chromium technology; support for more technologies is underway.

## Instructions for processing 10X Chromium 3' digital expression data

#### Getting started

The [getting started](http://pachterlab.github.io/kallisto/10xstarting.html) tutorial explains how to process the small example in the [example_dataset](https://github.com/lakigigar/scRNA-Seq-TCC-prep/tree/master/example_dataset) directory. This is a good starting point to make sure that the necessary programs are correctly installed. Note that you will need kallisto (≥ 0.43.0), python (≥ 2.7.10), scipy (≥ 0.16.0), scikit-learn (__= 0.16.1__) and Juypter Notebook (≥ 4.0.6) installed (the Jupyter requirement is not strictly necessary but highly recommended). __Please note that there appears to be a problem with tSNE in scikit-learn v0.17.1__.

#### Workflow organization

The processing workflow consists of four steps: 

0. Preparation of a configuration file that contains the parameters needed for the processing.
1. Identification of "true" cell barcodes according to read coverage followed by error correction when possible.
2. Creation of read/UMI files for each cell.
3. [Pseudoalignment](http://www.nature.com/nbt/journal/v34/n5/abs/nbt.3519.html) of reads associated with each cell using __kallisto__, deduplication according to UMIs, and generation of transcript compatibility counts (TCCs) for each cell.

Following the pre-processing, the transcript compatibility counts (TCC) matrix can be analyzed using a Jupyter Notebook. 

#### Creation of the configuration file

Parameters needed to run the processing require specification of a `config.json` file. The following parameters need to be specified:

- NUM_THREADS: the number of threads available for processing.
- WINDOW: this parameter contains a lower and upper threshold for the expected number of cells in the experiment. It is used in the determination of the number o cells in the experiment from reads coverage data.
- SOURCE_DIR: path to the source directory that contains the .py scripts
- BASE_DIR: this must contain the path to the (demultiplexed) FASTQ files from the sequencing. Note that our workflow does not currently demultiplex reads and you may have to do so with 10X's software; we plan to provide a demultiplexing script in the future.
- sample_idx: The sample index used for the run e.g., for SI-3A-A10 -> sample_idx: ["ACAGCAAC", "CGCAATTT", "GAGTTGCG", "TTTCGCGA"].
- SAVE_DIR: path to a directory where intermediate files will be saved.
- dmin: the minimum distance between barcodes needed for error correction to be performed.
- BARCODE_LENGTH: length of the barcodes.
- OUTPUT_DIR: directory in which to output results.
- kallisto: path to the binary for kallisto, location of the kallisto index file for the appropriate transcriptome and path where to save the TCC matrix.

#### Barcode analysis and selection

The workflow operates on demultiplexed Chromium-prepared sequencing samples (the raw barcode and read files can be converted to FASTQ using the [cellranger demux](http://software.10xgenomics.com/single-cell/pipelines/latest/what-is-cell-ranger) 10x software). Once the gzipped FASTQ files have been obtained, the first step in our workflow is to identify "true" barcodes, and to error correct barcodes that are close to true barcodes, yet associated with sufficiently low read coverage to be confidently identified as containing an error. The script `get_cell_barcodes.py` in the [source](https://github.com/lakigigar/scRNA-Seq-TCC-prep/tree/master/source) directory performs the identification and error correction and is called with `python get_cell_barcodes.py config.json`.

While `get_cell_barcodes.py` can be run from the command line, we strongly encourage users to instead perform this step using the Jupyter Notebook `10xGet_cell_barcodes.ipynb` in the [notebooks](https://github.com/lakigigar/scRNA-Seq-TCC-prep/tree/master/notebooks) directory. The interactive notebook produces summary statistics and figures that are useful for both quality control and for the setting of parameters for error correction. 

#### Cell file generation

Once barcodes have been identified and (some) erroneous barcodes corrected, the next step is to generate individual read and UMI files for each cell for processing by kallisto. This can be performed with the command `python error_correct_and_split.py config.json`. 

#### Pseudoalignment

The computation of transcript compatibility counts is performed using kallisto by running `python compute_TCCs.py config.json` followed by `python prep_TCC_matrix.py config.json`. The first script runs kallisto and the second step computes a pairwise distance matrix between cells that is essential for analysis. The result of running the two scripts is the generation of three files needed for analysis: `TCC_matrix.dat`, `pwise_dist_L1.dat` and `nonzero_ec.dat`. 

Note that the entire workflow can be run using the master script `10xDetect_and_Prep.py` although as explained above __we recommend examining the barcode data using the Jupyter Notebook `10xGet_cell_barcodes.ipynb`__. After the barcode analysis and selection step, the rest of the workflow can be completed by running `python 10xPrepData.py config.json`.
 
#### Analysis

The `TCC_matrix.dat` file contains a matrix that specifies, for each cell, a list of transcript sets with associated counts. Those counts, called transcript compatibility counts, are explained in [Ntranos _et al._ 2016](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0970-8). They are the starting point for downstream analysis of the data.

The analysis workflow for an experiment will depend on the specifics of the data and the questions associated with it. To help users get started, we have provided two examples based on datasets distributed by 10X: an experiment with both human and mouse cells and an analysis of peripheral blood mononuclear cells.

## Contributions

This 10X Chromium 3' digital expression processing workflow was designed and implemented by Vasilis Ntranos with some input from Lior Pachter. P&aacute;ll Melsted added the `--umi` option to kallisto which allows for deduplicating reads using associated unique molecular identifiers (UMIs).
