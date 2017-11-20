# Single-cell RNA-Seq TCC prep -- SC3Pv2

This repository contains scripts needed to generate transcript compatibility count (TCC) matrices from 10X Chromium 3' single-cell RNA-Seq data. Included is error-correction of barcodes, collapsing of UMIs and pseudoalignment of reads to a transcriptome to obtain transcript compatibility counts. The scripts utilize [kallisto](http://pachterlab.github.io/kallisto) for pseudoalignment.

The workflow operates on the *demultiplexed FASTQ files* that have been obtained from the Illumina BCL data using [cellranger demux](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/demultiplex).
The *output* is a matrix that specifies, for each cell, a list of transcript sets with associated counts. Those counts, called transcript compatibility counts, are explained in [Ntranos _et al._ 2016](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0970-8) and are the starting point for downstream analysis of the data.

Note: The scripts are modified to support data generated with 10X Chromium 3' SC3Pv2 chemistry (cellranger version ≥ 1.2.0). 

## Instructions for processing 10X Chromium 3' digital expression data

#### Getting started

The [getting started](http://) tutorial (TBD) explains how to process the small example in the [example_dataset](https://github.com/lakigigar/scRNA-Seq-TCC-prep/tree/10xTCCprep-SC3Pv2/example_dataset) directory. This is a good starting point to make sure that the necessary programs are correctly installed. Note that you will need kallisto (≥ 0.43.0), python (≥ 3.x.x) and Juypter Notebook (≥ 4.0.6) installed (the Jupyter requirement is not strictly necessary but highly recommended). 


#### Workflow organization

The processing workflow consists of four steps: 

0. Preparation of a configuration `config.json` file that contains the parameters needed for the processing.
1. Identification of "true" cell barcodes according to read coverage followed by error correction when possible.
2. Creation of read/UMI files for each cell.
3. [Pseudoalignment](http://www.nature.com/nbt/journal/v34/n5/abs/nbt.3519.html) of reads associated with each cell using __kallisto__, deduplication according to UMIs, and generation of transcript compatibility counts (TCCs) for each cell. 

The entire workflow can be completed using the master script `10xGet_TCCs.py`, by running `python 10xGet_TCCs.py config.json`.

#### Creation of the configuration file

Parameters needed to run the processing require specification of a `config.json` file. The following parameters need to be specified:

- NUM_THREADS: the number of threads available for processing.
- FASTQ_DIRS: this must contain the list of paths to the (demultiplexed) FASTQ files from the sequencing for all the samples that need to be jointly processed. Note that our workflow does not currently demultiplex reads and you may have to do so with 10X's software;
- SOURCE_DIR: path to the source directory that contains the .py scripts
- sample_names: The sample names used to get the FASTQ files. The fastqs that start with the same sample name will be considered to have been obtained from the same library and will processed together.
- EXP_CELLS: this parameter specifies the expected number of cells in each sample and is used in the determination of the number o cells in the experiment from reads coverage data. Note that EXP_CELLS should be given as a list whose entries correspond to sample_names.
- use_precomputed_barcodes: Binary parameter, default is zero. If set to 1, the workflow will skip cell detection and use a given set of cell barcodes (This may be useful if one wants to reanalyze a gene counts dataset and obtain TCCs with matching cell barcodes).
- precomputed_barcode_files: the list of paths to the precomputed barcode files for each sample. See below for file format requirements. 
- SAVE_DIR: path to a directory where intermediate files will be saved.
- dmin: the minimum distance between barcodes needed for error correction to be performed.
- BARCODE_LENGTH: length of the barcodes.
- OUTPUT_DIR: directory in which to output results.
- kallisto: path to the binary for kallisto, location of the kallisto index file for the appropriate transcriptome and path where to save the TCC matrix.

#### Barcode analysis and selection

Once the gzipped FASTQ files have been obtained, the first step in our workflow is to identify "true" barcodes, and to error correct barcodes that are close to true barcodes, yet associated with sufficiently low read coverage to be confidently identified as containing an error. The script `get_cell_barcodes_v2.py` in the [source](https://github.com/lakigigar/scRNA-Seq-TCC-prep/tree/10xTCCprep-SC3Pv2/source) directory performs the identification and error correction and is called with `python get_cell_barcodes_v2.py config.json`.

While `get_cell_barcodes_v2.py` can be run from the command line, we strongly encourage users to instead perform this step using the Jupyter Notebook `10xDetect_cell_barcodes.ipynb` in the [notebooks](https://github.com/lakigigar/scRNA-Seq-TCC-prep/tree/10xTCCprep-SC3Pv2/notebooks) directory. The interactive notebook produces summary statistics and figures that are useful for both quality control and for the setting of parameters for error correction. 

#### Cell file generation

Once barcodes have been identified and (some) erroneous barcodes corrected, the next step is to generate individual read and UMI files for each cell for processing by kallisto. This can be performed with the command `python error_correct_and_split_v2.py config.json`. 

#### Pseudoalignment

The computation of transcript compatibility counts is performed using kallisto by running `python compute_TCCs.py config.json`. 

Note that the entire workflow can be run using the master script `10xGet_TCCs.py` although as explained above __we recommend examining the barcode data first using the Jupyter Notebook `10xDetect_cell_barcodes.ipynb`__. After the barcode analysis and selection step, the rest of the workflow can still be completed by running `python 10xGet_TCCs.py config.json`, as the workflow will detect whether the `10xDetect_cell_barcodes.ipynb` has been executed and automatically skip the barcode detection step.
 

## Contributions

This 10X Chromium 3' digital expression processing workflow was designed and implemented by Vasilis Ntranos with some input from Lior Pachter. P&aacute;ll Melsted added the `--umi` option to kallisto which allows for deduplicating reads using associated unique molecular identifiers (UMIs). In [sc_read_kallisto_wrapper](https://github.com/vibbits/sc_read_kallisto_wrapper) Copyright (c) BITS VIB 2017, Alexander Botzki modified our original scripts to handle SC3Pv2 chemistry and replaced the cell number finding algorithm with one that is more similar to what Cell Ranger does (expected number of cells specified as input). In this version we integrated these modifications with further updates in our workflow (to add functionality and enable faster processing).
