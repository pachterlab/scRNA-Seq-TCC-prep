# single-cell RNA-Seq TCC prep

This repository contains scripts needed to generate transcript compatibility matrices from single-cell RNA-Seq data. Included is error-correction of barcodes, collapsing of UMIs and pseudoalignment of reads to a transcriptome to obtain transcript compatibility counts. The scripts utilize [kallisto](http://pachterlab.github.io/kallisto) for pseudoalignment.

We currently support the 10X Chromium technology; support for more technologies is underway.

## Instructions for processing 10X Chromium 3' digital expression data

### Workflow organization

The processing workflow consists of four steps: 

0. Preparation of a config.json file that contains the parameters needed for the processing.

1. Identification of "true" cell barcodes according to read coverage followed by error correction when possible.

2. Creation of read/UMI files for each cell.

3. Pseudoalignment of reads associated with each cell using __kallisto__, deduplication according to UMIs, and generation of transcript compatibility counts for each cell.

Following the pre-processing, the transcript compatibility counts (TCC) matrix can be analyzed using a Jupyter Notebook. 

## Creation of the config.json file

## Barcode analysis and selection

## Cell file generation

## Pseudoalignment

## Analysis

The [getting started](http://pachterlab.github.io/kallisto/10xstarting.html) tutorial explains how to process the small example in the [example_dataset](https://github.com/lakigigar/scRNA-Seq-TCC-prep/tree/master/example_dataset) directory. 