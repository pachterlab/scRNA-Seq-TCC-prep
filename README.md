# single cell RNA-Seq TCC prep

This repository contains scripts needed to generate transcript compatibility matrices from single-cell RNA-Seq data. This involves error-correction of barcodes, collapsing of UMIs and pseudoalignment of reads to a transcriptome to obtain transcript compatibility counts. The scripts utilize [kallisto](http://pachterlab.github.io/kallisto) for pseudoalignment.

We currently support the 10X Chromium technology; support for more technologies is underway.

