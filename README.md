# HetchromSilencing

This repository contains the computational pipeline of the paper "[Heterochromatin reduces RNA Polymerase II occupancy, transcriptional efficiency and RNA stability]. It contains the scripts and software necessary for reproducing the results in the paper.

## Overview

Current data suggest that heterochromatic silencing is a combination of transcriptional silencing and RNA degradation, however, the contribution of each pathway to silencing is not known. In this study we analyzed RNA Polymerase II (Pol II) occupancy, nascent RNA and steady state RNA levels to quantify the contribution of these pathways to heterochromatic silencing in fission yeast. We found that transcriptional silencing consists of two components, reduced RNA Pol II accessibility and reduced transcriptional efficiency. Using our data we determined and quantified which heterochromatic complexes contribute to reduced RNA Pol II occupancy, reduced transcriptional efficiency and reduced RNA stability.

In our paper we analyzed next-generation sequencing data for RNA Polymerase II (ChIP-seq), Pol II bound nascent RNA (RIP-seq), and steady state RNA levels (pA RNA-seq) in S. pombe fission yeast.

The [notebooks](https://github.com/canzarlab/HetchromSilencing/tree/master/notebooks) folders contains Jupyter notebooks for [preprocessing](https://github.com/canzarlab/HetchromSilencing/tree/master/notebooks/preprocess), [analysing](https://github.com/canzarlab/HetchromSilencing/tree/master/notebooks/analysis), and [visualizing the analysis results](https://github.com/canzarlab/HetchromSilencing/tree/master/notebooks/plots).


## Software requirements

The following programs are required to run the scripts in this repository.

* [STAR](https://github.com/alexdobin/STAR)
* [Samtools](http://www.htslib.org/download/)

To run the scripts in the Jupyter notebooks, the following Python modules are required.

* [numpy](http://docs.scipy.org/doc/numpy-1.10.1/user/install.html)
* [pandas](https://pandas.pydata.org/)
* [scikit-learn](http://scikit-learn.org/stable/install.html)
* [matplotlib](http://matplotlib.org/users/installing.html)
* [pysam](https://github.com/pysam-developers/pysam)