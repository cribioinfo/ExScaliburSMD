# ExScaliburSMD: Scalable and Accurate Detection of Somatic Mutations from Whole Exome Sequencing in the Cloud #

ExScaliburSMD is a somatic mutation pipeline for whole-exome sequencing data developed by the [Center for Research Informatics](cri.uchicago.edu) at the University of Chicago. It is designed to be flexible, robust, and to take advantage of multiple alignment and variant-calling algorithms.

The pipeline is designed to be

* **Scalable** *Manages analysis of tens to hundreds of samples*
* **Accurate** *Over 99% sensitivity and accuracy in variant detection*
* **Automate** *Full analysis pipeline from QC to annotated somatic variants with submission of one master script*
* **Comprehensive** *Implements multiple aligners and callers allowing for comparison and integration of different variant sets*
* **Flexible** *Fine control of analysis modules, parameters, and environmental settings*
* **Robust** *Extensive checkpoints and error detection features*
* **Reproducible** *Quick restart of analysis and sharing of results*
* **Real-time** *Easy access and monitor of the pipeline progress*
* **Transferable** *Runs on a desktop, work station, HPC and cloud with one single switch*

The pipeline is implemented in 

* [BigDataScript](http://pcingola.github.io/BigDataScript/) *Scripting language for data pipelines*
* [Python](https://www.python.org/) *Utilities*

The SMD pipeline is part of the ExScalibur suite (CRI, University of Chicago): https://exscalibur.cri.uchicago.edu.

# Documentation #

Please see the [Wiki](https://bitbucket.org/cribioinformatics/exscalibursmd/wiki) for full documentation.

# Communication #

For pipeline questions, please contact Kyle Hernandez (kmhernan at uchicago dot edu).
For other general questions, please contact the CRI bioinformatics team (bioinformatics at bsd dot uchicago dot edu).

# Release #

* **Version 0.5.0** *2015-03-31*

# Development #

* Add checkpoints
* Add hg38 genome

# More #

* [ExScaliburGMD](https://bitbucket.org/cribioinformatics/exscaliburgmd) *germline detection*
* [ExScaliburViz](https://bitbucket.org/cribioinformatics/exscaliburviz) *visualization of data analysis report*