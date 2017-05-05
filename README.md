# Genetic basis of sex determination in a parasitoid wasp species
Master project for the Molecular Life Sciences master in Bioinformatics, University of Lausanne

### Cyril Matthey-Doret
### Supervisor: Casper Van der Kooi
### Director: Tanja Schwander
---
In this project, we use restriction-site associated DNA-sequencing (RAD-seq) and build a custom pipeline to locate and identify the complementary sex determination (CSD) locus/loci in the parasitoid wasp species _Lysiphlebus fabarum_.

This repository contains a pipeline to map the reads using BWA and build loci using the different components of STACKS with optimal parameters. It was written to run on a cluster with LSF.

To run the pipeline:
* Replace paths accordingly in ```config.mk```
* Download the data (not available yet) into the main folder
* ```cd``` to the main folder
* Untar the data using ```tar -xzvf data```
* Type ```make``` to run the pipeline

The pipeline will later implement association mapping and more.
