* See `https://github.com/benjjneb/dada2/issues/1164` for discussion on losing reads in CCS
* Setting minQ=0 and maxEE to 4 can rescue a lot of reads
  * minQ=0 and maxEE=2 is similar to minQ=3 and maxEE=2 in MSA1003
* Using lima to demux first instead of using DADA2 can give a lot more reads. Downside is
  that we have to modify the `run_dada_ccs.R` script in QIIME plugin.
  * On this end, can we use custom script and export PATH during run time so we won't need to replace the default?
* Cluster ASV into OTUs using `https://docs.qiime2.org/2022.2/plugins/available/vsearch/cluster-features-open-reference/` 
* Implement quality filtering and filter for Q30 reads up-front
* Visualization of all results in a nice report
* Add Krona plot
* How to BLAST ASVs for strain level assignment
* Phylo tree and PCoA
* Dockerize. Remember to replace run_dada_ccs.R so it works with trimmed sequences
* Kraken 2 for reads-level assignment?
* Generate a metadata TSV if none is provided
* Use csvtk to dynamically filter for rarefaction depth (Covering >50% of samples)
* For Krona plot pip install, can we clone the repo locally and install from that
to freeze it (avoiding risk of repo going down)
* Output number of filtered reads in reads QC step
