* See `https://github.com/benjjneb/dada2/issues/1164` for discussion on losing reads in CCS
  * Setting minQ=0 and maxEE to 4 can rescue a lot of reads
    * minQ=0 and maxEE=2 is similar to minQ=3 and maxEE=2 in MSA1003
  * Expose minQ?
* DONE: Using lima to demux first instead of using DADA2 can give a lot more reads. Downside is
  that we have to modify the `run_dada_ccs.R` script in QIIME plugin.
  * DONE: On this end, can we use custom script and export PATH during run time so we won't need to replace the default?
* Cluster ASV into OTUs using `https://docs.qiime2.org/2022.2/plugins/available/vsearch/cluster-features-open-reference/` 
* DONE: Implement quality filtering and filter for Q30 reads up-front
* DONE: Visualization of all results in a nice report
* DONE: Add Krona plot
* How to BLAST ASVs for strain level assignment
* Phylo tree and PCoA and diversity analysis
* Dockerize
* Kraken 2 for reads-level assignment?
  * 2022-5-17: SILVA database does not have species in the taxonomy file. Will
  require effort to manually create a database suitable for species assignment.
  * Try NCBI database?
* Generate a metadata TSV if none is provided
* DONE: Use csvtk to dynamically filter for rarefaction depth (Covering >50% of samples)
* For Krona plot pip install, can we clone the repo locally and install from that
to freeze it (avoiding risk of repo going down)
* DONE: Output number of filtered reads in reads QC step
* Deal with uncultured bacterium and metagenome in SILVA?
* Document that cutadapt will orientate sequences making it usable for naive bayesian classifier and
  without that the results will be weird
* Clean up import qiime manifest do not need to put samples_demux_rate file
* Table too large with high amount of samples will cause HTML report to fail to load
* Count number of reads in top samples
* Make kraken2 reads classification optional
