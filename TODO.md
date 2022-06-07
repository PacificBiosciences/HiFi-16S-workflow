# TODO

* Document all outputs and clean up workflow names
* Cite all databases and tools and provide single Zenodo location for all
  databases.
* Split ASV and sample freq filtering from denoise step, so that
  changing filtering parameters won't cause denoise to restart
* For Krona plot pip install, can we clone the repo locally and install from that
* Add LotuS2 for OTU-based analysis? Or just add a guide on GitHub
* Allow using just VSEARCH result and make naive bayes optional
* Dockerize
* See `https://github.com/benjjneb/dada2/issues/1164` for discussion on losing reads in CCS
  * Setting minQ=0 and maxEE to 4 can rescue a lot of reads
    * minQ=0 and maxEE=2 is similar to minQ=3 and maxEE=2 in MSA1003
  * Expose minQ?
* Clean up import qiime manifest do not need to put samples_demux_rate file
* Parallelize classify step by splitting the ASV FASTA files into fixed number of ASV chunks
* Guide on how to BLAST ASVs for strain level assignment

# DONE

* DONE: Using lima to demux first instead of using DADA2 can give a lot more reads. Downside is
  that we have to modify the `run_dada_ccs.R` script in QIIME plugin.
  * DONE: On this end, can we use custom script and export PATH during run time so we won't need to replace the default?
* DONE: Implement quality filtering and filter for Q30 reads up-front
* DONE: Visualization of all results in a nice report
* DONE: Add Krona plot
* DONE: Use csvtk to dynamically filter for rarefaction depth (Covering >50% of samples)
to freeze it (avoiding risk of repo going down)
* DONE: Output number of filtered reads in reads QC step
* DONE: Make kraken2 reads classification optional
  * With higher confidence (0.5 and above), the assignment rate was found to be quite low.
  So deprecated 2022-5-31. Code still in there if needed to enable in the future
* ABANDONED: Kraken 2 for reads-level assignment?
  * 2022-5-17: SILVA database does not have species in the taxonomy file. Will
  require effort to manually create a database suitable for species assignment.
  * Try NCBI database?
* DONE: Table too large with high amount of samples will cause HTML report to fail to load
  * Don't show full tables. No point.
* WORKAROUND: Deal with uncultured bacterium and metagenome in SILVA?
  * Database curation is hard. Right now, if can be classified with GTDB or RefSeq, use those,
  otherwise keep them as is (the metagenome samples will likely help to "absorb" false-positive hits)
* DONE: Document that cutadapt will orientate sequences making it usable for naive bayesian classifier and
  without that the results will be weird
* DONE: Remove ASV in only one samples and less than 5 reads
* DONE: Handle one sample workflow
* DONE: Phylo tree and PCoA and diversity analysis
