# v0.4 changelog
* Fix some bugs wrongly displaying Naive-Bayes statistics as VSEARCH stats in HTML
* Allow manually specifying primers with `--forward_p` and `--adapter_p` for
  `cutadapt`
* Allow skipping Naive-Bayes classification via the `--skip_nb` parameter. This
  allows user to easily use any database with VSEARCH and generate report. See FAQ.
* Fixed bug preventing single sample heatmap from displaying
* Save versions and parameters in `parameters.txt` for reproducibility
* Updated QIIME 2 environment YAML.
* Fixed issue 17 reported by Stephane. Nextflow and Mamba do not need to be installed
  at base conda environment.

# v0.3 changelog
* Updated naive classifier approach to use GTDB r207 (updated from r202) and make GTDB r207 the first priority so it's more consistent with VSEARCH-based classification. I.e. in the naive classification approach, ASVs are classified firstly by GTDB before Silva and RefSeq (In v0.2, Silva was the first priority).
* New feature: Optional  `PICRUSt2` for pathway inference.
* Databases are now downloaded from a central Zenodo DOI to make process more robust: https://zenodo.org/record/6912512 

# v0.2 changelog
* Starts tracking version, no main changes compared to previous version except that the pipeline is now versioned.
