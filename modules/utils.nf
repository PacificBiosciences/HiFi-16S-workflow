// Helper function to display help message
def helpMessage() {
    log.info"""
    Usage:
    This pipeline takes in the standard sample manifest and metadata file used in
    QIIME 2 and produces QC summary, taxonomy classification results and visualization.

    For samples TSV, two columns named "sample-id" and "absolute-filepath" are
    required. For metadata TSV file, at least two columns named "sample_name" and
    "condition" to separate samples into different groups.

    nextflow run main.nf --input samples.tsv --metadata metadata.tsv \\
        --dada2_cpu 8 --vsearch_cpu 8

    By default, sequences are first trimmed with cutadapt. If adapters are already trimmed, you can skip 
    cutadapt by specifying "--skip_primer_trim".

    Other important options:
    --front_p    Forward primer sequence. Default to F27. (default: AGRGTTYGATYMTGGCTCAG)
    --adapter_p    Reverse primer sequence. Default to R1492. (default: AAGTCGTAACAAGGTARCY)
    --filterQ    Filter input reads above this Q value (default: 20).
    --downsample    Limit reads to a maximum of N reads if there are more than N reads (default: off)
    --max_ee    DADA2 max_EE parameter. Reads with number of expected errors higher than
                this value will be discarded (default: 2)
    --minQ    DADA2 minQ parameter. Reads with any base lower than this score 
                will be removed (default: 0)
    --min_len    Minimum length of sequences to keep (default: 1000)
    --max_len    Maximum length of sequences to keep (default: 1600)
    --pooling_method    QIIME 2 pooling method for DADA2 denoise see QIIME 2 
                        documentation for more details (default: "pseudo", alternative: "independent") 
    --maxreject    max-reject parameter for VSEARCH taxonomy classification method in QIIME 2
                    (default: 100)
    --maxaccept    max-accept parameter for VSEARCH taxonomy classification method in QIIME 2
                    (default: 100)
    --min_asv_totalfreq    Total frequency of any ASV must be above this threshold
                            across all samples to be retained. Set this to 0 to disable filtering
                            (default 5)
    --min_asv_sample    ASV must exist in at least min_asv_sample to be retained. 
                        Set this to 0 to disable. (default 1)
    --vsearch_identity    Minimum identity to be considered as hit (default 0.97)
    --rarefaction_depth    Rarefaction curve "max-depth" parameter. By default the pipeline
                            automatically select a cut-off above the minimum of the denoised 
                            reads for >80% of the samples. This cut-off is stored in a file called
                            "rarefaction_depth_suggested.txt" file in the results folder
                            (default: null)
    --omegac    OMEGA_C parameter for DADA2. (default: 1e-40)
    --learn_error_sample    Use this FASTQ to learn error model for DADA2 (e.g. low complexity
                            control in the same sequencing run with same library prep)
    --dada2_cpu    Number of threads for DADA2 denoising (default: 8)
    --vsearch_cpu    Number of threads for VSEARCH taxonomy classification (default: 8)
    --cutadapt_cpu    Number of threads for primer removal using cutadapt (default: 16)
    --outdir    Output directory name (default: "results")
    --vsearch_db	Location of VSEARCH database (e.g. silva-138-99-seqs.qza can be
                    downloaded from QIIME database)
    --vsearch_tax    Location of VSEARCH database taxonomy (e.g. silva-138-99-tax.qza can be
                    downloaded from QIIME database)
    --silva_db   Location of Silva 138 database for taxonomy classification 
    --gtdb_db    Location of GTDB r202 for taxonomy classification
    --gg2_db    Location of GreenGenes2 database for taxonomy classification
    --db_to_prioritize    Prioritize this database for taxonomy assignment. Must be one of the 
                          following: GG2, GTDB, Silva (default: GG2)
    --skip_primer_trim    Skip all primers trimming (switch off cutadapt and DADA2 primers
                            removal) (default: trim with cutadapt)
    --skip_nb    Skip Naive-Bayes classification (only uses VSEARCH) (default: false)
    --colorby    Columns in metadata TSV file to use for coloring the MDS plot
                in HTML report (default: condition)
    --run_picrust2    Run PICRUSt2 pipeline. Note that pathway inference with 16S using PICRUSt2
                        has not been tested systematically (default: false)
    --download_db    Download databases needed for taxonomy classification only. Will not
                    run the pipeline. Databases will be downloaded to a folder "databases"
                    in the Nextflow pipeline directory.
    --publish_dir_mode    Outputs mode based on Nextflow "publishDir" directive. Specify "copy"
                            if requires hard copies. (default: symlink)
    --version    Output version
    """
}

process write_log {
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest" 
    publishDir "$params.outdir", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    val(logs)

    output:
    path "parameters.txt"

    script:
    """
    echo '$logs' > parameters.txt
    """
}

process download_db {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$projectDir/databases", mode: "copy"
    label 'cpu_def'

    output:
    path "*"

    script:
    """
    echo "Downloading GG2 sequences for Naive Bayes classifier..."
    wget -N --content-disposition 'https://zenodo.org/records/14169078/files/gg2_2024_09_toSpecies_trainset.fa.gz?download=1'
    
    echo "Downloading SILVA sequences for Naive Bayes classifier..."
    wget -N --content-disposition 'https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1'
    
    echo "Downloading GTDB sequences for Naive Bayes classifier..."
    wget -N --content-disposition 'https://zenodo.org/records/13984843/files/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz?download=1'
    
    echo "Downloading SILVA sequences and taxonomies for VSEARCH..."
    wget -N --content-disposition 'https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza'
    wget -N --content-disposition 'https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza'

    echo "Downloading GTDB database for VSEARCH..."
    wget 'https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_all/ssu_all_r220.fna.gz'
    zgrep "^>" ssu_all_r220.fna.gz |  awk '{print gensub(/^>(.*)/, "\\\\1", "g", \$1),gensub(/^>.* (d__.*) \\[.*\\[.*\\[.*/, "\\\\1", "g", \$0)}' OFS=\$'\\t' > ssu_all_r220.taxonomy.tsv
    gunzip ssu_all_r220.fna.gz
    qiime tools import --type 'FeatureData[Sequence]' --input-path ssu_all_r220.fna --output-path GTDB_ssu_all_r220.qza
    qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ssu_all_r220.taxonomy.tsv --output-path GTDB_ssu_all_r220.taxonomy.qza
    rm ssu_all_r220.fna
    rm ssu_all_r220.taxonomy.tsv

    echo "Downloading GreenGenes2 database for VSEARCH..."
    wget -N --content-disposition 'http://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.full-length.fna.qza'
    wget -N --content-disposition 'http://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.tax.qza'
    """
}

process picrust2 {
    conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
    container "kpinpb/pb-16s-nf-tools:latest"
    label 'cpu32'
    publishDir "$params.outdir/results", mode: params.publish_dir_mode

    input:
    path dada2_asv
    path vsearch_biom

    output:
    path "picrust2"

    script:
    """
    picrust2_pipeline.py -s $dada2_asv -i $vsearch_biom \
        -o picrust2 \
        --in_traits EC,KO \
        --verbose -p $task.cpus
    """
}

// Function to validate input parameters
def validateParameters() {
    if (!params.input) {
        exit 1, "Input file not specified!"
    }
    if (!params.metadata) {
        exit 1, "Metadata file not specified!"
    }
    if (params.min_len >= params.max_len) {
        exit 1, "Minimum length must be less than maximum length!"
    }
}

// Function to print pipeline header
def printHeader() {
    log.info"""
    ===========================================
    PacBio 16S Analysis Pipeline
    ===========================================
    """
}

// Function to print pipeline completion
def printComplete() {
    log.info"""
    Pipeline completed!
    ===========================================
    """
} 