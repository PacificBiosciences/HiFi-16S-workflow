process qiime2_phylogeny_diversity {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results/phylogeny_diversity", mode: params.publish_dir_mode
    label 'cpu8' 

    input:
    path metadata
    path asv_seq
    path asv_freq
    val(rarefaction_depth)
    val(n_sample)

    output:
    path "*.qza"
    path "core-metrics-diversity"
    path "phylotree_mafft_rooted.nwk", emit:qiime2_tree
    path "core-metrics-diversity/bray_curtis_distance_matrix.tsv", emit:bray_mat
    path "core-metrics-diversity/weighted_unifrac_distance_matrix.tsv", emit:wunifrac_mat
    path "core-metrics-diversity/unweighted_unifrac_distance_matrix.tsv", emit:unifrac_mat

    script:
    if (n_sample > 1)
    """
    qiime phylogeny align-to-tree-mafft-fasttree \
        --p-n-threads $task.cpus \
        --i-sequences $asv_seq \
        --o-alignment mafft_alignment.qza \
        --o-masked-alignment mafft_alignment_masked.qza \
        --o-tree phylotree_mafft_unrooted.qza \
        --o-rooted-tree phylotree_mafft_rooted.qza

    qiime tools export --input-path phylotree_mafft_rooted.qza \
        --output-path ./
    mv tree.nwk phylotree_mafft_rooted.nwk

    qiime diversity core-metrics-phylogenetic \
        --p-n-jobs-or-threads $task.cpus \
        --i-phylogeny phylotree_mafft_rooted.qza \
        --i-table $asv_freq \
        --m-metadata-file $metadata \
        --p-sampling-depth $rarefaction_depth \
        --output-dir ./core-metrics-diversity

    # Export various matrix for plotting later
    qiime tools export --input-path ./core-metrics-diversity/bray_curtis_distance_matrix.qza \
        --output-path ./core-metrics-diversity
    mv ./core-metrics-diversity/distance-matrix.tsv \
        ./core-metrics-diversity/bray_curtis_distance_matrix.tsv
    qiime tools export --input-path ./core-metrics-diversity/weighted_unifrac_distance_matrix.qza \
        --output-path ./core-metrics-diversity
    mv ./core-metrics-diversity/distance-matrix.tsv \
        ./core-metrics-diversity/weighted_unifrac_distance_matrix.tsv
    qiime tools export --input-path ./core-metrics-diversity/unweighted_unifrac_distance_matrix.qza \
        --output-path ./core-metrics-diversity
    mv ./core-metrics-diversity/distance-matrix.tsv \
        ./core-metrics-diversity/unweighted_unifrac_distance_matrix.tsv
    """
    else
    """
    qiime phylogeny align-to-tree-mafft-fasttree \
        --p-n-threads $task.cpus \
        --i-sequences $asv_seq \
        --o-alignment mafft_alignment.qza \
        --o-masked-alignment mafft_alignment_masked.qza \
        --o-tree phylotree_mafft_unrooted.qza \
        --o-rooted-tree phylotree_mafft_rooted.qza

    qiime tools export --input-path phylotree_mafft_rooted.qza \
        --output-path ./
    mv tree.nwk phylotree_mafft_rooted.nwk
    mkdir ./core-metrics-diversity
    touch ./core-metrics-diversity/bray_curtis_distance_matrix.tsv \
        ./core-metrics-diversity/weighted_unifrac_distance_matrix.tsv \
        ./core-metrics-diversity/unweighted_unifrac_distance_matrix.tsv
    """
}

process barplot {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_tab
    path tax
    path metadata
    val(oname)

    output:
    path "${oname}"

    script:
    """
    qiime taxa barplot --i-table $asv_tab --i-taxonomy $tax \
        --m-metadata-file $metadata \
        --o-visualization $oname
    """
}

process barplot_nb {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_tab
    path tax
    path metadata
    val(oname)

    output:
    path "${oname}"

    script:
    """
    qiime taxa barplot --i-table $asv_tab --i-taxonomy $tax \
        --m-metadata-file $metadata \
        --o-visualization $oname
    """
}

process html_rep {
    conda (params.enable_conda ? "$projectDir/env/pb-16s-vis-conda.yml" : null)
    container "kpinpb/pb-16s-vis:latest"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path tax_freq_tab_tsv
    path metadata
    path sample_manifest
    path dada2_qc
    path reads_qc
    path summarised_reads_qc
    path cutadapt_summary_qc
    path vsearch_tax_tsv
    path bray_mat
    path unifrac_mat
    path wunifrac_mat
    val(colorby)
    path post_trim_readstats
    path rmd_vis_biom_script
    path rmd_helper

    output:
    path "visualize_biom.html", emit: html_report

    script:
    if (params.enable_container)
    """
    export R_LIBS_USER="/opt/conda/envs/pb-16S-vis/lib/R/library"
    cp $rmd_vis_biom_script vis_script.rmd
    Rscript -e 'rmarkdown::render("vis_script.rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat", post_trim_readstats="$post_trim_readstats"), output_dir="./")'
    mv vis_script.html visualize_biom.html
    """
    else
    """
    cp $rmd_vis_biom_script vis_script.rmd
    Rscript -e 'rmarkdown::render("vis_script.rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat", post_trim_readstats="$post_trim_readstats"), output_dir="./")'
    mv vis_script.html visualize_biom.html
    """
}

process html_rep_skip_cutadapt {
    conda (params.enable_conda ? "$projectDir/env/pb-16s-vis-conda.yml" : null)
    container "kpinpb/pb-16s-vis:latest"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path tax_freq_tab_tsv
    path metadata
    path sample_manifest
    path dada2_qc
    path reads_qc
    path summarised_reads_qc
    val(cutadapt_summary_qc)
    path vsearch_tax_tsv
    path bray_mat
    path unifrac_mat
    path wunifrac_mat
    val(colorby)
    val(post_trim_readstats)
    path rmd_vis_biom_script
    path rmd_helper

    output:
    path "visualize_biom.html", emit: html_report

    script:
    if (params.enable_container)
    """
    export R_LIBS_USER="/opt/conda/envs/pb-16S-vis/lib/R/library"
    cp $rmd_vis_biom_script vis_script.rmd
    Rscript -e 'rmarkdown::render("vis_script.rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat", post_trim_readstats="$post_trim_readstats"), output_dir="./")'
    mv vis_script.html visualize_biom.html
    """
    else
    """
    cp $rmd_vis_biom_script vis_script.rmd
    Rscript -e 'rmarkdown::render("vis_script.rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat", post_trim_readstats="$post_trim_readstats"), output_dir="./")'
    mv vis_script.html visualize_biom.html
    """
}

process krona_plot {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'
    errorStrategy = 'ignore'

    input:
    path asv_freq
    path taxonomy

    output:
    path "krona.qzv"
    path "krona_html"

    script:
    """
    pip install git+https://github.com/kaanb93/q2-krona.git
    conda install -c bioconda krona
    qiime krona collapse-and-plot --i-table $asv_freq --i-taxonomy $taxonomy --o-krona-plot krona.qzv
    qiime tools export --input-path krona.qzv --output-path krona_html
    """
} 