process learn_error {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/dada2", mode: params.publish_dir_mode
    cpus params.dada2_cpu

    input:
    path error_sample 
    path learnError_script

    output:
    path "errorfun.rds", emit: dada2_error_model
    path "plot_error_model.pdf"

    script:
    """
    Rscript --vanilla ${learnError_script} ${error_sample} $task.cpus
    """
}

process dada2_denoise {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/dada2", mode: params.publish_dir_mode
    cpus params.dada2_cpu

    input:
    path samples_qza
    path dada_ccs_script
    val(minQ)

    output:
    path "dada2-ccs_rep*.qza", emit: asv_seq
    path "dada2-ccs_table*.qza", emit: asv_freq
    path "dada2-ccs_stats*.qza", emit: asv_stats
    path "seqtab_nochim*.rds", emit: dada2_rds
    path "plot_error_model*.pdf"

    script:
    """
    mkdir -p dada2_custom_script
    cp $dada_ccs_script dada2_custom_script/run_dada.R
    sed -i 's/OMEGA_C=1e-40/OMEGA_C=$params.omegac/g' dada2_custom_script/run_dada.R
    chmod +x dada2_custom_script/run_dada.R
    export PATH="./dada2_custom_script:\$PATH"
    
    qiime dada2 denoise-ccs --i-demultiplexed-seqs $samples_qza \
        --o-table dada2-ccs_table.qza \
        --o-representative-sequences dada2-ccs_rep.qza \
        --o-denoising-stats dada2-ccs_stats.qza \
        --p-min-len $params.min_len --p-max-len $params.max_len \
        --p-max-ee $params.max_ee \
        --p-front 'none' \
        --p-adapter 'none' \
        --p-n-threads $task.cpus \
        --p-pooling-method \'$params.pooling_method\'
    """
}

process dada2_denoise_with_error_model {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/dada2", mode: params.publish_dir_mode
    cpus params.dada2_cpu

    input:
    path samples_qza
    path dada_ccs_script
    val(minQ)
    path error_model

    output:
    path "dada2-ccs_rep*.qza", emit: asv_seq
    path "dada2-ccs_table*.qza", emit: asv_freq
    path "dada2-ccs_stats*.qza", emit:asv_stats
    path "seqtab_nochim*.rds", emit: dada2_rds

    script:
    """
    mkdir -p dada2_custom_script
    cp $dada_ccs_script dada2_custom_script/run_dada.R
    sed -i 's/OMEGA_C=1e-40/OMEGA_C=$params.omegac/g' dada2_custom_script/run_dada.R
    sed -i '/### LEARN ERROR RATES ###/,/### PROCESS ALL SAMPLES ###/ s/^/#/' dada2_custom_script/run_dada.R
    sed -i '/#### LEARN ERROR RATES ###/a err <- readRDS("$error_model")' dada2_custom_script/run_dada.R
    chmod +x dada2_custom_script/run_dada.R
    export PATH="./dada2_custom_script:\$PATH"

    qiime dada2 denoise-ccs --i-demultiplexed-seqs $samples_qza \
        --o-table dada2-ccs_table.qza \
        --o-representative-sequences dada2-ccs_rep.qza \
        --o-denoising-stats dada2-ccs_stats.qza \
        --p-min-len $params.min_len --p-max-len $params.max_len \
        --p-max-ee $params.max_ee \
        --p-front 'none' \
        --p-adapter 'none' \
        --p-n-threads $task.cpus \
        --p-pooling-method \'$params.pooling_method\'
    """
}

process mergeASV {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/dada2", mode: params.publish_dir_mode
    cpus params.dada2_cpu

    input:
    path "dada2-ccs_rep*.qza"
    path "dada2-ccs_table*.qza"
    path "dada2-ccs_stats*.qza"

    output:
    path "dada2-ccs_rep_merged.qza", emit: asv_seq
    path "dada2-ccs_table_merged.qza", emit: asv_freq
    path "dada2-ccs_stats_merged.qza", emit: asv_stats

    script:
    """
    qiime feature-table merge \
        --i-tables dada2-ccs_table*.qza \
        --o-merged-table dada2-ccs_table_merged.qza
    
    qiime feature-table merge-seqs \
        --i-data dada2-ccs_rep*.qza \
        --o-merged-data dada2-ccs_rep_merged.qza

    for i in \$(ls dada2-ccs_stats*.qza);
    do
        qiime tools export --input-path \${i} \
        --output-path ./\${i%%.qza}_export/
    done

    cat ./*_export/stats.tsv | awk '{if (NR>2 && \$0~/^#|sample-id/) {next} else {print \$0}}' > merged_stats.tsv
    qiime tools import --input-path merged_stats.tsv --output-path dada2-ccs_stats_merged --type 'SampleData[DADA2Stats]'
    """
}

process filter_dada2 {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/dada2", mode: params.publish_dir_mode
    cpus params.dada2_cpu

    input:
    path asv_table
    path asv_seq
    val min_asv_totalfreq
    val min_asv_sample

    output:
    path "dada2-ccs_rep_filtered.qza", emit: asv_seq
    path "dada2-ccs_table_filtered.qza", emit: asv_freq
    path "dada2_ASV.fasta", emit:asv_seq_fasta

    script:
    """
    if [ $min_asv_sample -gt 0 ] && [ $min_asv_totalfreq -gt 0 ]
    then
        qiime feature-table filter-features \
            --i-table $asv_table \
            --p-min-frequency $min_asv_totalfreq \
            --p-min-samples $min_asv_sample \
            --p-filter-empty-samples \
            --o-filtered-table dada2-ccs_table_filtered.qza
    elif [ $min_asv_sample -gt 0 ] && [ $min_asv_totalfreq -eq 0 ]
    then
        qiime feature-table filter-features \
            --i-table $asv_table \
            --p-min-samples $min_asv_sample \
            --p-filter-empty-samples \
            --o-filtered-table dada2-ccs_table_filtered.qza
    elif [ $min_asv_sample -eq 0 ] && [ $min_asv_totalfreq -gt 0 ]
    then
        qiime feature-table filter-features \
            --i-table $asv_table \
            /para--p-min-frequency $min_asv_totalfreq \
            --p-filter-empty-samples \
            --o-filtered-table dada2-ccs_table_filtered.qza
    else
        mv $asv_table dada2-ccs_table_filtered.qza
    fi

    qiime feature-table filter-seqs \
        --i-data $asv_seq \
        --i-table dada2-ccs_table_filtered.qza \
        --o-filtered-data dada2-ccs_rep_filtered.qza

    qiime tools export --input-path dada2-ccs_rep_filtered.qza \
        --output-path .
    mv dna-sequences.fasta dada2_ASV.fasta
    """
}

process dada2_qc {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_stats
    path asv_freq
    path metadata

    output:
    path "dada2_stats.qzv", emit: dada2_stats
    path "dada2_table.qzv", emit: dada2_table
    path "stats.tsv", emit: dada2_stats_tsv
    path "dada2_qc.tsv", emit: dada2_qc_tsv
    env(int_rarefaction_d), emit: rarefaction_depth
    path "rarefaction_depth_suggested.txt"
    env(int_alpha_d), emit:alpha_depth

    script:
    """
    qiime metadata tabulate --m-input-file $asv_stats \
        --o-visualization dada2_stats.qzv

    qiime feature-table summarize --i-table $asv_freq \
        --o-visualization dada2_table.qzv \
        --m-sample-metadata-file $metadata

    qiime tools export --input-path dada2_table.qzv \
        --output-path dada2_table_summary

    sed -n '/<script id="table-data" type="application\\/json">/,/<\\/script>/p' dada2_table_summary/sample-frequency-detail.html | \
        sed '1d;\$d' | \
        jq -r '.Frequency | to_entries | .[] | [.key, .value] | @tsv' > \
        dada2_table_summary/sample-frequency-detail.tsv

    qiime tools export --input-path $asv_stats \
        --output-path ./

    number=`cat dada2_table_summary/sample-frequency-detail.tsv | wc -l`
    ninety=0.8
    if [ \${number} -le 2 ];
    then
        rarefaction_d=`cat dada2_table_summary/sample-frequency-detail.tsv | sort -k2 -nr | \
            tail -n1 | cut -f2`
        int_rarefaction_d=\${rarefaction_d%%.*}
        echo \${int_rarefaction_d} > rarefaction_depth_suggested.txt
        int_alpha_d=\${int_rarefaction_d}
    elif [ \${number} -gt 2 ] && [ \${number} -lt 5 ];
    then
        result=`awk -v n="\$number" -v p=0.8 'BEGIN{print int(n * p)}'`
        rarefaction_d=`cat dada2_table_summary/sample-frequency-detail.tsv | sort -k2 -nr | \
            head -n \${result} | tail -n1 | cut -f2`
        int_rarefaction_d=\${rarefaction_d%%.*}
        echo \${int_rarefaction_d} > rarefaction_depth_suggested.txt
        int_alpha_d=\${int_rarefaction_d}
    else
        result=`awk -v n="\$number" -v p=0.8 'BEGIN{print int(n * p)}'`
        rarefaction_d=`cat dada2_table_summary/sample-frequency-detail.tsv | sort -k2 -nr | \
            head -n \${result} | tail -n1 | cut -f2`
        int_rarefaction_d=\${rarefaction_d%%.*}
        echo \${int_rarefaction_d} > rarefaction_depth_suggested.txt
        alpha_res=`awk -v n="\$number" 'BEGIN{print int(n * 0.2)}'`
        alpha_d=`cat dada2_table_summary/sample-frequency-detail.tsv | sort -k2 -nr | \
            head -n \${alpha_res} | tail -n1 | cut -f2`
        int_alpha_d=\${alpha_d%%.*}
    fi

    qiime tools export --input-path dada2_stats.qzv \
        --output-path ./dada2_stats
    mv dada2_stats/metadata.tsv dada2_qc.tsv
    """
}

process dada2_rarefaction {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_freq
    path metadata
    val(alpha_d)

    output:
    path "alpha-rarefaction-curves.qzv", emit: alpha_rarefaction_curves

    script:
    """
    qiime diversity alpha-rarefaction --i-table $asv_freq \
        --m-metadata-file $metadata \
        --o-visualization alpha-rarefaction-curves.qzv \
        --p-min-depth 10 --p-max-depth $alpha_d
    """
}