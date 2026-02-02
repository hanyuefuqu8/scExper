#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process qualityControl {
    tag "qualityControl"
    container 'hanyuefuqu8/scexpr-qc:v1.0'
    
    input:
    path dataset_dir
    path sample_info
    val sct
    path mito_genes

    output:
    path "QC_results", emit: QC_results

    script:
    """
    mkdir -p QC_results
    Rscript "/scripts/qc.R" ${sample_info} ${dataset_dir} ${sct} ${mito_genes}
    Rscript "/scripts/soupx.R" -c ${sample_info} -d ${dataset_dir}
    """
}

process CHOIR_analysis {
    tag "CHOIR_analysis"
    container 'hanyuefuqu8/scexpr-choir:v1.0'
    
    input:
    path QC_results
    val sct

    output:
    path "CHOIR_results", emit: CHOIR_results

    script:
    """
    if [ "${sct}" -eq 1 ]; then
        Rscript "/scripts/choir.r" --input_rds ${QC_results}/*.rds -s -o CHOIR_results -n 8
    else
        Rscript "/scripts/choir.r" --input_rds ${QC_results}/*.rds -o CHOIR_results -n 8
    fi
    """
}

process Metaneighbor_analysis {
    tag "Metaneighbor_analysis"
    container 'seurat-v5-8:latest'
    
    input:
    path CHOIR_results
    val sct
    val cluster_threshold

    output:
    path "Metaneighbor_results_with_metaneighbor.rds", emit: Metaneighbor_results

    script:
    """
    if [ "${sct}" -eq 1 ]; then
        Rscript "/scripts/metaneighbor.r" -i ${CHOIR_results}/project_after*.rds -o Metaneighbor_results -b "sample_type" -c "CHOIR_clusters_0.05" -t ${cluster_threshold} -s
    else
        Rscript "/scripts/metaneighbor.r" -i ${CHOIR_results}/project_after*.rds -o Metaneighbor_results -b "sample_type" -c "CHOIR_clusters_0.05" -t ${cluster_threshold}
    fi  
    """
}

process Integrate_analysis{
    tag "Integrate_analysis"
    container 'hanyuefuqu8/scexpr-integrate:v1.0'

    publishDir "${params.outdir}/Integrate", mode: 'copy'


    input:
    path Metaneighbor_results
    val sct

    output:
    path "Integrate_results/unintegrated.h5ad", emit:unintegrated_h5ad
    path "Integrate_results/best*.h5ad", emit:best_h5ad

    script:
    """
    mkdir "Integrate_results"
    if [ "${sct}" -eq 1 ]; then
        Rscript "/scripts/harmony.r" -i ${Metaneighbor_results} -o "harmony_integrated.rds" -b "sample_type" -c "metaneighbor" -s
    #    Rscript "/scripts/scvi.r" -i ${Metaneighbor_results} -o "scvi_integrated.rds" -b "sample_type" -c "metaneighbor" -s
    else
        Rscript "/scripts/harmony.r" -i ${Metaneighbor_results} -o "harmony_integrated.rds" -b "sample_type" -c "metaneighbor"
    #    Rscript "/scripts/scvi.r" -i ${Metaneighbor_results} -o "scvi_integrated.rds" -b "sample_type" -c "metaneighbor"
    fi 
    
    Rscript "/scripts/rliger.r" -i ${Metaneighbor_results} -o "rliger_integrated.rds" -b "sample_type" -c "metaneighbor"

    find "Integrate_results" -name "*.rds" -print0 | while IFS= read -r -d '' file; do
        temp=\${file##*/}
        name=\${temp%.rds}
        Rscript "/scripts/Seurat2h5ad.R" -i \$file -o Integrate_results/integrated-\${name}.h5ad
    done
    
    Rscript "/scripts/Seurat2h5ad.R" -i ${Metaneighbor_results} -o Integrate_results/unintegrated.h5ad

    python /scripts/scib.py --unintegrated_h5ad  Integrate_results/unintegrated.h5ad\
                            --integrated_file Integrate_results/integrated-*.h5ad \
                            --batch_key sample_type \
                            --label_key metaneighbor \
                            --n_jobs 8
    """
}


process Cellrank{
    tag "Cellrank"
    container 'hanyuefuqu8/scexpr-integrate:v1.0'

    publishDir "${params.outdir}/Cellrank", mode: 'copy'

    input:
    path best_h5ad
    val specified_cluster

    output:
    path "Cellrank_results", emit: Cellrank_results
    script:
    """
    mkdir Cellrank_results
    cd Cellrank_results
    python /scripts/cellrank_realtime.py --dir ../${best_h5ad} \
                                         --celltype metaneighbor \
                                         --specified_key ${specified_cluster} \
    """
}


process Memento{
    tag "Memento"
    container 'hanyuefuqu8/scexpr-integrate:v1.0'
    publishDir "${params.outdir}/Memento", mode: 'copy'
    
    input:
    path unintegrated_h5ad
    val control_value

    output:
    path "Memento_results", emit: Memento_results
    path "Cytotrace_results", emit: ct_dir

    script:
    """
    mkdir Memento_results
    python /scripts/memento_analysis.py --input_h5ad ${unintegrated_h5ad} \
                              --group_key info \
                              --control_value ${control_value} \
                              --celltype_key metaneighbor \
                              --replication_key bio_rep
    mkdir Cytotrace_results
    python /scripts/cytotrace.py --input_h5ad ${unintegrated_h5ad}
    """
}

process genes2genes{
    tag "genes2genes"
    container 'hanyuefuqu8/scexpr-g2g:v1.0'
    
    publishDir "${params.outdir}/genes2genes", mode: 'copy'


    input:
    path ct_dir
    val control_value
    val treat_value 

    output:
    path "genes2genes_results", emit: genes2genes_results

    script:
    """
    mkdir genes2genes_results
    python /scripts/genes2genes_analysis.py --input_dir ${ct_dir} \
                                  --control_label ${control_value} \
                                  --treat_label ${treat_value} \
                                  --n_bins 10
    """
}



workflow {
    // 执行主要流程
    qualityControl(file(params.dataset_dir), file(params.sample_info), params.sct, file(params.mito_genes))
    CHOIR_analysis(qualityControl.out.QC_results, params.sct)
    Metaneighbor_analysis(CHOIR_analysis.out.CHOIR_results, params.sct, params.cluster_threshold)
    Integrate_analysis(Metaneighbor_analysis.out.Metaneighbor_results, params.sct)
    
    // 条件分支处理
    if (params.type == "time-series") {
        Cellrank(Integrate_analysis.out.best_h5ad, params.specified_cluster)
    } else if (params.type == "controlled-experiment") {
        Memento(Integrate_analysis.out.unintegrated_h5ad, params.control_value)
        genes2genes(Memento.out.ct_dir, params.control_value, params.treat_value)
    }
}  

