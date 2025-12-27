#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process qualityControl {
    tag "QualityControl"
    container 'seurat-v5-7:latest'
    
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
    container 'choir-docker-4:latest'
    
    input:
    path QC_results
    val sct

    output:
    path "CHOIR_results", emit: CHOIR_results

    script:
    """
    if [ "${sct}" -eq 1 ]; then
        Rscript "/scripts/choir.r" --input_rds ${QC_results}/*.rds -s -o CHOIR_results -n 16
    else
        Rscript "/scripts/choir.r" --input_rds ${QC_results}/*.rds -o CHOIR_results -n 16
    fi
    """
}

process Metaneighbor_analysis {
    tag "Metaneighbor_analysis"
    container 'seurat-v5-8:latest'
    
    input:
    path CHOIR_results
    val sct

    output:
    path "Metaneighbor_results", emit: Metaneighbor_results

    script:
    """
    if [ "${sct}" -eq 1 ]; then
        Rscript "/scripts/metaneighbor.r" -i ${CHOIR_results}/project_after*.rds -o Metaneighbor_results -b "sample_type" -c "CHOIR_clusters_0.05" -s
    else
        Rscript "/scripts/metaneighbor.r" -i ${CHOIR_results}/project_after*.rds -o Metaneighbor_results -b "sample_type" -c "CHOIR_clusters_0.05"
    fi  
    """
}

process Integrate_analysis{
    tag "Integrate_analysis"
    container 'choir-docker-14:latest'

    input:
    path Metaneighbor_results
    val sct

    output:
    path "Integrate_results/unintegrated.h5ad", emit:unintegrated_h5ad
    path "Integrate_results/best.h5ad", emit:best_h5ad

    script:
    """
    #!/bin/bash
    set -e
    mkdir Integrate_results
    if [ "${sct}" -eq 1 ]; then
    Rscript "/scripts/harmony.r" -i ${Metaneighbor_results} -o "harmony_integrated.rds" -b "sample_type" -c "metaneighbor" -s
    Rscript "/scripts/scvi.r" -i ${Metaneighbor_results} -o "scvi_integrated.rds" -b "sample_type" -c "metaneighbor" -s
    else
    Rscript "/scripts/harmony.r" -i ${Metaneighbor_results} -o "harmony_integrated.rds" -b "sample_type" -c "metaneighbor"
    Rscript "/scripts/scvi.r" -i ${Metaneighbor_results} -o "scvi_integrated.rds" -b "sample_type" -c "metaneighbor"
    fi  
    Rscript "/scripts/rliger.r" -i ${Metaneighbor_results} -o "rliger_integrated.rds" -b "sample_type" -c "metaneighbor"

    find Integrate_results -name "*.rds" -print0 | while IFS= read -r -d '' file; do
        temp=\${file##*/}
        name=\${temp%.rds}
        Rscript "/scripts/Seurat2h5ad.R" -i \$file -o Integrate_results/integrated-\${name}.h5ad
    done
    
    Rscript "/scripts/Seurat2h5ad.R" -i ${Metaneighbor_results} -o Integrate_results/unintegrated.h5ad

    python /scripts/scib.py --unintegrated_h5ad  Integrate_results/unintegrated.h5ad\
                            --integrated_file Integrate_results/integrated-*.h5ad \
                            --batch_key sample_type \
                            --label_key metaneighbor \
                            --njobs 16
    """
}


process Cellrank{
    tag "Cellrank_analysis"
    container 'choir-docker-13:latest'

    input:
    path best_h5ad

    output:
    path "Cellrank_results", emit: Cellrank_results
    script:
    """
    python /scripts/cellrank_realtime.py --dir ${best_h5ad} \
                                         --celltype metaneighbor
    """
}

process Memento{
    tag "Memento_analysis"
    container 'choir-docker-14:latest'
    
    input:
    file unintegrated_h5ad

    output:
    path "Memento_results", emit: Memento_results

    script:
    """
    python /scripts/memento.py --input_h5ad ${unintegrated_h5ad} \
                              --group_key sample_type \
                              --control_value ${control_value} \
                              --sample_key metaneighbor \
    """
}

process genes2genes{
    tag "genes2genes_analysis"
    container 'choir-docker-14:latest'
    
    input:
    file unintegrated_h5ad

    output:
    path "genes2genes_results", emit: genes2genes_results

    script:
    """
    python /scripts/cytotrace.py --input_h5ad ${unintegrated_h5ad}
    python /scripts/genes2genes.py --input_h5ad ${unintegrated_h5ad} \
                                  --group_key sample_type \
                                  --sample_key metaneighbor \
    """
}

workflow {
    qualityControl(file(params.dataset_dir), file(params.sample_info), params.sct, file(params.mito_genes))
    CHOIR_analysis(qualityControl.out.QC_results, params.sct)
    Metaneighbor_analysis(CHOIR_analysis.out.CHOIR_results, params.sct)
    Integrate_analysis(CHOIR_analysis.out.CHOIR_results, params.sct)
    if (params.type == "time-series") {
        Cellrank(Integrate_analysis.out.best_h5ad)
        emit: Cellrank_result for Cellrank.out.Cellrank_results
    } else if (params.type == "treatment") {
        Memento(Integrate_analysis.out.unintegrated_h5ad)
        genes2genes(Integrate_analysis.out.unintegrated_h5ad)
        emit: Memento_result for Memento.out.Memento_results
        emit: genes2genes_result for genes2genes.out.genes2genes_results
    }
}
