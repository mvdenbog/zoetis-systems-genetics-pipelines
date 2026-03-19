process METHYLATION_CONTROLS_QC {

input:
 tuple val(meta), path (mergecg_bed)

output:
    path ("ecoli_k12_mg1655/${meta.id}_ecoli_qc.bed"), emit: ecoli_bed

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir ecoli_k12_mg1655
    zcat ${mergecg_bed} | grep '^U00096.3'  > ecoli_k12_mg1655/${meta.id}_ecoli_qc.bed 
    """
}

process METHYLATION_CONTROLS_FIGURE {
    container "docker.io/mvdenbog/rockylinux93_r" 

input:
path (ecoli_beds)
output:
path "*.pdf"

script:
"""
test.R ${ecoli_beds}
"""
}
