process METASPADES {
  tag "${meta.id}"
  publishDir "${params.outdir}/02_assembly/02_1_primary_assembly", mode: 'copy', pattern: "metaspades_${meta.id}/*"
  label 'ASSEMBLY_SR'

  input:
  tuple val(meta), path(read1), path(read2)

  output:
  tuple val(meta), path("metaspades_${meta.id}/${meta.id}.contigs.fa"), emit: assembly
  tuple val(meta.id), path("metaspades_${meta.id}/${meta.id}.log"), path("metaspades_${meta.id}/${meta.id}.params.txt"), emit: report
  path "v_spades.txt", emit: v_metaspades

  script:
    (_,mem,unit) = (task.memory =~ /(\d+) ([A-Z]B)/)[0]
    if ( unit =~ /GB/ ) { 
      """
      echo "
[
    {
        orientation: \\"fr\\",
        type: \\"paired-end\\",
        right reads: [
            \\"${read1.join('\\",\n            \\"')}\\"
        ],
        left reads: [
            \\"${read2.join('\\",\n            \\"')}\\"
        ]
    }
]" > input.yaml

      metaspades.py -t ${task.cpus} -m $mem --dataset input.yaml -o metaspades_${meta.id}
      mv metaspades_${meta.id}/scaffolds.fasta metaspades_${meta.id}/${meta.id}.contigs.fa
      mv metaspades_${meta.id}/spades.log metaspades_${meta.id}/${meta.id}.log
      mv metaspades_${meta.id}/params.txt metaspades_${meta.id}/${meta.id}.params.txt 

      spades.py --version &> v_spades.txt
      """
    } else { 
      error "Memory setting for the ASSEMBLY process is in $unit, it must be in GB (check config files) "
    }
}


process MEGAHIT {
  tag "${meta.id}"
  publishDir "${params.outdir}/02_assembly/02_1_primary_assembly", mode: 'copy', pattern: "megahit_${meta.id}/*"
  label 'ASSEMBLY_SR'

  input:
  tuple val(meta), path(read1), path(read2)

  output:
  tuple val(meta), path("megahit_${meta.id}/${meta.id}.contigs.fa"), emit: assembly
  tuple val(meta.id), path("megahit_${meta.id}/${meta.id}.log"), path("megahit_${meta.id}/${meta.id}.params.txt"), emit: report
  path "v_megahit.txt", emit: v_megahit

  script:
    """
    megahit -t ${task.cpus} -1 ${read1.join(',')} -2 ${read2.join(',')} -o megahit_${meta.id} --out-prefix "${meta.id}"
    mv megahit_${meta.id}/options.json megahit_${meta.id}/${meta.id}.params.txt
    rm -r  megahit_${meta.id}/intermediate_contigs

    megahit --version &> v_megahit.txt
    """
}


process HIFIASM_META {
  tag "${meta.id}"
  publishDir "${params.outdir}/02_assembly/02_1_primary_assembly", mode: 'copy', pattern: "hifiasm_meta_${meta.id}/${meta.id}.*"
  label 'ASSEMBLY_HIFI'

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("hifiasm_meta_${meta.id}/${meta.id}.contigs.fa"), emit: assembly
  tuple val(meta.id), path("hifiasm_meta_${meta.id}/${meta.id}.log"), path("hifiasm_meta_${meta.id}/${meta.id}.params.txt"), emit: report
  path "v_hifiasm_meta.txt", emit: v_hifiasm_meta

  script:
  """
  mkdir hifiasm_meta_${meta.id}

  hifiasm_meta -t ${task.cpus} -o ${meta.id} ${reads.join(' ')} 2> hifiasm_meta_${meta.id}/${meta.id}.log

  # gfa to fasta format
  awk '/^S/{print ">"\$2"\\n"\$3}' ${meta.id}.p_ctg.gfa | fold > hifiasm_meta_${meta.id}/${meta.id}.contigs.fa

  mv ${meta.id}.cmd  hifiasm_meta_${meta.id}/${meta.id}.params.txt

  echo \$(hifiasm_meta --version 2>&1) > v_hifiasm_meta.txt
  """
}


process METAFLYE {
  tag "${meta.id}"
  publishDir "${params.outdir}/02_assembly/02_1_primary_assembly", mode: 'copy', pattern: "metaflye_${meta.id}/${meta.id}.*"
  label 'ASSEMBLY_HIFI'

  input:
  tuple val(meta), path(reads)

  output:
  tuple val(meta), path("metaflye_${meta.id}/${meta.id}.contigs.fa"), emit: assembly
  tuple val(meta.id), path("metaflye_${meta.id}/${meta.id}.log"), path("metaflye_${meta.id}/${meta.id}.params.json"), emit: report
  tuple val(meta), path("metaflye_${meta.id}/assembly_info.txt"), emit: infos
  path "v_metaflye.txt", emit: v_metaflye

  script:
  """
  mkdir metaflye_${meta.id}
  
  flye --pacbio-hifi ${reads.join(' ')} -o 'metaflye_${meta.id}' --meta -t ${task.cpus} 2> metaflye_${meta.id}/${meta.id}.log

  mv metaflye_${meta.id}/assembly.fasta  metaflye_${meta.id}/${meta.id}.contigs.fa
  mv metaflye_${meta.id}/params.json  metaflye_${meta.id}/${meta.id}.params.json

  flye --version &> v_metaflye.txt
  """
}