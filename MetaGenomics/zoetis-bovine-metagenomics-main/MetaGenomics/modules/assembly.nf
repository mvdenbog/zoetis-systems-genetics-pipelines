process SPADES{
errorStrategy 'ignore'
publishDir "${params.outdir}/${id}/Assembly/", mode: 'copy'
label "SPADES"
memory 60.GB
disk 80.GB
input:
tuple val(id), path(FQs)
output:
tuple val(id), path("${id}.*.fa"), emit: spades_fa_result
path "v_spades.txt", emit: v_spades
script:
"""

for i in Virus Bacteria Fungi Other unclassifieds ; do
mkdir ${id}_spades_output_\${i};
if ([ -s ${id}.\${i}.fq ]); then
#if ([ \$(cat ${id}.\${i}.fq | wc -l) -gt 8000 ]); then
#spades.py -s ${id}.\${i}.fq -o ${id}_spades_output_\${i};
#cp ${id}_spades_output_\${i}/contigs.fasta ${id}.\${i}.fa
#else
#seqkit fq2fa ${id}.\${i}.fq -o ${id}.\${i}.fa
#fi
 if ([ \$(cat ${id}.\${i}.fq | wc -l) -gt 8000 ]); then
   fqdeinterleave.pl ${id}.\${i}.fq
   if ([ -s ${id}.\${i}_sgl.fq ]); then 
      metaspades.py  -1 ${id}.\${i}_R1.fq -2 ${id}.\${i}_R2.fq  -s ${id}.\${i}_sgl.fq  -o ${id}_spades_output_\${i}
   else 
      metaspades.py -1 ${id}.\${i}_R1.fq -2 ${id}.\${i}_R2.fq  -o ${id}_spades_output_\${i}
   fi; 
   cp ${id}_spades_output_\${i}/contigs.fasta ${id}.\${i}.fa
 else
   seqkit fq2fa ${id}.\${i}.fq -o ${id}.\${i}.fa
 fi
fi
done

spades.py --version &> v_spades.txt
"""
}

process SPADES_KAIJU{
publishDir "${params.outdir}/${id}/Assembly/Classification", mode: 'copy'
label "SPADES"
input:
tuple val(id), path(FAs)
tuple path(nodes), path(fmi), path(names)
output:
path("${id}_kaiju_krona.html")
path "v_kaiju.txt", emit: v_kaiju
script:
"""
cat *.fa > ${id}.FA
kaiju -z ${task.cpus} -t ${nodes} -f ${fmi} -i ${id}.FA  -o ${id}_kaiju.out -v -a greedy -m 1 -e 5
kaiju-addTaxonNames -i ${id}_kaiju.out -o ${id}_kaiju_labels.out -t ${nodes} -n ${names} -u -p
kaiju2krona -i ${id}_kaiju.out -t ${nodes} -n ${names} -o ${id}.krona
ktImportText -o  ${id}_kaiju_krona.html ${id}.krona

echo \$(kaiju -h 2>&1) > v_kaiju.txt
"""
}

