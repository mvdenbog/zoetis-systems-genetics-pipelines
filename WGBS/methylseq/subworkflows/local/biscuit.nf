 include { BISCUIT_ALIGN } from '../../modules/nf-core/biscuit/align/main' 
 include { BISCUIT_BLASTER } from '../../modules/nf-core/biscuit/biscuitblaster/main'
 include { BISCUIT_BSCONV } from '../../modules/nf-core/biscuit/bsconv/main'
 include { BISCUIT_EPIREAD } from '../../modules/nf-core/biscuit/epiread/main'
 include { BISCUIT_INDEX } from '../../modules/nf-core/biscuit/index/main'
 include { BISCUIT_MERGECG } from '../../modules/nf-core/biscuit/mergecg/main'
 include { BISCUIT_PILEUP } from '../../modules/nf-core/biscuit/pileup/main'
 include { BISCUIT_QC } from '../../modules/nf-core/biscuit/qc/main'
 include { BISCUIT_VCF2BED } from '../../modules/nf-core/biscuit/vcf2bed/main'
 include { METHYLATION_CONTROLS_QC } from '../../modules/local/biscuit/methylation_controls_qc/main'
include { METHYLATION_CONTROLS_FIGURE } from '../../modules/local/biscuit/methylation_controls_qc/main'




include { SAMTOOLS_SORT as SAMTOOLS_SORT_DEDUPLICATED } from '../../modules/nf-core/samtools/sort/main'

process pileup {
container = "https://depot.galaxyproject.org/singularity/mulled-v2-13ecc02e46375c1530003669d682d88af478b062:8f5ef05bea4d307d50ffcf518735290a84db36a5-0"
// container = "https://depot.galaxyproject.org/singularity/biscuit:1.0.2.20220113--he272189_1"
        input:
        tuple val(id), path(bam), path(bai)
    path(ref)

        output:
        tuple val(id), path("*_pileup.vcf.gz")

        script:
    ref = params.fasta
        """
                biscuit pileup -@ 4 ${ref} ${bam} | bgzip -f > ${id}_pileup.vcf.gz
        """
}

process multiqc {
container = "https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0"
    input:
        path(files)
        val(tag)

    output:
        path("*.html")

    script:
    """
        multiqc . --title "${tag} Report" --filename "${tag}_multiqc.html"
    """
}

process index {
container = "https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0"
    input:
        tuple val(id), path(bam)
        val(tag)

    output:
        tuple val(id), path('*.bam', includeInputs: true), path("*.bai")

    script:
    """
        samtools index -@ 4 ${bam}
    """
}

process stats {
container = "https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0"
    input:
        tuple val(id), path(bam), path(index)
        path(ref)

    output:
        path("*.txt")

    script:
    ref = params.fasta
    """
        samtools stats -@ 2 -r ${ref} ${bam} > ${id}_stats.txt
    """
}

process filter_bam {
container = "https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0"
    input:
        tuple val(id), path(bam), path(bai)
        path(bed_file)

    output:
        tuple val(id), path("*.bam"), path("*.bai")

    script:
    """
        samtools view --write-index -hb -L ${bed_file} ${bam} -o ${id}_filtered.bam##idx##${id}_filtered.bam.bai
    """
}

process tabix_vcf {
 container = "https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0"
 
    input:
        tuple val(id), path(vcf)

    output:
        tuple val(id), path("*.vcf.gz", includeInputs: true), path("*.tbi")

    script:
    """
        tabix -p vcf ${vcf}
    """
}

process tabix_bed {
 container = "https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0"
    input:
        tuple val(id), path(bed)

    output:
        tuple val(id), path("*.bed.gz", includeInputs: true), path("*.tbi")

    script:
    """
        tabix -p bed ${bed}
    """
}

process biscuiteer {
container = "https://depot.galaxyproject.org/singularity/bioconductor-biscuiteer:1.8.0--r41hdfd78af_0"
    input:
        tuple val(id), path(bed), path(bed_tbi), path (vcf), path(vcf_tbi)

    output:
        tuple val(id), path("*.rds")

    script:
    """
        #!/usr/bin/env Rscript --vanilla
        library(biscuiteer)

        bisc <- readBiscuit(BEDfile = "${bed}", VCFfile = "${vcf}", merged = T, genome = "${params.genome}", sparse = T)
        saveRDS(bisc, file = "${id}.rds", compress = F)
    """
}

process bamtobed {
memory = 200.GB // 550.GB
// queue = 'all.q' // 'ram.q'
container = "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3"
    input:
        tuple val(id), path(bam), path(bai)

    output:
        tuple val(id), path("*.bed")

    script:
    """
        bedtools bamtobed -i ${bam} > ${id}.tmp.bed
        sortBed -i ${id}.tmp.bed > ${id}.tmp.sorted.bed
        cut -f1-3 ${id}.tmp.sorted.bed > ${id}.bed
        rm -f ${id}.tmp.bed ${id}.tmp.sorted.bed
    """
}

process preseq {
    label 'error_ignore'
container = "https://depot.galaxyproject.org/singularity/preseq:3.1.2--h2c25361_3"
    input:
        tuple val(id), path(bed)

    output:
        path("*.tsv")

    script:
    //defects = ["", "-defects", "-defects -Q"][task.attempt - 1]
     defects = "-defects -Q"
    """
        preseq gc_extrap \
            -o preseq_gc-extrap_${id}.tsv \
            -e 3E9 \
            -s 10E6 \
            -bed \
            ${defects} ${bed}
    """
}

process biscuit_qc {
memory = 200.GB
// queue = 'all.q'
container = "https://depot.galaxyproject.org/singularity/mulled-v2-13ecc02e46375c1530003669d682d88af478b062:8f5ef05bea4d307d50ffcf518735290a84db36a5-0"
// container = "https://depot.galaxyproject.org/singularity/biscuit:1.0.2.20220113--he272189_1"
        input:
        tuple val(id), path(bam), path(bai), path(vcf)
    path(ref)
    path(assets)

        output:
        path("*.txt")

        script:
    ref = params.fasta
    assets = params.bqc_assets
        """
        ${params.bqc_path} -o . -v ${vcf} ${assets} ${ref} ${id} ${bam}
        """
}


workflow BISCUIT {
    take:
    reads              // channel: [ val(meta), [ reads ] ]
    biscuit_index      // channel: /path/to/BiscuitIndex/
    skip_deduplication // boolean: whether to deduplicate alignments

    main:
    versions = Channel.empty()
    assets = Channel.fromPath(params.bqc_assets).first()
    target_bed = Channel.fromPath(params.target_bed)
    skip_deduplication = false


    BISCUIT_ALIGN( reads, biscuit_index )
    versions = versions.mix(BISCUIT_ALIGN.out.versions)

bams_with_index = BISCUIT_ALIGN.out.bam.join(BISCUIT_ALIGN.out.bai, by: 0)
//bams_with_index.view()

bams_with_index
    .map { item, bam, bai -> 
        def meta = [:] 
        meta.id = item.id
        return [item.id,bam,bai]
      }
    .set{bams_with_index}
// bams_with_index.view()


reference = params.fasta
pileup(bams_with_index, reference)
 if (!params.skip_biscuit_qc) {
        biscuit_qc(bams_with_index.join(pileup.out, by: 0), reference, assets)
        bamtobed(bams_with_index) | preseq
        stats(bams_with_index, reference)

        biscuit_qc.out
                .mix(preseq.out)
                .mix(stats.out)
                .collect()
                .set{final_multiqc_files}
} else {
        bamtobed(bams_with_index) | preseq
        stats(bams_with_index, reference)
                
        preseq.out
                .mix(stats.out)
                .collect()
                .set{final_multiqc_files}


}
        multiqc(final_multiqc_files, "post-align")

input = BISCUIT_ALIGN.out.bam.join(BISCUIT_ALIGN.out.bai, by: 0)
input
    .map { item, bam, bai -> 
        def meta = [:] 
        meta.id = item.id
        return [meta,bam,bai,[],[]]
      }
    .set{input}
// input.view()




//    BISCUIT_PILEUP(BISCUIT_ALIGN.out.bam.join(BISCUIT_ALIGN.out.bai, by: 0).join(Channel.empty(), by: 0).join(Channel.empty(), by: 0), biscuit_index)
BISCUIT_PILEUP(input, biscuit_index)

// Works, but not used further :
// BISCUIT_BLASTER(reads, biscuit_index)

bams_with_index = BISCUIT_ALIGN.out.bam.join(BISCUIT_ALIGN.out.bai, by: 0)
bams_with_index
    .map { item, bam, bai -> 
        def meta = [:] 
        meta.id = "${item.id}.bsconv"
        return [meta,bam,bai]
      }
    .set{bams_with_index}
// bams_with_index.view()

BISCUIT_BSCONV(bams_with_index, biscuit_index)

BISCUIT_VCF2BED(BISCUIT_PILEUP.out.vcf)

input = BISCUIT_ALIGN.out.bam.join(BISCUIT_ALIGN.out.bai, by: 0)
// input.view()
input
    .map { item, bam, bai -> 
        def meta = [:] 
        meta.id = item.id   // "id:${item.id}"
        return [meta,bam,bai]
      }
    .set{input}
// input.view()

input_bed = BISCUIT_VCF2BED.out.bed
// input_bed.view()
input_bed
 .map { item, bed ->
        def meta = [:] 
        meta.id = item.id
        meta.id -= ~/.bsconv/
        return [meta,bed]
      }
    .set{input_bed}

// input.view()
input = input.join(input_bed, by: 0)
// input.view()
input
    .map { item, bam, bai, bed -> 
        def meta = [:] 
        meta.id = item.id
        meta.id -= ~/.bsconv/
        return [meta,bam,bai,bed]
      }
    .set{input}
// input.view()


// for those 2, works, but output is empty ... see if works with real data?
//BISCUIT_EPIREAD(input, biscuit_index)
BISCUIT_MERGECG(BISCUIT_VCF2BED.out.bed, biscuit_index)

if (!params.skip_methylation_controls) {
all_beds = METHYLATION_CONTROLS_QC(BISCUIT_MERGECG.out.mergecg_bed).collect()

METHYLATION_CONTROLS_FIGURE(all_beds)
}

//works:
    BISCUIT_QC(BISCUIT_ALIGN.out.bam, biscuit_index)

 // BISCUIT_QC.out.biscuit_qc_reports.view()


//    filter_bam(BISCUIT_ALIGN.out, params.target_bed)
//    pileup(filter_bam.out, biscuit_index) | tabix_vcf | vcf2bed | tabix_b1

//    BISCUIT_MERGECG(tabix_b1.out, biscuit_index) | tabix_b2
//    biscuiteer(tabix_b2.out.join(tabix_vcf.out))

    emit:
    bam        = BISCUIT_ALIGN.out.bam        // channel: [ val(meta), [ bam ] ] ## sorted, non-deduplicated (raw) BAM from aligner
    mqc        = multiqc.out                        // path: *{html,txt}
    dedup      = BISCUIT_ALIGN.out.bam // SAMTOOLS_SORT_DEDUPLICATED.out.bam 
    versions                                       // path: *.version.txt


}

