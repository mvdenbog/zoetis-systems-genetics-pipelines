cp /sc/kzd/proj/SysGen/MetagenomeReseach/Pipelines/ref/GCF_002263795.3_ARS-UCD2.0_genomic.fna .
 mv GCF_002263795.3_ARS-UCD2.0_genomic.fna GCF_002263795.3_ARS-UCD2.0_genomic.fa

module purge
module load  sge/2011.11  nextflow/21.10 singularity/3.8.1

export NXF_VER=23.04.0
export NXF_HOME=$PWD/.nextflow

# maybe necessary to:
export PATH=~/jdk-11.0.22/bin/:$PATH

nextflow run main.nf -c nextflow.config  --input ../../BioLizard/zoetis-wgbs-pipeline/BosTaurus_WGBS_samplesheet.csv --outdir BosTaurus_WGBS_Output --aligner bismark --fasta $PWD/GCF_002263795.3_ARS-UCD2.0_genomic.fa   -profile singularity,cluster

nextflow run main.nf -c nextflow.config  --input $PWD/zoetisdata_samplesheet.csv --outdir BosTaurus_WGBS_Output --aligner bismark --fasta /sc/kzd/proj/SysGen/MetagenomeReseach/Pipelines/IndexHost/BosTaurus/Bismark/GCF_002263795.3_ARS-UCD2.0_genomic.fa --fasta_index /sc/kzd/proj/SysGen/MetagenomeReseach/Pipelines/IndexHost/BosTaurus/Faidx/GCF_002263795.3_ARS-UCD2.0_genomic.fna.fai --bismark_index /sc/kzd/proj/SysGen/MetagenomeReseach/Pipelines/IndexHost/BosTaurus/Bismark/   -profile singularity,cluster

nextflow run main.nf -c nextflow.config  --input $PWD/zoetisdata_samplesheet.csv --outdir BosTaurus_WGBS_Output --aligner bismark --fasta /sc/kzd/proj/SysGen/MetagenomeReseach/Pipelines/IndexHost/BosTaurus/Bismark/GCF_002263795.3_ARS-UCD2.0_genomic.fa --fasta_index /sc/kzd/proj/SysGen/MetagenomeReseach/Pipelines/IndexHost/BosTaurus/Faidx/GCF_002263795.3_ARS-UCD2.0_genomic.fna.fai --bismark_index /sc/kzd/proj/SysGen/MetagenomeReseach/Pipelines/IndexHost/BosTaurus/Bismark/   -profile singularity,cluster -resume

