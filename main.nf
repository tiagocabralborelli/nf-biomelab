process QUALITY_CONTROL {
    publishDir 'results/fastp', mode:'symlink'
    conda '${params.envs.index}'
    input:
        tuple val(sample_name), path(reads)
    output:
        tuple val(sample_name), path('*trim_{1,2}.fastq.gz')
    script:
    """
    fastp --in1 ${reads[0]}\
    --in2 ${reads[1]}\
    --out1 ${sample_name}_trim_1.fastq.gz\
    --out2 ${sample_name}_trim_2.fastq.gz\
    --trim_front1 20\
    --trim_tail1 20\
    --trim_front2 20\
    --trim_tail2 20\
    --qualified_quality_phred 20\
    --length_required 50\
    -h ${sample_name}.html\
    -j ${sample_name}.json
    """
}
process DECONTAMINATION {
    publishDir 'results/bowtie2', mode:'symlink'
    conda '${params.envs.training}'

    input:
        tuple val(sample_name), path(reads)
    output:
        tuple val(sample_name), path('*bowtie_{1,2}.fastq.gz')
    script:
    """
    bowtie2 -x ${params.index.bowtie} \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -S ${sample_name}_bowtie.sam \
    --threads 12 

    samtools view -bS -f 13 ${sample_name}_bowtie.sam > ${sample_name}_bowtie.bam
    samtools fastq -1 ${sample_name}_bowtie_1.fastq.gz -2 ${sample_name}_bowtie_2.fastq.gz ${sample_name}_bowtie.bam
    rm ${sample_name}_bowtie.bam
    rm ${sample_name}_bowtie.sam
    """
}
process ASSEMBLY {
    publishDir 'results/spades', mode:'symlink'
    conda '${params.envs.training}'
    input:        
        tuple val(sample_name), path(reads)
    output:
    tuple val(sample_name), path("${sample_name}_contigs.fasta")
    script:
    """
    spades.py \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o ${sample_name} \
    --threads 8 \
    --memory 64
    cp ./${sample_name}/contigs.fasta ./${sample_name}_contigs.fasta
    """
}
process ORF_FINDER {
    publishDir 'results/prodigal', mode:'symlink'
    conda '${params.envs.training}'
    input:        
        tuple val(sample_name), path(contigs)
    output:
        tuple val(sample_name), path("${sample_name}.fna"), path("${sample_name}.faa")
    script:
    """
    prodigal -q -i ${contigs} -p single -o /dev/null -d ${sample_name}.fna -a ${sample_name}.faa -f gff
    """        
}

include {ANNOTAION as ANNOT_ARGS} from './modules/hmmsearch.nf'
include {ANNOTAION as ANNOT_CAZ} from './modules/hmmsearch.nf'
include {ANNOTAION as ANNOT_PLT} from './modules/hmmsearch.nf'

workflow {
    reads_ch = Channel.fromFilePairs("${params.input.reads}/*_{1,2}.fastq")
    contigs = ASSEMBLY(reads_ch)
    orfs = ORF_FINDER(contigs)
    ANNOT_ARGS(orfs,params.index.hammer_arg,"args")
    ANNOT_PLT(orfs,params.index.hammer_plastic,"plastic")
}