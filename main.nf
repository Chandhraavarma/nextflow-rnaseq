nextflow.enable.dsl=2

params.input  = null
params.fasta  = null
params.gtf    = null
params.outdir = 'results'
params.skip_deseq2 = true

process FASTQC_RAW {
    tag "${meta.id}"
    label 'rnaseq'
    publishDir "${params.outdir}/qc/raw_fastqc", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path('*.html'), path('*.zip'), emit: reports
    script:
    def read_args = reads.collect { it.toString() }.join(' ')
    """
    fastqc --threads ${task.cpus} ${read_args}
    """
}

process FASTP {
    tag "${meta.id}"
    label 'rnaseq'
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.id}_R1.trim.fastq.gz"), path("${meta.id}_R2.trim.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.html"), path("${meta.id}.fastp.json"), emit: reports
    script:
    """
    fastp \
      --thread ${task.cpus} \
      --in1 ${reads[0]} --in2 ${reads[1]} \
      --out1 ${meta.id}_R1.trim.fastq.gz \
      --out2 ${meta.id}_R2.trim.fastq.gz \
      --html ${meta.id}.fastp.html \
      --json ${meta.id}.fastp.json
    """
}

process STAR_INDEX {
    tag "STAR index"
    label 'star_index'
    publishDir "${params.outdir}/star_index", mode: 'copy'
    input:
    path fasta
    path gtf
    output:
    path 'star_index', emit: index
    script:
    """
    mkdir star_index
    STAR --runThreadN ${task.cpus} \
      --runMode genomeGenerate \
      --genomeDir star_index \
      --genomeFastaFiles ${fasta} \
      --sjdbGTFfile ${gtf}
    """
}

process STAR_ALIGN {
    tag "${meta.id}"
    label 'star_align'
    publishDir "${params.outdir}/alignment", mode: 'copy'
    input:
    tuple val(meta), path(r1), path(r2)
    path index
    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}.Log.final.out"), emit: log
    script:
    """
    STAR --runThreadN ${task.cpus} \
      --genomeDir ${index} \
      --readFilesIn ${r1} ${r2} \
      --readFilesCommand zcat \
      --outFileNamePrefix ${meta.id}. \
      --outSAMtype BAM SortedByCoordinate
    """
}

process FEATURECOUNTS {
    tag "${meta.id}"
    label 'rnaseq'
    publishDir "${params.outdir}/counts", mode: 'copy'
    input:
    tuple val(meta), path(bam)
    path gtf
    output:
    tuple val(meta), path("${meta.id}.counts.txt"), emit: counts
    path "${meta.id}.counts.txt.summary", emit: summaries
    script:
    """
    featureCounts -T ${task.cpus} -p --countReadPairs \
      -a ${gtf} -o ${meta.id}.counts.txt ${bam}
    """
}

process DESEQ2 {
    tag "DESeq2"
    label 'rnaseq'
    publishDir "${params.outdir}/deseq2", mode: 'copy'
    input:
    path counts
    path samplesheet
    path script_file
    output:
    path "deseq2_results.csv", emit: results
    path "*.pdf", optional: true, emit: plots
    script:
    """
    Rscript ${script_file} \
      --counts ${counts.join(',')} \
      --samplesheet ${samplesheet} \
      --outdir .
    """
}

process MULTIQC {
    tag "MultiQC"
    label 'rnaseq'
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    input:
    path qc_files
    output:
    path 'multiqc_report.html'
    path 'multiqc_data'
    script:
    """
    multiqc . --force --outdir .
    """
}

workflow {
    if (!params.input || !params.fasta || !params.gtf) {
        error "Required: --input samplesheet.csv --fasta genome.fa --gtf annotation.gtf"
    }

    samplesheet = file(params.input, checkIfExists: true)
    fasta = file(params.fasta, checkIfExists: true)
    gtf = file(params.gtf, checkIfExists: true)

    reads_ch = Channel
      .fromPath(params.input, checkIfExists: true)
      .splitCsv(header: true)
      .map { row ->
          if (!row.sample || !row.condition || !row.fastq_1 || !row.fastq_2)
              error "Samplesheet columns required: sample,condition,fastq_1,fastq_2"
          def meta = [id: row.sample.toString(), condition: row.condition.toString()]
          tuple(meta,
                [file(row.fastq_1.toString(), checkIfExists: true),
                 file(row.fastq_2.toString(), checkIfExists: true)])
      }

    FASTQC_RAW(reads_ch)
    FASTP(reads_ch)
    STAR_INDEX(fasta, gtf)
    STAR_ALIGN(FASTP.out.reads, STAR_INDEX.out.index)
    FEATURECOUNTS(STAR_ALIGN.out.bam, gtf)

    if (!params.skip_deseq2) {
        count_paths = FEATURECOUNTS.out.counts.map { meta, f -> f }.collect()
        DESEQ2(count_paths, samplesheet, file("${projectDir}/bin/deseq2.R"))
    }

    multiqc_inputs = FASTQC_RAW.out.reports
        .flatMap { meta, html, zip -> [html, zip] }
        .mix(FASTP.out.reports.flatMap { meta, html, json -> [html, json] })
        .mix(STAR_ALIGN.out.log.map { meta, log -> log })
        .mix(FEATURECOUNTS.out.summaries)
        .collect()
    MULTIQC(multiqc_inputs)
}
