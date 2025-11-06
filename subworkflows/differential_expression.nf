process build_minimap_index_transcriptome {
    label "isoforms"
    cpus params.threads
    memory "31 GB"
    input:
        path reference
    output:
        tuple path("genome_index.mmi"), path(reference), emit: index
    script:
    """
    minimap2 -t "${task.cpus}" ${params.minimap2_index_opts} -I 1000G -d genome_index.mmi "${reference}"
    """
}

process map_transcriptome {
    label "isoforms"
    cpus params.threads
    memory "16 GB"
    input:
        tuple val(meta), path(fastq_reads), path(index)
    output:
        tuple val(meta), path("${meta.alias}_reads_aln_sorted.bam"), emit: bam
        path("${meta.alias}.flagstat.stats"), emit: align_stats
    script:
    """
    minimap2 -t ${task.cpus} -ax splice -uf -p 1.0 "${index}" "${fastq_reads}" \
        | samtools view -Sb - > output.bam
    samtools sort -@ ${task.cpus} output.bam -o "${meta.alias}_reads_aln_sorted.bam"
    samtools flagstat -O json "${meta.alias}_reads_aln_sorted.bam" > "${meta.alias}.flagstat.stats"
    """
}

process count_transcripts {
    label "isoforms"
    cpus params.threads
    memory "31 GB"
    input:
        tuple val(meta), path(bam), path(ref_transcriptome)
    output:
        path "*transcript_counts.tsv", emit: counts
    script:
    """
    echo ">>> Running Salmon quantification for \${meta.alias}"
    salmon quant --noErrorModel -p "${task.cpus}" -t "${ref_transcriptome}" -l SF -a "${bam}" -o counts
    mv counts/quant.sf "${meta.alias}.transcript_counts.tsv"
    """
}

process mergeCounts {
    label "isoforms"
    cpus 1
    memory "2 GB"
    input:
        path counts
    output:
        path "unfiltered_transcript_counts.tsv"
    script:
    """
    workflow-glue merge_count_tsvs -z -o unfiltered_transcript_counts.tsv -tsvs ${counts}
    """
}

process mergeTPM {
    label "isoforms"
    cpus 1
    memory "2 GB"
    input:
        path counts
    output:
        path "unfiltered_tpm_transcript_counts.tsv"
    script:
    """
    workflow-glue merge_count_tsvs -o unfiltered_tpm_transcript_counts.tsv -z -tpm True -tsvs $counts
    """
}

workflow differential_expression {
    take:
        ref_transcriptome
        full_len_reads
        sample_sheet
        ref_annotation

    main:
        log.info "▶ Running simplified quantification-only workflow (Salmon quantification only)."

        t_index = build_minimap_index_transcriptome(ref_transcriptome)
        mapped = map_transcriptome(
            full_len_reads.combine(t_index)
                .map { meta, fastq, reference, transcriptome -> tuple(meta, fastq, reference) }
        )

        counted = count_transcripts(
            mapped.bam.combine(t_index.map { mmi, reference -> reference })
        )

        merged = mergeCounts(counted.out.counts.collect())
        merged_TPM = mergeTPM(counted.out.counts.collect())

        de_report = merged_TPM.collect()
        de_outputs_concat = merged_TPM.collect()
        collected_de_alignment_stats = mapped.align_stats.collect()

    emit:
        all_de = de_report
        de_alignment_stats = collected_de_alignment_stats
        de_outputs = de_outputs_concat
}
