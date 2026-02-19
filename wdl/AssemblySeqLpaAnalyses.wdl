version 1.0

import "Structs.wdl"


workflow AssemblySeqLpaAnalyses {

    input {
        File bam1
        File bam2
        File bai1
        File bai2
        File fasta1
        File fasta2
        File annotation_file

        File script_run_blast_from_table
        File script_extract_spanned_regions
        String chrom
        Int start
        Int end
        String prefix
        String docker_image
        RuntimeAttr? runtime_attr_extract_region
        RuntimeAttr? runtime_attr_get_unique_reads
        RuntimeAttr? runtime_attr_extract_fasta_reads
        RuntimeAttr? runtime_attr_annotate_sequences
        RuntimeAttr? runtime_attr_cut_sequences_to_regions
        RuntimeAttr? runtime_attr_combine_sequences
    }

    call ExtractRegion as extract_region1 {
        input:
            bam = bam1,
            bai = bai1,
            chrom = chrom,
            start = start,
            end = end,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_extract_region
    }

    call ExtractRegion as extract_region2 {
        input:
            bam = bam2,
            bai = bai2,
            chrom = chrom,
            start = start,
            end = end,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_extract_region
    }

    call GetUniqueReads as unique_reads1 {
        input:
            bam = extract_region1.region_bam,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_get_unique_reads
    }

    call GetUniqueReads as unique_reads2 {
        input:
            bam = extract_region2.region_bam,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_get_unique_reads
    }

    call ExtractFastaReads as fasta_reads1 {
        input:
            fasta = fasta1,
            read_list = unique_reads1.reads,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_extract_fasta_reads
    }

    call ExtractFastaReads as fasta_reads2 {
        input:
            fasta = fasta2,
            read_list = unique_reads2.reads,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_extract_fasta_reads
    }

    call AnnotateSequences as annotate_seq1{
        input:
            target_fa = fasta_reads1.target_fasta,
            target_fai = fasta_reads1.target_fasta_fai,
            annotation_file = annotation_file,
            run_blast_from_table = script_run_blast_from_table,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_annotate_sequences
    }

    call AnnotateSequences as annotate_seq2{
        input:
            target_fa = fasta_reads2.target_fasta,
            target_fai = fasta_reads2.target_fasta_fai,
            annotation_file = annotation_file,
            run_blast_from_table = script_run_blast_from_table,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_annotate_sequences
    }

    call CutSequencesToRegions as cut_to_region1 {
        input:
            target_fa = fasta_reads1.target_fasta,
            target_fai = fasta_reads1.target_fasta_fai,
            annotations = annotate_seq1.annotated,
            extract_spanned_regions = script_extract_spanned_regions,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_cut_sequences_to_regions
    }

    call CutSequencesToRegions as cut_to_region2 {
        input:
            target_fa = fasta_reads2.target_fasta,
            target_fai = fasta_reads2.target_fasta_fai,
            annotations = annotate_seq2.annotated,
            extract_spanned_regions = script_extract_spanned_regions,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_cut_sequences_to_regions
    }

    call CombineSequences {
        input:
            prefix = prefix,
            sequences = [cut_to_region1.cut_sequences, cut_to_region2.cut_sequences],
            hap_labels = ["hap1", "hap2"],
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_combine_sequences
    }

    output {
        File combined_fasta = CombineSequences.combined_fasta
        File combined_fai = CombineSequences.combined_fai
    }
}

# -------------------
# Tasks
# -------------------

task ExtractRegion {
    input {
        File bam
        File bai
        String chrom
        Int start
        Int end
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String output_prefix = "region_" + chrom + "_" + start + "_" + end

    command <<<
        set -euo pipefail

        gatk PrintReads \
            -I ~{bam} \
            -L ~{chrom}:~{start}-~{end} \
            --interval-merging-rule ALL \
            -O ~{output_prefix}.bam
    >>>

    output {
        File region_bam = "~{output_prefix}.bam"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(bam, "GiB")),
        disk_gb: 20 + ceil(size(bam, "GiB")),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GetUniqueReads {
    input {
        File bam
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String output_prefix = basename(bam, ".bam")

    command <<<
        set -euo pipefail

        samtools view -H ~{bam} > /dev/null
        samtools view ~{bam} | cut -f1 | sort | uniq > ~{output_prefix}_reads.txt
    >>>

    output {
        File reads = "~{output_prefix}_reads.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8 + ceil(size(bam, "GiB")),
        disk_gb: 10 + ceil(size(bam, "GiB")),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExtractFastaReads {
    input {
        File fasta
        File read_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String output_prefix = basename(fasta, ".fa")

    command <<<
        set -euo pipefail

        samtools faidx ~{fasta} $(cat ~{read_list}) > ~{output_prefix}_reads.fa
        samtools faidx ~{output_prefix}_reads.fa

    >>>

    output {
        File target_fasta = "~{output_prefix}_reads.fa"
        File target_fasta_fai = "~{output_prefix}_reads.fa.fai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8 + ceil(size(fasta, "GiB")),
        disk_gb: 10 + ceil(size(fasta, "GiB")),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AnnotateSequences {
    input {
        File target_fa
        File target_fai
        File annotation_file
        File run_blast_from_table
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(target_fa, ".fa")
    command <<<
        set -euo pipefail

        bash ~{run_blast_from_table} ~{annotation_file} ~{target_fa} ~{prefix}.tsv

    >>>

    output {
        File annotated = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CutSequencesToRegions {
    input {
        File target_fa
        File target_fai
        File annotations
        File extract_spanned_regions
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(target_fa,".fa")
    command <<<
        set -euo pipefail

        python ~{extract_spanned_regions} ~{target_fa} ~{annotations}  ~{prefix}.target.fa

    >>>

    output {
        File cut_sequences = "~{prefix}.target.fa"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CombineSequences {
    input {
        Array[File] sequences
        Array[String] hap_labels
        String prefix
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=" " sequences} > ~{prefix}.fa
        samtools faidx ~{prefix}.fa
    >>>

    output {
        File combined_fasta = "~{prefix}.fa"
        File combined_fai = "~{prefix}.fa.fai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}