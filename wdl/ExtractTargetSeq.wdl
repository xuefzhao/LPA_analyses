version 1.0

import "Structs.wdl"

workflow ExtractTargetSeq {
    input {
        File input_bam_h1
        File input_bai_h1
        File input_bam_h2
        File input_bai_h2

        String prefix
        String midfix

        File extract_target_sequence_py

        String chrom
        Int start
        Int end

        String sv_pipeline_base_docker

        File? monitoring_script

        RuntimeAttr? runtime_attr_align_h1
        RuntimeAttr? runtime_attr_align_h2
        RuntimeAttr? runtime_attr_compress_index_h1
        RuntimeAttr? runtime_attr_compress_index_h2
    }


    call ExtractSeq as extract_seq_from_h1 {
        input:
            bam = input_bam_h1,
            bai = input_bai_h1,
            extract_target_sequence_py = extract_target_sequence_py, 
            chrom = chrom,
            start = start,
            end = end,
            docker_file = sv_pipeline_base_docker,
            monitoring_script = monitoring_script,
            runtime_attr_override = runtime_attr_align_h1
    }

    call ExtractSeq as extract_seq_from_h2 {
        input:
            bam = input_bam_h2,
            bai = input_bai_h2, 
            extract_target_sequence_py = extract_target_sequence_py, 
            chrom = chrom,
            start = start,
            end = end,
            docker_file = sv_pipeline_base_docker,
            monitoring_script = monitoring_script,
            runtime_attr_override = runtime_attr_align_h1
    }

    call ConcatFasta{
        input:
            fasta_list = [extract_seq_from_h1.seq, extract_seq_from_h2.seq],
            output_prefix = "~{prefix}.~{midfix}",
            docker_image = sv_pipeline_base_docker
    }

    output {
        File seq = ConcatFasta.merged_fasta
        File fai = ConcatFasta.fasta_index    
    }
}

 
task ExtractSeq {
    input {
        File bam
        File bai

        File extract_target_sequence_py

        String chrom
        Int start
        Int end

        File? monitoring_script

        String docker_file
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, '.bam')

    command <<<

        set -euxo pipefail

        python ~{extract_target_sequence_py} \
          -b ~{bam} \
          -r ~{chrom}:~{start}-~{end} \
          -o ~{prefix}.fa

    >>>


    output {
        File seq = "~{prefix}.fa"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 20 + ceil(size(bam, "GiB")*3),
        disk_gb: 40 + ceil(size(bam, "GiB")*3),
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
        docker: docker_file
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatFasta {
  input {
    Array[File] fasta_list
    String output_prefix
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(10 + size(fasta_list, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    out_fa="~{output_prefix}.fa"

    for fa in ~{sep=' ' fasta_list}; do
      name=$(basename "${fa}" .fa)

      # Extract sequence, remove headers/whitespace, uppercase
      seq=$(grep -v '^>' "${fa}" \
            | tr -d ' \t\r\n' \
            | tr 'acgtn' 'ACGTN')

      echo ">${name}" >> "${out_fa}"
      echo "${seq}" >> "${out_fa}"
    done 

    samtools faidx "~{output_prefix}.fa"

    >>>

  output {
    File merged_fasta = "~{output_prefix}.fa"
    File fasta_index  = "~{output_prefix}.fa.fai"
  }

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
