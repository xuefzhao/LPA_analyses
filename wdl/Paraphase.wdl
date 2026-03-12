version 1.0

import "Structs.wdl"

workflow paraphase_workflow {
    input {
        File bam
        File bai
        File fasta
        File fasta_fai

        String sample_name = "paraphase_output"

        String docker_image

        RuntimeAttr? runtime_attr_paraphase
    }

    call RunParaphase {
        input:
            bam = bam,
            bai = bai,
            fasta = fasta,
            fasta_fai = fasta_fai,
            docker_image = docker_image,
            sample_name = sample_name,
            runtime_attr_override = runtime_attr_paraphase

    }

    output {
        File paraphase_tar_gz = RunParaphase.tarball
    }
}

task RunParaphase {
    input {
        File bam
        File bai
        File fasta
        File fasta_fai
        String  docker_image
        String sample_name
        RuntimeAttr? runtime_attr_override

    }

    command <<<
        set -euo pipefail

        mkdir output_dir

        paraphase \
            -b ~{bam} \
            -o output_dir \
            -r ~{fasta}

        tar -cvf ~{sample_name}.tar output_dir
        bgzip ~{sample_name}.tar
    >>>

    output {
        File tarball = "~{sample_name}.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(bam, "GB"))*2 + 10,
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

