version 1.0

import "Structs.wdl"

workflow paraphase_workflow {
    input {
        File bam
        File bai
        String prefix
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

    call ProcessTarGz {
        input:
            input_tar = RunParaphase.tarball,
            prefix = prefix,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_process_tar_gz
    }

    output {
        File paraphase_bam = ProcessTarGz.bam
        File paraphase_bai = ProcessTarGz.bam_bai
        File paraphase_json = ProcessTarGz.json
        File paraphase_vcfs = ProcessTarGz.vcfs_tar
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

task ProcessTarGz {

    input {
        File input_tar
        String prefix
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # unpack input
        tar -xzf ~{input_tar}

        # locate files
        BAM=$(ls output_dir/*.bam)
        BAI=$(ls output_dir/*.bam.bai)
        JSON=$(ls output_dir/*.json)

        # find folder containing "_vcfs"
        VCF_DIR=$(ls -d output_dir/*_vcfs*)

        # compress the vcfs folder
        tar -czf ~{prefix}.vcfs.tar.gz ${VCF_DIR}
    >>>

    output {
        File bam = glob("output_dir/*.bam")[0]
        File bam_bai = glob("output_dir/*.bam.bai")[0]
        File json = glob("output_dir/*.json")[0]
        File vcfs_tar = "~{prefix}.vcfs.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(input_tar, "GB"))*4 + 10,
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
