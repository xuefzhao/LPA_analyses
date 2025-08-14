version 1.0

import "Structs.wdl"

workflow AlignAsm {
    input {
        String input_sample_name

        File input_asm_h1
        File input_asm_h2
        File ref

        String minimap_docker
        String sv_pipeline_base_docker

        File? monitoring_script

        RuntimeAttr? runtime_attr_align_h1
        RuntimeAttr? runtime_attr_align_h2
        RuntimeAttr? runtime_attr_compress_index_h1
        RuntimeAttr? runtime_attr_compress_index_h2
    }
    meta {
        workflow_description: "Creates callerset for a single sample"
    }

    call alignToRef as alignToRefH1 {
        input:
            asmIn = input_asm_h1,
            sample = input_sample_name,
            hap = "h1",
            ref = ref,
            threads = 4,
            docker_file = minimap_docker,
            monitoring_script = monitoring_script,
            runtime_attr_override = runtime_attr_align_h1
    }
    call compressAndIndex as compressAndIndexH1 {
        input:
            samIn = alignToRefH1.samOut,
            sample = input_sample_name,
            hap = "h1",
            docker_file = sv_pipeline_base_docker,
            monitoring_script = monitoring_script,
            runtime_attr_override = runtime_attr_compress_index_h1
    }

    call alignToRef as alignToRefH2 {
        input:
            asmIn = input_asm_h2,
            sample = input_sample_name,
            hap = "h2",
            ref = ref,
            threads = 4,
            docker_file = minimap_docker,
            monitoring_script = monitoring_script,
            runtime_attr_override = runtime_attr_align_h2
    }
    call compressAndIndex as compressAndIndexH2 {
        input:
            samIn = alignToRefH2.samOut,
            sample = input_sample_name,
            hap = "h2",
            docker_file = sv_pipeline_base_docker,
            monitoring_script = monitoring_script,
            runtime_attr_override = runtime_attr_compress_index_h2
   }

    output {
        File sample_h1_bam = compressAndIndexH1.bamOut
        File sample_h1_idx = compressAndIndexH1.indexOut
        File sample_h2_bam = compressAndIndexH2.bamOut
        File sample_h2_idx = compressAndIndexH2.indexOut
    }
}

 
task alignToRef {
    input {
        File asmIn
        String sample
        File ref
        String hap
        Int threads
        File? monitoring_script
        String docker_file
        RuntimeAttr? runtime_attr_override
    }

    command <<<

        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        minimap2 -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 --secondary=no -a -t ~{threads} --eqx -Y ~{ref} ~{asmIn} > ~{sample + "-asm_" + hap + ".minimap2.sam"}
    >>>


    output {
        File samOut = "~{sample}-asm_~{hap}.minimap2.sam"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(asmIn, "GiB")*5),
        disk_gb: 20 + ceil(size(asmIn, "GiB")*5),
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

task compressAndIndex {
    input {
        File samIn
        String sample
        String hap
        String docker_file
        File? monitoring_script
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        samtools view -b ~{samIn} | samtools sort -O bam -o ~{sample + "-asm_" + hap + ".minimap2.bam"} -
        samtools index ~{sample + "-asm_" + hap + ".minimap2.bam"}

    >>>


    output {
        File bamOut = "~{sample}-asm_~{hap}.minimap2.bam"
        File indexOut = "~{sample}-asm_~{hap}.minimap2.bam.bai"
    }


    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(samIn, "GiB")*2),
        disk_gb: 20 + ceil(size(samIn, "GiB")*2),
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

