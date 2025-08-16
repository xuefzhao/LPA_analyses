version 1.0

import "Structs.wdl"

workflow AlignAsm {
    input {
        String input_sample_name

        File input_asm_h1
        File input_asm_h2
        File ref

        Int reads_per_file = 50

        String minimap_docker
        String sv_pipeline_base_docker

        File? monitoring_script

        RuntimeAttr? runtime_attr_split_fasta
        RuntimeAttr? runtime_attr_align
        RuntimeAttr? runtime_attr_compress_index
        RuntimeAttr? runtime_attr_merge_bam
    }
    meta {
        workflow_description: "Creates callerset for a single sample"
    }

    call SplitFasta as split_fasta_h1{
        input:
            fasta_file = input_asm_h1,
            reads_per_file = reads_per_file,
            docker_file = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_split_fasta
    }

    call SplitFasta as split_fasta_h2{
        input:
            fasta_file = input_asm_h2,
            reads_per_file = reads_per_file,
            docker_file = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_split_fasta
    }

    scatter (fa_h1 in split_fasta_h1.output_fastas){
        call alignToRef as alignToRefH1 {
            input:
                asmIn = fa_h1,
                sample = input_sample_name,
                hap = "h1",
                ref = ref,
                threads = 4,
                docker_file = minimap_docker,
                monitoring_script = monitoring_script,
                runtime_attr_override = runtime_attr_align
        }

        call compressAndIndex as compressAndIndexH1 {
            input:
                samIn = alignToRefH1.samOut,
                sample = input_sample_name,
                hap = "h1",
                docker_file = sv_pipeline_base_docker,
                monitoring_script = monitoring_script,
                runtime_attr_override = runtime_attr_compress_index
        }
    }

    scatter (fa_h2 in split_fasta_h2.output_fastas){
        call alignToRef as alignToRefH2 {
            input:
                asmIn = fa_h2,
                sample = input_sample_name,
                hap = "h2",
                ref = ref,
                threads = 4,
                docker_file = minimap_docker,
                monitoring_script = monitoring_script,
                runtime_attr_override = runtime_attr_align
        }

        call compressAndIndex as compressAndIndexH2 {
            input:
                samIn = alignToRefH2.samOut,
                sample = input_sample_name,
                hap = "h2",
                docker_file = sv_pipeline_base_docker,
                monitoring_script = monitoring_script,
                runtime_attr_override = runtime_attr_compress_index
        }
    }

    call MergeBams as merge_bam_h1{
        input:
            bam_files = compressAndIndexH1.bamOut,
            bai_files = compressAndIndexH1.baiOut,
            prefix = "~{input_sample_name}_h1",
            docker_file = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_merge_bam
    }

    call MergeBams as merge_bam_h2{
        input:
            bam_files = compressAndIndexH2.bamOut,
            bai_files = compressAndIndexH2.baiOut,
            prefix = "~{input_sample_name}_h2",
            docker_file = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_merge_bam
    }


    output {
        File sample_h1_bam = merge_bam_h1.merged_bam
        File sample_h1_bai = merge_bam_h1.merged_bai
        File sample_h2_bam = merge_bam_h2.merged_bam
        File sample_h2_bai = merge_bam_h1.merged_bai
    }
}

 
task SplitFasta {
    input {
        File fasta_file
        Int reads_per_file
        String docker_file
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(fasta_file)

    command <<<
        set -euxo pipefail

        # Index fasta
        samtools faidx ~{fasta_file}

        # Split fai into groups of N reads
        mkdir split_output
        split -l ~{reads_per_file} ~{fasta_file}.fai split_output/fai_chunk_

        # Extract reads for each chunk
        for f in split_output/fai_chunk_*; do
            chunk_name=$(basename $f)
            out="split_output/${chunk_name}.fasta"
            cut -f1 $f | xargs samtools faidx ~{fasta_file} > $out
        done



    >>>

    output {
        Array[File] output_fastas = glob("split_output/*.fasta")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 20 + ceil(size(fasta_file, "GiB")*3),
        disk_gb: 40 + ceil(size(fasta_file, "GiB")*3),
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
        mem_gb: 100 + ceil(size(asmIn, "GiB")*50),
        disk_gb: 200 + ceil(size(asmIn, "GiB")*50),
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
        File baiOut = "~{sample}-asm_~{hap}.minimap2.bam.bai"
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

task MergeBams {
    input {
        Array[File] bam_files
        Array[File] bai_files
        String prefix
        String docker_file
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        samtools merge -@ 4 merged.bam ~{sep=' ' bam_files}
        mv merged.bam "~{prefix}.bam"
        samtools index "~{prefix}.bam"
    >>>

    output {
        File merged_bam = "~{prefix}.bam"
        File merged_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 20 + ceil(size(bam_files, "GiB")*5),
        disk_gb: 40 + ceil(size(bam_files, "GiB")*5),
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