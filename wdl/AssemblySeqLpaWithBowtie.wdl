version 1.0

import "Structs.wdl"

workflow LPA_alignment_pipeline {

    input {
        File genome_fasta                # HG01573.fa
        String genome_prefix             # HG01573
        File exon_fasta                  # LPA_seq.cds.fa
        File polish_bam_script           # polish_bam.py
        File extract_exon_script         # extract_exon_seq.from_bam.py
        File recognize_pattern_script    # recognize_gene_pattern.py
        String bowtie_docker
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    call BuildBowtieIndex {
        input:
            fasta = genome_fasta,
            prefix = genome_prefix,
            docker_image = bowtie_docker,
            runtime_attr_override = runtime_attr_override
    }

    call AlignExonSeq {
        input:
            prefix = genome_prefix,
            exon_fasta = exon_fasta,
            docker_image = bowtie_docker,
            runtime_attr_override = runtime_attr_override
    }

    call PolishBam {
        input:
            bam_file = AlignExonSeq.bam,
            script = polish_bam_script,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_override
    }

    call ExtractGeneStructure {
        input:
            polished_bam = PolishBam.polished_bam,
            genome_fasta = genome_fasta,
            script = extract_exon_script,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_override
    }

    call RecognizeGenePattern {
        input:
            tsv_file = ExtractGeneStructure.tsv,
            script = recognize_pattern_script,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File polished_bam = PolishBam.polished_bam
        File polished_vcf = PolishBam.vcf
        File polished_tsv = ExtractGeneStructure.tsv
        File gene_structure = RecognizeGenePattern.output_struc
    }
}

# ----------------------------
# Task 1: Build Bowtie2 Index
# ----------------------------
task BuildBowtieIndex {
    input {
        File fasta
        String prefix
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        bowtie2-build ~{fasta} ~{prefix}
    >>>

    output {
        File index_prefix = "~{prefix}.1.bt2"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(fasta,"GiB")),
        disk_gb: 20 + ceil(size(fasta,"GiB")),
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

# ----------------------------
# Task 2: Align Exon Sequences
# ----------------------------
task AlignExonSeq {
    input {
        String prefix
        File exon_fasta
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        bowtie2 -x ~{prefix} -f -U ~{exon_fasta} -k 100 | samtools sort -o LPA_seq.cds.vs.~{prefix}.bam
    >>>

    output {
        File bam = "LPA_seq.cds.vs.~{prefix}.bam"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 32 + ceil(size(exon_fasta,"GiB")),
        disk_gb: 50 + ceil(size(exon_fasta,"GiB")),
        boot_disk_gb: 20,
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

# ----------------------------
# Task 3: Polish BAM
# ----------------------------
task PolishBam {
    input {
        File bam_file
        File script
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 ~{script} ~{bam_file} LPA_seq.cds.vs.HG01573.polished.bam LPA_seq.cds.vs.HG01573.vcf
    >>>

    output {
        File polished_bam = "LPA_seq.cds.vs.HG01573.polished.bam"
        File vcf = "LPA_seq.cds.vs.HG01573.vcf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(bam_file,"GiB")),
        disk_gb: 20 + ceil(size(bam_file,"GiB")),
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

# ----------------------------
# Task 4: Extract Gene Structure
# ----------------------------
task ExtractGeneStructure {
    input {
        File polished_bam
        File genome_fasta
        File script
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 ~{script} -b ~{polished_bam} -f ~{genome_fasta} -o LPA_seq.cds.vs.HG01573.polished.tsv
    >>>

    output {
        File tsv = "LPA_seq.cds.vs.HG01573.polished.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(polished_bam,"GiB")),
        disk_gb: 20 + ceil(size(polished_bam,"GiB")),
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

# ----------------------------
# Task 5: Recognize Gene Pattern
# ----------------------------
task RecognizeGenePattern {
    input {
        File tsv_file
        File script
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 ~{script} ~{tsv_file} LPA_seq.cds.vs.HG01573.polished.struc
    >>>

    output {
        File output_struc = "LPA_seq.cds.vs.HG01573.polished.struc"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8 + ceil(size(tsv_file,"GiB")),
        disk_gb: 10 + ceil(size(tsv_file,"GiB")),
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