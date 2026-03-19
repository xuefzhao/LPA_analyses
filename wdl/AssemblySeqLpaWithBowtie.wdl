version 1.0

import "Structs.wdl"

import "AssemblySeqLpaAnalyses.wdl" as AssemblySeqLpaAnalyses

workflow LPA_alignment_pipeline {

    input {

        String bam1
        String bam2
        String bai1
        String bai2
        File fasta1
        File fasta2

        String chrom
        Int start
        Int end
        Int flank_length
        String genome_prefix             # HG01573
        String midfix_gene_name            # eg. LPA, SMN1_SMN2

        File annotation_file
        File exon_fasta                  # LPA_seq.cds.fa

        File polish_bam_script           # polish_bam.py
        File extract_exon_script         # extract_exon_seq.from_bam.py
        File recognize_pattern_Rscript    # recognize_gene_pattern.py

        String bowtie_docker
        String python_docker
        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_extract_region
        RuntimeAttr? runtime_attr_get_unique_reads
        RuntimeAttr? runtime_attr_extract_fasta_reads
        RuntimeAttr? runtime_attr_annotate_sequences
        RuntimeAttr? runtime_attr_cut_sequences_to_regions
        RuntimeAttr? runtime_attr_combine_sequences
        RuntimeAttr? runtime_attr_add_prefix_to_read
        RuntimeAttr? runtime_attr_override
    }


    call AssemblySeqLpaAnalyses.ExtractRegion as extract_region1 {
        input:
            bam = bam1,
            bai = bai1,
            chrom = chrom,
            start = start,
            end = end,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_extract_region
    }

    call AssemblySeqLpaAnalyses.ExtractRegion as extract_region2 {
        input:
            bam = bam2,
            bai = bai2,
            chrom = chrom,
            start = start,
            end = end,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_extract_region
    }

    call AssemblySeqLpaAnalyses.GetUniqueReads as unique_reads1 {
        input:
            bam = extract_region1.region_bam,
            docker_image = sv_pipeline_base_docker,
            chrom = chrom, 
            start = start,
            end = end,
            runtime_attr_override = runtime_attr_get_unique_reads
    }

    call AssemblySeqLpaAnalyses.GetUniqueReads as unique_reads2 {
        input:
            bam = extract_region2.region_bam,
            docker_image = sv_pipeline_base_docker,
            chrom = chrom,
            start = start,
            end = end,
            runtime_attr_override = runtime_attr_get_unique_reads
    }

    call AssemblySeqLpaAnalyses.ExtractFastaReads as fasta_reads1 {
        input:
            fasta = fasta1,
            read_list = unique_reads1.reads,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_extract_fasta_reads
    }

    call AssemblySeqLpaAnalyses.ExtractFastaReads as fasta_reads2 {
        input:
            fasta = fasta2,
            read_list = unique_reads2.reads,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_extract_fasta_reads
    }


    call AssemblySeqLpaAnalyses.CombineSequences {
        input:
            prefix = genome_prefix,
            sequences = [fasta_reads1.target_fasta, fasta_reads2.target_fasta],
            hap_labels = ["hap1", "hap2"],
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_combine_sequences
    }

    call AssemblySeqLpaAnalyses.AddPrefixIfMissing{
        input:
            fasta = CombineSequences.combined_fasta,
            prefix_string = genome_prefix,
            docker_image = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_add_prefix_to_read
    }


    call BuildBowtieIndex {
        input:
            fasta = AddPrefixIfMissing.updated_fasta,
            prefix = genome_prefix,
            docker_image = bowtie_docker,
            runtime_attr_override = runtime_attr_override
    }

    call AlignExonSeq {
        input:
            prefix = genome_prefix,
            midfix = midfix_gene_name,
            index_tar = BuildBowtieIndex.index_tar,
            exon_fasta = exon_fasta,
            docker_image = bowtie_docker,
            runtime_attr_override = runtime_attr_override
    }

    call PolishBam {
        input:
            bam_file = AlignExonSeq.bam,
            script = polish_bam_script,
            docker_image = python_docker,
            runtime_attr_override = runtime_attr_override
    }

    call ExtractGeneStructure {
        input:
            polished_bam = PolishBam.polished_bam,
            genome_fasta = AddPrefixIfMissing.updated_fasta,
            script = extract_exon_script,
            docker_image = python_docker,
            runtime_attr_override = runtime_attr_override
    }

    call RecognizeGenePattern {
        input:
            tsv_file = ExtractGeneStructure.tsv,
            script = recognize_pattern_Rscript,
            docker_image = python_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File target_fa = AddPrefixIfMissing.updated_fasta
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
        mkdir -p index_dir
        bowtie2-build ~{fasta} index_dir/~{prefix}
        tar czf ~{prefix}_index.tar.gz -C index_dir .

    >>>

    output {
        File index_tar = "~{prefix}_index.tar.gz"
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
        String midfix
        File index_tar
        File exon_fasta
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        mkdir -p index_dir
        tar xzf ~{index_tar} -C index_dir
        cd index_dir
        bowtie2 -x ~{prefix} -f -U ~{exon_fasta} -k 100 | samtools sort -o ../~{midfix}.cds.vs.~{prefix}.bam
    >>>

    output {
        File bam = "~{midfix}.cds.vs.~{prefix}.bam"
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

    String prefix = basename(bam_file, ".bam")
    command <<<
        set -euo pipefail
        python3 ~{script} ~{bam_file} ~{prefix}.polished.bam ~{prefix}.vcf
    >>>

    output {
        File polished_bam = "~{prefix}.polished.bam"
        File vcf = "~{prefix}.vcf"
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

    String prefix = basename(polished_bam, ".bam")
    command <<<
        set -euo pipefail
        python3 ~{script} -b ~{polished_bam} -f ~{genome_fasta} -o ~{prefix}.tsv
    >>>

    output {
        File tsv = "~{prefix}.tsv"
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

    String prefix = basename(tsv_file, ".tsv")
    command <<<
        set -euo pipefail
        Rscript ~{script} -i ~{tsv_file} -o ~{prefix}.struc
    >>>

    output {
        File output_struc = "~{prefix}.struc"
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