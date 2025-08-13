ersion 1.0

import "Structs.wdl"

workflow SplitHaplotypeReads{
    input{
        File input_fa
        String sample 

        String sv_base_mini_docker
        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_split_hap
    }


    Call ExtractHaplotypeSeq{
        input:
            input_fa = input_fa,
            sample = sample,
            RuntimeAttr = runtime_attr_split_hap
    } 

    output{
        File hap1_fasta = ExtractHaplotypeSeq.hap1_fa
        File hap2_fasta = ExtractHaplotypeSeq.hap2_fa
    }
}




task ExtractHaplotypeSeq {
    input {
        File input_fa
        String sample
        RuntimeAttr? runtime_attr_override
    }


    String prefix=basename(input_fa, ".fna.gz")

    command <<<
        set -euo pipefail

        # Get prefix without .fna.gz

        # Uncompress and recompress with bgzip
        gunzip -c "~{input_fa}" > "${prefix}.fna"
        bgzip -f "${prefix}.fna"

        # Index fasta
        samtools faidx "${prefix}.fna.gz"

        # Extract sequence IDs
        zcat "${prefix}.fna.gz" | grep '^>' | grep 'h1tg' > seq_ID.h1.tsv
        zcat "${prefix}.fna.gz" | grep '^>' | grep 'h2tg' > seq_ID.h2.tsv

        # Extract sequences for hap1
        samtools faidx "${prefix}.fna.gz" $(awk '{print $1}' seq_ID.h1.tsv | sed -e 's/>//') > "~{sample}.h1.hifiasm.fa"

        # Extract sequences for hap2
        samtools faidx "${prefix}.fna.gz" $(awk '{print $1}' seq_ID.h2.tsv | sed -e 's/>//') > "~{sample}.h2.hifiasm.fa"
    >>>

    output {
        File hap1_fa = "~{sample}.h1.hifiasm.fa"
        File hap2_fa = "~{sample}.h2.hifiasm.fa"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10 + ceil(size(input_fa, "GiB")*4),
        disk_gb: 15 + ceil(size(input_fa, "GiB")*4),
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
}