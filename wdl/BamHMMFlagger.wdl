version 1.0

import "Structs.wdl"

workflow BamHMMFlagger {
    input {
        File bam
        File bai
        File fasta
        File fai

        String flagger_docker
        RuntimeAttr? runtime_attr_run_hmm_flagger

    }

    call RunHMMFlagger {
        input:
            bam = bam,
            bai = bai,
            fasta = fasta,
            fai = fai,

            docker_image = flagger_docker,
            runtime_attr_override = runtime_attr_run_hmm_flagger
    }

    output {
        File coverage_file = RunHMMFlagger.coverage_file
        File coverage_index = RunHMMFlagger.coverage_index
        File hmm_outputs = RunHMMFlagger.hmm_outputs
    }
}


task RunHMMFlagger {

    input {
        File bam
        File bai
        File fasta
        File fai

        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Generate whole genome BED from FASTA index
        awk '{print $1"\t0\t"$2}' ~{fai} > whole_genome.bed

        # Create annotation JSON
        cat > annotations_path.json <<EOF
{
  "whole_genome": "whole_genome.bed"
}
EOF

        # Run bam2cov
        bam2cov \
            --bam ~{bam} \
            --output coverage_file.cov.gz \
            --annotationJson annotations_path.json \
            --threads 16 \
            --baselineAnnotation whole_genome

        # Run hmm_flagger
        mkdir -p hmm_flagger_outputs

        hmm_flagger \
            --input coverage_file.cov.gz \
            --outputDir hmm_flagger_outputs \
            --alphaTsv /home/programs/config/alpha_optimum_trunc_exp_gaussian_w_4000_n_50.tsv \
            --labelNames Err,Dup,Hap,Col \
            --threads 16

        tar czvf hmm_flagger.tar.gz hmm_flagger_outputs/
    >>>

    output {
        File coverage_file = "coverage_file.cov.gz"
        File coverage_index = "coverage_file.cov.gz.tbi"
        File hmm_outputs = "hmm_flagger.tar.gz"
    }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 20 + ceil(size(bam, "GiB")*3),
    disk_gb: 25 + ceil(size(bam, "GiB")*3),
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
