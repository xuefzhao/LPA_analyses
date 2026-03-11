version 1.0

workflow paraphase_workflow {
    input {
        File bam
        File bai
        File fasta
        File fasta_fai
        String docker_image
        String sample_name = "paraphase_output"
    }

    call RunParaphase {
        input:
            bam = bam,
            bai = bai,
            fasta = fasta,
            fasta_fai = fasta_fai,
            docker_image = docker_image,
            sample_name = sample_name
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

    runtime {
        docker: docker_image
        cpu: 4
        memory: "8G"
    }
}

