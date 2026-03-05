version 1.0

import "Structs.wdl"

workflow MosDepthLocal {
    input {
        File bam
        File bai
        String chrom
        int start
        int end
        int flank
        String prefix
        String midfix
        Boolean quantize_mode
        String sv_pipeline_base_docker
    }

    call ExtractRegionWithFlank {
        input:
            bam = bam,
            bai = bai,
            chrom = chrom,
            start = start, 
            end = end,
            flank = flank,
            prefix = prefix,
            midfix = midfix,
            docker_image = sv_pipeline_base_docker
    }

    call RunMosDepth {
        input:
            bam = ExtractRegionWithFlank.region_bam,
            bai = ExtractRegionWithFlank.region_bai,
            contig = chrom,  # pass the original region
            prefix = prefix,
            quantize_mode = quantize_mode
    }


    output {
        File mosdepth_dist = RunMosDepth.dist
        File mosdepth_summary = RunMosDepth.summary
        File mosdepth_per_base = RunMosDepth.per_base
        File mosdepth_per_base_csi = RunMosDepth.per_base_csi
        File? mosdepth_quantized_bed = RunMosDepth.quantized_bed
        File? mosdepth_quantized_bed_csi = RunMosDepth.quantized_bed_csi
    }
}




task ExtractRegionWithFlank {
    input {
        File bam
        File bai
        String chrom  # format: chr:start-end
        Int start
        Int end
        Int flank
        String prefix
        String midfix  # optional string to include in output BAM name
        Int flank = 10000  # 10Kb flank
        String docker_image = "quay.io/biocontainers/samtools:1.17--h8ee4bcc_0"
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # add flanks
        start_flank=$((start - flank))
        if [ $start_flank -lt 0 ]; then start_flank=0; fi
        end_flank=$((end + flank))

        output_bam="~{prefix}.~{midfix}.region.bam"

        samtools view -b ~{bam} ${chr}:${start_flank}-${end_flank} > $output_bam
        samtools index $output_bam
    >>>

    output {
        File region_bam = "~{prefix}.~{midfix}.region.bam"
        File region_bai = "~{prefix}.~{midfix}.region.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4 + ceil(size(bam,"GiB")),
        disk_gb: 10 + ceil(size(bam,"GiB")),
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

task RunMosDepth {
    input {
        File bam
        File bai
        String contig
        String prefix
        Boolean quantize_mode
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ "~{quantize_mode}" == "true" ]; then
            export MOSDEPTH_Q0=NO_COVERAGE
            export MOSDEPTH_Q1=LOW_COVERAGE
            export MOSDEPTH_Q2=CALLABLE
            export MOSDEPTH_Q3=HIGH_COVERAGE

            mosdepth \
                -t 4 \
                -c "~{contig}" \
                -Q 1 \
                -x \
                --quantize 0:1:5:150: \
                ~{prefix} \
                ~{bam}
        else
            mosdepth \
                -t 2 \
                -c "~{contig}" \
                -Q 1 \
                -x \
                ~{prefix} \
                ~{bam}
        fi
    >>>

    output {
        File dist = "~{prefix}.mosdepth.global.dist.txt"
        File summary = "~{prefix}.mosdepth.summary.txt"
        File per_base = "~{prefix}.per-base.bed.gz"
        File per_base_csi = "~{prefix}.per-base.bed.gz.csi"
        File? quantized_bed = "~{prefix}.quantized.bed.gz"
        File? quantized_bed_csi = "~{prefix}.quantized.bed.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(bam, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: "us.gcr.io/broad-dsp-lrma/lr-mosdepth:0.3.1"
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}