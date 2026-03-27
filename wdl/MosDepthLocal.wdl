version 1.0

import "Structs.wdl"
import "RemoteTabixWithGatk.wdl" as RemoteTabixWithGatk
import "ExtractRegionFromBam.wdl" as ExtractRegionFromBam

workflow MosDepthLocal {
    input {
        File bam
        File bai
        String prefix   
        Int start
        Int end
        String chrom
        String midfix

        # If regional_bam (and regional_bai) are already extracted, skip ExtractRegion
        File? regional_bam
        File? regional_bai

        File Rscript_calculate_kiv2_depth

        String gatk_docker
        RuntimeAttr? runtime_attr_extract_region

        Boolean? quantize_mode = false
    }

    if (!defined(regional_bam)) {
        call ExtractRegionFromBam.ExtractRegion as ExtractRegion {
            input:
                bam = bam,
                bai = bai,
                start = start,
                end = end,
                chrom = chrom,
                mid_fix = midfix,
                docker_image = gatk_docker,
                runtime_attr_override = runtime_attr_extract_region
        }
    }

    File bam_for_mosdepth = select_first([regional_bam, ExtractRegion.regional_bam])
    File bai_for_mosdepth = select_first([regional_bai, ExtractRegion.regional_bai])

    call RunMosDepth {
            input:
                bam = bam_for_mosdepth,
                bai = bai_for_mosdepth,
                contig = chrom,  # pass the original region
                prefix = "~{prefix}.~{midfix}",
                quantize_mode = quantize_mode
    }

    call EstimateKIV2Depth{
        input:
            depth_file = RunMosDepth.per_base,
            Rscript_calculate_kiv2_depth = Rscript_calculate_kiv2_depth,
            pos1 = 160531482, #start of LPA gene
            pos2 = 160611722, #start of KIV2 units
            pos3 = 160650498, #end of KIV1 units
            pos4 = 160664275, #end of LPA gene
            sample_id = prefix 
    }

    output {
        File mosdepth_dist = RunMosDepth.dist
        File mosdepth_summary = RunMosDepth.summary
        File mosdepth_per_base = RunMosDepth.per_base
        File mosdepth_per_base_csi = RunMosDepth.per_base_csi
        File? mosdepth_quantized_bed = RunMosDepth.quantized_bed
        File? mosdepth_quantized_bed_csi = RunMosDepth.quantized_bed_csi
        File lpa_depth_plot = EstimateKIV2Depth.depth_plot
        File lpa_depth_stats = EstimateKIV2Depth.depth_stats
        Float lpa_kiv2_copy = EstimateKIV2Depth.estimated_kiv2

    }
}



task RunMosDepth {
    input {
        File bam
        File bai
        String contig
        String prefix
        Boolean? quantize_mode
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


task EstimateKIV2Depth {
    input {
        File depth_file
        File Rscript_calculate_kiv2_depth
        Int pos1
        Int pos2
        Int pos3
        Int pos4
        String sample_id

        # Runtime parameters
        String docker_image = "rocker/tidyverse:latest"
        Int cpu = 2
        Int memory_gb = 4
    }

    command {
        Rscript ~{Rscript_calculate_kiv2_depth} \
            --input ~{depth_file} \
            --output ~{sample_id}_depth_plot.pdf \
            --pos1 ~{pos1} \
            --pos2 ~{pos2} \
            --pos3 ~{pos3} \
            --pos4 ~{pos4} \
            --stats ~{sample_id}_kiv2_results.txt
    }

    output {
        File depth_plot = "~{sample_id}_depth_plot.pdf"
        File depth_stats = "~{sample_id}_kiv2_results.txt"
        Float estimated_kiv2 = read_map("~{sample_id}_kiv2_results.txt")["estimated_KIV2"]
    }

    runtime {
        docker: docker_image
        cpu: cpu
        memory: "~{memory_gb} GB"
    }
}