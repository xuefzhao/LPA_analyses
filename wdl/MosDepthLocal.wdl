version 1.0

import "Structs.wdl"
import "RemoteTabixWithGatk.wdl" as RemoteTabixWithGatk
import "ExtractRegionFromBam.wdl" as ExtractRegionFromBam

workflow MosDepthLocal {
    input {

        Array[File] bam_list          # List of BAMs
        Array[File] bai_list          # List of BAIs (same order as BAM)
        Int start
        Int end
        String chrom
        String mid_fix
        String gatk_docker
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_tabix_bam
        RuntimeAttr? runtime_attr_extract_region
        RuntimeAttr? runtime_attr_bam_to_fastq

        Boolean? quantize_mode = false
        String sv_pipeline_base_docker
    }

    call ExtractRegionFromBam.ExtractRegionFromBam as ExtractRegionFromBam{
        input:
            bam_list = bam_list,
            bai_list = bai_list,
            start = start,
            end = end, 
            chrom = chrom, 
            mid_fix  = mid_fix,
            gatk_docker = gatk_docker,
            sv_pipeline_base_docker = sv_pipeline_base_docker,
            runtime_attr_tabix_bam = runtime_attr_tabix_bam,
            runtime_attr_extract_region = runtime_attr_extract_region,
            runtime_attr_bam_to_fastq = runtime_attr_bam_to_fastq

    }

    scatter(i in range(length(ExtractRegionFromBam.regional_bams))) {
        call RunMosDepth {
            input:
                bam = ExtractRegionFromBam.regional_bams[i],
                bai = ExtractRegionFromBam.regional_bais[i],
                contig = chrom,  # pass the original region
                prefix = prefix,
                quantize_mode = quantize_mode
        }
    }

    output {
        Array[File] mosdepth_dist = RunMosDepth.dist
        Array[File] mosdepth_summary = RunMosDepth.summary
        Array[File] mosdepth_per_base = RunMosDepth.per_base
        Array[File] mosdepth_per_base_csi = RunMosDepth.per_base_csi
        Array[File]? mosdepth_quantized_bed = RunMosDepth.quantized_bed
        Array[File]? mosdepth_quantized_bed_csi = RunMosDepth.quantized_bed_csi
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