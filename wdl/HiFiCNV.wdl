version 1.0

import "Structs.wdl"

workflow HiFiCNV {
    meta {
        desciption:
        "Runs the PacBio HiFiCNV tool on a single (human) HiFi bam."
    }
    input {
        File bam
        File bai
        File ref_fasta
        File ref_fasta_fai
        File exclude_bed
        File sex_specific_cn
        String gcs_out_root_dir
        String docker_pb_hifi_cnv
        String docker_finalize_log
    }
    parameter_meta {
        exclude_bed:      "BED holding regions that are known to cause artifacts during HiFiCNV data processing (e.g. centromeres)."
        sex_specific_cn:  "Sex-specific files annotating the PAR regions on and expected copy numbers of sex chromosomes."
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String workflow_name = 'HiFiCNV'
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}"


    call InferSampleName { 
        input: 
            bam = bam, 
            bai = bai
        }

    call PacBioHiFiCNV { 
        input:
            bam = bam, bai = bai,
            sample_name = InferSampleName.sample_name,
            output_prefix = InferSampleName.sample_name,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            exclude_bed = exclude_bed,
            docker_image = docker_pb_hifi_cnv,
            sex_specific_cn = sex_specific_cn
        }

    call FinalizeToFile as FinalizeLog { 
        input: 
            outdir = outdir + '/~{InferSampleName.sample_name}', 
            file = PacBioHiFiCNV.log,
            docker_image = docker_finalize_log
            }

    call FinalizeToFile as FinalizeVCF { 
        input: 
            outdir = outdir + '/~{InferSampleName.sample_name}', 
            file = PacBioHiFiCNV.vcf,
            docker_image = docker_finalize_log
        }

    call FinalizeToFile as FinalizeBedGraph { 
        input: 
            outdir = outdir + '/~{InferSampleName.sample_name}', 
            file = PacBioHiFiCNV.bedgraph ,
            docker_image = docker_finalize_log
        }

    call FinalizeToFile as FinalizeBigWig { 
        input: 
            outdir = outdir + '/~{InferSampleName.sample_name}', 
            file = PacBioHiFiCNV.depth_bw ,
            docker_image = docker_finalize_log
        }

    output {
        Map[String, String] hificnv_outs = {'vcf': FinalizeVCF.gcs_path,
                                            'bedgraph': FinalizeLog.gcs_path,
                                            'log': FinalizeBedGraph.gcs_path,
                                            'depth_bw': FinalizeBigWig.gcs_path
                                            }
    }
}

task PacBioHiFiCNV {
    input {
        File bam
        File bai
        String sample_name
        String output_prefix
        File ref_fasta
        File ref_fasta_fai

        String docker_image

        File exclude_bed
        File sex_specific_cn
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eux

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        hificnv \
            --bam ~{bam} \
            --ref ~{ref_fasta} \
            --exclude ~{exclude_bed} \
            --expected-cn ~{sex_specific_cn} \
            --threads "${num_core}" \
            --output-prefix ~{output_prefix}

        tree
    >>>
    output {
        File vcf = "~{output_prefix}.${sample_name}.vcf.gz"
        File bedgraph = "~{output_prefix}.${sample_name}.copynum.bedgraph"
        File log = "~{output_prefix}.log"
        File depth_bw = "~{output_prefix}.${sample_name}.depth.bw"
    }

    #########################
    Int min_disk = 40
    Float disk_multiplier = 1
    Int disk_size = ceil(disk_multiplier * size(bam, "GiB")) + 20
    Int use_this_disk_sz = if (min_disk>disk_size) then min_disk else disk_size

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             6,
        disk_gb:            use_this_disk_sz,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/hificnv:1.0.0"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 docker_image
    }
}

task InferSampleName {
    meta {
        description: "Infer sample name encoded on the @RG line of the header section. Fails if multiple values found, or if SM ~= unnamedsample."
    }

    parameter_meta {
        bam: {
            localization_optional: true,
            description: "BAM file"
        }
    }

    input {
        File bam
        File bai
    }



    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} > header.txt
        if ! grep -q '^@RG' header.txt; then echo "No read group line found!" && exit 1; fi

        grep '^@RG' header.txt | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g' | sort | uniq > sample.names.txt
        if [[ $(wc -l sample.names.txt) -gt 1 ]]; then echo "Multiple sample names found!" && exit 1; fi
        if grep -iq "unnamedsample" sample.names.txt; then echo "Sample name found to be unnamedsample!" && exit 1; fi
    >>>

    output {
        String sample_name = read_string("sample.names.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task FinalizeToFile {

    meta{
        description: "Copies the given file to the specified bucket."
    }

    parameter_meta {
        file: {
            description: "file to finalize",
            localization_optional: true
        }
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        outdir: "directory to which files should be uploaded"
        name:   "name to set for uploaded file"
    }

    input {
        File file
        String outdir
        String? name

        String docker_image
        File? keyfile

        RuntimeAttr? runtime_attr_override
    }



    String gcs_output_dir = sub(outdir, "/+$", "")
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, basename(file)])

    command <<<
        set -euxo pipefail

        gsutil -m cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 docker_image
    }
}


