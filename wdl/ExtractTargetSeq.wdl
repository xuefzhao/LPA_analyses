version 1.0

import "Structs.wdl"

workflow ExtractTargetSeq {
    input {
        File input_bam_h1
        File input_bai_h1
        File input_bam_h2
        File input_bai_h2

        String chrom
        Int start
        Int end

        String minimap_docker
        String sv_pipeline_base_docker

        File? monitoring_script

        RuntimeAttr? runtime_attr_align_h1
        RuntimeAttr? runtime_attr_align_h2
        RuntimeAttr? runtime_attr_compress_index_h1
        RuntimeAttr? runtime_attr_compress_index_h2
    }


    call ExtractSeq as extract_seq_from_h1 {
        input:
            bam = input_bam_h1,
            bai = input_bai_h1, 
            chrom = chrom,
            start = start,
            end = end,
            docker_file = sv_pipeline_base_docker,
            monitoring_script = monitoring_script,
            runtime_attr_override = runtime_attr_align_h1
    }

    call ExtractSeq as extract_seq_from_h2 {
        input:
            bam = input_bam_h2,
            bai = input_bai_h2, 
            chrom = chrom,
            start = start,
            end = end,
            docker_file = sv_pipeline_base_docker,
            monitoring_script = monitoring_script,
            runtime_attr_override = runtime_attr_align_h1
    }



    output {
        File seq_h1 = extract_seq_from_h1.seq
        File seq_h2 = extract_seq_from_h2.seq
    }
}

 
task ExtractSeq {
    input {
        File bam
        File bai

        String chrom
        Int start
        Int end

        File? monitoring_script

        String docker_file
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, '.bam')

    command <<<

        set -euxo pipefail


        python3 <<CODE

        import pysam

        def stitch_reads_to_region(bam_file, chromosome, start, end):
            """
            Stitch together all reads covering a target region into a single reference-aligned sequence.
            Gaps not covered by any read are filled with 'N'.
            Parameters:
                bam_file (str): Path to BAM file.
                chromosome (str): Chromosome name (e.g., 'chr1').
                start (int): 0-based start position of the region.
                end (int): End position (exclusive).
            Returns:
                str: Stitched sequence of length (end-start), with uncovered positions as 'N'.
            """

            region_length = end - start
            stitched_seq = ['N'] * region_length  # initialize with N
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam.fetch(chromosome, start, end):
                    if read.is_unmapped:
                        continue
                    ref_pos = read.reference_start
                    read_pos = 0
                    for (cigar_op, length) in read.cigartuples:
                        if cigar_op in (0, 7, 8):  # M, =, X
                            for i in range(length):
                                if start <= ref_pos < end:
                                    idx = ref_pos - start
                                    stitched_seq[idx] = read.query_sequence[read_pos]
                                ref_pos += 1
                                read_pos += 1
                        elif cigar_op in (1, 4, 5):  # I, S, H
                            read_pos += length
                        elif cigar_op in (2, 3):  # D, N
                            ref_pos += length
                        elif cigar_op == 6:  # P
                            continue
                        else:
                            raise ValueError(f"Unknown CIGAR operation: {cigar_op}")
            return ''.join(stitched_seq)

        def write_stitched_sequence_to_fasta(stitched_sequence, prefix, chromosome, start, end, output_file):
            """
            Write stitched sequence to a FASTA file.
            """
            header = f">~{prefix}_{chromosome}:{start}-{end}"
            with open(output_file, "w") as f:
                f.write(f"{header}\n")
                # Write sequence in 60 bp lines (FASTA standard)
                f.write(stitched_sequence)

        start = int("~{start}")
        end = int("~{end}")

        stitched_sequence = stitch_reads_to_region("~{bam}", "~{chrom}", start, end)

        write_stitched_sequence_to_fasta(stitched_sequence, "~{prefix}", "~{chrom}", start, end, "~{prefix}.seq" )


        CODE

    >>>


    output {
        File seq = "~{prefix}.seq"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 20 + ceil(size(bam, "GiB")*3),
        disk_gb: 40 + ceil(size(bam, "GiB")*3),
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

