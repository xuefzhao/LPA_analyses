version 1.0

import "Structs.wdl"


workflow AssemblySeqLpaAnalyses {

    input {
        String bam1
        String bam2
        String bai1
        String bai2
        File fasta1
        File fasta2
        File annotation_file

        File script_run_blast_from_table
        File script_extract_spanned_regions
        String chrom
        Int start
        Int end
        String prefix
        String docker_image
        RuntimeAttr? runtime_attr_extract_region
        RuntimeAttr? runtime_attr_get_unique_reads
        RuntimeAttr? runtime_attr_extract_fasta_reads
        RuntimeAttr? runtime_attr_annotate_sequences
        RuntimeAttr? runtime_attr_cut_sequences_to_regions
        RuntimeAttr? runtime_attr_combine_sequences
        RuntimeAttr? runtime_attr_add_prefix_to_read
    }

    call ExtractRegion as extract_region1 {
        input:
            bam = bam1,
            bai = bai1,
            chrom = chrom,
            start = start,
            end = end,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_extract_region
    }

    call ExtractRegion as extract_region2 {
        input:
            bam = bam2,
            bai = bai2,
            chrom = chrom,
            start = start,
            end = end,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_extract_region
    }

    call GetUniqueReads as unique_reads1 {
        input:
            bam = extract_region1.region_bam,
            docker_image = docker_image,
            chrom = chrom, 
            start = start,
            end = end,
            runtime_attr_override = runtime_attr_get_unique_reads
    }

    call GetUniqueReads as unique_reads2 {
        input:
            bam = extract_region2.region_bam,
            docker_image = docker_image,
            chrom = chrom,
            start = start,
            end = end,
            runtime_attr_override = runtime_attr_get_unique_reads
    }

    call ExtractFastaReads as fasta_reads1 {
        input:
            fasta = fasta1,
            read_list = unique_reads1.reads,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_extract_fasta_reads
    }

    call ExtractFastaReads as fasta_reads2 {
        input:
            fasta = fasta2,
            read_list = unique_reads2.reads,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_extract_fasta_reads
    }

    call AnnotateSequences as annotate_seq1{
        input:
            target_fa = fasta_reads1.target_fasta,
            target_fai = fasta_reads1.target_fasta_fai,
            annotation_file = annotation_file,
            run_blast_from_table = script_run_blast_from_table,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_annotate_sequences
    }

    call AnnotateSequences as annotate_seq2{
        input:
            target_fa = fasta_reads2.target_fasta,
            target_fai = fasta_reads2.target_fasta_fai,
            annotation_file = annotation_file,
            run_blast_from_table = script_run_blast_from_table,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_annotate_sequences
    }

    call CutSequencesToRegions as cut_to_region1 {
        input:
            target_fa = fasta_reads1.target_fasta,
            target_fai = fasta_reads1.target_fasta_fai,
            annotations = annotate_seq1.annotated,
            extract_spanned_regions = script_extract_spanned_regions,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_cut_sequences_to_regions
    }

    call CutSequencesToRegions as cut_to_region2 {
        input:
            target_fa = fasta_reads2.target_fasta,
            target_fai = fasta_reads2.target_fasta_fai,
            annotations = annotate_seq2.annotated,
            extract_spanned_regions = script_extract_spanned_regions,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_cut_sequences_to_regions
    }

    call CombineSequences {
        input:
            prefix = prefix,
            sequences = [cut_to_region1.cut_sequences, cut_to_region2.cut_sequences],
            hap_labels = ["hap1", "hap2"],
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_combine_sequences
    }

    call AddPrefixIfMissing{
        input:
            fasta = CombineSequences.combined_fasta,
            prefix_string = prefix,
            docker_image = docker_image,
            runtime_attr_override = runtime_attr_add_prefix_to_read
    }

    output {
        File combined_fasta = CombineSequences.combined_fasta
        File combined_fai = CombineSequences.combined_fai
    }
}

# -------------------
# Tasks
# -------------------

task ExtractRegion {
    input {
        String bam
        String bai
        String chrom
        Int start
        Int end
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String output_prefix = "region_" + chrom + "_" + start + "_" + end
    String prefix = basename(bam,".bam")
    command <<<
        set -euo pipefail

        gsutil cp ~{bam} ./
        gsutil cp ~{bai} ./
        samtools view -b ~{prefix}.bam ~{chrom}:~{start}-~{end} > ~{prefix}.~{output_prefix}.bam
        samtools index ~{prefix}.~{output_prefix}.bam
    >>>

    output {
        File region_bam = "~{prefix}.~{output_prefix}.bam"
        File region_bai = "~{prefix}.~{output_prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(bam, "GiB")),
        disk_gb: 20 + ceil(size(bam, "GiB")),
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


task GetUniqueReads {

  input {
    File bam               # Already restricted to target region
    String chrom
    Int start              # 1-based inclusive
    Int end                # 1-based inclusive
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(bam , ".bam")
  command <<<
    set -euo pipefail

    REGION_START=~{start}
    REGION_END=~{end}

    # Extract alignments (exclude unmapped)
    samtools view -F 4 ~{bam} > reads.sam

    # Always initialize alignment table with header
    echo -e "read_id\treference_region\tread_aligned_region\taligned_read_length\taligned_reference_length\tread_total_length" > alignment_table.tsv

    # Count unique read IDs
    cut -f1 reads.sam | sort | uniq > unique_reads.txt
    COUNT=$(wc -l < unique_reads.txt)

    if [ "$COUNT" -eq 0 ]; then
            echo "No reads found in region" > result.tsv
            # alignment_table.tsv already exists (header only)
            mv result.tsv ~{prefix}.reads
            mv alignment_table.tsv ~{prefix}.alignment_table.tsv
            exit 0
    fi

    ############################################
    # Build alignment table (for >=1 reads)
    ############################################

    awk -v chrom="~{chrom}" '
    function cigar_parse(cigar,        len,type,aligned_read,aligned_ref) {
            aligned_read=0
            aligned_ref=0
            while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {
                    len=substr(cigar, RSTART, RLENGTH-1)
                    type=substr(cigar, RSTART+RLENGTH-1, 1)

                    if (type=="M" || type=="=" || type=="X") {
                            aligned_read+=len
                            aligned_ref+=len
                    }
                    else if (type=="I") {
                            aligned_read+=len
                    }
                    else if (type=="D" || type=="N") {
                            aligned_ref+=len
                    }
                    cigar=substr(cigar, RSTART+RLENGTH)
            }
            return aligned_read "|" aligned_ref
    }

    {
            read_id=$1
            pos=$4
            cigar=$6
            read_len=length($10)

            split(cigar_parse(cigar), arr, "|")
            aligned_read=arr[1]
            aligned_ref=arr[2]

            ref_start=pos
            ref_end=pos + aligned_ref - 1

            read_aln_start=1
            read_aln_end=aligned_read

            print read_id "\t" chrom ":" ref_start "-" ref_end "\t" \
                        read_aln_start "-" read_aln_end "\t" \
                        aligned_read "\t" aligned_ref "\t" read_len >> "alignment_table.tsv"

            print read_id "\t" ref_start "\t" ref_end >> "coords.tmp"
    }
    ' reads.sam

    ############################################
    # If only one unique read
    ############################################

    if [ "$COUNT" -eq 1 ]; then
            cat unique_reads.txt > result.tsv
            mv result.tsv ~{prefix}.reads
            mv alignment_table.tsv ~{prefix}.alignment_table.tsv
            exit 0
    fi

    ############################################
    # Collapse multiple alignments per read
    ############################################

    awk '
    {
        id=$1; s=$2; e=$3;
        if (!(id in min) || s < min[id]) min[id]=s;
        if (!(id in max) || e > max[id]) max[id]=e;
    }
    END {
        for (id in min)
            print id "\t" min[id] "\t" max[id];
    }
    ' coords.tmp > merged_coords.tmp

    ############################################
    # Check single-read full coverage
    ############################################

    awk -v rstart="$REGION_START" -v rend="$REGION_END" '
    ($2 <= rstart && $3 >= rend) { print $1 > "result.tsv"; exit 0 }
    ' merged_coords.tmp

    if [ -f result.tsv ]; then
            mv result.tsv ~{prefix}.reads
            mv alignment_table.tsv ~{prefix}.alignment_table.tsv
            exit 0
    fi

    ############################################
    # Test two-read union
    ############################################

    FOUND=0
    while read id1 s1 e1; do
        while read id2 s2 e2; do
            if [ "$id1" != "$id2" ]; then
                MIN_START=$(( s1 < s2 ? s1 : s2 ))
                MAX_END=$(( e1 > e2 ? e1 : e2 ))

                if [ "$MIN_START" -le "$REGION_START" ] && [ "$MAX_END" -ge "$REGION_END" ]; then
                    echo -e "$id1\n$id2" > result.tsv
                    FOUND=1
                    break 2
                fi
            fi
        done < merged_coords.tmp
    done < merged_coords.tmp

    if [ "$FOUND" -eq 0 ]; then
            echo "No single or pair of reads fully cover region" > result.tsv
    fi

    mv result.tsv ~{prefix}.reads
    mv alignment_table.tsv ~{prefix}.alignment_table.tsv
  >>>


    output {
        File reads = "~{prefix}.reads"
        File alignment_table = "~{prefix}.alignment_table.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8 + ceil(size(bam, "GiB")),
        disk_gb: 10 + ceil(size(bam, "GiB")),
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

task ExtractFastaReads {
    input {
        File fasta
        File read_list
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String output_prefix = basename(fasta, ".fa")

    command <<<
        set -euo pipefail

        samtools faidx ~{fasta} $(cat ~{read_list}) > ~{output_prefix}_reads.fa
        samtools faidx ~{output_prefix}_reads.fa

    >>>

    output {
        File target_fasta = "~{output_prefix}_reads.fa"
        File target_fasta_fai = "~{output_prefix}_reads.fa.fai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8 + ceil(size(fasta, "GiB")),
        disk_gb: 10 + ceil(size(fasta, "GiB")),
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

task AnnotateSequences {
    input {
        File target_fa
        File target_fai
        File annotation_file
        File run_blast_from_table
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(target_fa, ".fa")
    command <<<
        set -euo pipefail
        wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
        tar zxvf ncbi-blast-2.17.0+-x64-linux.tar.gz
        bash ~{run_blast_from_table} ~{annotation_file} ~{target_fa} ~{prefix}.tsv

    >>>

    output {
        File annotated = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10,
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

task CutSequencesToRegions {
    input {
        File target_fa
        File target_fai
        File annotations
        File extract_spanned_regions
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(target_fa,".fa")
    command <<<
        set -euo pipefail

        python ~{extract_spanned_regions} ~{target_fa} ~{annotations}  ~{prefix}.target.fa

    >>>

    output {
        File cut_sequences = "~{prefix}.target.fa"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10,
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

task CombineSequences {
    input {
        Array[File] sequences
        Array[String] hap_labels
        String prefix
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=" " sequences} > ~{prefix}.fa
        samtools faidx ~{prefix}.fa
    >>>

    output {
        File combined_fasta = "~{prefix}.fa"
        File combined_fai = "~{prefix}.fa.fai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10,
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

task AddPrefixIfMissing {

    input {
        File fasta
        String prefix_string
        String docker_image
        RuntimeAttr? runtime_attr_override

    }

    String file_prefix = basename(fasta, '.fa')
    command <<<
        set -euo pipefail

        awk -v prefix="~{prefix_string}" '
        /^>/ {
            header = substr($0, 2)
            split(header, parts, " ")
            readname = parts[1]

            if (index(readname, prefix) == 0) {
                parts[1] = prefix "_" readname
            }

            printf(">%s", parts[1])

            for (i = 2; i <= length(parts); i++) {
                printf(" %s", parts[i])
            }
            printf("\n")
            next
        }
        { print }
        ' ~{fasta} > ~{file_prefix}.read_id_fix.fa
    >>>

    output {
        File updated_fasta = " ~{file_prefix}.read_id_fix.fa"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: 10,
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
