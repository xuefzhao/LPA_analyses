version 1.0


workflow AlignAsm {
    input {
      String input_sample_name
      File input_asm_h1
      File input_asm_h2
      File ref
    }
    meta {
        workflow_description: "Creates callerset for a single sample"
    }

    call alignToRef as alignToRefH1 {
        input:
            asmIn = input_asm_h1,
            sample = input_sample_name,
            hap = "h1",
            ref = ref,
            threads = 4,
    }
    call compressAndIndex as compressAndIndexH1 {
        input:
            samIn = alignToRefH1.samOut,
            sample = input_sample_name,
            hap = "h1",
    }

    call alignToRef as alignToRefH2 {
        input:
            asmIn = input_asm_h2,
            sample = input_sample_name,
            hap = "h2",
            ref = ref,
            threads = 4,
    }
    call compressAndIndex as compressAndIndexH2 {
        input:
            samIn = alignToRefH2.samOut,
            sample = input_sample_name,
            hap = "h2",
    }

    output {
        File sample_h1_bam = compressAndIndexH1.bamOut
        File sample_h1_paf = compressAndIndexH1.pafOut
        File sample_h1_idx = compressAndIndexH1.indexOut
        File sample_h2_bam = compressAndIndexH2.bamOut
        File sample_h2_paf = compressAndIndexH2.pafOut
        File sample_h2_idx = compressAndIndexH2.indexOut
    }
}

 
task alignToRef {
  input {
    File asmIn
    String sample
    File ref
    String hap
    Int threads

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    minimap2 -x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5 --secondary=no -a -t ~{threads} --eqx -Y ~{ref} ~{asmIn} > ~{sample + "-asm_" + hap + ".minimap2.sam"}
  >>>


  output {
    File samOut = "~{sample}-asm_~{hap}.minimap2.sam"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          threads,
      mem_gb:             32,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  1,
      max_retries:        1,
      docker:             "biocontainers/minimap2:v2.15dfsg-1-deb_cv1"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
      cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
      memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
      disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
      preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
      docker:                 select_first([runtime_attr.docker,            default_attr.docker])
  }
}


task compressAndIndex {
  input {
    File samIn
    String sample
    String hap

    RuntimeAttr? runtime_attr_override
  }

  command <<<
    samtools view -b ~{samIn} | samtools sort -O bam -o ~{sample + "-asm_" + hap + ".minimap2.bam"} -
    samtools index ~{sample + "-asm_" + hap + ".minimap2.bam"}

    MM2_VERSION="2.24"

    /minimap2-${MM2_VERSION}_x64-linux/k8 \
    /minimap2-${MM2_VERSION}_x64-linux/paftools.js \
    sam2paf \
    -L \
    ~{samIn} \
    > ~{sample + "-asm_" + hap + ".minimap2.paf"}
  >>>


  output {
    File bamOut = "~{sample}-asm_~{hap}.minimap2.bam"
    File pafOut = "~{sample}-asm_~{hap}.minimap2.paf"
    File indexOut = "~{sample}-asm_~{hap}.minimap2.bam.bai"
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             4,
      disk_gb:            10,
      boot_disk_gb:       10,
      preemptible_tries:  1,
      max_retries:        1,
      docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  runtime {
      cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
      memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
      disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
      bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
      preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
      maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
      docker:                 select_first([runtime_attr.docker,            default_attr.docker])
  }
} 

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}