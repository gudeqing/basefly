version development

workflow pipeline {
    input {
        Directory index_dir
    }

    call getFastqInfo{}
    scatter (each in keys(getFastqInfo.fastq_info)) { 
        String sample = each
        File read1 = getFastqInfo.fastq_info[each][0][0]
        File read2 = getFastqInfo.fastq_info[each][1][0]
        call fastp {
            input: 
            read1 = read1,
            read2 = read2,
            out1 = "~{sample}.clean.R1.fq",
            out2 = "~{sample}.clean.R2.fq"
        }

        call salmon {
            input: 
            indexDir = index_dir,
            read1 = fastp.out1,
            read2 = fastp.out2,
            outDir = sample
        }

    }

    call MergeTranscriptTPM {
        input: 
        quants = salmon.outDir
    }

    call MergeTranscriptCount {
        input: 
        quants = salmon.outDir
    }

    meta {
        name: "PipelineExample"
        desc: "This is a simple pipeline for fast gene/transcript quantification. workflow = [fastq -> Fastp -> Salmon]"
        author: "unknown"
        source: "source URL for the tool"
        version: "unknown"
    }

    output{
        Array[File] fastp_out1 = fastp.out1
        Array[File] fastp_out2 = fastp.out2
        Array[File] salmon_transcript = salmon.transcript
        Array[Directory] salmon_outDir = salmon.outDir
        File MergeTranscriptTPM_result = MergeTranscriptTPM.result
        File MergeTranscriptCount_result = MergeTranscriptCount.result
    }

}



task getFastqInfo{
    input {
        Array[Directory]? fastq_dirs
        Array[File]? fastq_files
        String r1_name = '(.*).read1.fastq.gz'
        String r2_name = '(.*).read2.fastq.gz'
        String docker = 'gudeqing/getfastqinfo:1.0'
    }

    command <<<
        set -e
        python /get_fastq_info.py \
            ~{if defined(fastq_dirs) then "-fastq_dirs " else ""}~{sep=" " fastq_dirs} \
            ~{if defined(fastq_files) then "-fastq_files " else ""}~{sep=" " fastq_files} \
            -r1_name '~{r1_name}' \
            -r2_name '~{r2_name}' \
            -out fastq.info.json
    >>>

    output {
        Map[String, Array[Array[File]]] fastq_info = read_json("fastq.info.json")
        File fastq_info_json = "fastq.info.json"
    }

    runtime {
        docker: docker
    }

    parameter_meta {
        fastq_dirs: {desc: "directory list, target fastq files should be in these directories. All target files in 'fastq_files' or 'fastq_dirs' will be used", level: "optional", type: "indir", range: "", default: ""}
        fastq_files: {desc: "target fastq file list. 'fastq_files' or 'fastq_dirs' must be provided.", level: "optional", type: "infile", range: "", default: ""}
        r1_name: {desc: "python regExp that describes the full name of read1 fastq file name. It requires at least one pair of small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'", level: "required", type: "str", range: "", default: ""}
        r2_name: {desc: "python regExp that describes the full name of read2 fastq file name. It requires at least one pair of small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'", level: "required", type: "str", range: "", default: ""}
    }
}
    
task fastp{
    input {
        File read1
        File read2
        String out1
        String out2
        # for runtime
        String memory = "1000"
        Int cpu = 2
        String max_memory = "0"
        String max_cpu = "0"
        String docker = "gudeqing/fastp:0.21.0"
    }

    command <<<
        set -e 
        fastp \
        ~{"-i " + read1} \
        ~{"-I " + read2} \
        ~{"-o " + out1} \
        ~{"-O " + out2} 
    >>>

    output {
        File out1 = "~{out1}"
        File out2 = "~{out2}"
    }

    runtime {
        memory: memory
        cpu: cpu
        max_memory: max_memory
        max_cpu: max_cpu
        docker: docker
    }

    meta {
        name: "fastp"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
        version: "unknown"
    }

    parameter_meta {
        read1: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "read1 fastq file"}
        read2: {prefix: "-I ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "read2 fastq file"}
        out1: {prefix: "-o ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "clean read1 output fastq file"}
        out2: {prefix: "-O ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "clean read2 output fastq file"}
    }

}

task salmon{
    input {
        String libType = "A"
        Directory indexDir
        File read1
        File read2
        String outDir = "quant"
        Boolean gcBias = true
        # for runtime
        String memory = "2147483648"
        Int cpu = 2
        String max_memory = "0"
        String max_cpu = "0"
        String docker = "combinelab/salmon:latest"
    }

    command <<<
        set -e 
        salmon quant \
        ~{"--libType " + libType} \
        ~{"-i " + indexDir} \
        ~{"-1 " + read1} \
        ~{"-2 " + read2} \
        ~{"-o " + outDir} \
        ~{if gcBias then "--gcBias  " else ""} 
    >>>

    output {
        File transcript = "~{outDir}/quant.sf"
        Directory outDir = "~{outDir}"
    }

    runtime {
        memory: memory
        cpu: cpu
        max_memory: max_memory
        max_cpu: max_cpu
        docker: docker
    }

    meta {
        name: "salmon"
        desc: "transcript expression quantification"
        author: "unknown"
        source: "source URL for the tool"
        version: "unknown"
    }

    parameter_meta {
        libType: {prefix: "--libType ", type: "str", level: "required", default: "A", range: "None", array: "False", desc: "This is description of the argument."}
        indexDir: {prefix: "-i ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "transcript fasta index directory"}
        read1: {prefix: "-1 ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "read1 fastq file"}
        read2: {prefix: "-2 ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "read2 fastq file"}
        outDir: {prefix: "-o ", type: "str", level: "required", default: "quant", range: "None", array: "False", desc: "output directory"}
        gcBias: {prefix: "--gcBias ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "perform gc Bias correction"}
    }

}

task MergeTranscriptTPM{
    input {
        Array[Directory] quants
        Array[String]? names
        String out = "merged.TPM.txt"
        # for runtime
        String memory = "1000"
        Int cpu = 2
        String max_memory = "0"
        String max_cpu = "0"
        String docker = "combinelab/salmon:latest"
    }

    command <<<
        set -e 
        salmon quantmerge \
        ~{if defined(quants) then "--quants  " else ""}~{sep=" " quants} \
        ~{if defined(names) then "--names  " else ""}~{sep=" " names} \
        --column TPM \
        ~{"--output " + out} 
    >>>

    output {
        File result = "~{out}"
    }

    runtime {
        memory: memory
        cpu: cpu
        max_memory: max_memory
        max_cpu: max_cpu
        docker: docker
    }

    meta {
        name: "MergeTranscriptTPM"
        desc: "Merge multiple quantification results into a single file"
        author: "unknown"
        source: "source URL for the tool"
        version: "unknown"
    }

    parameter_meta {
        quants: {prefix: "--quants ", type: "indir", level: "required", default: "None", range: "None", array: "True", desc: "salmon quant dir list"}
        names: {prefix: "--names ", type: "str", level: "optional", default: "None", range: "None", array: "True", desc: "This is description of the argument."}
        out: {prefix: "--output ", type: "str", level: "required", default: "merged.TPM.txt", range: "None", array: "False", desc: "This is description of the argument."}
    }

}

task MergeTranscriptCount{
    input {
        Array[Directory] quants
        Array[String]? names
        String out = "merged.NumReads.txt"
        # for runtime
        String memory = "1000"
        Int cpu = 2
        String max_memory = "0"
        String max_cpu = "0"
        String docker = "combinelab/salmon:latest"
    }

    command <<<
        set -e 
        salmon quantmerge \
        ~{if defined(quants) then "--quants  " else ""}~{sep=" " quants} \
        ~{if defined(names) then "--names  " else ""}~{sep=" " names} \
        --column NumReads \
        ~{"--output " + out} 
    >>>

    output {
        File result = "~{out}"
    }

    runtime {
        memory: memory
        cpu: cpu
        max_memory: max_memory
        max_cpu: max_cpu
        docker: docker
    }

    meta {
        name: "MergeTranscriptCount"
        desc: "Merge multiple quantification results into a single file"
        author: "unknown"
        source: "source URL for the tool"
        version: "unknown"
    }

    parameter_meta {
        quants: {prefix: "--quants ", type: "indir", level: "required", default: "None", range: "None", array: "True", desc: "salmon quant dir list"}
        names: {prefix: "--names ", type: "str", level: "optional", default: "None", range: "None", array: "True", desc: "This is description of the argument."}
        out: {prefix: "--output ", type: "str", level: "required", default: "merged.NumReads.txt", range: "None", array: "False", desc: "This is description of the argument."}
    }

}

