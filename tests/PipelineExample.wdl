version development

workflow pipeline {
    input {
        Array[File] read1
        Array[File] read2
        Array[String] names
    }

    Array[Pair[File, File]] reads = zip(read1, read2)
    Array[Pair[String, Pair[File, File]]] init_array = zip(names, reads)

    scatter (each in init_array) { 
        call fastp {
            input: 
            read1 = each.right.left,
            read2 = each.right.right,
            out1 = "~{each.left}.clean.R1.fq",
            out2 = "~{each.left}.clean.R2.fq"
        }

        call salmon {
            input: 
            read1 = fastp.out1,
            read2 = fastp.out2,
            outDir = each.left
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
    }

    output{
        Array[File] out1 = fastp.out1
        Array[File] out2 = fastp.out2
        Array[File] transcript = salmon.transcript
#        Array[Directory] outDir = salmon.outDir
        File TPM = MergeTranscriptTPM.result
        File Count = MergeTranscriptCount.result
    }

}

task fastp{
    input {
        File read1
        File read2
        String out1
        String out2
        # for runtime
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
        docker: docker
    }

    meta {
        name: "fastp"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
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
        docker: docker
    }

    meta {
        name: "salmon"
        desc: "transcript expression quantification"
        author: "unknown"
        source: "source URL for the tool"
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
        docker: docker
    }

    meta {
        name: "MergeTranscriptTPM"
        desc: "Merge multiple quantification results into a single file"
        author: "unknown"
        source: "source URL for the tool"
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
        docker: docker
    }

    meta {
        name: "MergeTranscriptCount"
        desc: "Merge multiple quantification results into a single file"
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        quants: {prefix: "--quants ", type: "indir", level: "required", default: "None", range: "None", array: "True", desc: "salmon quant dir list"}
        names: {prefix: "--names ", type: "str", level: "optional", default: "None", range: "None", array: "True", desc: "This is description of the argument."}
        out: {prefix: "--output ", type: "str", level: "required", default: "merged.NumReads.txt", range: "None", array: "False", desc: "This is description of the argument."}
    }

}

