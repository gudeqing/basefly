version development
# 由tumor_normal.py生成草稿后修改完成
# refer https://github.com/Sentieon/sentieon-scripts/blob/master/example_pipelines/somatic/TNseq/

workflow pipeline {
    input {
        File ref
        Array[File] ref_idxes
        Array[File] known_dbsnp
        # known indels example [1000G_phase1.indels.b37.vcf.gz, Mills_and_1000G_gold_standard.indels.b37.vcf.gz]
        Array[File] known_indels
#        File known_mills
        File? pon
        File? germline_vcf
#        String tumor_sample
#        String normal_sample
        Array[File] snpeff_databse
        Int thread_number = 15
        String platform = "ILLUMINA"
        # tumor normal pair info, two-column txt file, first column is tumor sample name
        File pair_info
        Array[File] intervals = []
        Boolean skip_fastp = false
    }

    call getFastqInfo{}

    Array[String] sample_array = keys(getFastqInfo.fastq_info)
    scatter (each in sample_array) {
        String sample = each
        File read1 = getFastqInfo.fastq_info[each][0][0]
        File read2 = getFastqInfo.fastq_info[each][1][0]

        if (! skip_fastp) {
            call fastp {
                input:
                    read1 = read1,
                    read2 = read2,
                    sample_name = sample
            }
        }

        call bwa_mem {
            input: 
            readgroup = "@RG\\tID:~{sample}\\tSM:~{sample}\\tPL:~{platform}",
            t = thread_number,
            ref = ref,
            ref_idxes  = ref_idxes,
            read1 = select_first([fastp.out_read1_file, read1]),
            read2 = select_first([fastp.out_read2_file, read2]),
            ref2 = ref,
            out = "~{sample}.sorted.bam",
            t2 = thread_number
        }

        call get_metrics {
            input: 
            t = thread_number,
            ref = ref,
            ref_idxes  = ref_idxes,
            intervals = intervals,
            bam = bwa_mem.out, bam_bai = bwa_mem.out_bai,
            mq_metrics = "~{sample}.mq_metrics.txt",
            qd_metrics = "~{sample}.qd_metrics.txt",
            gc_summary = "~{sample}.gc_summary.txt",
            gc_metrics = "~{sample}.gc_metrics.txt",
            aln_metrics = "~{sample}.aln_metrics.txt",
            insert_metrics = "~{sample}.insert_metrics.txt"
        }

        call plotGCBias {
            input: 
            out = "~{sample}.GCBias.pdf",
            i = get_metrics.gc_metrics
        }

        call plotMeanQualityByCycle {
            input: 
            out = "~{sample}.MeanQualityByCycle.pdf",
            i = get_metrics.mq_metrics
        }

        call plotQualDistribution {
            input: 
            out = "~{sample}.QualDistribution.pdf",
            i = get_metrics.qd_metrics
        }

        call plotInsertSize {
            input: 
            out = "~{sample}.InsertSizeMetricAlgo.pdf",
            i = get_metrics.insert_metrics
        }

        call LocusCollector {
            input:
            t = thread_number,
            intervals = intervals,
            bam = bwa_mem.out, bam_bai = bwa_mem.out_bai,
            score = "~{sample}.score.txt"
        }

        call DeDup {
            input:
            t = thread_number,
            intervals = intervals,
            bam = bwa_mem.out, bam_bai = bwa_mem.out_bai,
            score = LocusCollector.score,
            score_idx = LocusCollector.score_idx,
            dedup_metrics = "~{sample}.dedup.metrics.txt",
            deduped_bam = "~{sample}.deduped.bam"
        }

        call CoverageMetrics {
            input:
            t = thread_number,
            intervals = intervals,
            ref = ref,
            ref_idxes  = ref_idxes,
            bam = DeDup.deduped_bam, bam_bai = DeDup.deduped_bam_bai,
            coverage_metrics = "~{sample}.cov.metrics.txt"
        }

        call realign {
            input:
            t = thread_number,
            intervals = intervals,
            ref = ref,
            ref_idxes  = ref_idxes,
            bam = DeDup.deduped_bam, bam_bai = DeDup.deduped_bam_bai,
            database = known_indels,
            realigned_bam = "~{sample}.realigned.bam"
        }

        call recalibration {
            input: 
            ref = ref,
            intervals = intervals,
            ref_idxes  = ref_idxes,
            t = thread_number,
            bam = realign.realigned_bam, bam_bai = realign.realigned_bam_bai,
            database = flatten([known_dbsnp, known_indels]),
            recal_data = "~{sample}.recal_data.table"
        }

    }

    Map[String, File] bam_dict = as_map(zip(sample_array, realign.realigned_bam))
    Map[String, File] bai_dict = as_map(zip(sample_array, realign.realigned_bam_bai))
    Map[String, File] recal_dict = as_map(zip(sample_array, recalibration.recal_data))
    Array[Array[String]] pair_array = read_tsv(pair_info)

    scatter (each in pair_array) {
        String tumor_sample = each[0]
        String normal_sample = each[1]

        call TNhaplotyper2 {
            input:
            ref = ref,
            ref_idxes  = ref_idxes,
            t = thread_number,
            intervals = intervals,
            bams = [bam_dict[tumor_sample], bam_dict[normal_sample]],
            bam_bais = [bai_dict[tumor_sample], bai_dict[normal_sample]],
            recal_datas = [recal_dict[tumor_sample], recal_dict[normal_sample]],
            tumor_sample = "~{tumor_sample}",
            normal_sample = "~{normal_sample}",
            germline_vcf = germline_vcf,
            pon = pon,
            out_vcf = "~{tumor_sample}.TNhaplotyper2.vcf.gz",
            orientation_sample = "~{tumor_sample}",
            orientation_data = "~{tumor_sample}.orientation.data",
            contamination_tumor = "~{tumor_sample}",
            contamination_normal = "~{normal_sample}",
            germline_vcf2 = germline_vcf,
            tumor_segments = "~{tumor_sample}.contamination.segments",
            contamination_data = "~{tumor_sample}.contamination.data"
        }

        call TNfilter {
            input:
            ref = ref,
            ref_idxes  = ref_idxes,
            tumor_sample = "~{tumor_sample}",
            normal_sample = "~{normal_sample}",
            tmp_vcf = TNhaplotyper2.out_vcf,
            contamination = TNhaplotyper2.contamination_data,
            tumor_segments = TNhaplotyper2.tumor_segments,
            orientation_data = TNhaplotyper2.orientation_data,
            out_vcf = "~{tumor_sample}.final.vcf.gz"
        }

        call snpEff {
            input:
            database_files = snpeff_databse,
            in_vcf = TNfilter.out_vcf,
            out_vcf = "~{tumor_sample}.final.annot.vcf"
        }
    }

    meta {
        name: "TN_pipeline"
        desc: "typical bioinformatics pipeline using sentieon TNSeq and snpEff"
        author: "unknown"
        source: "source URL for the tool"
    }

    output{
        Array[File] bwa_mem_out = bwa_mem.out
        Array[File] get_metrics_mq_metrics = get_metrics.mq_metrics
        Array[File] get_metrics_qd_metrics = get_metrics.qd_metrics
        Array[File] get_metrics_gc_summary = get_metrics.gc_summary
        Array[File] get_metrics_gc_metrics = get_metrics.gc_metrics
        Array[File] get_metrics_aln_metrics = get_metrics.aln_metrics
        Array[File] get_metrics_insert_metrics = get_metrics.insert_metrics
        Array[File] plotGCBias_out = plotGCBias.out
        Array[File] plotMeanQualityByCycle_out = plotMeanQualityByCycle.out
        Array[File] plotQualDistribution_out = plotQualDistribution.out
        Array[File] plotInsertSize_out = plotInsertSize.out
        Array[File] LocusCollector_score = LocusCollector.score
        Array[File] DeDup_dedup_metrics = DeDup.dedup_metrics
        Array[File] DeDup_deduped_bam = DeDup.deduped_bam
        Array[File] DeDup_deduped_bam_bai = DeDup.deduped_bam_bai
        Array[File] CoverageMetrics_coverage_metrics = CoverageMetrics.coverage_metrics
        Array[File] realign_realigned_bam = realign.realigned_bam
        Array[File] realign_realigned_bam_bai = realign.realigned_bam_bai
        Array[File] recalibration_recal_data = recalibration.recal_data
        Array[File] TNhaplotyper2_out_vcf = TNhaplotyper2.out_vcf
        Array[File] TNhaplotyper2_orientation_data = TNhaplotyper2.orientation_data
        Array[File] TNhaplotyper2_tumor_segments = TNhaplotyper2.tumor_segments
        Array[File] TNhaplotyper2_contamination_data = TNhaplotyper2.contamination_data
        Array[File] TNfilter_out_vcf = TNfilter.out_vcf
        Array[File] snpEff_out_vcf = snpEff.out_vcf
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
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/getfastqinfo:1.0"
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
        String? other_parameters
        Int threads = 4
        File read1
        File? read2
        String sample_name
        String? adapter_r1
        String? adapter_r2
        # for runtime
        String docker = "gudeqing/fastp:0.21.0"
        String memory = "5 GiB"
        Int cpu = 4
        String disks = "5 GiB"
        Int time_minutes = 10080
    }

    command <<<
        set -e
        fastp \
        ~{other_parameters} \
        ~{"--thread " + threads} \
        ~{"-i " + read1} \
        ~{"-I " + read2} \
        ~{"--adapter_sequence " + adapter_r1} \
        ~{"--adapter_sequence_r2 " + adapter_r2} \
        ~{"-o " + sample_name + ".clean.R1.fastq.gz"} \
        ~{if defined(read2) then "-O " + sample_name + ".clean.R2.fastq.gz" else ""} \
        ~{"-h " + sample_name + ".report.html"} \
        ~{"-j " + sample_name + ".report.json"}
    >>>

    output {
        File out_read1_file = sample_name + ".clean.R1.fastq.gz"
        File? out_read2_file = sample_name + ".clean.R2.fastq.gz"
        File html_report_file = sample_name + ".report.html"
        File json_report_file = sample_name + ".report.json"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/fastp:0.21.0"
        memory: memory
        cpu: cpu
        disks: disks
        time_minutes: time_minutes
    }

    meta {
        name: "fastp"
        desc: "A tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance. For more detail please refer to https://github.com/OpenGene/fastp"
        logo: "see fastp.png"
        version: "0.21.0"
        basecmd: "./fastp"
    }

    parameter_meta {
        other_parameters: {desc: "其他参数，你可以通过该参数输入一个或多个任何其他当前软件支持的参数，例如'-i x -j y'", level: "optional", type: "str", range: "", default: ""}
        threads: {desc: "Number of threads to use", level: "required", type: "int", range: "", default: "2"}
        read1: {desc: "read1 fastq file", level: "required", type: "infile", range: "", default: "sample.R1.fastq.gz"}
        read2: {desc: "read2 fastq file", level: "optional", type: "infile", range: "", default: "sample.R2.fastq.gz"}
        sample_name: {desc: "sample name, will be used in output files' name", level: "required", type: "str", range: "", default: "sample_name"}
        adapter_r1: {desc: "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])", level: "optional", type: "str", range: "", default: "auto"}
        adapter_r2: {desc: "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=auto])", level: "optional", type: "str", range: "", default: "auto"}
    }

}

task bwa_mem{
    input {
        String readgroup
        Int t = 16
        Int k = 10000000
        File ref
        Array[File] ref_idxes
        File read1
        File read2
        File ref2
        String out
        Int t2 = 16
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon bwa mem -M \
        -R  '~{readgroup}' \
        ~{"-t " + t} \
        ~{"-K " + k} \
        ~{ref} \
        ~{read1} \
        ~{read2} \
         | sentieon util sort \
        ~{"-r " + ref2} \
        ~{"-o " + out} \
        ~{"-t " + t2} \
        --sam2bam -i - 
    >>>

    output {
        File out = "~{out}"
        File out_bai = "~{out}.bai"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "bwa_mem"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        readgroup: {prefix: "-R ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "read group info"}
        t: {prefix: "-t ", type: "int", level: "required", default: "16", range: "None", array: "False", desc: "number of threads to use in computation, set to number of cores in the server"}
        k: {prefix: "-K ", type: "int", level: "required", default: "10000000", range: "None", array: "False", desc: "This is description of the argument."}
        ref: {prefix: "", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        read1: {prefix: "", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "read1 fastq file"}
        read2: {prefix: "", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "read2 fastq file"}
        ref2: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        out: {prefix: "-o ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output bam file"}
        t2: {prefix: "-t ", type: "int", level: "required", default: "16", range: "None", array: "False", desc: "number of threads to use"}
    }

}

task get_metrics{
    input {
        Int t = 16
        File ref
        Array[File] ref_idxes
        Array[File] intervals
        File bam
        File bam_bai
        String mq_metrics
        String qd_metrics
        String gc_summary
        String gc_metrics
        String aln_metrics
        String insert_metrics
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{prefix("--interval ", intervals)} \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        ~{"--algo MeanQualityByCycle " + mq_metrics} \
        ~{"--algo QualDistribution " + qd_metrics} \
        ~{"--algo GCBias --summary " + gc_summary} \
        ~{gc_metrics} \
        ~{"--algo AlignmentStat " + aln_metrics} \
        ~{"--algo InsertSizeMetricAlgo " + insert_metrics} 
    >>>

    output {
        File mq_metrics = "~{mq_metrics}"
        File qd_metrics = "~{qd_metrics}"
        File gc_summary = "~{gc_summary}"
        File gc_metrics = "~{gc_metrics}"
        File aln_metrics = "~{aln_metrics}"
        File insert_metrics = "~{insert_metrics}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "get_metrics"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        t: {prefix: "-t ", type: "int", level: "required", default: "16", range: "None", array: "False", desc: "number of threads to use in computation, set to number of cores in the server"}
        ref: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        bam: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input bam file"}
        mq_metrics: {prefix: "--algo MeanQualityByCycle ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "metric file of MeanQualityByCycle"}
        qd_metrics: {prefix: "--algo QualDistribution ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "metric file of QualDistribution"}
        gc_summary: {prefix: "--algo GCBias --summary ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "summary file of GCBias"}
        gc_metrics: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "metrics file of GCBias"}
        aln_metrics: {prefix: "--algo AlignmentStat ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "aln_metrics file of AlignmentStat"}
        insert_metrics: {prefix: "--algo InsertSizeMetricAlgo ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "insert_metrics file of InsertSizeMetricAlgo"}
    }

}

task plotGCBias{
    input {
        String method = "GCBias"
        String out
        File i
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon plot \
        ~{method} \
        ~{"-o " + out} \
        ~{i} 
    >>>

    output {
        File out = "~{out}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "plotGCBias"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        method: {prefix: "", type: "str", level: "required", default: "GCBias", range: "None", array: "False", desc: "method of plot"}
        out: {prefix: "-o ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "plot file"}
        i: {prefix: "", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input metrics file for plot"}
    }

}

task plotMeanQualityByCycle{
    input {
        String method = "MeanQualityByCycle"
        String out
        File i
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon plot \
        ~{method} \
        ~{"-o " + out} \
        ~{i} 
    >>>

    output {
        File out = "~{out}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "plotMeanQualityByCycle"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        method: {prefix: "", type: "str", level: "required", default: "MeanQualityByCycle", range: "None", array: "False", desc: "method of plot"}
        out: {prefix: "-o ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "plot file"}
        i: {prefix: "", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input metrics file for plot"}
    }

}

task plotQualDistribution{
    input {
        String method = "QualDistribution"
        String out
        File i
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon plot \
        ~{method} \
        ~{"-o " + out} \
        ~{i} 
    >>>

    output {
        File out = "~{out}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "plotQualDistribution"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        method: {prefix: "", type: "str", level: "required", default: "QualDistribution", range: "None", array: "False", desc: "method of plot"}
        out: {prefix: "-o ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "plot file"}
        i: {prefix: "", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input metrics file for plot"}
    }

}

task plotInsertSize{
    input {
        String method = "InsertSizeMetricAlgo"
        String out
        File i
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon plot \
        ~{method} \
        ~{"-o " + out} \
        ~{i} 
    >>>

    output {
        File out = "~{out}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "plotInsertSize"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        method: {prefix: "", type: "str", level: "required", default: "InsertSizeMetricAlgo", range: "None", array: "False", desc: "method of plot"}
        out: {prefix: "-o ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "plot file"}
        i: {prefix: "", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input metrics file for plot"}
    }

}

task LocusCollector{
    input {
        Int t = 16
        File bam
        File bam_bai
        Array[File] intervals
        String score
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{prefix("--interval ", intervals)} \
        ~{"-t " + t} \
        ~{"-i " + bam} \
        ~{"--algo LocusCollector --fun score_info " + score} 
    >>>

    output {
        File score = "~{score}"
        File score_idx = "~{score}.idx"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "LocusCollector"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        t: {prefix: "-t ", type: "int", level: "required", default: "16", range: "None", array: "False", desc: "number of threads to use in computation, set to number of cores in the server"}
        bam: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input bam file"}
        score: {prefix: "--algo LocusCollector --fun score_info ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output score file"}
    }

}

task DeDup{
    input {
        Int t = 16
        File bam
        File bam_bai
        Array[File] intervals
        File score
        File score_idx
        String dedup_metrics
        String deduped_bam
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{prefix("--interval ", intervals)} \
        ~{"-t " + t} \
        ~{"-i " + bam} \
        --algo Dedup \
        ~{"--score_info " + score} \
        ~{"--metrics " + dedup_metrics} \
        ~{deduped_bam} 
    >>>

    output {
        File dedup_metrics = "~{dedup_metrics}"
        File deduped_bam = "~{deduped_bam}"
        File deduped_bam_bai = "~{deduped_bam}.bai"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "DeDup"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        t: {prefix: "-t ", type: "int", level: "required", default: "16", range: "None", array: "False", desc: "number of threads to use in computation, set to number of cores in the server"}
        bam: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input bam file"}
        score: {prefix: "--score_info ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "score info file"}
        dedup_metrics: {prefix: "--metrics ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output metrics info file"}
        deduped_bam: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output metrics info file"}
    }

}

task CoverageMetrics{
    input {
        Int t = 16
        File ref
        Array[File] ref_idxes
        Array[File] intervals
        File bam
        File bam_bai
        String coverage_metrics
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{prefix("--interval ", intervals)} \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        ~{"--algo CoverageMetrics " + coverage_metrics} 
    >>>

    output {
        File coverage_metrics = "~{coverage_metrics}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "CoverageMetrics"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        t: {prefix: "-t ", type: "int", level: "required", default: "16", range: "None", array: "False", desc: "number of threads to use in computation, set to number of cores in the server"}
        ref: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        bam: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input bam file"}
        coverage_metrics: {prefix: "--algo CoverageMetrics ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output coverage metrics file"}
    }

}

task realign{
    input {
        Int t = 16
        File ref
        Array[File] ref_idxes
        Array[File] intervals
        File bam
        File bam_bai
        Array[File] database
        String realigned_bam
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{prefix("--interval ", intervals)} \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        --algo Realigner \
        ~{prefix("-k ", database)} \
        ~{realigned_bam} 
    >>>

    output {
        File realigned_bam = "~{realigned_bam}"
        File realigned_bam_bai = "~{realigned_bam}.bai"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "realign"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        t: {prefix: "-t ", type: "int", level: "required", default: "16", range: "None", array: "False", desc: "number of threads to use in computation, set to number of cores in the server"}
        ref: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        bam: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input bam file"}
        database: {prefix: "-k ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "known indel vcf file"}
        realigned_bam: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output realigned bam file"}
    }

}

task recalibration{
    input {
        Int t = 16
        File ref
        Array[File] ref_idxes
        Array[File] intervals
        File bam
        File bam_bai
        Array[File] database
        String recal_data
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{prefix("--interval ", intervals)} \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        --algo QualCal \
        ~{prefix("-k  ", database)} \
        ~{recal_data} 
    >>>

    output {
        File recal_data = "~{recal_data}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "recalibration"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        t: {prefix: "-t ", type: "int", level: "required", default: "16", range: "None", array: "False", desc: "number of threads to use in computation, set to number of cores in the server"}
        ref: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        bam: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input bam file"}
        database: {prefix: "-k ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "known indel vcf file"}
        recal_data: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output recal_data.table"}
    }

}

task TNhaplotyper2{
    input {
        Int t = 16
        File ref
        Array[File] ref_idxes
        Array[File] intervals
        Array[File] bams
        Array[File] bam_bais
        Array[File] recal_datas
        String tumor_sample = "tumor"
        String? normal_sample
        File? germline_vcf
        File? pon
        String out_vcf
        String? orientation_sample
        String? orientation_data
        String? contamination_tumor = "tumor"
        String? contamination_normal
        File? germline_vcf2
        String? tumor_segments
        String? contamination_data
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{prefix("--interval ", intervals)} \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{prefix("-i  ", bams)} \
        ~{prefix("-q  ", recal_datas)} \
        --algo TNhaplotyper2 \
        ~{"--tumor_sample " + tumor_sample} \
        ~{"--normal_sample " + normal_sample} \
        ~{"--germline_vcf " + germline_vcf} \
        ~{"--pon " + pon} \
        ~{out_vcf} \
        ~{"--algo OrientationBias --tumor_sample " + orientation_sample} \
        ~{orientation_data} \
        ~{"--algo ContaminationModel --tumor_sample " + contamination_tumor} \
        ~{"--normal_sample " + contamination_normal} \
        ~{"--vcf " + germline_vcf2} \
        ~{"--tumor_segments " + tumor_segments} \
        ~{contamination_data} 
    >>>

    output {
        File out_vcf = "~{out_vcf}"
        File orientation_data = "~{orientation_data}"
        File tumor_segments = "~{tumor_segments}"
        File contamination_data = "~{contamination_data}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "TNhaplotyper2"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        t: {prefix: "-t ", type: "int", level: "required", default: "16", range: "None", array: "False", desc: "number of threads to use in computation, set to number of cores in the server"}
        ref: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        bams: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reccaled tumor and normal bam list"}
        recal_datas: {prefix: "-q ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "tumor and normal recal data list"}
        tumor_sample: {prefix: "--tumor_sample ", type: "str", level: "required", default: "tumor", range: "None", array: "False", desc: "tumor sample name"}
        normal_sample: {prefix: "--normal_sample ", type: "str", level: "optional", default: "None", range: "None", array: "False", desc: "normal sample name"}
        germline_vcf: {prefix: "--germline_vcf ", type: "infile", level: "optional", default: "None", range: "None", array: "False", desc: "the location of the population germline resource"}
        pon: {prefix: "--pon ", type: "infile", level: "optional", default: "None", range: "None", array: "False", desc: "the location and name of panel of normal VCF file"}
        out_vcf: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output vcf file of TNhaplotyper2, this will be used later for filtering"}
        orientation_sample: {prefix: "--algo OrientationBias --tumor_sample ", type: "str", level: "optional", default: "None", range: "None", array: "False", desc: "tumor sample name"}
        orientation_data: {prefix: "", type: "str", level: "optional", default: "None", range: "None", array: "False", desc: "output orientation bias result file"}
        contamination_tumor: {prefix: "--algo ContaminationModel --tumor_sample ", type: "str", level: "optional", default: "tumor", range: "None", array: "False", desc: "tumor sample name"}
        contamination_normal: {prefix: "--normal_sample ", type: "str", level: "optional", default: "None", range: "None", array: "False", desc: "normal sample name"}
        germline_vcf2: {prefix: "--vcf ", type: "infile", level: "optional", default: "None", range: "None", array: "False", desc: "the location of the population germline resource"}
        tumor_segments: {prefix: "--tumor_segments ", type: "str", level: "optional", default: "None", range: "None", array: "False", desc: "output file name of the file containing the tumor segments information produced by ContaminationModel"}
        contamination_data: {prefix: "", type: "str", level: "optional", default: "None", range: "None", array: "False", desc: "output file containing the contamination information produced by ContaminationModel"}
    }

}

task TNfilter{
    input {
        File ref
        Array[File] ref_idxes
        String tumor_sample = "tumor"
        String? normal_sample
        File tmp_vcf
        File? contamination
        File? tumor_segments
        File? orientation_data
        String out_vcf
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{"-r " + ref} \
        --algo TNfilter \
        ~{"--tumor_sample " + tumor_sample} \
        ~{"--normal_sample " + normal_sample} \
        ~{"-v " + tmp_vcf} \
        ~{"--contamination " + contamination} \
        ~{"--tumor_segments " + tumor_segments} \
        ~{"--orientation_priors " + orientation_data} \
        ~{out_vcf} 
    >>>

    output {
        File out_vcf = "~{out_vcf}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon-joint-call:latest"
    }

    meta {
        name: "TNfilter"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        ref: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        tumor_sample: {prefix: "--tumor_sample ", type: "str", level: "required", default: "tumor", range: "None", array: "False", desc: "tumor sample name"}
        normal_sample: {prefix: "--normal_sample ", type: "str", level: "optional", default: "None", range: "None", array: "False", desc: "normal sample name"}
        tmp_vcf: {prefix: "-v ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "vcf file from TNhaplotyper2"}
        contamination: {prefix: "--contamination ", type: "infile", level: "optional", default: "None", range: "None", array: "False", desc: "file containing the contamination information produced by ContaminationModel"}
        tumor_segments: {prefix: "--tumor_segments ", type: "infile", level: "optional", default: "None", range: "None", array: "False", desc: "file containing the tumor segments information produced by ContaminationModel"}
        orientation_data: {prefix: "--orientation_priors ", type: "infile", level: "optional", default: "None", range: "None", array: "False", desc: "file containing the orientation bias information produced by OrientationBias"}
        out_vcf: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "final output vcf"}
    }

}

task snpEff{
    input {
#        String genome_version = "hg19"
        Array[File] database_files
        Boolean cancer = true
        File? cancerSamples
        Boolean canon = false
        String? other_args
        File in_vcf
        String out_vcf
        # for runtime
        String docker = "?"
    }
    String database_dir = sub(database_files[0], basename(database_files[0]), "")
    String genome_version = basename(database_dir)
    String database_parent_dir = sub(database_dir, basename(database_dir), "")
    command <<<
        set -e 
        java -Xmx9g snpEff.jar ann \
        ~{genome_version} \
        -dataDir  ~{database_parent_dir} \
        ~{if cancer then "-cancer  " else ""} \
        ~{"-cancerSamples " + cancerSamples} \
        ~{if canon then "-canon  " else ""} \
        ~{other_args} \
        ~{in_vcf} \
        > \
        ~{out_vcf} 
    >>>

    output {
        File out_vcf = "~{out_vcf}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/snpeff:5.0e"
    }

    meta {
        name: "snpEff"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
#        genome_version: {prefix: "", type: "str", level: "required", default: "hg19", range: "None", array: "False", desc: "human genome version"}
        database_files: {prefix: "-dataDir ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Override data_dir parameter from config file"}
        cancer: {prefix: "-cancer ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Perform cancer comparisons (Somatic vs Germline)"}
        cancerSamples: {prefix: "-cancerSamples ", type: "infile", level: "optional", default: "None", range: "None", array: "False", desc: "Two column TXT file defining 'oringinal derived' samples. If '-cancer' used and the file is missing, then the last sample will be assumed as tumor sample."}
        canon: {prefix: "-canon ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Only use canonical transcripts"}
        other_args: {prefix: "", type: "str", level: "optional", default: "None", range: "None", array: "False", desc: "other arguments that you want to input for the program, such as '-motif'"}
        in_vcf: {prefix: "", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input variant file"}
        out_vcf: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output annotated file"}
    }

}

