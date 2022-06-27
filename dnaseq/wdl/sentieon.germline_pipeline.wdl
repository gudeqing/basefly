version development
# The following steps compose the typical bioinformatics pipeline:
# 1. Map reads to reference: This step aligns the reads contained in the FASTQ files to map to a reference genome contained in the FASTA file. This step ensures that the data can be placed in context.
# 2. Calculate data metrics: This step produces a statistical summary of the data quality and the pipeline data analysis quality.
# 3. Remove or mark duplicates: This step detects reads indicative that the same DNA molecules were sequenced several times. These duplicates are not informative and should not be counted as additional evidence.
# 4. (optional) Indel realignment: This step performs a local realignment around indels. This step is necessary as reads mapped on the edges of indels often get mapped with mismatching bases that are mapping artifacts. However, when using haplotype based callers such as Haplotyper or DNAscope, this step is not required, as the variant caller performs a local reassembly that provides most of the accuracy increase that results from Indel Realignment.
# 5. Base quality score recalibration (BQSR): This step modifies the quality scores assigned to individual read bases of the sequence read data. This action removes experimental biases caused by the sequencing methodology.
# 6. Variant calling: This step identifies the sites where your data displays variation relative to the reference genome, and calculates genotypes for each sample at that site.
# 7. variant annotation with VEP
# 8. HLA-I typing with OptiType

workflow pipeline {
    input {
        File ref
        Array[File] ref_idxes
        Array[File] bwa_idxes
        Array[File] known_dbsnp
        Array[File] known_dbsnp_idx
        # known indels example [1000G_phase1.indels.b37.vcf.gz, Mills_and_1000G_gold_standard.indels.b37.vcf.gz]
        Array[File] known_indels
        Array[File] known_indels_idx

        # annotation database
        File vep_cache
        File vep_plugins_zip
        Array[String] vep_plugin_names = ['Frameshift', 'Wildtype']
        Int thread_number = 15
        String platform = "ILLUMINA"

        Array[File] intervals = []

        Boolean skip_fastp = false
        Boolean skip_vep = false
        Boolean skip_optiType = true
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

        if (!skip_optiType){
            call OptiType {
                input:
                reads = [select_first([fastp.out_read1_file, read1]), select_first([fastp.out_read2_file, read2])],
                prefix = "~{sample}",
                threads = thread_number - 1
            }
        }

        call bwa_mem {
            input:
            readgroup = "@RG\\tID:~{sample}\\tSM:~{sample}\\tPL:~{platform}",
            t = thread_number,
            ref = ref,
            ref_idxes = ref_idxes,
            bwa_idxes = bwa_idxes,
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
            bam = bwa_mem.out, bam_bai = bwa_mem.out_bai,
            score = "~{sample}.score.txt"
        }

        call DeDup {
            input:
            t = thread_number,
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

        call recalibration {
            input:
            ref = ref,
            ref_idxes = ref_idxes,
            t = thread_number,
            intervals = intervals,
            bam = DeDup.deduped_bam,
            bam_bai = DeDup.deduped_bam_bai,
            database = flatten([known_dbsnp, known_indels]),
            database_idx = flatten([known_dbsnp_idx, known_indels_idx]),
            recal_data = "~{sample}.recal_data.table"
        }

        call Haplotyper {
            input:
            ref = ref,
            ref_idxes = ref_idxes,
            intervals = intervals,
            bam = DeDup.deduped_bam,
            bam_bai = DeDup.deduped_bam_bai,
            recal_data = recalibration.recal_data,
            out_vcf = "~{sample}.haplotyper.gvcf.gz"
        }

        call GVCFtyper {
            input:
            ref = ref,
            ref_idxes = ref_idxes,
            in_gvcf = [Haplotyper.out_vcf],
            in_gvcf_idx = [Haplotyper.out_vcf_idx],
            known_dbsnp = known_dbsnp[0],
            known_dbsnp_idx = known_dbsnp_idx[0],
            out_vcf = "~{sample}.germline.vcf.gz"
        }

        if (! skip_vep) {
            call VEP {
                input:
                input_file = GVCFtyper.out_vcf,
                input_file_idx = GVCFtyper.out_vcf_idx,
                cache_targz = vep_cache,
                plugins_zip= vep_plugins_zip,
                plugin_names = vep_plugin_names,
                output_file = "~{sample}.somatic.vep.vcf.gz",
                stats_file = "~{sample}.somatic.vep.summary.html",
                fasta = ref,
                fasta_idx = ref_idxes,
                fork = thread_number
            }
        }
    }

    meta {
        name: "senteion_germline_pipeline"
        desc: "typical bioinformatics pipeline using sentieon haplotype and vep"
        author: "unknown"
        source: "source URL for the tool"
    }

    output{
        File fastq_info_json = getFastqInfo.fastq_info_json
#        Array[File] bwa_mem_out = bwa_mem.out
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
        Array[File] recalibration_recal_data = recalibration.recal_data

        Array[File] Haplotyper_gvcf = Haplotyper.out_vcf
        Array[File] Haplotyper_gvcf_idx = Haplotyper.out_vcf_idx

        Array[File] GVCFtyper_vcf = GVCFtyper.out_vcf
        Array[File] GVCFtyper_vcf_idx = GVCFtyper.out_vcf_idx

        Array[File?] vep_out_vcf = VEP.out_vcf
        Array[File?] vep_out_vcf_idx = VEP.out_vcf_idx
        Array[File?] vep_out_vcf_stats = VEP.stats_file

        Array[File?] HLA_ABC = OptiType.result_tsv
        Array[File?] HLA_ABC_coverage = OptiType.result_pdf
    }

}

task getFastqInfo{
    input {
        Array[Directory]? fastq_dirs
        Array[File]? fastq_files
        String r1_name = '(.*?)_.*R1_001.fastq.gz'
        String r2_name = '(.*?)_.*R2_001.fastq.gz'
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
        String? other_parameters = ' '
        Int threads = 4
        File read1
        File? read2
        String sample_name
        String? adapter_r1 = "auto"
        String? adapter_r2 = "auto"
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
        Array[File] bwa_idxes
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
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        ~{sep=' ' if length(intervals) > 0 then prefix("--interval ", intervals) else [' ']} \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        ~{"--algo MeanQualityByCycle " + mq_metrics} \
        ~{"--algo QualDistribution " + qd_metrics} \
        ~{"--algo GCBias --summary " + gc_summary} ~{gc_metrics} \
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
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        String score
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e
        sentieon driver \
        ~{"-t " + t} \
        ~{"-i " + bam} \
        ~{"--algo LocusCollector --fun score_info " + score}
    >>>

    output {
        File score = "~{score}"
        File score_idx = "~{score}.idx"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        ~{sep=' ' if length(intervals) > 0 then prefix("--interval ", intervals) else [' ']} \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        ~{"--algo CoverageMetrics " + coverage_metrics}
    >>>

    output {
        File coverage_metrics = "~{coverage_metrics}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        File bam
        File bam_bai
        Array[File] database
        Array[File] database_idx
        String realigned_bam
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e
        sentieon driver \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        --algo Realigner \
        ~{sep=' ' prefix("-k ", database)} \
        ~{realigned_bam}
    >>>

    output {
        File realigned_bam = "~{realigned_bam}"
        File realigned_bam_bai = "~{realigned_bam}.bai"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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
        Array[File] database_idx
        String recal_data
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e
        sentieon driver \
        ~{sep=' ' if length(intervals) > 0 then prefix("--interval ", intervals) else [' ']} \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        --algo QualCal \
        ~{sep=' ' prefix("-k  ", database)} \
        ~{recal_data}
    >>>

    output {
        File recal_data = "~{recal_data}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
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

task Haplotyper{
    input {
        Array[File] intervals
        File bam
        File bam_bai
        File recal_data
        File ref
        Array[File] ref_idxes
        String emit_mode = "gvcf"
        Int ploidy = 2
        String out_vcf
        # for runtime
        String docker = "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
    }

    command <<<
        set -e
        sentieon driver \
        ~{sep=' ' if length(intervals) > 0 then prefix("--interval ", intervals) else [' ']} \
        ~{"-i " + bam} \
        ~{"-q " + recal_data} \
        ~{"-r " + ref} \
        --algo Haplotyper \
        ~{"--emit_mode " + emit_mode} \
        ~{"--ploidy " + ploidy} \
        ~{out_vcf}
    >>>

    output {
        File out_vcf = "~{out_vcf}"
        File out_vcf_idx = "~{out_vcf}.tbi"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
    }

    meta {
        name: "Haplotyper"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        intervals: {prefix: "--interval ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "interval file, support bed file or picard interval or vcf format"}
        bam: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reccaled bam file"}
        recal_data: {prefix: "-q ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "recal data"}
        ref: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        emit_mode: {prefix: "--emit_mode ", type: "str", level: "required", default: "gvcf", range: "None", array: "False", desc: "determines what calls will be emitted. possible values:variant,confident,all,gvcf"}
        ploidy: {prefix: "--ploidy ", type: "int", level: "required", default: "2", range: "None", array: "False", desc: "determines the ploidy number of the sample being processed. The default value is 2."}
        out_vcf: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output vcf file"}
    }

}

task GVCFtyper{
    input {
        File ref
        Array[File] ref_idxes
        Array[File] in_gvcf
        Array[File] in_gvcf_idx
        File known_dbsnp
        File known_dbsnp_idx
        Int call_conf = 30
        String genotype_model = "multinomial"
        String out_vcf
        # for runtime
        String docker = "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
    }

    command <<<
        set -e
        sentieon driver \
        ~{"-r " + ref} \
        --algo GVCFtyper \
        ~{sep=" " prefix("-v ", in_gvcf)} \
        ~{"-d " + known_dbsnp} \
        ~{"--call_conf " + call_conf} \
        ~{"--genotype_model " + genotype_model} \
        ~{out_vcf}
    >>>

    output {
        File out_vcf = "~{out_vcf}"
        File out_vcf_idx = "~{out_vcf}.tbi"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
    }

    meta {
        name: "GVCFtyper"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        ref: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        in_gvcf: {prefix: "-v ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input gvcf file"}
        known_dbsnp: {prefix: "-d ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "dbsnp file"}
        call_conf: {prefix: "--call_conf ", type: "int", level: "required", default: "30", range: "None", array: "False", desc: "determine the threshold of variant quality to emit a variant. Variants with quality less than CONFIDENCE will be not be added to the output VCF file."}
        genotype_model: {prefix: "--genotype_model ", type: "str", level: "required", default: "multinomial", range: "{'multinomial', 'coalescent'}", array: "False", desc: "determines which model to use for genotyping and QUAL calculation"}
        out_vcf: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output vcf file"}
    }

}

task VEP{
    input {
        File input_file
        File? input_file_idx
        File fasta
        Array[File] fasta_idx
        String output_file = "tumor.vep.vcf.gz"
        String output_format = "vcf"
        String compress_output = "bgzip"
        Boolean force_overwrite = true
        Int fork = 4
        # 解压缩需要7min左右的时间，但没有更好的办法
        File cache_targz
        File plugins_zip
        Array[String] plugin_names = ['Frameshift', 'Wildtype']
        String species = "homo_sapiens"
        String assembly_version = "GRCh37"
        String stats_file = "tumor.vep.summary.html"
        Boolean cache = true
        Boolean offline = true
        Boolean merged = false
        Boolean variant_class = true
        String sift = "b"
        String polyphen = "b"
        String nearest = "transcript"
        Boolean gene_phenotype = true
        Boolean regulatory = true
        Boolean phased = true
        Boolean numbers = true
        Boolean hgvs = true
        Boolean transcript_version = true
        Boolean symbol = true
        Boolean tsl = true
        Boolean canonical = true
        Boolean biotype = true
        Boolean mane = true
        Boolean max_af = true
        Boolean af_1kg = true
        Boolean af_gnomad = true
        Boolean af_esp = false
        Boolean coding_only = false
        Boolean pick = false
        Boolean flag_pick = true
        Boolean filter_common = true
        # for runtime
        String docker = "ensemblorg/ensembl-vep:release_104.3"
    }

    command <<<
        set -e
        tar zxf ~{cache_targz}
        unzip ~{plugins_zip} -d Plugins
        vep \
        ~{"-i " + input_file} \
        ~{"--fasta " + fasta} \
        ~{"-o " + output_file} \
        ~{"--" + output_format} \
        ~{"--compress_output " + compress_output} \
        ~{if force_overwrite then "--force_overwrite  " else ""} \
        ~{"--fork " + fork} \
        ~{"--species " + species} \
        ~{"--assembly " + assembly_version} \
        --dir_cache ./ \
        --dir_plugins ./Plugins \
        ~{"--stats_file " + stats_file} \
        ~{if cache then "--cache  " else ""} \
        ~{if offline then "--offline  " else ""} \
        ~{if merged then "--merged  " else ""} \
        ~{sep=" " prefix("--plugin ", plugin_names)} \
        ~{if variant_class then "--variant_class  " else ""} \
        ~{"--sift " + sift} \
        ~{"--polyphen " + polyphen} \
        ~{"--nearest " + nearest} \
        ~{if gene_phenotype then "--gene_phenotype  " else ""} \
        ~{if regulatory then "--regulatory  " else ""} \
        ~{if phased then "--phased  " else ""} \
        ~{if numbers then "--numbers  " else ""} \
        ~{if hgvs then "--hgvs  " else ""} \
        ~{if transcript_version then "--transcript_version  " else ""} \
        ~{if symbol then "--symbol  " else ""} \
        ~{if tsl then "--tsl  " else ""} \
        ~{if canonical then "--canonical  " else ""} \
        ~{if biotype then "--biotype  " else ""} \
        ~{if mane then "--mane " else ""} \
        ~{if max_af then "--max_af  " else ""} \
        ~{if af_1kg then "--af_1kg  " else ""} \
        ~{if af_gnomad then "--af_gnomad  " else ""} \
        ~{if af_esp then "--af_esp  " else ""} \
        ~{if coding_only then "--af_esp  " else ""} \
        ~{if pick then "--pick " else ""} \
        ~{if flag_pick then "--flag_pick  " else ""} \
        ~{if filter_common then "--filter_common  " else ""}
        rm -fr homo_sapiens Plugins
    >>>

    output {
        File out_vcf = "~{output_file}"
        File out_vcf_idx = "~{output_file}.tbi"
        File stats_file = "~{stats_file}"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/vep:release_104.3"
    }

    meta {
        name: "VEP"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        input_file: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input file"}
        fasta: {prefix: "--fasta ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence. The first time you run VEP with this parameter an index will be built which can take a few minutes. This is required if fetching HGVS annotations (--hgvs) or checking reference sequences (--check_ref) in offline mode (--offline), and optional with some performance increase in cache mode (--cache)."}
        output_file: {prefix: "-o ", type: "str", level: "required", default: "tumor.vep.vcf.gz", range: "None", array: "False", desc: "output file"}
        output_format: {prefix: "--", type: "str", level: "required", default: "vcf", range: "{'vcf', 'json', 'tab'}", array: "False", desc: "If we choose to write output in VCF format. Consequences are added in the INFO field of the VCF file, using the key 'CSQ'. Data fields are encoded separated by '|'; the order of fields is written in the VCF header. Output fields in the 'CSQ' INFO field can be selected by using --fields."}
        compress_output: {prefix: "--compress_output ", type: "str", level: "required", default: "bgzip", range: "None", array: "False", desc: "Writes output compressed using either gzip or bgzip"}
        force_overwrite: {prefix: "--force_overwrite ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Force overwriting of output file"}
        fork: {prefix: "--fork ", type: "int", level: "required", default: "4", range: "None", array: "False", desc: "Use forking(multi-cpu/threads) to improve script runtime"}
        species: {prefix: "--species ", type: "str", level: "required", default: "homo_sapiens", range: "None", array: "False", desc: "Species for your data. This can be the latin name e.g. homo_sapiens or any Ensembl alias e.g. mouse."}
        assembly_version: {prefix: "--assembly ", type: "str", level: "required", default: "GRCh37", range: "None", array: "False", desc: "Select the assembly version to use if more than one available."}
        cache_targz: {prefix: "--dir_cache ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Specify the cache files of VEP to use"}
        plugins_zip: {prefix: "--dir_plugins ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Specify the plugin files to use"}
        stats_file: {prefix: "--stats_file ", type: "str", level: "required", default: "tumor.vep.summary.html", range: "None", array: "False", desc: "Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end with <.html>."}
        cache: {prefix: "--cache ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Enables use of cache"}
        offline: {prefix: "--offline ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Enables offline mode. No database connections, and a cache file or GFF/GTF file is required for annotation"}
        merged: {prefix: "--merged ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used."}
        plugin_names: {prefix: "--plugin ", type: "str", level: "required", default: "['Frameshift', 'Wildtype']", range: "None", array: "False", desc: "Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory.Multiple plugins can be used by supplying the --plugin flag multiple times"}
        variant_class: {prefix: "--variant_class ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Output the Sequence Ontology variant class."}
        sift: {prefix: "--sift ", type: "str", level: "required", default: "b", range: "{'b', 's', 'p'}", array: "False", desc: "Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both."}
        polyphen: {prefix: "--polyphen ", type: "str", level: "required", default: "b", range: "{'b', 's', 'p'}", array: "False", desc: "Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both."}
        nearest: {prefix: "--nearest ", type: "str", level: "required", default: "transcript", range: "{'gene', 'transcript', 'symbol'}", array: "False", desc: "Retrieve the transcript or gene with the nearest protein-coding transcription start site (TSS) to each input variant. Use transcript to retrieve the transcript stable ID, gene to retrieve the gene stable ID, or symbol to retrieve the gene symbol. Note that the nearest TSS may not belong to a transcript that overlaps the input variant, and more than one may be reported in the case where two are equidistant from the input coordinates."}
        gene_phenotype: {prefix: "--gene_phenotype ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Indicates if the overlapped gene is associated with a phenotype, disease or trait."}
        regulatory: {prefix: "--regulatory ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature."}
        phased: {prefix: "--phased ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data."}
        numbers: {prefix: "--numbers ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds affected exon and intron numbering to to output. Format is Number/Total"}
        hgvs: {prefix: "--hgvs ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate."}
        transcript_version: {prefix: "--transcript_version ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add version numbers to Ensembl transcript identifiers"}
        symbol: {prefix: "--symbol ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the gene symbol (e.g. HGNC) (where available) to the output."}
        tsl: {prefix: "--tsl ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the transcript support level for this transcript to the output."}
        canonical: {prefix: "--canonical ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds a flag indicating if the transcript is the canonical transcript for the gene."}
        biotype: {prefix: "--biotype ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the biotype of the transcript or regulatory feature."}
        mane: {prefix: "--mane ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds a flag indicating if the transcript is the MANE Select or MANE Plus Clinical transcript for the gene."}
        max_af: {prefix: "--max_af ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD."}
        af_1kg: {prefix: "--af_1kg ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output."}
        af_gnomad: {prefix: "--af_gnomad ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included"}
        af_esp: {prefix: "--af_esp ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Include allele frequency from NHLBI-ESP populations."}
        coding_only: {prefix: "--af_esp ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Only return consequences that fall in the coding regions of transcripts. Not used by default"}
        pick: {prefix: "--pick", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Pick one line or block of consequence data per variant, including transcript-specific columns. This is the best method to use if you are interested only in one consequence per variant"}
        flag_pick: {prefix: "--flag_pick ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others."}
        filter_common: {prefix: "--filter_common ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters."}
    }

}

task OptiType{
    input {
        Array[File] reads
        Boolean is_dna = true
        Boolean is_rna = false
        Int threads = 5
        Int enumerate = 1
        String outdir = "."
        String prefix
        # for runtime
        String docker = "fred2/optitype:1.3.1"
    }

    command <<<
        set -e
        cp /usr/local/bin/OptiType/config.ini .
        sed -i 's/threads=1/threads=~{threads}/g' config.ini
        OptiTypePipeline.py \
        ~{if defined(reads) then "--input  " else ""}~{sep=" " reads} \
        ~{if is_dna then "--dna " else ""} \
        ~{if is_rna then "--rna " else ""} \
        ~{"--enumerate " + enumerate} \
        ~{"--outdir " + outdir} \
        ~{"--prefix " + prefix} \
        --config config.ini
    >>>

    output {
        File result_tsv = "~{prefix}_result.tsv"
        File result_pdf = "~{prefix}_coverage_plot.pdf"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/optitype:1.3.1"
    }

    meta {
        name: "OptiType"
        desc: "OptiType: 4-digit HLA typer"
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        reads: {prefix: "--input ", type: "infile", level: "required", default: "None", range: "None", array: "True", desc: "fastq file(s) (fished or raw) or .bam files stored for re-use, generated by an earlier OptiType run."}
        is_dna: {prefix: "--dna", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "use with DNA sequencing data"}
        is_rna: {prefix: "--rna", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "use with RNA sequencing data"}
        enumerate: {prefix: "--enumerate ", type: "int", level: "required", default: "1", range: "None", array: "False", desc: "Number of enumerations. OptiType will output the optimal solution and the top N-1 suboptimal solutions in the results CSV."}
        outdir: {prefix: "--outdir ", type: "str", level: "required", default: ".", range: "None", array: "False", desc: "Specifies the out directory to which all files should be written."}
        prefix: {prefix: "--prefix ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "prefix of output files"}
    }

}