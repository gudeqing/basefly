version development

workflow pipeline {
    input {
        Int thread_number = 16
        File ref = "human_g1k_v37_decoy.fasta"
        File known_dbsnp = "dbsnp_138.b37.vcf.gz"
        File known_indel = "1000G_phase1.indels.b37.vcf.gz"
        File known_mills = "Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
        File pon = "PanelOfNormal.vcf"
        File germline_vcf = "germline.vcf"
        Directory snpeff_databse = "data/"
        Directory vep_cache_dir = "vep/"
        Directory vep_plugin_dir = "vep/"
        File intervals = "['interval.1', 'interval.2']"
    }

    call getFastqInfo{}
    scatter (each in keys(getFastqInfo.fastq_info)) { 
        String sample = each
        File read1 = getFastqInfo.fastq_info[each][0][0]
        File read2 = getFastqInfo.fastq_info[each][1][0]
        call OptiType {
            input: 
            reads = [],
            prefix = "~{sample}"
        }

        call bwa_mem {
            input: 
            readgroup = "@RG\tID:~{sample}\tSM:~{sample}\tPL:ILLUMINA",
            t = thread_number,
            ref = ref,
            read1 = "~{read1}",
            read2 = "~{read2}",
            ref2 = ref,
            out = "~{sample}.sorted.bam",
            t2 = thread_number
        }

        call get_metrics {
            input: 
            t = thread_number,
            ref = ref,
            bam = bwa_mem.out,
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
            bam = bwa_mem.out,
            score = "~{sample}.score.txt"
        }

        call DeDup {
            input: 
            t = thread_number,
            bam = bwa_mem.out,
            score = LocusCollector.score,
            dedup_metrics = "~{sample}.dedup.metrics.txt",
            deduped_bam = "~{sample}.deduped.bam"
        }

        call CoverageMetrics {
            input: 
            t = thread_number,
            ref = ref,
            bam = DeDup.deduped_bam,
            coverage_metrics = "~{sample}.cov.metrics.txt"
        }

        call realign {
            input: 
            t = thread_number,
            ref = ref,
            bam = DeDup.deduped_bam,
            database = ['known_indel', 'known_mills'],
            realigned_bam = "~{sample}.realigned.bam"
        }

        call recalibration {
            input: 
            t = thread_number,
            ref = ref,
            bam = realign.realigned_bam,
            database = ['known_indel', 'known_dbsnp', 'known_mills'],
            recal_data = "~{sample}.recal_data.table"
        }

    }

    call Haplotyper {
        input: 
        intervals = intervals,
        bam = realign.realigned_bam,
        recal_data = recalibration.recal_data,
        ref = ref,
        out_vcf = "~{normal_sample}.g.vcf.gz"
    }

    call GVCFtyper {
        input: 
        ref = ref,
        in_gvcf = Haplotyper.out_vcf,
        known_dbsnp = known_dbsnp,
        out_vcf = "~{normal_sample}.vcf.gz"
    }

    call TNhaplotyper2 {
        input: 
        t = thread_number,
        ref = ref,
        bams = realign.realigned_bam,
        recal_datas = recalibration.recal_data,
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
        data_dir = snpeff_databse,
        in_vcf = TNfilter.out_vcf,
        out_vcf = "~{tumor_sample}.final.annot.vcf"
    }

    call VEP {
        input: 
        input_file = TNfilter.out_vcf,
        fasta = ref,
        dir_cache = vep_cache_dir,
        dir_plugins = vep_plugin_dir
    }

    call CombineVariants {
        input: 
        ref = ref,
        variant = ['TNfilter.out_vcf', 'GVCFtyper.out_vcf']
    }

    call SortVcf {
        input: 
        in_vcf = CombineVariants.combined_vcf
    }

    call ReadBackedPhasing {
        input: 
        ref = ref,
        bam = realign.realigned_bam,
        variant = SortVcf.sorted_vcf,
        interval = SortVcf.sorted_vcf
    }

    call vep_phased {
        input: 
        input_file = ReadBackedPhasing.phased_vcf,
        fasta = ref,
        dir_cache = vep_cache_dir,
        dir_plugins = vep_plugin_dir
    }

    meta {
        name: "TN_pipeline"
        desc: "typical bioinformatics pipeline using sentieon TNSeq and snpEff"
        author: "unknown"
        source: "source URL for the tool"
    }

    output{
        Array[File] OptiType_result_tsv = OptiType.result_tsv
        Array[File] OptiType_result_pdf = OptiType.result_pdf
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
        Array[File] CoverageMetrics_coverage_metrics = CoverageMetrics.coverage_metrics
        Array[File] realign_realigned_bam = realign.realigned_bam
        Array[File] recalibration_recal_data = recalibration.recal_data
        File Haplotyper_out_vcf = Haplotyper.out_vcf
        File Haplotyper_out_vcf_idx = Haplotyper.out_vcf_idx
        File GVCFtyper_out_vcf = GVCFtyper.out_vcf
        File GVCFtyper_out_vcf_idx = GVCFtyper.out_vcf_idx
        File TNhaplotyper2_out_vcf = TNhaplotyper2.out_vcf
        File TNhaplotyper2_orientation_data = TNhaplotyper2.orientation_data
        File TNhaplotyper2_tumor_segments = TNhaplotyper2.tumor_segments
        File TNhaplotyper2_contamination_data = TNhaplotyper2.contamination_data
        File TNfilter_out_vcf = TNfilter.out_vcf
        File snpEff_out_vcf = snpEff.out_vcf
        File VEP_out_vcf = VEP.out_vcf
        File VEP_out_vcf_idx = VEP.out_vcf_idx
        File CombineVariants_combined_vcf = CombineVariants.combined_vcf
        File SortVcf_sorted_vcf = SortVcf.sorted_vcf
        File ReadBackedPhasing_phased_vcf = ReadBackedPhasing.phased_vcf
        File vep_phased_out_vcf = vep_phased.out_vcf
        File vep_phased_out_vcf_idx = vep_phased.out_vcf_idx
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
        python /get_fastq_info.py             ~{if defined(fastq_dirs) then "-fastq_dirs " else ""}~{sep=" " fastq_dirs}             ~{if defined(fastq_files) then "-fastq_files " else ""}~{sep=" " fastq_files}             -r1_name '~{r1_name}'             -r2_name '~{r2_name}'             -out fastq.info.json
    >>>

    output {
        Map[String, Array[Array[File]]] fastq_info = read_json("fastq.info.json")
    }

    runtime {
        docker: "registry-xdp-v3-pre-yifang.basebit.me/basebitai/getfastqinfo:1.0"
    }

    parameter_meta {
        fastq_dirs: {desc: "directory list, target fastq files should be in these directories. All target files in 'fastq_files' or 'fastq_dirs' will be used", level: "optional", type: "indir", range: "", default: ""}
        fastq_files: {desc: "target fastq file list. 'fastq_files' or 'fastq_dirs' must be provided.", level: "optional", type: "infile", range: "", default: ""}
        r1_name: {desc: "python regExp that describes the full name of read1 fastq file name. It requires at least one pair of small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'", level: "required", type: "str", range: "", default: ""}
        r2_name: {desc: "python regExp that describes the full name of read2 fastq file name. It requires at least one pair of small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'", level: "required", type: "str", range: "", default: ""}
    }
}
    
task OptiType{
    input {
        Array[File] reads
        Boolean is_dna = true
        Boolean is_rna = false
        Int enumerate = 1
        String outdir = "."
        String prefix
        String config = "config.ini"
        # for runtime
        String docker = "fred2/optitype:1.3.1"
    }

    command <<<
        set -e 
        OptiTypePipeline.py \
        ~{if defined(reads) then "--input  " else ""}~{sep=" " reads} \
        ~{if is_dna then "--dna " else ""} \
        ~{if is_rna then "--rna " else ""} \
        ~{"--enumerate " + enumerate} \
        ~{"--outdir " + outdir} \
        ~{"--prefix " + prefix} \
        ~{"--config " + config} 
    >>>

    output {
        File result_tsv = "~{prefix}_result.tsv"
        File result_pdf = "~{prefix}_coverage_plot.pdf"
    }

    runtime {
        docker: docker
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
        config: {prefix: "--config ", type: "str", level: "required", default: "config.ini", range: "None", array: "False", desc: "config.ini file"}
    }

}

task bwa_mem{
    input {
        String readgroup
        Int t = 16
        Int k = 10000000
        File ref
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
        ~{"-R " + readgroup} \
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
    }

    runtime {
        docker: docker
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
        File bam
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
        docker: docker
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
        docker: docker
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
        docker: docker
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
        docker: docker
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
        docker: docker
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
    }

    runtime {
        docker: docker
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
        File score
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
    }

    runtime {
        docker: docker
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
        File bam
        String coverage_metrics
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        ~{"--algo CoverageMetrics " + coverage_metrics} 
    >>>

    output {
        File coverage_metrics = "~{coverage_metrics}"
    }

    runtime {
        docker: docker
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
        File bam
        Array[File] database
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
        ~{sep=" " prefix("-k ", database)} \
        ~{realigned_bam} 
    >>>

    output {
        File realigned_bam = "~{realigned_bam}"
    }

    runtime {
        docker: docker
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
        File bam
        Array[File] database
        String recal_data
        # for runtime
        String docker = "docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{"-i " + bam} \
        --algo QualCal \
        ~{sep=" " prefix("-k ", database)} \
        ~{recal_data} 
    >>>

    output {
        File recal_data = "~{recal_data}"
    }

    runtime {
        docker: docker
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
        File recal_data
        File ref
        String emit_mode = "gvcf"
        Int ploidy = 2
        String out_vcf
        # for runtime
        String docker = "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02"
    }

    command <<<
        set -e 
        sentieon driver \
        ~{sep=" " prefix("--interval ", intervals)} \
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
        docker: docker
    }

    meta {
        name: "Haplotyper"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        intervals: {prefix: "--interval ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "interval file, support bed file or picard interval or vcf format"}
        bam: {prefix: "-i ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reccaled tumor and normal bam list"}
        recal_data: {prefix: "-q ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "tumor and normal recal data list"}
        ref: {prefix: "-r ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        emit_mode: {prefix: "--emit_mode ", type: "str", level: "required", default: "gvcf", range: "None", array: "False", desc: "determines what calls will be emitted. possible values:variant,confident,all,gvcf"}
        ploidy: {prefix: "--ploidy ", type: "int", level: "required", default: "2", range: "None", array: "False", desc: "determines the ploidy number of the sample being processed. The default value is 2."}
        out_vcf: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output vcf file"}
    }

}

task GVCFtyper{
    input {
        File ref
        Array[File] in_gvcf
        File known_dbsnp
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
        docker: docker
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
        genotype_model: {prefix: "--genotype_model ", type: "str", level: "required", default: "multinomial", range: "{'coalescent', 'multinomial'}", array: "False", desc: "determines which model to use for genotyping and QUAL calculation"}
        out_vcf: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output vcf file"}
    }

}

task TNhaplotyper2{
    input {
        Int t = 16
        File ref
        Array[File] bams
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
        ~{"-t " + t} \
        ~{"-r " + ref} \
        ~{sep=" " prefix("-i ", bams)} \
        ~{sep=" " prefix("-q ", recal_datas)} \
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
        docker: docker
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
        docker: docker
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
        String genome_version = "hg19"
        Directory data_dir
        Boolean cancer = true
        File? cancerSamples
        Boolean canon = false
        Array[File]? interval
        String? other_args
        File in_vcf
        String out_vcf
        # for runtime
        String docker = "?"
    }

    command <<<
        set -e 
        java -Xmx9g snpEff.jar ann \
        ~{genome_version} \
        ~{"-dataDir " + data_dir} \
        ~{if cancer then "-cancer  " else ""} \
        ~{"-cancerSamples " + cancerSamples} \
        ~{if canon then "-canon  " else ""} \
        ~{sep=" " prefix("-interval ", interval)} \
        ~{other_args} \
        ~{in_vcf} \
        > \
        ~{out_vcf} 
    >>>

    output {
        File out_vcf = "~{out_vcf}"
    }

    runtime {
        docker: docker
    }

    meta {
        name: "snpEff"
        desc: "This is description of the tool/workflow."
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        genome_version: {prefix: "", type: "str", level: "required", default: "hg19", range: "None", array: "False", desc: "human genome version"}
        data_dir: {prefix: "-dataDir ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Override data_dir parameter from config file"}
        cancer: {prefix: "-cancer ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Perform cancer comparisons (Somatic vs Germline)"}
        cancerSamples: {prefix: "-cancerSamples ", type: "infile", level: "optional", default: "None", range: "None", array: "False", desc: "Two column TXT file defining 'oringinal derived' samples. If '-cancer' used and the file is missing, then the last sample will be assumed as tumor sample."}
        canon: {prefix: "-canon ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Only use canonical transcripts"}
        interval: {prefix: "-interval ", type: "infile", level: "optional", default: "None", range: "None", array: "False", desc: "Use a custom intervals in TXT/BED/BigBed/VCF/GFF file (you may use this option many times)"}
        other_args: {prefix: "", type: "str", level: "optional", default: "None", range: "None", array: "False", desc: "other arguments that you want to input for the program, such as '-motif'"}
        in_vcf: {prefix: "", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input variant file"}
        out_vcf: {prefix: "", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output annotated file"}
    }

}

task VEP{
    input {
        File input_file
        File fasta
        String output_file = "tumor.vep.vcf.gz"
        String output_format = "vcf"
        String compress_output = "bgzip"
        Boolean force_overwrite = true
        Int fork = 4
        String species = "homo_sapiens"
        String assembly_version = "GRCh37"
        Directory dir_cache
        Directory dir_plugins
        String stats_file = "tumor.vep.summary.html"
        Boolean cache = true
        Boolean offline = true
        Boolean merged = false
        Array[String] plugins = "['Frameshift', 'Wildtype']"
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
        Boolean max_af = true
        Boolean af_1kg = true
        Boolean af_gnomad = true
        Boolean af_esp = false
        Boolean coding_only = false
        Boolean pick = false
        Boolean flag_pick = true
        Boolean filter_common = true
        String other_args = ""
        # for runtime
        String docker = "ensemblorg/ensembl-vep:2.0.3"
    }

    command <<<
        set -e 
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
        ~{"--dir_cache " + dir_cache} \
        ~{"--dir_plugins " + dir_plugins} \
        ~{"--stats_file " + stats_file} \
        ~{if cache then "--cache  " else ""} \
        ~{if offline then "--offline  " else ""} \
        ~{if merged then "--merged  " else ""} \
        ~{sep=" " prefix("--plugin ", plugins)} \
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
        ~{if max_af then "--max_af  " else ""} \
        ~{if af_1kg then "--af_1kg  " else ""} \
        ~{if af_gnomad then "--af_gnomad  " else ""} \
        ~{if af_esp then "--af_esp  " else ""} \
        ~{if coding_only then "--af_esp  " else ""} \
        ~{if pick then "--pick " else ""} \
        ~{if flag_pick then "--flag_pick  " else ""} \
        ~{if filter_common then "--filter_common  " else ""} \
        ~{other_args} 
    >>>

    output {
        File out_vcf = "~{output_file}"
        File out_vcf_idx = "~{output_file}.tbi"
    }

    runtime {
        docker: docker
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
        dir_cache: {prefix: "--dir_cache ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Specify the cache directory to use"}
        dir_plugins: {prefix: "--dir_plugins ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Specify the plugin directory to use"}
        stats_file: {prefix: "--stats_file ", type: "str", level: "required", default: "tumor.vep.summary.html", range: "None", array: "False", desc: "Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end with <.html>."}
        cache: {prefix: "--cache ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Enables use of cache"}
        offline: {prefix: "--offline ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Enables offline mode. No database connections, and a cache file or GFF/GTF file is required for annotation"}
        merged: {prefix: "--merged ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used."}
        plugins: {prefix: "--plugin ", type: "str", level: "required", default: "['Frameshift', 'Wildtype']", range: "None", array: "False", desc: "Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory.Multiple plugins can be used by supplying the --plugin flag multiple times"}
        variant_class: {prefix: "--variant_class ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Output the Sequence Ontology variant class."}
        sift: {prefix: "--sift ", type: "str", level: "required", default: "b", range: "{'p', 'b', 's'}", array: "False", desc: "Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both."}
        polyphen: {prefix: "--polyphen ", type: "str", level: "required", default: "b", range: "{'p', 'b', 's'}", array: "False", desc: "Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both."}
        nearest: {prefix: "--nearest ", type: "str", level: "required", default: "transcript", range: "{'gene', 'symbol', 'transcript'}", array: "False", desc: "Retrieve the transcript or gene with the nearest protein-coding transcription start site (TSS) to each input variant. Use transcript to retrieve the transcript stable ID, gene to retrieve the gene stable ID, or symbol to retrieve the gene symbol. Note that the nearest TSS may not belong to a transcript that overlaps the input variant, and more than one may be reported in the case where two are equidistant from the input coordinates."}
        gene_phenotype: {prefix: "--gene_phenotype ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Indicates if the overlapped gene is associated with a phenotype, disease or trait."}
        regulatory: {prefix: "--regulatory ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature."}
        phased: {prefix: "--phased ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data."}
        numbers: {prefix: "--numbers ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds affected exon and intron numbering to to output. Format is Number/Total"}
        hgvs: {prefix: "--hgvs ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate."}
        transcript_version: {prefix: "--transcript_version ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add version numbers to Ensembl transcript identifiers"}
        symbol: {prefix: "--symbol ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the gene symbol (e.g. HGNC) (where available) to the output."}
        tsl: {prefix: "--tsl ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the transcript support level for this transcript to the output."}
        canonical: {prefix: "--canonical ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds a flag indicating if the transcript is the canonical transcript for the gene"}
        biotype: {prefix: "--biotype ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the biotype of the transcript or regulatory feature."}
        max_af: {prefix: "--max_af ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD"}
        af_1kg: {prefix: "--af_1kg ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output."}
        af_gnomad: {prefix: "--af_gnomad ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included"}
        af_esp: {prefix: "--af_esp ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Include allele frequency from NHLBI-ESP populations."}
        coding_only: {prefix: "--af_esp ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Only return consequences that fall in the coding regions of transcripts. Not used by default"}
        pick: {prefix: "--pick", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Pick one line or block of consequence data per variant, including transcript-specific columns. This is the best method to use if you are interested only in one consequence per variant"}
        flag_pick: {prefix: "--flag_pick ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others."}
        filter_common: {prefix: "--filter_common ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters."}
        other_args: {prefix: "", type: "str", level: "required", default: "", range: "None", array: "False", desc: "specify other arguments that you want to append to the command"}
    }

}

task CombineVariants{
    input {
        File ref
        Array[File] variant
        String out_vcf
        Boolean assumeIdenticalSamples = false
        # for runtime
        String docker = "broadinstitute/gatk3:3.8-1"
    }

    command <<<
        set -e 
        java -Xmx10g -jar GenomeAnalysisTK.jar -T CombineVariants \
        ~{"-R " + ref} \
        ~{sep=" " prefix("--variant ", variant)} \
        ~{"-o " + out_vcf} \
        ~{if assumeIdenticalSamples then "--assumeIdenticalSamples " else ""} 
    >>>

    output {
        File combined_vcf = "~{out_vcf}"
    }

    runtime {
        docker: docker
    }

    meta {
        name: "CombineVariants"
        desc: "Combine variants"
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        ref: {prefix: "-R ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        variant: {prefix: "--variant ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "variant vcf file array"}
        out_vcf: {prefix: "-o ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "This is description of the argument."}
        assumeIdenticalSamples: {prefix: "--assumeIdenticalSamples", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "If true, assume input VCFs have identical sample sets and disjoint calls. This option allows the user to perform a simple merge (concatenation) to combine the VCFs."}
    }

}

task SortVcf{
    input {
        File in_vcf
        File out_vcf
        # for runtime
        String docker = "broadinstitute/picard:latest"
    }

    command <<<
        set -e 
        java -jar /usr/picard/picard.jar SortVcf \
        ~{"I=" + in_vcf} \
        ~{"O=" + out_vcf} 
    >>>

    output {
        File sorted_vcf = "~{out_vcf}"
    }

    runtime {
        docker: docker
    }

    meta {
        name: "SortVcf"
        desc: "sort vcf"
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        in_vcf: {prefix: "I=", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input vcf to sort"}
        out_vcf: {prefix: "O=", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "output sorted vcf"}
    }

}

task ReadBackedPhasing{
    input {
        File ref
        File bam
        File variant
        File interval
        String out_vcf
        # for runtime
        String docker = "broadinstitute/gatk3:3.8-1"
    }

    command <<<
        set -e 
        java -Xmx10g -jar GenomeAnalysisTK.jar -T ReadBackedPhasing \
        ~{"-R " + ref} \
        ~{"-I " + bam} \
        ~{"--variant " + variant} \
        ~{"-L " + interval} \
        ~{"-o " + out_vcf} 
    >>>

    output {
        File phased_vcf = "~{out_vcf}"
    }

    runtime {
        docker: docker
    }

    meta {
        name: "ReadBackedPhasing"
        desc: "ReadBackedPhasing"
        author: "unknown"
        source: "source URL for the tool"
    }

    parameter_meta {
        ref: {prefix: "-R ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "reference fasta file"}
        bam: {prefix: "-I ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "tumor bam file"}
        variant: {prefix: "--variant ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input vcf file"}
        interval: {prefix: "-L ", type: "infile", level: "required", default: "None", range: "None", array: "False", desc: "input vcf file"}
        out_vcf: {prefix: "-o ", type: "str", level: "required", default: "None", range: "None", array: "False", desc: "output vcf file"}
    }

}

task vep_phased{
    input {
        File input_file
        File fasta
        String output_file = "tumor.vep.vcf.gz"
        String output_format = "vcf"
        String compress_output = "bgzip"
        Boolean force_overwrite = true
        Int fork = 4
        String species = "homo_sapiens"
        String assembly_version = "GRCh37"
        Directory dir_cache
        Directory dir_plugins
        String stats_file = "tumor.vep.summary.html"
        Boolean cache = true
        Boolean offline = true
        Boolean merged = false
        Array[String] plugins = "['Frameshift', 'Wildtype']"
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
        Boolean max_af = true
        Boolean af_1kg = true
        Boolean af_gnomad = true
        Boolean af_esp = false
        Boolean coding_only = false
        Boolean pick = false
        Boolean flag_pick = true
        Boolean filter_common = true
        String other_args = ""
        # for runtime
        String docker = "ensemblorg/ensembl-vep:2.0.3"
    }

    command <<<
        set -e 
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
        ~{"--dir_cache " + dir_cache} \
        ~{"--dir_plugins " + dir_plugins} \
        ~{"--stats_file " + stats_file} \
        ~{if cache then "--cache  " else ""} \
        ~{if offline then "--offline  " else ""} \
        ~{if merged then "--merged  " else ""} \
        ~{sep=" " prefix("--plugin ", plugins)} \
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
        ~{if max_af then "--max_af  " else ""} \
        ~{if af_1kg then "--af_1kg  " else ""} \
        ~{if af_gnomad then "--af_gnomad  " else ""} \
        ~{if af_esp then "--af_esp  " else ""} \
        ~{if coding_only then "--af_esp  " else ""} \
        ~{if pick then "--pick " else ""} \
        ~{if flag_pick then "--flag_pick  " else ""} \
        ~{if filter_common then "--filter_common  " else ""} \
        ~{other_args} 
    >>>

    output {
        File out_vcf = "~{output_file}"
        File out_vcf_idx = "~{output_file}.tbi"
    }

    runtime {
        docker: docker
    }

    meta {
        name: "vep_phased"
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
        dir_cache: {prefix: "--dir_cache ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Specify the cache directory to use"}
        dir_plugins: {prefix: "--dir_plugins ", type: "indir", level: "required", default: "None", range: "None", array: "False", desc: "Specify the plugin directory to use"}
        stats_file: {prefix: "--stats_file ", type: "str", level: "required", default: "tumor.vep.summary.html", range: "None", array: "False", desc: "Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end with <.html>."}
        cache: {prefix: "--cache ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Enables use of cache"}
        offline: {prefix: "--offline ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Enables offline mode. No database connections, and a cache file or GFF/GTF file is required for annotation"}
        merged: {prefix: "--merged ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used."}
        plugins: {prefix: "--plugin ", type: "str", level: "required", default: "['Frameshift', 'Wildtype']", range: "None", array: "False", desc: "Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory.Multiple plugins can be used by supplying the --plugin flag multiple times"}
        variant_class: {prefix: "--variant_class ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Output the Sequence Ontology variant class."}
        sift: {prefix: "--sift ", type: "str", level: "required", default: "b", range: "{'p', 'b', 's'}", array: "False", desc: "Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both."}
        polyphen: {prefix: "--polyphen ", type: "str", level: "required", default: "b", range: "{'p', 'b', 's'}", array: "False", desc: "Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both."}
        nearest: {prefix: "--nearest ", type: "str", level: "required", default: "transcript", range: "{'gene', 'symbol', 'transcript'}", array: "False", desc: "Retrieve the transcript or gene with the nearest protein-coding transcription start site (TSS) to each input variant. Use transcript to retrieve the transcript stable ID, gene to retrieve the gene stable ID, or symbol to retrieve the gene symbol. Note that the nearest TSS may not belong to a transcript that overlaps the input variant, and more than one may be reported in the case where two are equidistant from the input coordinates."}
        gene_phenotype: {prefix: "--gene_phenotype ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Indicates if the overlapped gene is associated with a phenotype, disease or trait."}
        regulatory: {prefix: "--regulatory ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature."}
        phased: {prefix: "--phased ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data."}
        numbers: {prefix: "--numbers ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds affected exon and intron numbering to to output. Format is Number/Total"}
        hgvs: {prefix: "--hgvs ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate."}
        transcript_version: {prefix: "--transcript_version ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add version numbers to Ensembl transcript identifiers"}
        symbol: {prefix: "--symbol ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the gene symbol (e.g. HGNC) (where available) to the output."}
        tsl: {prefix: "--tsl ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the transcript support level for this transcript to the output."}
        canonical: {prefix: "--canonical ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds a flag indicating if the transcript is the canonical transcript for the gene"}
        biotype: {prefix: "--biotype ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Adds the biotype of the transcript or regulatory feature."}
        max_af: {prefix: "--max_af ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD"}
        af_1kg: {prefix: "--af_1kg ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output."}
        af_gnomad: {prefix: "--af_gnomad ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included"}
        af_esp: {prefix: "--af_esp ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Include allele frequency from NHLBI-ESP populations."}
        coding_only: {prefix: "--af_esp ", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Only return consequences that fall in the coding regions of transcripts. Not used by default"}
        pick: {prefix: "--pick", type: "bool", level: "required", default: "False", range: "{False, True}", array: "False", desc: "Pick one line or block of consequence data per variant, including transcript-specific columns. This is the best method to use if you are interested only in one consequence per variant"}
        flag_pick: {prefix: "--flag_pick ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others."}
        filter_common: {prefix: "--filter_common ", type: "bool", level: "required", default: "True", range: "{False, True}", array: "False", desc: "Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters."}
        other_args: {prefix: "", type: "str", level: "required", default: "", range: "None", array: "False", desc: "specify other arguments that you want to append to the command"}
    }

}

