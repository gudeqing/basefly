## 系统概述
### 系统名称
* GATK-DNAseq-Workflow
### 系统说明
* 系统版本号：1.0
* 适用范围：
* 参考来源：
### 系统主要功能

    当前流程是参考博得研究所最新的GATK-Best-Practice流程构建的DNAseq突变检测分析流程改写而成
    流程包含的主要功能：
    * 使用fastp进行接头自动去除
    * 使用BWA进行比对分析
    * 使用GATK检测small SNP/Indel, 同时支持tumor-only和tumor-normal配对模式
    * 使用GATK检测germline突变，如果输入多个normal样本，则直接进行joint-calling
    * 基于VEP进行突变注释
    
## 系统总体结构
### 模块列表
| name | desc | source | vesion |
| :--- | :---: | :---: | :---: |
|fastp|This is description of the tool/workflow.|source URL for the tool|unknown|
|FastqToSam|convert fastq to sam|source URL for the tool|unknown|
|uBam2FastqBwaMem|ubam to fastq and then mapping|source URL for the tool|unknown|
|MergeBamAlignment|merge bam alignment|source URL for the tool|unknown|
|MarkDuplicates|merge bam alignment|source URL for the tool|unknown|
|SortAndFixTags|sort bam and fix tags|source URL for the tool|unknown|
|BaseRecalibrator|This is description of the tool/workflow.|source URL for the tool|unknown|
|GatherBQSRReports|Gathers scattered BQSR recalibration reports into a single file|source URL for the tool|unknown|
|ApplyBQSR|Apply Base Quality Score Recalibration (BQSR) model|source URL for the tool|unknown|
|GatherBamFiles|Concatenate efficiently BAM files that resulted from a scattered parallel analysis.|source URL for the tool|unknown|
|SplitIntervals|This tool takes in intervals via the standard arguments of IntervalArgumentCollection and splits them into interval files for scattering. The resulting files contain equal number of bases.|source URL for the tool|unknown|
|Mutect2|Call somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations.|https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2|unknown|
|GetPileupSummaries|Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination|https://gatk.broadinstitute.org/hc/en-us/articles/360037593451-GetPileupSummaries|unknown|
|LearnReadOrientationModel|Learn the prior probability of read orientation artifact from the output of CollectF1R2Counts of Mutect2|source URL for the tool|unknown|
|MergeVcfs|Combines multiple variant files into a single variant file.|https://gatk.broadinstitute.org/hc/en-us/articles/360056969852-MergeVcfs-Picard-|unknown|
|MergeMutectStats|This is description of the tool/workflow.|source URL for the tool|unknown|
|FilterMutectCalls|FilterMutectCalls applies filters to the raw output of Mutect2|source URL for the tool|unknown|
|FilterAlignmentArtifacts|Alignment artifacts can occur whenever there is sufficient sequence similarity between two or more regions in the genome to confuse the alignment algorithm. This can occur when the aligner for whatever reason overestimate how uniquely a read maps, thereby assigning it too high of a mapping quality. It can also occur through no fault of the aligner due to gaps in the reference, which can also hide the true position to which a read should map. By using a good alignment algorithm (the GATK wrapper of BWA-MEM), giving it sensitive settings (which may have been impractically slow for the original bam alignment) and mapping to the best available reference we can avoid these pitfalls. The last point is especially important: one can (and should) use a BWA-MEM index image corresponding to the best reference, regardless of the reference to which the bam was aligned.|https://gatk.broadinstitute.org/hc/en-us/articles/4418051467035-FilterAlignmentArtifacts-EXPERIMENTAL-|unknown|
|VcfLeftNorm| Left-align and normalize indels; check if REF alleles match the reference; split multiallelic sites into multiple rows; recover multiallelics from multiple rows|source URL for the tool|unknown|
|Haplotyper|This is description of the tool/workflow.|source URL for the tool|unknown|
|IndelsRecalibrator|This is description of the tool/workflow.|source URL for the tool|unknown|
|SNPsRecalibrator|This is description of the tool/workflow.|source URL for the tool|unknown|
### 系统总体结构图
本系统为数据分析流程，流程结构图可能随参数有所变化，一下仅为一个典型示例图
![分析流程的有向无循环(DAG)图](./state.svg "GATK-DNAseq-Workflow")
## 系统运行环境
该系统主要内容是分析流程，需要在Linux操作系统下运行，需安装python > 3.7的环境，环境中需要安装graphviz和xcmds等包. 一般情况下，该环境完全由docker容器提供，包括流程运行每个环节所需要的软件或工具，具体所需软件工具可参考前面章节的”模块列表“.
## 系统运行控制
系统跟随容器一起运行，理论上没有时间限制，主要取决于容器运行周期，容器运行周期又取决于流程本身的具体任务执行情况，分析任务完成后, 容器自动停止。针对系统所包含的分析流程而言，有如下参数控制流程的运行，以方便用户能够根据具体情况分析数据

* help: show this help message and exit
* outdir: 结果目录
* run: 默认不运行流程, 仅生成流程; 如果outdir目录已经存在cmd_state.txt文件，则自动续跑
* docker: 指示是否使用docker镜像创造运行环境
* plot: 该参数用于通知是否对分析流程进行实时绘图，生成的流程图可以显示任务之间的依赖关系,也可以显示每个任务的参数信息，同时还能根据每个节点的颜色判断任务的运行状态以及任务完成需要的时间'
* threads: 允许的最大并行的任务数量, 默认3
* update_args: 输入参数配置文件, json格式. 流程每次运行时都会输出wf.args.json, 由于其中包含一些样本信息参数，不可直接使用，但可以根据此模板修改
* skip: 指定要跳过的步骤或具体task,空格分隔,默认程序会自动跳过依赖他们的步骤, 使用--list_cmd or --list_task可查看候选
* no_skip_depend: 当使用skip参数时, 如果同时指定该参数，则不会自动跳过依赖的步骤
* rerun_steps: 指定需要重跑的步骤，不论其是否已经成功完成，空格分隔, 这样做的可能原因可能是: 重新设置了命令参数. 使用--list_task可查看候选,可使用task的前缀，并且以'*'结尾，将自动匹配符合前缀条件的所有task
* assume_success_steps: 假定哪些步骤已经成功运行，不论其是否真的已经成功完成，空格分隔, 这样做的可能原因: 利用之前已经成功运行的结果(需要把之前的运行结果放到当前结果目录). 使用--list_task可查看候选, 可使用task的前缀，并且以'*'结尾，将自动匹配符合前缀条件的所有task,也可以使用cmd.meata.name
* retry: 某步骤运行失败后再尝试运行的次数, 默认1次. 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini
* list_cmd: 仅仅显示当前流程包含的主步骤, 不会显示指定跳过的步骤
* show_cmd: 提供一个cmd名称,输出该cmd的样例
* list_task: 仅仅显示当前流程包含的详细步骤, 且已经排除指定跳过的步骤
* monitor_resource: 是否监控每一步运行时的资源消耗, 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini
* wait_resource_time: 等待资源的时间上限, 默认每次等待时间为900秒, 等待时间超过这个时间且资源不足时判定任务失败
* no_check_resource_before_run: 指示运行某步骤前检测指定的资源是否足够, 如不足, 则该步骤失败; 如果设置该参数, 则运行前不检查资源. 如需对某一步设置不同的值,可运行前修改pipeline.ini. 如需更改指定的资源, 可在运行流程前修改pipeline.ini
* only_write_docs: 仅仅输出markdown格式的流程说明文档
* fastq_info: A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json]
* r1_name: python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'
* r2_name: python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'
* exclude_samples: samples to exclude from analysis
* pair_info: tumor normal pair info, two-column txt file, first column is tumor sample name. sample not in pair info will be skipped
* ref: reference fasta file
* scatter: scatter number used for interval splitting of variant calling steps
* dbsnp: dbsnp vcf file
* axiomPoly: high confidence known indel vcf file, such as Axiom_Exome_Plus.genotypes.all_populations.poly.vcf.gz. 1,249 individuals are included, drawn from the International HapMap Project and 1000 Genomes Project sample collection. https://www.thermofisher.cn/cn/zh/home/life-science/microarray-analysis/microarray-data-analysis/microarray-analysis-sample-data/axiom-exome-sample-data-set.html
* mills: high confidence known indel vcf file. such as Mills_and_1000G_gold_standard.indels.hg38.vcf
* hapmap: high confidence known snp vcf file. 来自国际人类单倍体型图计划, 这个数据集包含了大量家系数据，并且有非常严格的质控和严密的实验验证，因此它的准确性是目前公认最高的
* omni: high confidence known snp vcf file. 这个数据源自Illumina的Omni基因型芯片，大概2.5百万个位点，它的验证结果常常作为基因型的金标准
* G1000: high confidence known snp vcf file. source from 1000 genomes project
* pon: panel of normal vcf file for germline variant filtering, this will be required for tumor only analysis
* germline_vcf: germline vcf, will be used for germline variant filtering and contamination analysis
* alleles: The set of alleles to force-call regardless of evidence
* contamination_vcf: germline vcf such as small_exac_common_3_b37.vcf, will be used for contamination analysis
* bwaMemIndexImage: bwa-mem-index-mage for artifact alignment filtering. you may created it with tool BwaMemIndexImageCreator with only fasta as input
* vep_cache_dir: VEP cache directory
* vep_plugin_dir: VEP plugin directory
* intervals: interval file, support bed file or picard interval file.
## 模块 fastp
### 模块说明
* 简介: This is description of the tool/workflow.
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/fastp:0.21.0
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### read1
+ 类型: infile
+ 默认值: None
+ 描述: read1 fastq file
#### read2
+ 类型: infile
+ 默认值: None
+ 描述: read2 fastq file
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: other arguments you want to use, such as '-x val'
#### threads
+ 类型: int
+ 默认值: 7
+ 描述: thread number
#### out1
+ 类型: str
+ 默认值: None
+ 描述: clean read1 output fastq file
#### out2
+ 类型: str
+ 默认值: None
+ 描述: clean read2 output fastq file
#### html
+ 类型: str
+ 默认值: None
+ 描述: html report file
#### json
+ 类型: str
+ 默认值: None
+ 描述: json report file
### 模块输出
#### out1
+ type: outfile
+ value: ${{mode:outdir}}/fastp-0113LQCL230110-0/0113LQCL230110-0.clean.R1.fq.gz
+ desc: None
#### out2
+ type: outfile
+ value: ${{mode:outdir}}/fastp-0113LQCL230110-0/0113LQCL230110-0.clean.R2.fq.gz
+ desc: None
#### html
+ type: outfile
+ value: 0113LQCL230110-0.fastp.html
+ desc: None
#### json
+ type: outfile
+ value: 0113LQCL230110-0.fastp.json
+ desc: None
## 模块 FastqToSam
### 模块说明
* 简介: convert fastq to sam
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 5368709120
### 模块输入文件参数
#### read1
+ 类型: infile
+ 默认值: None
+ 描述: read1 fastq file
#### read2
+ 类型: infile
+ 默认值: None
+ 描述: read2 fastq file
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### out
+ 类型: str
+ 默认值: 0113LQCL230110.unmapped.bam
+ 描述: output sam file
#### read_group_name
+ 类型: str
+ 默认值: 0113LQCL230110
+ 描述: read group name
#### sample_name
+ 类型: str
+ 默认值: 0113LQCL230110
+ 描述: sample name
#### library_name
+ 类型: str
+ 默认值: 0113LQCL230110
+ 描述: library name
#### platform
+ 类型: str
+ 默认值: illumina
+ 描述: sequencing platform name
#### tmpdir
+ 类型: str
+ 默认值: .
+ 描述: directorie with space available to be used by this program for temporary storage of working files
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/FastqToSam-0113LQCL230110-0/0113LQCL230110-0.unmapped.bam
+ desc: None
## 模块 uBam2FastqBwaMem
### 模块说明
* 简介: ubam to fastq and then mapping
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 4
* Memory: 5368709120
### 模块输入文件参数
#### ubam
+ 类型: infile
+ 默认值: None
+ 描述: input ubam file
#### ref
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### paired
+ 类型: str
+ 默认值: true
+ 描述: if input is paired fastq, set it be true, else set it be false
#### k
+ 类型: int
+ 默认值: 10000000
+ 描述: This is description of the argument.
#### t
+ 类型: int
+ 默认值: 4
+ 描述: number of threads to use in mapping step
#### out
+ 类型: str
+ 默认值: None
+ 描述: output bam file
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/uBam2FastqBwaMem-0113LQCL230110-0/0113LQCL230110-0.unmerged.bam
+ desc: None
## 模块 MergeBamAlignment
### 模块说明
* 简介: merge bam alignment
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 5368709120
### 模块输入文件参数
#### ALIGNED_BAM
+ 类型: infile
+ 默认值: None
+ 描述: SAM or BAM file
#### UNMAPPED_BAM
+ 类型: infile
+ 默认值: None
+ 描述: unmapped bam file
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### VALIDATION_STRINGENCY
+ 类型: str
+ 默认值: SILENT
+ 描述: This is description of the argument.
#### EXPECTED_ORIENTATIONS
+ 类型: str
+ 默认值: None
+ 描述: This is description of the argument.
#### ATTRIBUTES_TO_RETAIN
+ 类型: str
+ 默认值: X0
+ 描述: This is description of the argument.
#### OUTPUT
+ 类型: str
+ 默认值: 0113LQCL230110-0.merged.unsorted.bam
+ 描述: output bam file
#### SORT_ORDER
+ 类型: str
+ 默认值: "unsorted"
+ 描述: This is description of the argument.
#### CLIP_ADAPTERS
+ 类型: str
+ 默认值: false
+ 描述: Whether to clip adapters where identified.
#### MAX_RECORDS_IN_RAM
+ 类型: int
+ 默认值: 2000000
+ 描述: When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed.
#### MAX_INSERTIONS_OR_DELETIONS
+ 类型: int
+ 默认值: -1
+ 描述: The maximum number of insertions or deletions permitted for an alignment to be included.
#### PRIMARY_ALIGNMENT_STRATEGY
+ 类型: str
+ 默认值: MostDistant
+ 描述: This is description of the argument.
#### UNMAPPED_READ_STRATEGY
+ 类型: str
+ 默认值: COPY_TO_TAG
+ 描述: How to deal with alignment information in reads that are being unmapped (e.g. due to cross-species contamination.) Currently ignored unless UNMAP_CONTAMINANT_READS = true. Note that the DO_NOT_CHANGE strategy will actually reset the cigar and set the mapping quality on unmapped reads since otherwisethe result will be an invalid record. To force no change use the DO_NOT_CHANGE_INVALID strategy.
#### ALIGNER_PROPER_PAIR_FLAGS
+ 类型: str
+ 默认值: true
+ 描述: Use the aligner's idea of what a proper pair is rather than computing in this program.
#### UNMAP_CONTAMINANT_READS
+ 类型: str
+ 默认值: true
+ 描述: Detect reads originating from foreign organisms (e.g. bacterial DNA in a non-bacterial sample),and unmap + label those reads accordingly.
#### PROGRAM_RECORD_ID
+ 类型: str
+ 默认值: "bwamem"
+ 描述: This is description of the argument.
#### PROGRAM_GROUP_VERSION
+ 类型: str
+ 默认值: bwa-mem2-2.2.1_x64-linux
+ 描述: This is description of the argument.
#### PROGRAM_GROUP_COMMAND_LINE
+ 类型: str
+ 默认值: "bwa-mem2 mem -M -Y -p -v 3 -K 10000000 -t 4 ref.fa"
+ 描述: This is description of the argument.
#### PROGRAM_GROUP_NAME
+ 类型: str
+ 默认值: "bwamem"
+ 描述: This is description of the argument.
#### tmpdir
+ 类型: str
+ 默认值: .
+ 描述: directorie with space available to be used by this program for temporary storage of working files
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/MergeBamAlignment-0113LQCL230110-0/0113LQCL230110-0.merged.unsorted.bam
+ desc: None
## 模块 MarkDuplicates
### 模块说明
* 简介: merge bam alignment
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 5368709120
### 模块输入文件参数
#### INPUT
+ 类型: infile
+ 默认值: None
+ 描述: input bam file list
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### OUTPUT
+ 类型: str
+ 默认值: 0113LQCL230110.unsorted.dup_marked.bam
+ 描述: output bam file
#### METRICS_FILE
+ 类型: str
+ 默认值: 0113LQCL230110.dup_metrics.txt
+ 描述: This is description of the argument.
#### VALIDATION_STRINGENCY
+ 类型: str
+ 默认值: SILENT
+ 描述: This is description of the argument.
#### OPTICAL_DUPLICATE_PIXEL_DISTANCE
+ 类型: int
+ 默认值: 2500
+ 描述: The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is moreappropriate. For other platforms and models, users should experiment to find what works best.
#### ASSUME_SORT_ORDER
+ 类型: str
+ 默认值: "queryname"
+ 描述: This is description of the argument.
#### tmpdir
+ 类型: str
+ 默认值: .
+ 描述: directorie with space available to be used by this program for temporary storage of working files
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/MarkDuplicates-0113LQCL230110/0113LQCL230110.unsorted.dup_marked.bam
+ desc: None
## 模块 SortAndFixTags
### 模块说明
* 简介: sort bam and fix tags
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 5368709120
### 模块输入文件参数
#### INPUT
+ 类型: infile
+ 默认值: None
+ 描述: input bam file list
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### SORT_ORDER
+ 类型: str
+ 默认值: "coordinate"
+ 描述: This is description of the argument.
#### CREATE_INDEX
+ 类型: str
+ 默认值: false
+ 描述: This is description of the argument.
#### tmpdir
+ 类型: str
+ 默认值: .
+ 描述: directorie with space available to be used by this program for temporary storage of working files
#### OUTPUT
+ 类型: str
+ 默认值: 0113LQCL230110.sorted.dup_marked.bam
+ 描述: output bam file
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/SortAndFixTags-0113LQCL230110/0113LQCL230110.sorted.dup_marked.bam
+ desc: None
#### out_idx
+ type: outfile
+ value: 0113LQCL230110.sorted.dup_marked.bam.bai
+ desc: None
## 模块 BaseRecalibrator
### 模块说明
* 简介: This is description of the tool/workflow.
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 5368709120
### 模块输入文件参数
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
#### INPUT
+ 类型: infile
+ 默认值: None
+ 描述: input bam file
#### known-sites
+ 类型: infile
+ 默认值: None
+ 描述: One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### use-original-qualities
+ 类型: bool
+ 默认值: True
+ 描述: This is description of the argument.
#### OUTPUT
+ 类型: str
+ 默认值: 0113LQCL230110-0.recal_table
+ 描述: The output recalibration table file to create
#### intervals
+ 类型: str
+ 默认值: None
+ 描述: One or more genomic intervals over which to operate
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/BaseRecalibrator/BaseRecalibrator-0113LQCL230110-0/0113LQCL230110-0.recal_table
+ desc: None
## 模块 GatherBQSRReports
### 模块说明
* 简介: Gathers scattered BQSR recalibration reports into a single file
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 3221225472
### 模块输入文件参数
#### INPUT
+ 类型: infile
+ 默认值: None
+ 描述: input bqsr file list
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### OUTPUT
+ 类型: str
+ 默认值: 0113LQCL230110.recal_table
+ 描述: The output recalibration table file to create
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/GatherBQSRReports-0113LQCL230110/0113LQCL230110.recal_table
+ desc: None
## 模块 ApplyBQSR
### 模块说明
* 简介: Apply Base Quality Score Recalibration (BQSR) model
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 3221225472
### 模块输入文件参数
#### INPUT
+ 类型: infile
+ 默认值: None
+ 描述: input bam file
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
#### bqsr
+ 类型: infile
+ 默认值: None
+ 描述: Input recalibration table for BQSR
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### OUTPUT
+ 类型: str
+ 默认值: 0113LQCL230110-0.recalibrated.bam
+ 描述: The output recalibration table file to create
#### intervals
+ 类型: str
+ 默认值: None
+ 描述: One or more genomic intervals over which to operate
#### static-quantized-quals
+ 类型: int
+ 默认值: [10, 20, 30]
+ 描述: This is description of the argument.
#### add-output-sam-program-record
+ 类型: bool
+ 默认值: True
+ 描述: This is description of the argument.
#### use-original-qualities
+ 类型: bool
+ 默认值: True
+ 描述: This is description of the argument.
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/ApplyBQSR/ApplyBQSR-0113LQCL230110-0/0113LQCL230110-0.recalibrated.bam
+ desc: None
## 模块 GatherBamFiles
### 模块说明
* 简介: Concatenate efficiently BAM files that resulted from a scattered parallel analysis.
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 3221225472
### 模块输入文件参数
#### INPUT
+ 类型: infile
+ 默认值: None
+ 描述: input bam file
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### OUTPUT
+ 类型: str
+ 默认值: 0113LQCL230110.recalibrated.bam
+ 描述: output bam file
#### CREATE_INDEX
+ 类型: str
+ 默认值: true
+ 描述: This is description of the argument.
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/GatherBamFiles-0113LQCL230110/0113LQCL230110.recalibrated.bam
+ desc: None
#### out_idx
+ type: outfile
+ value: 0113LQCL230110.recalibrated.bam.bai
+ desc: None
## 模块 SplitIntervals
### 模块说明
* 简介: This tool takes in intervals via the standard arguments of IntervalArgumentCollection and splits them into interval files for scattering. The resulting files contain equal number of bases.
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 3221225472
### 模块输入文件参数
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
#### intervals
+ 类型: infile
+ 默认值: None
+ 描述: One or more genomic intervals over which to operate
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### scatter
+ 类型: fix
+ 默认值: 10
+ 描述: number of output interval files to split into
#### mode
+ 类型: str
+ 默认值: None
+ 描述: How to divide intervals.
#### interval-merging-rule
+ 类型: str
+ 默认值: OVERLAPPING_ONLY
+ 描述: This is description of the argument.
#### outdir
+ 类型: str
+ 默认值: intervals_folder
+ 描述: The directory into which to write the scattered interval sub-directories.
### 模块输出
#### out0
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0000-scattered.interval_list
+ desc: None
#### out1
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0001-scattered.interval_list
+ desc: None
#### out2
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0002-scattered.interval_list
+ desc: None
#### out3
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0003-scattered.interval_list
+ desc: None
#### out4
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0004-scattered.interval_list
+ desc: None
#### out5
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0005-scattered.interval_list
+ desc: None
#### out6
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0006-scattered.interval_list
+ desc: None
#### out7
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0007-scattered.interval_list
+ desc: None
#### out8
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0008-scattered.interval_list
+ desc: None
#### out9
+ type: outfile
+ value: ${{mode:outdir}}/SplitIntervals-ForCaller/./0009-scattered.interval_list
+ desc: None
## 模块 Mutect2
### 模块说明
* 简介: Call somatic short mutations via local assembly of haplotypes. Short mutations include single nucleotide (SNA) and insertion and deletion (indel) alterations.
* 参考：https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 3
* Memory: 5368709120
### 模块输入文件参数
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
#### tumor_bam
+ 类型: infile
+ 默认值: None
+ 描述: tumor bam
#### normal_bam
+ 类型: infile
+ 默认值: None
+ 描述: normal bam
#### germline-resource
+ 类型: infile
+ 默认值: None
+ 描述: optional database of known germline variants (and its index) (see http://gnomad.broadinstitute.org/downloads)
#### pon
+ 类型: infile
+ 默认值: None
+ 描述: 
#### intervals
+ 类型: infile
+ 默认值: None
+ 描述: One or more genomic intervals over which to operate
#### alleles
+ 类型: infile
+ 默认值: None
+ 描述: The set of alleles to force-call regardless of evidence
#### f1r2-tar-gz
+ 类型: infile
+ 默认值: 0113LQCL230110-0.f1r2.tar.gz
+ 描述: If specified, collect F1R2 counts and output files into this tar.gz file
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### tumor_name
+ 类型: str
+ 默认值: None
+ 描述: tumor sample name
#### normal_name
+ 类型: str
+ 默认值: None
+ 描述: normal sample name
#### string_intervals
+ 类型: str
+ 默认值: None
+ 描述: interval string such as "chrM"
#### out
+ 类型: str
+ 默认值: 0113LQCL230110-0.vcf.gz
+ 描述: output vcf
#### bam_output
+ 类型: str
+ 默认值: None
+ 描述: output bam file
#### mitochondria
+ 类型: bool
+ 默认值: False
+ 描述: if to turn on mitochondria mode. Specifically, the mode sets --initial-tumor-lod to 0, --tumor-lod-to-emit to 0, --af-of-alleles-not-in-resource to 4e-3, and the advanced parameter --pruning-lod-thres
#### tmpdir
+ 类型: str
+ 默认值: .
+ 描述: directorie with space available to be used by this program for temporary storage of working files
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/Mutect2/Mutect2-0113LQCL230110-0/0113LQCL230110-0.vcf.gz
+ desc: None
#### f1r2
+ type: outfile
+ value: ${{mode:outdir}}/Mutect2/Mutect2-0113LQCL230110-0/0113LQCL230110-0.f1r2.tar.gz
+ desc: None
#### stats
+ type: outfile
+ value: ${{mode:outdir}}/Mutect2/Mutect2-0113LQCL230110-0/0113LQCL230110-0.vcf.gz.stats
+ desc: None
## 模块 GetPileupSummaries
### 模块说明
* 简介: Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination
* 参考：https://gatk.broadinstitute.org/hc/en-us/articles/360037593451-GetPileupSummaries
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
#### bam
+ 类型: infile
+ 默认值: None
+ 描述: BAM/SAM/CRAM file containing reads
#### intervals
+ 类型: infile
+ 默认值: None
+ 描述: One or more genomic intervals over which to operate
#### variants_for_contamination
+ 类型: infile
+ 默认值: None
+ 描述: A VCF file containing variants and allele frequencies
#### out
+ 类型: infile
+ 默认值: 0113LQCL230110-0.tumor-pileups.table
+ 描述: The output table
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### interval-set-rule
+ 类型: str
+ 默认值: INTERSECTION
+ 描述: Set merging approach to use for combining interval inputs
### 模块输出
#### out
+ type: outfile
+ value: 0113LQCL230110-0.tumor-pileups.table
+ desc: None
## 模块 LearnReadOrientationModel
### 模块说明
* 简介: Learn the prior probability of read orientation artifact from the output of CollectF1R2Counts of Mutect2
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### inputs
+ 类型: infile
+ 默认值: None
+ 描述: One or more .tar.gz containing outputs of CollectF1R2Counts
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### out
+ 类型: str
+ 默认值: 0113LQCL230110.artifact-priors.tar.gz
+ 描述: tar.gz of artifact prior tables
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/LearnReadOrientationModel-0113LQCL230110/0113LQCL230110.artifact-priors.tar.gz
+ desc: None
## 模块 MergeVcfs
### 模块说明
* 简介: Combines multiple variant files into a single variant file.
* 参考：https://gatk.broadinstitute.org/hc/en-us/articles/360056969852-MergeVcfs-Picard-
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### inputs
+ 类型: infile
+ 默认值: None
+ 描述: input vcf list
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### out
+ 类型: str
+ 默认值: 0113LQCL230110.vcf.gz
+ 描述: The merged VCF or BCF file. File format is determined by file extension.
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/MergeVcfs-0113LQCL230110/0113LQCL230110.vcf.gz
+ desc: None
## 模块 MergeMutectStats
### 模块说明
* 简介: This is description of the tool/workflow.
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### stats
+ 类型: infile
+ 默认值: None
+ 描述: This is description of the argument.
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### out
+ 类型: str
+ 默认值: 0113LQCL230110.vcf.stats
+ 描述: output merged stat files
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/MergeMutectStats-0113LQCL230110/0113LQCL230110.vcf.stats
+ desc: None
## 模块 FilterMutectCalls
### 模块说明
* 简介: FilterMutectCalls applies filters to the raw output of Mutect2
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### vcf
+ 类型: infile
+ 默认值: None
+ 描述: A VCF file containing variants
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
#### contamination-table
+ 类型: infile
+ 默认值: None
+ 描述: This is description of the argument.
#### tumor-segmentation
+ 类型: infile
+ 默认值: None
+ 描述: This is description of the argument.
#### ob-priors
+ 类型: infile
+ 默认值: None
+ 描述: This is description of the argument.
#### stats
+ 类型: infile
+ 默认值: None
+ 描述: This is description of the argument.
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### out
+ 类型: str
+ 默认值: 0113LQCL230110.filtered.vcf.gz
+ 描述: output vcf file
#### filtering-stats
+ 类型: str
+ 默认值: 0113LQCL230110.filtering.stats
+ 描述: output filtering stat file
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/FilterMutectCalls-0113LQCL230110/0113LQCL230110.filtered.vcf.gz
+ desc: None
#### filtering-stats
+ type: outfile
+ value: 0113LQCL230110.filtering.stats
+ desc: None
## 模块 FilterAlignmentArtifacts
### 模块说明
* 简介: Alignment artifacts can occur whenever there is sufficient sequence similarity between two or more regions in the genome to confuse the alignment algorithm. This can occur when the aligner for whatever reason overestimate how uniquely a read maps, thereby assigning it too high of a mapping quality. It can also occur through no fault of the aligner due to gaps in the reference, which can also hide the true position to which a read should map. By using a good alignment algorithm (the GATK wrapper of BWA-MEM), giving it sensitive settings (which may have been impractically slow for the original bam alignment) and mapping to the best available reference we can avoid these pitfalls. The last point is especially important: one can (and should) use a BWA-MEM index image corresponding to the best reference, regardless of the reference to which the bam was aligned.
* 参考：https://gatk.broadinstitute.org/hc/en-us/articles/4418051467035-FilterAlignmentArtifacts-EXPERIMENTAL-
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### vcf
+ 类型: infile
+ 默认值: None
+ 描述: A VCF file containing variants
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
#### bam
+ 类型: infile
+ 默认值: None
+ 描述: input bam file
#### bwa-mem-index-image
+ 类型: infile
+ 默认值: None
+ 描述: BWA-mem index image
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### out
+ 类型: str
+ 默认值: 0113LQCL230110.align_artifacts_filtered.vcf.gz
+ 描述: output vcf file
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/FilterAlignmentArtifacts-0113LQCL230110/0113LQCL230110.align_artifacts_filtered.vcf.gz
+ desc: None
## 模块 VcfLeftNorm
### 模块说明
* 简介:  Left-align and normalize indels; check if REF alleles match the reference; split multiallelic sites into multiple rows; recover multiallelics from multiple rows
* 参考：source URL for the tool
### 模块运行环境
* 镜像：gudeqing/dnaseq:1.0
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### fasta-ref
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
#### vcf
+ 类型: infile
+ 默认值: None
+ 描述: input vcf file
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### multiallelics
+ 类型: str
+ 默认值: None
+ 描述: Split multiallelics (-) or join biallelics (+), type: snps|indels|both|any [both]
#### out
+ 类型: str
+ 默认值: None
+ 描述: Write output to a file [standard output]
#### output-type
+ 类型: str
+ 默认值: v
+ 描述: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]
#### threads
+ 类型: int
+ 默认值: 2
+ 描述: Use multithreading with <int> worker threads [2]
### 模块输出
#### out
+ type: outfile
+ value: 0113LQCL230110.somatic.raw.vcf
+ desc: None
## 模块 Haplotyper
### 模块说明
* 简介: This is description of the tool/workflow.
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### intervals
+ 类型: infile
+ 默认值: None
+ 描述: interval file, support bed file or picard interval or vcf format
#### bam
+ 类型: infile
+ 默认值: None
+ 描述: reccaled tumor and normal bam list
#### REFERENCE_SEQUENCE
+ 类型: infile
+ 默认值: None
+ 描述: reference fasta file
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### emit_mode
+ 类型: str
+ 默认值: GVCF
+ 描述: The reference confidence mode makes it possible to emit a per-bp or summarized confidence estimate for a site being strictly homozygous-reference.
#### ploidy
+ 类型: int
+ 默认值: 2
+ 描述: determines the ploidy number of the sample being processed. The default value is 2.
#### annotation-group
+ 类型: str
+ 默认值: ['StandardAnnotation', 'StandardHCAnnotation', 'AS_StandardAnnotation']
+ 描述: This is description of the argument.
#### gvcf-gq-bands
+ 类型: str
+ 默认值: [10, 20, 30, 40, 50, 60, 70, 80, 90]
+ 描述: Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)
#### tmpdir
+ 类型: str
+ 默认值: .
+ 描述: directorie with space available to be used by this program for temporary storage of working files
#### out_vcf
+ 类型: str
+ 默认值: LQCL230117-0.g.vcf.gz
+ 描述: output vcf file
### 模块输出
#### out
+ type: outfile
+ value: ${{mode:outdir}}/Haplotyper/Haplotyper-LQCL230117-0/LQCL230117-0.g.vcf.gz
+ desc: None
#### out_idx
+ type: outfile
+ value: LQCL230117-0.g.vcf.gz.tbi
+ desc: None
## 模块 IndelsRecalibrator
### 模块说明
* 简介: This is description of the tool/workflow.
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### vcf
+ 类型: infile
+ 默认值: None
+ 描述: site only variant filtered input vcf
#### mills
+ 类型: infile
+ 默认值: None
+ 描述: mills resource vcf
#### axiomPloly
+ 类型: infile
+ 默认值: None
+ 描述: axiomPoly resource vcf
#### dbsnp
+ 类型: infile
+ 默认值: None
+ 描述: dbsnp resource vcf
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### out
+ 类型: str
+ 默认值: Joint.indel.recal
+ 描述: The output recal file used by ApplyVQSR
#### out-tranches-file
+ 类型: str
+ 默认值: Joint.indel.tranches
+ 描述: The output tranches file used by ApplyVQSR
#### tranche
+ 类型: str
+ 默认值: ['100.0', '99.95', '99.9', '99.5', '99.0', '97.0', '96.0', '95.0', '94.0', '93.5', '93.0', '92.0', '91.0', '90.0']
+ 描述: recalibration tranche values
#### use-annotation
+ 类型: str
+ 默认值: ['AS_FS', 'AS_ReadPosRankSum', 'AS_MQRankSum', 'AS_QD', 'AS_SOR']
+ 描述: The names of the annotations which should used for calculations
#### trust-all-polymorphic
+ 类型: bool
+ 默认值: True
+ 描述: This is description of the argument.
#### use-allele-specific-annotations
+ 类型: bool
+ 默认值: True
+ 描述: This is description of the argument.
#### mode
+ 类型: fix
+ 默认值: INDEL
+ 描述: Recalibration mode to employ
#### max-gaussians
+ 类型: int
+ 默认值: 4
+ 描述: Max number of Gaussians for the positive model
### 模块输出
#### out
+ type: outfile
+ value: Joint.indel.recal
+ desc: None
#### tranches
+ type: outfile
+ value: Joint.indel.tranches
+ desc: None
## 模块 SNPsRecalibrator
### 模块说明
* 简介: This is description of the tool/workflow.
* 参考：source URL for the tool
### 模块运行环境
* 镜像：registry-xdp-v3-yifang.xdp.basebit.me/basebitai/gatk:4.2.6.1
* CPU: 2
* Memory: 1024
### 模块输入文件参数
#### vcf
+ 类型: infile
+ 默认值: None
+ 描述: site only variant filtered input vcf
#### hapmap
+ 类型: infile
+ 默认值: None
+ 描述: hapmap resource vcf
#### omni
+ 类型: infile
+ 默认值: None
+ 描述: omni resource vcf
#### 1000G
+ 类型: infile
+ 默认值: None
+ 描述: one thousand genomes resource vcf
#### dbsnp
+ 类型: infile
+ 默认值: None
+ 描述: dbsnp resource vcf
### 模块普通参数
#### other_args
+ 类型: str
+ 默认值: 
+ 描述: This argument is designed to provide any arguments that are not wrapped in Command
#### out
+ 类型: str
+ 默认值: Joint.snp.recal
+ 描述: The output recal file used by ApplyVQSR
#### out-tranches-file
+ 类型: str
+ 默认值: Joint.snp.tranches
+ 描述: The output tranches file used by ApplyVQSR
#### tranche
+ 类型: str
+ 默认值: ['100.0', '99.95', '99.9', '99.8', '99.7', '99.6', '99.5', '99.4', '99.3', '99.0', '98.0', '97.0', '90.0']
+ 描述: recalibration tranche values
#### use-annotation
+ 类型: str
+ 默认值: ['AS_QD', 'AS_MQRankSum', 'AS_ReadPosRankSum', 'AS_FS', 'AS_MQ', 'AS_SOR']
+ 描述: The names of the annotations which should used for calculations
#### trust-all-polymorphic
+ 类型: bool
+ 默认值: True
+ 描述: This is description of the argument.
#### use-allele-specific-annotations
+ 类型: bool
+ 默认值: True
+ 描述: This is description of the argument.
#### mode
+ 类型: fix
+ 默认值: SNP
+ 描述: Recalibration mode to employ
#### max-gaussians
+ 类型: int
+ 默认值: 6
+ 描述: Max number of Gaussians for the positive model
#### output-model
+ 类型: str
+ 默认值: None
+ 描述: If specified, the variant recalibrator will output the VQSR model to this file path.
### 模块输出
#### out
+ type: outfile
+ value: Joint.snp.recal
+ desc: None
#### tranches
+ type: outfile
+ value: Joint.snp.tranches
+ desc: None
#### model
+ type: outfile
+ value: None
+ desc: None
