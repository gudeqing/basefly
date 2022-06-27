## 胚系突变检测流程
* 流程WDL文件：sentieon.germline_pipeline.wdl
* 流程参数说明表：germline.args.detail.xlsx
* 该流程使用到的软件：sentieon，VEP

## 体细胞突变检测流程
* 流程WDL文件：sentieon.TN_pipeline.wdl
* 流程参数说明表：TN.args.detail.xlsx
* 该流程使用到的软件：sentieon，VEP

## 流程的输入文件
### 关于fastq文件的输入，主要由下面四个参数完成，因为fastq_dirs和fastq_file均提供了原始数据信息，两个参数必须至少入一个。
pipeline.getFastqInfo.fastq_dirs    原始数据所在目录	TNpipelineTestData/fastqdir/
pipeline.getFastqInfo.fastq_files	具体原始数据文件    测试时可以空着
pipeline.getFastqInfo.r1_name   python正则表达式，第一个小括号匹配到的内容将作为样本名称	= (.*).R1.fq.gz
pipeline.getFastqInfo.r2_name   python正则表达式，第一个小括号匹配到的内容将作为样本名称	= (.*).R2.fq.gz
### 基于配对数据进行体细胞
pipeline.pair_info  文件，第一列为肿瘤样本名，第二列为对照样本名，必须和上述python正则表达式匹配到的名称一致	TNpipelineTestData/fastqdir/pair.info
### 参考基因组信息
pipeline.ref    参考基因组fasta文件	TNpipelineTestData/references/chr17.fa
pipeline.ref_idxes  参考基因组fasta文件对应的索引文件，有2个文件	[TNpipelineTestData/references/chr17.dict	TNpipelineTestData/references/chr17.fa.fai]
pipeline.bwa_idxes  参考基因组对应的bwa索引，有至少5个文件	TNpipelineTestData/references/*
### 捕获区间文件，如果无则不填，如果不是WGS测序，建议必须填写
pipeline.intervals	捕获区间文件  测试时可以空着
### 数据库文件，包括dbsnp数据库和其他已知indel突变数据库，将用于突变质量矫正和注释，下面测试举例
pipeline.known_dbsnp	TNpipelineTestData/chr17.dbsnp_146.hg38.vcf.gz
pipeline.known_dbsnp_idx	TNpipelineTestData/chr17.dbsnp_146.hg38.vcf.gz.tbi
pipeline.known_indels   可以输入多个文件	TNpipelineTestData/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
pipeline.known_indels_idx   可以输入多个文件	TNpipelineTestData/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
pipeline.pon    测试配对样本时，建议空着	TNpipelineTestData/chr17_m2pon.vcf.gz
pipeline.pon_idx    测试配对样本时，建议空着	TNpipelineTestData/chr17_m2pon.vcf.gz.tbi
pipeline.germline_vcf	TNpipelineTestData/chr17_small_exac_common_3_grch38.vcf.gz
pipeline.germline_vcf_idx	TNpipelineTestData/chr17_small_exac_common_3_grch38.vcf.gz.tbi

### 突变注释软件VEP的输入，因为WDL流程语法问题，下面的输入必须要写，即使选择skip_vep = true
pipeline.vep_cache	设置pipeline.skip_vep = true， 然后随便选择一个文件去测试
pipeline.vep_plugins_zip	设置pipeline.skip_vep = true， 然后随便选择一个文件去测试

### HLA基因定型，可选步骤
pipeline.skip_optiType	建议设置为true，因为测试数据无法进行该项分析
