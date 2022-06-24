## 胚系突变检测流程
* 流程WDL文件：sentieon.germline_pipeline.wdl
* 流程参数数目：
  * 该流程基于sentieon进行胚系突变检测
   



pipeline.getFastqInfo.fastq_dirs	TNpipelineTestData/fastqdir/
pipeline.getFastqInfo.fastq_files	测试时可以空着
pipeline.getFastqInfo.r1_name	= (.*).R1.fq.gz
pipeline.getFastqInfo.r2_name	= (.*).R2.fq.gz
pipeline.pair_info	TNpipelineTestData/fastqdir/pair.info
pipeline.ref	TNpipelineTestData/references/chr17.fa
pipeline.ref_idxes	[TNpipelineTestData/references/chr17.dict	TNpipelineTestData/references/chr17.fa.fai]
pipeline.bwa_idxes	TNpipelineTestData/references/*
pipeline.known_dbsnp	TNpipelineTestData/chr17.dbsnp_146.hg38.vcf.gz
pipeline.known_dbsnp_idx	TNpipelineTestData/chr17.dbsnp_146.hg38.vcf.gz.tbi
pipeline.known_indels	TNpipelineTestData/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
pipeline.known_indels_idx	TNpipelineTestData/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
pipeline.intervals	无，可以空着
pipeline.pon	TNpipelineTestData/chr17_m2pon.vcf.gz
pipeline.pon_idx	TNpipelineTestData/chr17_m2pon.vcf.gz.tbi
pipeline.germline_vcf	TNpipelineTestData/chr17_small_exac_common_3_grch38.vcf.gz
pipeline.germline_vcf_idx	TNpipelineTestData/chr17_small_exac_common_3_grch38.vcf.gz.tbi
pipeline.vep_cache	设置pipeline.skip_vep = true， 然后随便选择一个文件去测试
pipeline.vep_plugins_zip	设置pipeline.skip_vep = true， 然后随便选择一个文件去测试
pipeline.skip_optiType	建议设置为true，因为测试数据无法进行该项分析
