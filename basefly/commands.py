from .basefly import Argument, Output, Command, TmpVar, TopVar

"""
Command开发形式或步骤：
1. 初始化command对象：cmd = Command()
2. 添加meta信息：cmd.meta.{name, desc, source, version, author, any_other_attr}
3. 添加runtime信息： cmd.runtime.{tool,tool_dir, image, memory, cpu, timeout, max_memory, max_cpu} 
4. 定义参数：通过Argument对象来定义参数，需按照使用顺序定义; 用13属性描述一个参数，每个属性都有默认值, 如cmd.args[arg_name] = Argument(...)
5. 定义outputs：需通过Ouput对象，使用value和type属性描述，value中可以使用‘{}'引用之前的参数，如value = ’{prefix}.txt'
6. return cmd

注意:
pycharm里设置code completion 允许“suggest variable and parameter name”, 可以极大方便流程编写
0. 写workflow时，参数赋值规范建议：args[X].value = TopVar[?] | task.Outputs[?] | TmpVar()
*. 如果不是为了写wdl流程，可以不使用TmpVar，直接赋值就ok
1. 一定要正确定义参数的类型, 尤其是输入文件对应的参数类型。 type is one of ['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix']
    其中‘fix'可以用于表示命令行中的固定字符串或固定参数, 如 “bwa xxx | samtools sort -" 中的‘| samtools sort -’ 可以用fix固定
2. 参数的添加顺序对于命令行的正确形成很重要，这里字典的有序性得到利用
3. 定义output时，value(或path）属性对应的值可以直接用{}引用cmd.args的key，

关于runtime:
memory和cpu是定义最小计算资源需求
max_memory和max_cpu定义计算资源上限
image: 定义docker镜像
tool：工具命令，即命令行的第一个参数
tool_dir: 定义tool所在路径

关于meta：
name: 定义命令行的名称，会参与具体task的name的形成，建议组成：[数字，字母，’-‘], 下划线会自动被替换为中划线’-‘
其他字段都是描述工具的开发作者(author)，链接(source)，版本号(version)，简介（desc)

* 下面是已经写好的常用生信command, 写生信分析流程时可以直接调用, 供参考，欢迎补充
"""


def fastp(sample):
    cmd = Command()
    cmd.meta.name = 'fastp'
    # cmd.runtime.image = 'gudeqing/fastp:0.21.0'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.tool = 'fastp'
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', desc='read2 fastq file')
    cmd.args['threads'] = Argument(prefix='-w ', default=7, desc='thread number')
    cmd.args['other_args'] = Argument(prefix='', default='', desc="other arguments you want to use, such as '-x val'")
    cmd.args['out1'] = Argument(prefix='-o ', value=TmpVar(value=f'{sample}.clean.R1.fq.gz', name='~{sample}.clean.R1.fq.gz'), type='str', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', value=TmpVar(value=f'{sample}.clean.R2.fq.gz', name='~{sample}.clean.R2.fq.gz'), type='str', desc='clean read2 output fastq file')
    cmd.args['html'] = Argument(prefix='-h ', value=TmpVar(value=f'{sample}.fastp.html', name='~{sample}.fastp.html'), type='str', desc='html report file')
    cmd.args['json'] = Argument(prefix='-j ', value=TmpVar(value=f'{sample}.fastp.json', name='~{sample}.fastp.json') , type='str', desc='html report file')
    # 下面的outputs设置起初是为了能够生成wdl设置,
    cmd.outputs['out1'] = Output(value="{out1}", type='outfile')  # 这里使用”{}“引用其他Argument对象作为输入
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}")
    cmd.outputs['json'] = Output(value="{json}")
    return cmd


def bwa_mem(sample, platform):
    cmd = Command()
    cmd.meta.name = 'bwa_mem'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon bwa mem -M'
    cmd.args['readgroup'] = Argument(prefix='-R ', desc='read group info', value=f'"@RG\\tID:{sample}\\tSM:{sample}\\tPL:{platform}"')
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['k'] = Argument(prefix='-K ', default=10000000)
    cmd.args['ref'] = Argument(type='infile', desc='reference fasta file')
    cmd.args['read1'] = Argument(type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(type='infile', desc='read2 fastq file')
    cmd.args['_x2'] = Argument(type='fix', value=' | sentieon util sort')
    cmd.args['ref2'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['out'] = Argument(prefix='-o ', desc='output bam file', value=f'{sample}.sorted.bam')
    cmd.args['t2'] = Argument(prefix='-t ', default=16, desc='number of threads to use')
    cmd.args['_x3'] = Argument(type='fix', value='--sam2bam -i -')
    cmd.outputs['out'] = Output(value="{out}", type='outfile')
    return cmd


def get_metrics(sample):
    cmd = Command()
    cmd.meta.name = 'get_metrics'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['intervals'] = Argument(prefix='--interval ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile',  desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile',  desc='input bam file')
    cmd.args['mq_metrics'] = Argument(prefix='--algo MeanQualityByCycle ', value=f'{sample}.mq_metrics.txt', desc='metric file of MeanQualityByCycle')
    cmd.args['qd_metrics'] = Argument(prefix='--algo QualDistribution ', value=f'{sample}.qd_metrics.txt', desc='metric file of QualDistribution')
    cmd.args['gc_summary'] = Argument(prefix='--algo GCBias --summary ', value=f'{sample}.gc_summary.txt', desc='summary file of GCBias')
    cmd.args['gc_metrics'] = Argument(desc='metrics file of GCBias', value=f'{sample}.gc_metrics.txt')
    cmd.args['aln_metrics'] = Argument(prefix='--algo AlignmentStat ', value=f'{sample}.aln_metrics.txt', desc='aln_metrics file of AlignmentStat')
    cmd.args['insert_metrics'] = Argument(prefix='--algo InsertSizeMetricAlgo ', value=f'{sample}.insert_metrics.txt', desc='insert_metrics file of InsertSizeMetricAlgo')
    cmd.outputs['mq_metrics'] = Output(value='{mq_metrics}')
    cmd.outputs['qd_metrics'] = Output(value='{qd_metrics}')
    cmd.outputs['gc_summary'] = Output(value='{gc_summary}')
    cmd.outputs['gc_metrics'] = Output(value='{gc_metrics}')
    cmd.outputs['aln_metrics'] = Output(value='{aln_metrics}')
    cmd.outputs['insert_metrics'] = Output(value='{insert_metrics}')
    return cmd


def plot_metrics(sample, method='GCBias'):
    cmd = Command()
    cmd.meta.name = f'plot{method}'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon plot'
    cmd.args['method'] = Argument(desc='method of plot', default=method)
    cmd.args['out'] = Argument(prefix='-o ', desc='plot file', value=f'{sample}.{method}.pdf')
    cmd.args['i'] = Argument(type='infile', desc='input metrics file for plot')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def locus_collector(sample):
    cmd = Command()
    cmd.meta.name = 'LocusCollector'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args['score'] = Argument(prefix='--algo LocusCollector --fun score_info ', desc='output score file', value=f'{sample}.score.txt')
    cmd.outputs['score'] = Output(value='{score}')
    return cmd


def dedup(sample):
    cmd = Command()
    cmd.meta.name = 'DeDup'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args['_x'] = Argument(type='fix', value='--algo Dedup')
    cmd.args['score'] = Argument(prefix='--score_info ', type='infile', desc='score info file')
    cmd.args['dedup_metrics'] = Argument(prefix='--metrics ', desc='output metrics info file', value=f'{sample}.dedup.metrics.txt')
    cmd.args['deduped_bam'] = Argument(desc='output metrics info file', value=f'{sample}.deduped.bam')
    cmd.outputs['dedup_metrics'] = Output(value='{dedup_metrics}')
    cmd.outputs['out_bam'] = Output(value='{deduped_bam}')
    return cmd


def coverage_metrics(sample):
    cmd = Command()
    cmd.meta.name = 'CoverageMetrics'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['intervals'] = Argument(prefix='--interval ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args['coverage_metrics'] = Argument(prefix='--algo CoverageMetrics ', value=f'{sample}.cov.metrics.txt', desc='output coverage metrics file')
    cmd.outputs['coverage_metrics'] = Output(value="{coverage_metrics}")
    return cmd


def realign(sample):
    cmd = Command()
    cmd.meta.name = 'realign'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args._x = Argument(type='fix', value='--algo Realigner')
    cmd.args['database'] = Argument(prefix='-k ', type='infile', multi_times=True, desc='known indel vcf file')
    cmd.args['realigned_bam'] = Argument(desc='output realigned bam file', value=f'{sample}.realigned.bam')
    cmd.outputs['out_bam'] = Output(value='{realigned_bam}')
    return cmd


def recalibration(sample):
    cmd = Command()
    cmd.meta.name = 'recalibration'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['intervals'] = Argument(prefix='--interval ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args['_x'] = Argument(type='fix', value='--algo QualCal')
    cmd.args['database'] = Argument(prefix='-k ', type='infile', multi_times=True, desc='known indel vcf file')
    cmd.args['recal_data'] = Argument(desc="output recal_data.table", value=f'{sample}.recal_data.table')
    cmd.outputs['recal_data'] = Output(value='{recal_data}')
    return cmd


def TNhaplotyper2(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'TNhaplotyper2'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=16, desc='number of threads to use in computation, set to number of cores in the server')
    cmd.args['intervals'] = Argument(prefix='--interval ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    # basic inputs
    cmd.args['bams'] = Argument(prefix='-i ', type='infile', multi_times=True, desc='reccaled tumor and normal bam list')
    cmd.args['recal_datas'] = Argument(prefix='-q ', type='infile', multi_times=True, desc='tumor and normal recal data list')
    cmd.args['method'] = Argument(prefix='--algo ', value='TNhaplotyper2', type='fix')
    cmd.args['tumor_sample'] = Argument(prefix='--tumor_sample ', desc='tumor sample name', default=tumor_sample)
    cmd.args['normal_sample'] = Argument(prefix='--normal_sample ', desc='normal sample name', level='optional')
    # optional inputs
    cmd.args['germline_vcf'] = Argument(prefix='--germline_vcf ', type='infile', level='optional', desc='the location of the population germline resource')
    cmd.args['pon'] = Argument(prefix='--pon ', type='infile', level='optional', desc='the location and name of panel of normal VCF file')
    cmd.args['out_vcf'] = Argument(value=f'{tumor_sample}.TNhaplotyper2.vcf.gz', desc='output vcf file of TNhaplotyper2, this will be used later for filtering')
    # orientation
    cmd.args['orientation_sample'] = Argument(prefix='--algo OrientationBias --tumor_sample ', level='optional', desc='tumor sample name')
    cmd.args['orientation_data'] = Argument(level='optional', default=f'{tumor_sample}.orientation.data', desc='output orientation bias result file')
    # contamination, 如果无对照样本或者germline vcf，则无该项分析
    cmd.args['contamination_tumor'] = Argument(prefix="--algo ContaminationModel --tumor_sample ", level='optional', desc='tumor sample name', default=tumor_sample)
    cmd.args['contamination_normal'] = Argument(prefix="--normal_sample ", level='optional', desc='normal sample name')
    cmd.args['germline_vcf2'] = Argument(prefix='--vcf ', type='infile', level='optional', desc='the location of the population germline resource')
    cmd.args['tumor_segments'] = Argument(prefix='--tumor_segments ', level='optional', default=f'{tumor_sample}.contamination.segments', desc='output file name of the file containing the tumor segments information produced by ContaminationModel')
    cmd.args['contamination_data'] = Argument(level='optional', default=f'{tumor_sample}.contamination.data', desc='output file containing the contamination information produced by ContaminationModel')
    cmd.outputs['out_vcf'] = Output(value='{out_vcf}')
    cmd.outputs['orientation_data'] = Output(value='{orientation_data}')
    cmd.outputs['tumor_segments'] = Output(value='{tumor_segments}')
    cmd.outputs['contamination_data'] = Output(value='{contamination_data}')
    return cmd


def TNfilter(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'TNfilter'
    cmd.runtime.image = 'docker-reg.basebit.me:5000/pipelines/sentieon-joint-call:2019.11'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['method'] = Argument(type='fix', value='--algo TNfilter')
    cmd.args['tumor_sample'] = Argument(prefix='--tumor_sample ', desc='tumor sample name', default=tumor_sample)
    cmd.args['normal_sample'] = Argument(prefix='--normal_sample ', level='optional', desc='normal sample name')
    cmd.args['tmp_vcf'] = Argument(prefix='-v ', type='infile', desc='vcf file from TNhaplotyper2')
    cmd.args['contamination'] = Argument(prefix='--contamination ', type='infile', level='optional', desc='file containing the contamination information produced by ContaminationModel')
    cmd.args['tumor_segments'] = Argument(prefix='--tumor_segments ', type='infile', level='optional', desc='file containing the tumor segments information produced by ContaminationModel')
    cmd.args['orientation_data'] = Argument(prefix='--orientation_priors ', type='infile', level='optional', desc='file containing the orientation bias information produced by OrientationBias')
    cmd.args['out_vcf'] = Argument(desc='final output vcf', value=f'{tumor_sample}.somatic.vcf.gz')
    cmd.outputs['out_vcf'] = Output(value='{out_vcf}')
    return cmd


def Haplotyper(sample):
    cmd = Command()
    cmd.meta.name = 'Haplotyper'
    cmd.runtime.image = 'registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['intervals'] = Argument(prefix='--interval ', level='optional', type='infile', multi_times=True, desc="interval file, support bed file or picard interval or vcf format")
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='reccaled tumor and normal bam list')
    cmd.args['recal_data'] = Argument(prefix='-q ', type='infile', desc='tumor and normal recal data list')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['method'] = Argument(type='fix', value='--algo Haplotyper')
    cmd.args['emit_mode'] = Argument(prefix='--emit_mode ', level='optional', desc='determines what calls will be emitted. possible values:variant,confident,all,gvcf')
    cmd.args['ploidy'] = Argument(prefix='--ploidy ', type='int', default=2, desc='determines the ploidy number of the sample being processed. The default value is 2.')
    cmd.args['trim_soft_clip'] = Argument(prefix='--trim_soft_clip', type='bool', default=False, desc='used for rnaseq variant calling')
    cmd.args['call_conf'] = Argument(prefix='--call_conf ', default=30, desc='determine the threshold of variant quality to call a variant. Recommend 20 for rnaseq')
    cmd.args['emit_conf'] = Argument(prefix='--emit_conf ', default=30, desc='determine the threshold of variant quality to emit a variant.')
    cmd.args['out_vcf'] = Argument(value=f'{sample}.vcf.gz', desc='output vcf file')
    cmd.outputs['out_vcf'] = Output(value='{out_vcf}')
    cmd.outputs['out_vcf_idx'] = Output(value='{out_vcf}.tbi')
    return cmd


def GVCFtyper(normal_sample):
    cmd = Command()
    cmd.meta.name = 'GVCFtyper'
    cmd.runtime.image = 'registry-xdp-v3-yifang.xdp.basebit.me/basebitai/sentieon:202010.02'
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['method'] = Argument(type='fix', value='--algo GVCFtyper')
    cmd.args['in_gvcf'] = Argument(prefix='-v ', type='infile', multi_times=True,  desc='input gvcf file')
    cmd.args['known_dbsnp'] = Argument(prefix='-d ', type='infile', desc='dbsnp file')
    cmd.args['call_conf'] = Argument(prefix='--call_conf ', type='int', default=30, desc="determine the threshold of variant quality to emit a variant. Variants with quality less than CONFIDENCE will be not be added to the output VCF file.")
    cmd.args['genotype_model'] = Argument(prefix='--genotype_model ', range={"coalescent", "multinomial"}, default='multinomial', desc="determines which model to use for genotyping and QUAL calculation")
    cmd.args['out_vcf'] = Argument(value=f'{normal_sample}.germline.vcf.gz', desc='output vcf file')
    cmd.outputs['out_vcf'] = Output(value='{out_vcf}')
    cmd.outputs['out_vcf_idx'] = Output(value='{out_vcf}.tbi')
    return cmd


def vep(sample):
    cmd = Command()
    cmd.meta.name = 'VEP'
    cmd.runtime.image = 'ensemblorg/ensembl-vep:2.0.3'
    cmd.runtime.tool = 'vep'
    cmd.args['input_file'] = Argument(prefix='-i ', type='infile', desc='input file')
    cmd.args['fasta'] = Argument(prefix='--fasta ', type='infile', desc="Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence. The first time you run VEP with this parameter an index will be built which can take a few minutes. This is required if fetching HGVS annotations (--hgvs) or checking reference sequences (--check_ref) in offline mode (--offline), and optional with some performance increase in cache mode (--cache).")
    cmd.args['output_file'] = Argument(prefix='-o ', default=f'{sample}.vep.vcf.gz', desc='output file')
    cmd.args['output_format'] = Argument(prefix='--', range={'vcf', 'json', 'tab'}, default='vcf', desc="If we choose to write output in VCF format. Consequences are added in the INFO field of the VCF file, using the key 'CSQ'. Data fields are encoded separated by '|'; the order of fields is written in the VCF header. Output fields in the 'CSQ' INFO field can be selected by using --fields.")
    cmd.args['compress_output'] = Argument(prefix='--compress_output ', default='bgzip', desc="Writes output compressed using either gzip or bgzip")
    cmd.args['force_overwrite'] = Argument(prefix="--force_overwrite ", type='bool', default=True, desc="Force overwriting of output file")
    cmd.args['fork'] = Argument(prefix='--fork ', type='int', default=7, desc='Use forking(multi-cpu/threads) to improve script runtime')
    cmd.args['species'] = Argument(prefix='--species ', default='homo_sapiens', desc='Species for your data. This can be the latin name e.g. homo_sapiens or any Ensembl alias e.g. mouse.')
    cmd.args['assembly_version'] = Argument(prefix='--assembly ', default='GRCh37', desc='Select the assembly version to use if more than one available.')
    cmd.args['dir_cache'] = Argument(prefix='--dir_cache ', type='indir', desc='Specify the cache directory to use')
    cmd.args['dir_plugins'] = Argument(prefix='--dir_plugins ', type='indir', desc='Specify the plugin directory to use')
    cmd.args['stats_file'] = Argument(prefix='--stats_file ', default=f'{sample}.vep.summary.html', desc='Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end with <.html>.')
    cmd.args['cache'] = Argument(prefix='--cache ', type='bool', default=True, desc='Enables use of cache')
    cmd.args['offline'] = Argument(prefix='--offline ', type='bool', default=True, desc='Enables offline mode. No database connections, and a cache file or GFF/GTF file is required for annotation')
    cmd.args['merged'] = Argument(prefix='--merged ', type='bool', default=False, desc='Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used.')
    cmd.args['plugins'] = Argument(prefix='--plugin ', multi_times=True, default=['Frameshift', 'Wildtype'], desc='Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory.Multiple plugins can be used by supplying the --plugin flag multiple times')
    cmd.args['variant_class'] = Argument(prefix='--variant_class ', type='bool', default=True, desc='Output the Sequence Ontology variant class.')
    cmd.args['sift'] = Argument(prefix='--sift ', default='b', range={'p', 's', 'b'}, desc="Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both.")
    cmd.args['polyphen'] = Argument(prefix='--polyphen ', default='b', range={'p', 's', 'b'}, desc="Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both.")
    cmd.args['nearest'] = Argument(prefix='--nearest ', default='transcript', range={'transcript', 'gene', 'symbol'}, desc='Retrieve the transcript or gene with the nearest protein-coding transcription start site (TSS) to each input variant. Use transcript to retrieve the transcript stable ID, gene to retrieve the gene stable ID, or symbol to retrieve the gene symbol. Note that the nearest TSS may not belong to a transcript that overlaps the input variant, and more than one may be reported in the case where two are equidistant from the input coordinates.')
    cmd.args['gene_phenotype'] = Argument(prefix='--gene_phenotype ', type='bool', default=True, desc='Indicates if the overlapped gene is associated with a phenotype, disease or trait.')
    cmd.args['regulatory'] = Argument(prefix='--regulatory ', type='bool', default=True, desc="Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature.")
    cmd.args['phased'] = Argument(prefix='--phased ', type='bool', default=True, desc="Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data.")
    cmd.args['numbers'] = Argument(prefix='--numbers ', type='bool', default=True, desc="Adds affected exon and intron numbering to to output. Format is Number/Total")
    cmd.args['hgvs'] = Argument(prefix='--hgvs ',  type='bool', default=True, desc="Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate.")
    cmd.args['transcript_version'] = Argument(prefix='--transcript_version ', type='bool', default=True, desc="Add version numbers to Ensembl transcript identifiers")
    cmd.args['gencode_basic'] = Argument(prefix='--gencode_basic ', type='bool', default=True, desc="Limit your analysis to transcripts belonging to the GENCODE basic set. This set has fragmented or problematic transcripts removed. ")
    cmd.args['symbol'] = Argument(prefix='--symbol ', type='bool', default=True, desc="Adds the gene symbol (e.g. HGNC) (where available) to the output.")
    cmd.args['tsl'] = Argument(prefix='--tsl ', type='bool', default=True, desc="Adds the transcript support level for this transcript to the output.")
    cmd.args['canonical'] = Argument(prefix='--canonical ', type='bool', default=True, desc="Adds a flag indicating if the transcript is the canonical transcript for the gene")
    cmd.args['biotype'] = Argument(prefix='--biotype ', type='bool', default=True, desc="Adds the biotype of the transcript or regulatory feature.")
    cmd.args['max_af'] = Argument(prefix='--max_af ', type='bool', default=True, desc="Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD")
    cmd.args['af_1kg'] = Argument(prefix='--af_1kg ', type='bool', default=True, desc="Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output.")
    cmd.args['af_gnomad'] = Argument(prefix='--af_gnomad ', type='bool', default=True, desc="Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included")
    cmd.args['af_esp'] = Argument(prefix='--af_esp ', type='bool', default=False, desc="Include allele frequency from NHLBI-ESP populations.")
    cmd.args['coding_only'] = Argument(prefix='--af_esp ', type='bool', default=False, desc="Only return consequences that fall in the coding regions of transcripts. Not used by default")
    cmd.args['pick'] = Argument(prefix='--pick', type='bool', default=False, desc="Pick one line or block of consequence data per variant, including transcript-specific columns. This is the best method to use if you are interested only in one consequence per variant")
    cmd.args['flag_pick'] = Argument(prefix='--flag_pick ', type='bool', default=True, desc="As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others.")
    cmd.args['filter_common'] = Argument(prefix='--filter_common ', type='bool', default=False, desc="Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters.")
    cmd.args['other_args'] = Argument(default='', desc='specify other arguments that you want to append to the command')
    cmd.args['_create_index'] = Argument(value='&& tabix *vcf.gz', type='fix')
    cmd.outputs['out_vcf'] = Output(value='{output_file}')
    cmd.outputs['out_vcf_idx'] = Output(value='{output_file}.tbi')
    return cmd


def CombineVariants(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'CombineVariants'
    cmd.meta.desc = "Combine variants"
    cmd.runtime.image = 'broadinstitute/gatk3:3.8-1'
    cmd.runtime.tool = 'java -Xmx10g -jar GenomeAnalysisTK.jar -T CombineVariants'
    cmd.args['ref'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['variant'] = Argument(prefix='--variant ', multi_times=True, type='infile', desc='variant vcf file array')
    cmd.args['out_vcf'] = Argument(prefix='-o ', value=f'{tumor_sample}.combined_germline.vcf')
    cmd.args['assumeIdenticalSamples'] = Argument(prefix='--assumeIdenticalSamples', type='bool', desc='If true, assume input VCFs have identical sample sets and disjoint calls. This option allows the user to perform a simple merge (concatenation) to combine the VCFs.')
    cmd.outputs['combined_vcf'] = Output(value='{out_vcf}')
    return cmd


def SortVcf(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'SortVcf'
    cmd.meta.desc = "sort vcf"
    cmd.runtime.image = 'broadinstitute/picard:latest'
    cmd.runtime.tool = 'java -jar /usr/picard/picard.jar SortVcf'
    cmd.args['in_vcf'] = Argument(prefix='I=', type='infile', desc='input vcf to sort')
    cmd.args['out_vcf'] = Argument(prefix='O=', value=f'{tumor_sample}.combined_germline.sorted.vcf', type='infile', desc='output sorted vcf')
    cmd.outputs['sorted_vcf'] = Output(value='{out_vcf}')
    return cmd


def ReadBackedPhasing(tumor_sample):
    cmd = Command()
    cmd.meta.name = 'ReadBackedPhasing'
    cmd.meta.desc = "ReadBackedPhasing"
    cmd.runtime.image = 'broadinstitute/gatk3:3.8-1'
    cmd.runtime.tool = 'java -Xmx10g -jar GenomeAnalysisTK.jar -T ReadBackedPhasing'
    cmd.args['ref'] = Argument(prefix='-R ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-I ', type='infile', desc='tumor bam file')
    cmd.args['variant'] = Argument(prefix='--variant ', type='infile', desc='input vcf file')
    cmd.args['interval'] = Argument(prefix='-L ', type='infile', desc='input vcf file')
    cmd.args['out_vcf'] = Argument(prefix='-o ', value=f'{tumor_sample}.phased.vcf', desc='output vcf file')
    cmd.outputs['phased_vcf'] = Output(value='{out_vcf}')
    return cmd


def HLA_ABC_typer(sample):
    cmd = Command()
    cmd.meta.name = 'OptiType'
    cmd.meta.desc = "OptiType:4-digit HLA typer"
    cmd.runtime.image = 'fred2/optitype:1.3.1'
    cmd.runtime.tool = 'OptiTypePipeline.py'
    cmd.args['reads'] = Argument(prefix='--input ', type='infile', array=True, desc='fastq file(s) (fished or raw) or .bam files stored for re-use, generated by an earlier OptiType run.')
    cmd.args['is_dna'] = Argument(prefix='--dna', type='bool', default=True, desc='use with DNA sequencing data')
    cmd.args['is_rna'] = Argument(prefix='--rna', type='bool', default=False, desc='use with RNA sequencing data')
    cmd.args['enumerate'] = Argument(prefix='--enumerate ', type='int', default=1, desc='Number of enumerations. OptiType will output the optimal solution and the top N-1 suboptimal solutions in the results CSV.')
    cmd.args['outdir'] = Argument(prefix='--outdir ', default='.', desc='Specifies the out directory to which all files should be written.')
    cmd.args['prefix'] = Argument(prefix='--prefix ', value=sample, desc='prefix of output files')
    cmd.args['config'] = Argument(prefix='--config ', default='config.ini', desc='config.ini file')
    cmd.outputs['result_tsv'] = Output(value='{prefix}_result.tsv')
    cmd.outputs['result_pdf'] = Output(value='{prefix}_coverage_plot.pdf')
    return cmd


def hisat_genotype():
    """
    hisatgenotype --base hla --locus-list A,B,C,DRB1,DQA1 -1 ILMN/NA12892.extracted.1.fq.gz -2 ILMN/NA12892.extracted.2.fq.gz
    hisatgenotype_toolkit parse-results --csv --in-dir hisatgenotype_out
    """
    cmd = Command()
    cmd.meta.name = 'HisatGenotype'
    cmd.meta.desc = " HLA-typing using hisat"
    cmd.runtime.image = ''
    cmd.runtime.tool = 'hisatgenotype'
    cmd.runtime.cpu = 6
    cmd.runtime.memory = 8*1024**3
    cmd.args['base'] = Argument(prefix='--base ', default='hla', desc='Base file name for index, variants, haplotypes, etc. (e.g. hla, rbg, codis)')
    cmd.args['locus'] = Argument(prefix='--locus-list ', level='optional', array=True, delimiter=',', desc='A comma-separated list of gene names (default: empty, all genes)')
    cmd.args['read1'] = Argument(prefix='-1 ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-2 ', type='infile', desc='read2 fastq file')
    cmd.args['_read_dir'] = Argument(prefix='--in-dir ', value='/', type='fix')
    cmd.args['threads'] = Argument(prefix='--threads ', default=5, desc='Number of threads')
    cmd.args['hisat_threads'] = Argument(prefix='--pp ', default=7, desc='Number of threads')
    cmd.args['indicies'] = Argument(prefix='--index_dir ', level='optional', type='indir', desc="Set location to use for indicies")
    cmd.args['_outdir'] = Argument(prefix='--out-dir ', value='./', type='fix')
    cmd.args['_parse_result'] = Argument(value='&& hisatgenotype_toolkit parse-results --csv --in-dir .', type='fix')
    cmd.args['level'] = Argument(prefix='-t ', default=3, desc='Trim allele to specific field level (example : A*01:01:01:01 trim 2 A*01:01)')
    cmd.args['out'] = Argument(prefix='--output-file ', desc='output of csv file')
    cmd.outputs['out'] = Output(value='{out}', type='outfile')
    return cmd


def mutscan():
    cmd = Command()
    cmd.meta.name = 'mutscan'
    cmd.runtime.image = 'gudeqing/neoantigen:1.0'
    cmd.runtime.tool = '/opt/mutscan'
    cmd.args['read1'] = Argument(prefix='-1 ', type='infile', desc='read1 file name')
    cmd.args['read2'] = Argument(prefix='-2 ', type='infile', desc='read2 file name')
    cmd.args['vcf'] = Argument(prefix='-m ', type='infile', desc='mutation file name, can be a CSV format or a VCF format')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file name (only needed when mutation file is a VCF)')
    cmd.args['html'] = Argument(prefix='-h ', desc='html report')
    cmd.args['json'] = Argument(prefix='-j ', desc='json report')
    cmd.args['threads'] = Argument(prefix='-t ', default=4, desc='worker thread number, default is 4')
    cmd.args['support'] = Argument(prefix='-S ', default=3, desc='min read support required to report a mutation')
    cmd.args['mark'] = Argument(prefix='-k ', level='optional', desc='when mutation file is a vcf file, --mark means only process the records with FILTER column is M')
    cmd.outputs['html'] = Output(value='{html}')
    cmd.outputs['json'] = Output(value='{json}')
    return cmd


def add_exp_to_vcf():
    cmd = Command()
    cmd.meta.name = 'vcf-expression-annotator'
    cmd.runtime.image = 'gudeqing/neoantigen:1.0'
    cmd.runtime.tool = 'vcf-expression-annotator'
    cmd.args['id-column'] = Argument(prefix='--id-column ', default='Name', desc='The column header in the expression_file for the column containing gene names/transcript ids.')
    cmd.args['expression-column'] = Argument(prefix='--expression-column ',  default='TPM', desc='The column header in the expression_file for the column containing expression data.')
    cmd.args['sample-name'] = Argument(prefix='--sample-name ',  desc='If the input_vcf contains multiple samples, the name of the sample to annotate')
    cmd.args['output-vcf'] = Argument(prefix='--output-vcf ', desc='Path to write the output VCF file.')
    cmd.args['ignore-ensembl-id-version'] = Argument(prefix='--ignore-ensembl-id-version', type='bool', default=True, desc='Assumes that the final period and number denotes the Ensembl ID version and ignores it ')
    cmd.args['input-vcf'] = Argument(prefix='', type='infile', desc='A VEP-annotated VCF file')
    cmd.args['expression-file'] = Argument(prefix='', type='infile', desc='A TSV file containing expression estimates')
    cmd.args['_x'] = Argument(value='custom', type='fix')
    cmd.args['exp-type'] = Argument(prefix='', default='gene', range=['gene', 'transcript'], desc='The type of expression data in the expression_file')
    cmd.outputs['output-vcf'] = Output(value='{output-vcf}', type='outfile')
    return cmd


def bam_read_count():
    cmd = Command()
    cmd.meta.name = 'bamReadCount'
    cmd.runtime.image = 'gudeqing/neoantigen:1.0'
    cmd.runtime.tool = 'python /opt/bam_readcount_helper.py'
    cmd.args['vcf'] = Argument(prefix='-vcf_file ', type='infile', desc='somatic vcf file')
    cmd.args['bam'] = Argument(prefix='-bam_file ', type='infile', desc='rnaseq bam file')
    cmd.args['sample'] = Argument(prefix='-sample ', desc='sample name, will be used as prefix of output file')
    cmd.args['ref_fasta'] = Argument(prefix='-ref_fasta ', type='infile', desc='reference sequence in the fasta format.')
    cmd.args['min-base-quality'] = Argument(prefix='-min_base_qual ', default=20, desc='minimum base quality at a position to use the read for counting.')
    cmd.args['output_dir'] = Argument(prefix='-output_dir ', default='./', desc='output dir')
    cmd.args['not_only_pass'] = Argument(prefix='-not_only_pass', type='bool', default=False, desc='indicate if to only include filter="PASS" variant')
    cmd.outputs['indel_readcount'] = Output(value='{sample}.bam_readcount.indel.tsv')
    cmd.outputs['snv_readcount'] = Output(value='{sample}.bam_readcount.snv.tsv')
    return cmd


def add_read_count_to_vcf():
    cmd = Command()
    cmd.meta.name = 'vcf-readcount-annotator'
    cmd.runtime.image = 'gudeqing/neoantigen:1.0'
    cmd.runtime.tool = 'vcf-readcount-annotator'
    cmd.args['sample'] = Argument(prefix='-s ', desc='sample to annotate')
    cmd.args['out_vcf'] = Argument(prefix='-o ', desc='output vcf')
    cmd.args['variant_type'] = Argument(prefix='-t ', default='all', desc='The type of variant to process.')
    cmd.args['in_vcf'] = Argument(type='infile', desc='input vcf file')
    cmd.args['read_count_file'] = Argument(type='infile', desc='input bam_readcount file')
    cmd.args['seq_type'] = Argument(default='RNA', desc='The type of data in the bam_readcount_file. If `DNA` is chosen, the readcounts will be written to the AD, AF, and DP fields. If `RNA` is chosen, the readcounts will be written to the RAD, RAF, and RDP fields.')
    cmd.outputs['output-vcf'] = Output(value='{out_vcf}')
    return cmd


def pvacseq():
    """
    usage: pvacseq run [-h] [-e1 CLASS_I_EPITOPE_LENGTH]
           [-e2 CLASS_II_EPITOPE_LENGTH]
           [--iedb-install-directory IEDB_INSTALL_DIRECTORY]
           [-b BINDING_THRESHOLD]
           [--percentile-threshold PERCENTILE_THRESHOLD]
           [--allele-specific-binding-thresholds] [-m {lowest,median}]
           [-r IEDB_RETRIES] [-k] [-t N_THREADS]
           [--net-chop-method {cterm,20s}] [--netmhc-stab]
           [--net-chop-threshold NET_CHOP_THRESHOLD]
           [--run-reference-proteome-similarity] [-a {sample_name}]
           [-s FASTA_SIZE] [--exclude-NAs]
           [-d DOWNSTREAM_SEQUENCE_LENGTH]
           [--normal-sample-name NORMAL_SAMPLE_NAME]
           [-p PHASED_PROXIMAL_VARIANTS_VCF] [-c MINIMUM_FOLD_CHANGE]
           [--normal-cov NORMAL_COV] [--tdna-cov TDNA_COV]
           [--trna-cov TRNA_COV] [--normal-vaf NORMAL_VAF]
           [--tdna-vaf TDNA_VAF] [--trna-vaf TRNA_VAF]
           [--expn-val EXPN_VAL]
           [--maximum-transcript-support-level {1,2,3,4,5}]
           [--pass-only]
           input_file sample_name allele
           {MHCflurry,MHCnuggetsI,MHCnuggetsII,NNalign,NetMHC,NetMHCIIpan,NetMHCcons,NetMHCpan,PickPocket,SMM,SMMPMBEC,SMMalign,all,all_class_i,all_class_ii}
           [{MHCflurry,MHCnuggetsI,MHCnuggetsII,NNalign,NetMHC,NetMHCIIpan,NetMHCcons,NetMHCpan,PickPocket,SMM,SMMPMBEC,SMMalign,all,all_class_i,all_class_ii} ...]
           output_dir
    :return:
    """
    cmd = Command()
    cmd.meta.name = 'pvacseq'
    cmd.meta.version = '2.0.4'
    cmd.runtime.image = 'griffithlab/pvactools:latest'
    # cmd.runtime.image = 'gudeqing/neoantigen:1.0'
    cmd.runtime.tool = 'pvacseq run'
    cmd.args['input_file'] = Argument(prefix='', type='infile', desc="A VEP-annotated single- or multi-sample VCF containing genotype, transcript, Wildtype protein sequence, and Downstream protein sequence information.The VCF may be gzipped (requires tabix index)")
    cmd.args['tumor-sample-name'] = Argument(prefix='', desc='The name of the tumor sample being processed. When processing a multi-sample VCF the sample name must be a sample ID in the input VCF #CHROM header line.')
    cmd.args['allele'] = Argument(prefix='', desc="Name of the allele to use for epitope prediction. Multiple alleles can be specified using a comma-separated list.")
    cmd.args['algorithms'] = Argument(prefix='', array=True, default=['MHCflurry', 'MHCnuggetsII'], desc="The epitope prediction algorithms to use. Multiple prediction algorithms can be specified, separated by spaces. ['MHCflurry', 'MHCnuggetsI', 'MHCnuggetsII', 'NNalign', 'NetMHC', 'NetMHCIIpan', 'NetMHCcons', 'NetMHCpan', 'PickPocket', 'SMM', 'SMMPMBEC', 'SMMalign']")
    cmd.args['output_dir'] = Argument(prefix='', default='.', desc='The directory for writing all result files.')
    cmd.args['e1'] = Argument(prefix='-e1 ', type='int', array=True, delimiter=',', default=[8, 9, 10, 11], desc="CLASS_I_EPITOPE_LENGTH Length of MHC Class I subpeptides (neoepitopes) to predict. Multiple epitope lengths can be specified using a comma-separated list. Typical epitope lengths vary between 8-15.")
    cmd.args['e2'] = Argument(prefix='-e2 ', type='int', array=True, delimiter=',', default=[12, 13, 14, 15, 16, 17, 18], desc="CLASS_II_EPITOPE_LENGTH Length of MHC Class II subpeptides (neoepitopes) to predict. Multiple epitope lengths can be specified using a comma-separated list. Typical epitope lengths vary between 11-30. Required for Class II prediction algorithms.")
    cmd.args['iedb-install-directory'] = Argument(prefix='--iedb-install-directory ', default='/opt/iedb/', type='indir', level='optional', desc="Directory that contains the local installation of IEDB MHC I and/or MHC II")
    cmd.args['binding-threshold'] = Argument(prefix='--binding-threshold ', default=500, desc='Report only epitopes where the mutant allele has ic50 binding scores below this value')
    cmd.args['percentile-threshold'] = Argument(prefix='--percentile-threshold ', level='optional', desc="Report only epitopes where the mutant allele has a percentile rank below this value")
    cmd.args['allele-specific-binding-thresholds'] = Argument(prefix='--allele-specific-binding-thresholds', type='bool', default=False, desc="Use allele-specific binding thresholds. To print the allele-specific binding thresholds run `pvacseq allele_specific_cutoffs`. If an allele does not have a special threshold value, the `--binding-threshold` value will be used.")
    cmd.args['top-score-metric'] = Argument(prefix='--top-score-metric ', default='median', desc="The ic50 scoring metric to use when filtering epitopes by binding-threshold or minimum fold change. lowest:Use the best MT Score and Corresponding Fold Change (i.e. the lowest MT ic50 binding score and corresponding fold change of all chosen prediction methods). median: Use the median MT Score and Median Fold Change (i.e. the median MT ic50 binding score and fold change of all chosen prediction methods).")
    cmd.args['keep-tmp-files'] = Argument(prefix='--keep-tmp-files', type='bool', default=False, desc='Keep intermediate output files')
    cmd.args['threads'] = Argument(prefix='--n-threads ', default=5, desc="Number of threads to use for parallelizing peptide-MHC binding prediction calls")
    cmd.args['netmhc-stab'] = Argument(prefix='--netmhc-stab', type='bool', default=False, desc='Run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes.')
    cmd.args['run-reference-proteome-similarity'] = Argument(prefix='--run-reference-proteome-similarity', default=True, desc="Blast peptides against the reference proteome.")
    cmd.args['additional-report-columns'] = Argument(prefix='--additional-report-columns ', level='optional', default='sample_name', desc='Additional columns to output in the final report. If sample_name is chosen, this will add a column with the sample name in every row of the output. This can be useful if you later want to concatenate results from multiple individuals into a single file.')
    cmd.args['fasta-size'] = Argument(prefix='--fasta-size ', default=200, desc="Number of FASTA entries per IEDB request. For some resource-intensive prediction algorithms like Pickpocket and NetMHCpan it might be helpful to reduce this number. Needs to be an even number.")
    cmd.args['exclude-NAs'] = Argument(prefix='--exclude-NAs', default=False, desc="Exclude NA values from the filtered output.")
    cmd.args['downstream-sequence-length'] = Argument(prefix='--downstream-sequence-length ', type='str', default='1000', desc="DOWNSTREAM_SEQUENCE_LENGTH Cap to limit the downstream sequence length for frameshifts when creating the FASTA file. Use 'full' to include the full downstream sequence")
    cmd.args['normal-sample-name'] = Argument(prefix='--normal-sample-name ', level='optional', desc="In a multi-sample VCF, the name of the matched normal sample")
    cmd.args['phased-proximal-variants-vcf'] = Argument(prefix='--phased-proximal-variants-vcf ', type='infile', level='optional', desc='A VCF with phased proximal variant information. Must be gzipped and tabix indexed.')
    cmd.args['minimum-fold-change'] = Argument(prefix='--minimum-fold-change ', default=1.0, desc='Minimum fold change between mutant (MT) binding score and wild-type (WT) score (fold change = WT/MT). The default is 0, which filters no results, but 1 is often a sensible choice (requiring that binding is better to the MT than WT peptide). This fold change is sometimes referred to as a differential agretopicity index.')
    # 如果使用normal cov去筛选突变，意味normal样本必需在该位点也能检测到，考虑到分析结果中，新抗原的结果并不多，将该值设为0
    cmd.args['normal-cov'] = Argument(prefix='--normal-cov ', default=0, desc='Normal Coverage Cutoff')
    cmd.args['tdna-cov'] = Argument(prefix='--tdna-cov ', default=15, desc='Tumor DNA Coverage Cutoff.')
    cmd.args['trna-cov'] = Argument(prefix='--trna-cov ', default=5, desc='Tumor RNA Coverage Cutoff. Only sites above this read depth cutoff will be considered')
    cmd.args['tdna-vaf'] = Argument(prefix='--tdna-vaf ', default=0.05, desc='Tumor DNA VAF Cutoff. Only sites above this cutoff will be considered')
    cmd.args['trna-vaf'] = Argument(prefix='--trna-vaf ', default=0.0, desc='Tumor RNA VAF Cutoff. Only sites above this cutoff will be considered')
    # 下面这个条件放松，如果正常样本存在一定的肿瘤样本污染，还是有可能出现较高vaf
    cmd.args['normal-vaf'] = Argument(prefix='--normal-vaf ', default=0.2, desc='Normal VAF Cutoff. Only sites BELOW this cutoff in normal will be considered')
    cmd.args['expn-val'] = Argument(prefix='--expn-val ', default=1.0, desc='Gene and Transcript Expression cutoff. Only sites above this cutoff will be considered.')
    cmd.args['maximum-transcript-support-level'] = Argument(prefix='--maximum-transcript-support-level ', default=1, range=[1,2,3,4,5], desc="The threshold to use for filtering epitopes on the Ensembl transcript support level (TSL). Keep all epitopes with a transcript support level <= to this cutoff")
    cmd.args['pass-only'] = Argument(prefix='--pass-only', type='bool', default=True, desc='Only process VCF entries with a PASS status.')
    cmd.outputs['outdir'] = Output(value='{output_dir}', type='outdir')
    cmd.outputs['all_epitopes'] = Output(value='{output_dir}/*.all_epitopes.tsv', type='outfile')
    cmd.outputs['aggregated_epitopes'] = Output(value='{output_dir}/*.aggregated.tsv', type='outfile')
    return cmd


def star(sample, platform='illumina', sentieon=False):
    """
    star alignment
    """
    cmd = Command()
    cmd.meta.name = 'star'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 25*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'sentieon STAR' if sentieon else 'STAR'
    cmd.args['threads'] = Argument(prefix='--runThreadN ', default=4, desc='threads to use')
    cmd.args['genomeDir'] = Argument(prefix='--genomeDir ', type='indir', desc='genome index directory')
    cmd.args['read1'] = Argument(prefix='--readFilesIn ', type='infile', array=True, desc='input read1 fastq file list', delimiter=',')
    cmd.args['read2'] = Argument(prefix=' ', type='infile', array=True, level='optional', desc='input read2 fastq file list', delimiter=',')
    cmd.args['outFileNamePrefix'] = Argument(prefix='--outFileNamePrefix ', value=TmpVar(name='~{sample}.', value=f'{sample}.'), desc='output file name prefix')
    cmd.args["outSAMattrRGline"] = Argument(prefix="--outSAMattrRGline ", value=TmpVar(name="ID:~{sample} SM:~{sample} PL:~{platform}", value=f"ID:{sample} SM:{sample} PL:{platform}"), desc="SAM/BAM read group line. The first word contains the read group identifier and must start with 'ID:")
    cmd.args['outSAMtype'] = Argument(prefix='--outSAMtype ', default="BAM SortedByCoordinate", desc='out bam type')
    cmd.args["outSAMunmapped"] = Argument(prefix="--outSAMunmapped ", default="Within", desc='if to include unmapped reads in sam')
    cmd.args["readFilesCommand"] = Argument(prefix="--readFilesCommand ", default="zcat", desc='method to unzip fastq')
    cmd.args["twopassMode"] = Argument(prefix="--twopassMode ", default="Basic", desc='')
    cmd.args["outFilterMultimapNmax"] = Argument(prefix="--outFilterMultimapNmax ", default=20, desc="max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped")
    cmd.args["alignSJoverhangMin"] = Argument(prefix="--alignSJoverhangMin ", default=8, desc="minimum overhang for unannotated junctions")
    cmd.args["alignSJDBoverhangMin"] = Argument(prefix="--alignSJDBoverhangMin ", default=10, desc="minimum overhang for annotated junctions")
    cmd.args["chimSegmentMin"] = Argument(prefix="--chimSegmentMin ", default=12, desc="parameter controls the minimum mapped length of the two segments that is allowed. For example, if you have 2x75 reads and used --chimSegmentMin 20, a chimeric alignment with 130b on one chromosome and 20b on the other will be output, while 135 + 15 won’t be")
    cmd.args["outFilterMismatchNoverLmax"] = Argument(prefix="--outFilterMismatchNoverLmax ", default=0.04, desc="alignment will be output only if its ratio of mismatches to mapped length is less than or equal to this value.")
    cmd.args["outFilterType"] = Argument(prefix="--outFilterType ", default="BySJout", desc="type of filtering")
    cmd.args["outSAMstrandField"] = Argument(prefix="--outSAMstrandField ", default="intronMotif", desc="include for potential use with StringTie for assembly")
    cmd.args["quantMode"] = Argument(prefix="--quantMode ", default="TranscriptomeSAM", desc='output bam file for downstream gene/transcript quantification')
    cmd.args["limitBAMsortRAM"] = Argument(prefix="--limitBAMsortRAM ", default=35000000000, desc="int>=0: maximum available RAM (bytes) for sorting BAM. If =0, it will be set to the genome index size. 0 value can only be used with –genomeLoad NoSharedMemory option.")
    cmd.args["limitIObufferSize"] = Argument(prefix="--limitIObufferSize ", default=200000000, desc="int>0: max available buffers size (bytes) for input/output, per thread")
    cmd.args["outSAMattrIHstart"] = Argument(prefix="--outSAMattrIHstart ", default=0, desc="start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.")
    cmd.args["alignMatesGapMax"] = Argument(prefix="--alignMatesGapMax ", default=100000, desc="maximum gap between two mates")
    cmd.args["alignIntronMax"] = Argument(prefix="--alignIntronMax ", default=100000, desc="maximum intron size")
    cmd.args["alignSJstitchMismatchNmax"] = Argument(prefix="--alignSJstitchMismatchNmax ", default="5 -1 5 5", desc="maximum number of mismatches for stitching of the splice junctions.(1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.")
    cmd.args["chimJunctionOverhangMin"] = Argument(prefix="--chimJunctionOverhangMin ", default=8, desc="minimum overhang for a chimeric junction")
    cmd.args["chimMultimapScoreRange"] = Argument(prefix="--chimMultimapScoreRange ", default=3, desc="the score range for multi-mapping chimeras below the best chimeric score. Only works with –chimMultimapNmax > 1")
    cmd.args["chimMultimapNmax"] = Argument(prefix="--chimMultimapNmax ", default=20, desc="maximum number of chimeric multi-alignments")
    cmd.args["chimNonchimScoreDropMin"] = Argument(prefix="--chimNonchimScoreDropMin ", default=10, desc="to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be smaller than this value")
    cmd.args["chimOutJunctionFormat"] = Argument(prefix="--chimOutJunctionFormat ", default=1, desc="if 1: add comment lines at the end of the file: command line and Nreads: total, unique, multi")
    cmd.args["peOverlapNbasesMin"] = Argument(prefix="--peOverlapNbasesMin ", default=12, desc="minimum number of overlap bases to trigger mates merging and realignment")
    cmd.args["peOverlapMMp"] = Argument(prefix="--peOverlapMMp ", default=0.1, desc="maximum proportion of mismatched bases in the overlap area")
    cmd.args["chimOutType"] = Argument(prefix="--chimOutType ", default="Junctions WithinBAM", desc="type of chimeric output")
    cmd.args["alignInsertionFlush"] = Argument(prefix="--alignInsertionFlush ", default="Right", desc="how to flush ambiguous insertion positions")
    cmd.args["chimScoreJunctionNonGTAG"] = Argument(prefix="--chimScoreJunctionNonGTAG ", default=-4, desc="penalty for a non-GT/AG chimeric junction")
    cmd.args["alignSplicedMateMapLminOverLmate"] = Argument(prefix="--alignSplicedMateMapLminOverLmate ", default=0, desc="alignSplicedMateMapLmin normalized to mate length")
    cmd.args["alignSplicedMateMapLmin"] = Argument(prefix="--alignSplicedMateMapLmin ", default=30, desc="minimum mapped length for a read mate that is spliced")
    cmd.args["quantTranscriptomeBan"] = Argument(prefix="--quantTranscriptomeBan ", default="Singleend", range=["IndelSoftclipSingleend", "Singleend"], desc= "prohibit various alignment type. Singleend is allowed for salmon.")
    # index, flagstats,idxstats
    cmd.args["_1"] = Argument(value=' && ', type='fix')
    cmd.args['_index'] = Argument(value=f"samtools index {sample}.Aligned.sortedByCoord.out.bam", type='fix')
    cmd.args["_2"] = Argument(value=' && ', type='fix')
    cmd.args['_flagstats'] = Argument(value=f'samtools flagstat {sample}.Aligned.sortedByCoord.out.bam > {sample}.align.flagstat.txt', type='fix')
    cmd.args["_3"] = Argument(value=' && ', type='fix')
    cmd.args['_idxstats'] = Argument(value=f"samtools idxstats {sample}.Aligned.sortedByCoord.out.bam > {sample}.align.idxstats.txt", type='fix')
    cmd.args["_4"] = Argument(value=' && ', type='fix')
    # 增加一行命令检查是否有chimeric输出，如果没有，则写一个空文件出来，这样不会导致后面的star-fusion的步骤失败
    cmd.args['_5'] = Argument(value=f"if [ ! -f {sample}.Chimeric.out.junction ]; then echo -e '# No output of chimeric\\n# Nreads 1\tNreadsUnique 1\tNreadsMulti 0' >  {sample}.Chimeric.out.junction ; fi", type='fix')
    cmd.outputs['bam'] = Output(value=f"{sample}.Aligned.sortedByCoord.out.bam")
    cmd.outputs['bam_bai'] = Output(value=f"{sample}.Aligned.sortedByCoord.out.bam.bai")
    cmd.outputs['transcript_bam'] = Output(value=f"{sample}.Aligned.toTranscriptome.out.bam")
    cmd.outputs['align_log'] = Output(value=f"{sample}.Log.final.out")
    cmd.outputs['chimeric'] = Output(value=f"{sample}.Chimeric.out.junction")
    cmd.outputs['sj'] = Output(value=f"{sample}.SJ.out.tab")
    cmd.outputs['flagstat'] = Output(value=f"{sample}.align.flagstat.txt")
    cmd.outputs['idxstats'] = Output(value=f"{sample}.align.idxstats.txt")
    return cmd


def salmon():
    """
    salmon quant in alignment-based mode
    ./bin/salmon quant -t transcripts.fa -l <LIBTYPE> -a aln.bam -o salmon_quant
    """
    cmd = Command()
    cmd.meta.name = 'salmon'
    cmd.meta.desc = 'gene/transcript expression quantification'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 2*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'salmon quant'
    cmd.args['threads'] = Argument(prefix='--threads ', default=4, desc='The number of threads that will be used for quasi-mapping, quantification, and bootstrapping / posterior sampling (if enabled).')
    cmd.args['libType'] = Argument(prefix='--libType ', default='A', desc="Salmon also has the ability to automatically infer (i.e. guess) the library type based on how the first few thousand reads map to the transcriptome. ")
    cmd.args['transcripts'] = Argument(prefix='-t ', type='infile', desc='transcript fasta file')
    cmd.args['geneMap'] = Argument(prefix='-g ', type='infile', desc="The transcript to gene mapping should be provided as either a GTF file, or a in a simple tab-delimited format where each line contains the name of a transcript and the gene to which it belongs separated by a tab.")
    cmd.args['gencode'] = Argument(prefix='--gencode ', type='bool', default=True, desc='expect the input transcript fasta to be in GENCODE format, and will split the transcript name at the first | character.')
    cmd.args['bam'] = Argument(prefix='-a ', type='infile', array=True, delimiter=' ',  desc='transcriptome bam files')
    cmd.args['outDir'] = Argument(prefix='-o ', type='str', default='.', desc='output directory')
    cmd.args['no-version-check'] = Argument(prefix='--no-version-check ', type='bool', default=True, desc='check version')
    cmd.args['gcBias'] = Argument(prefix='--gcBias ', type='bool', default=True, desc='perform gc Bias correction')
    cmd.args['seqBias'] = Argument(prefix='--seqBias ', type='bool', default=True, desc='perform seq Bias correction')
    cmd.args['posBias'] = Argument(prefix='--posBias ', type='bool', default=True, desc='This is meant to model non-uniform coverage biases that are sometimes present in RNA-seq data (e.g. 5’ or 3’ positional bias).')
    cmd.outputs['trans_exp'] = Output(path="{outDir}" + "/quant.sf")
    cmd.outputs['gene_exp'] = Output(path="{outDir}" + "/quant.genes.sf")
    cmd.outputs['outDir'] = Output(path="{outDir}", type='outdir')
    return cmd


def star_fusion():
    """
    star-fusion
    """
    cmd = Command()
    cmd.meta.name = 'star-fusion'
    cmd.meta.source = "https://github.com/STAR-Fusion/STAR-Fusion"
    cmd.meta.version = 'v1.10.0'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 10*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = ' STAR-Fusion'
    cmd.args['threads'] = Argument(prefix='--CPU ', default=4, desc='The number of threads')
    cmd.args['read1'] = Argument(prefix='--left_fq ', level='optional', type='infile', array=True, delimiter=',', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='--right_fq ', level='optional', type='infile',  array=True, delimiter=',', desc='read2 fastq file')
    cmd.args['chimeric_junction'] = Argument(prefix='-J ', type='infile', desc="generated file called 'Chimeric.out.junction' by STAR alignment")
    cmd.args['genomeLibDir'] = Argument(prefix='--genome_lib_dir ', type='indir', desc='ctat_genome_lib_build_dir which contains fusion database')
    cmd.args['outdir'] = Argument(prefix='--output_dir ', default='.', desc='output dir')
    cmd.args['tmpdir'] = Argument(prefix='--tmpdir ', default='.', desc='temp file dir')
    cmd.args['FusionInspector'] = Argument(prefix='--FusionInspector ', level='optional', range=["inspect", "validate"], desc="FusionInspector that provides a more in-depth view of the evidence supporting the predicted fusions.")
    cmd.args['examine_coding_effect'] = Argument(prefix='--examine_coding_effect', type='bool', desc="explore impact of fusions on coding sequences")
    cmd.args['denovo_reconstruct'] = Argument(prefix='--denovo_reconstruct', type='bool', desc="attempt to reconstruct fusion transcripts using Trinity de novo assembly (requires --FusionInspector)")
    cmd.outputs['fusion_predictions'] = Output(value="{outdir}/star-fusion.fusion_predictions.tsv")
    cmd.outputs['fusion_predictions_abridged'] = Output(value="{outdir}/star-fusion.fusion_predictions.abridged.tsv")
    cmd.outputs['fusion_inspector_validate_fusions_abridged'] = Output(value="{outdir}/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv")
    cmd.outputs['fusion_inspector_validate_web'] = Output(value="{outdir}/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv")
    cmd.outputs['fusion_inspector_inspect_web'] = Output(value="{outdir}/FusionInspector-validate/finspector.fusion_inspector_web.html")
    cmd.outputs['fusion_inspector_inspect_fusions_abridged'] = Output(value="{outdir}/FusionInspector-inspect/finspector.FusionInspector.fusions.abridged.tsv")
    return cmd


def collect_metrics(sample):
    """
    picard to collect rnaseq metrics
    """
    cmd = Command()
    cmd.meta.name = 'CollectRnaSeqMetrics'
    # cmd.runtime.image = 'broadinstitute/picard:latest'
    # cmd.runtime.tool = 'java -jar /usr/picard/picard.jar CollectRnaSeqMetrics'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 10*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'java -jar /usr/local/src/picard.jar CollectRnaSeqMetrics'
    cmd.args['bam'] = Argument(prefix='I=', type='infile', desc='input bam file')
    cmd.args['out_metrics'] = Argument(prefix='O=', value=f'{sample}.RnaSeqMetrics.txt', desc='output metric file')
    cmd.args['strandness'] = Argument(prefix='STRAND_SPECIFICITY=', default='NONE', desc='strand-specificity')
    cmd.args['ref_flat'] = Argument(prefix='REF_FLAT=', type='infile', desc='Gene annotations in refFlat form.  Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat')
    cmd.args['ribosomal_interval'] = Argument(prefix='RIBOSOMAL_INTERVALS=', type='infile', desc="Location of rRNA sequences in genome, in interval_list format.  If not specified no bases will be identified as being ribosomal.")
    # cmd.args['cov_pdf'] = Argument(prefix='CHART_OUTPUT=', value=f'{sample}.coverage.pdf', desc='coverage output file')
    cmd.outputs['metrics'] = Output(value="{out_metrics}")
    # cmd.outputs['cov_pdf'] = Output(value="{cov_pdf}")
    return cmd


def arcas_hla(threads=4):
    """
    set -e
    arcasHLA extract --unmapped -t ~{threads} -o ./ ~{bam}
    arcasHLA genotype --min_count 75 -t ~{threads} -o ./ *.1.fq.gz *.2.fq.gz
    arcasHLA merge
    """
    cmd = Command()
    cmd.meta.name = 'arcasHLA'
    cmd.meta.version = '0.2.5'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 5*1024**3
    cmd.runtime.cpu = 2
    # 软链接数据库
    cmd.args['database'] = Argument(prefix='rm -r /home/arcasHLA-master/dat && ln -s ', type='indir', level='optional', desc='database of arcas_software')
    cmd.args['link'] = Argument(prefix='/home/arcasHLA-master/dat &', type='bool', default=False)
    # run software
    cmd.args['_1'] = Argument(value=f'arcasHLA extract --temp ./ --unmapped -t {threads} -o .', type='fix')
    cmd.args['bam'] = Argument(value='', type='infile', desc='input bam file')
    cmd.args['_2'] = Argument(value=' && arcasHLA genotype', type='fix')
    cmd.args['_3'] = Argument(value=f'--min_count 75 -t {threads} -o ./ *.1.fq.gz *.2.fq.gz &&', type='fix')
    cmd.args['_4'] = Argument(value='arcasHLA merge', type='fix')
    cmd.outputs['hla_genotype'] = Output(value='*.genotype.json')
    return cmd


def quant_merge():
    cmd = Command()
    cmd.meta.name = 'quantMerge'
    cmd.meta.desc = 'Merge multiple quantification results into a single file'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 2*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'salmon quantmerge'
    # 下面的quants参数对应的是目录，所以type='indir'
    cmd.args['quants'] = Argument(prefix="--quants ", array=True, type='indir', desc='salmon quant dir list')
    cmd.args['names'] = Argument(prefix='--names ', array=True, level='optional', desc='sample name to assign for each input')
    cmd.args['column'] = Argument(prefix='--column ', default='TPM', range=['TPM', 'NumReads'], desc='indicate which column to merge')
    cmd.args['genes'] = Argument(prefix='--genes ', type='bool', default=False, desc='indicate if to merge gene data, default to merge transcript data')
    cmd.args['out'] = Argument(prefix='--output ', desc='merged result file')
    cmd.outputs['out'] = Output(path="{out}")
    return cmd


def RNASplitReadsAtJunction(sample):
    cmd = Command()
    cmd.meta.name = 'RNASplitReadsAtJunction'
    cmd.meta.desc = 'Perform the splitting of reads into exon segments and reassigning the mapping qualities from STAR'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.0'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'sentieon driver'
    cmd.args['t'] = Argument(prefix='-t ', default=4, desc='number of threads to use in computation')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', desc='reference fasta file')
    cmd.args['bam'] = Argument(prefix='-i ', type='infile', desc='input bam file')
    cmd.args['_x'] = Argument(type='fix', value='--algo RNASplitReadsAtJunction')
    cmd.args['reassign_mapq'] = Argument(prefix='--reassign_mapq ', default='255:60', desc='reassign mapq expression')
    cmd.args['out_bam'] = Argument(value=f'{sample}.junctionSplit.bam', desc='output bam file')
    cmd.outputs['out_bam'] = Output(value='{out_bam}')
    return cmd


def fastv():
    """"""
    cmd = Command()
    cmd.meta.name = 'fastv'
    cmd.runtime.image = 'viurs-fastv-tool:002'
    cmd.runtime.tool = '/data_analysis/software/fastv-0.8.1/fastv'
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', level='optional', type='infile', desc='read2 fastq file')
    cmd.args['threads'] = Argument(prefix='-w ', default=4, desc='thread number')
    cmd.args['other_args'] = Argument(prefix='', default='', desc="other arguments you want to use, such as '-x val'")
    cmd.args['kmer'] = Argument(prefix='-k ', type='infile', default='/data_analysis/software/fastv-0.8.1/data/SARS-CoV-2.kmer.fa', desc='the unique k-mer file of the detection target in fasta format. data/SARS-CoV-2.kmer.fa will be used if none of k-mer/Genomes/k-mer_Collection file is specified ')
    cmd.args['genomes'] = Argument(prefix='-g ', type='infile', default='/data_analysis/software/fastv-0.8.1/data/SARS-CoV-2.genomes.fa', desc='the genomes file of the detection target in fasta format. data/SARS-CoV-2.genomes.fa will be used if none of k-mer/Genomes/k-mer_Collection file is specified')
    cmd.args['out1'] = Argument(prefix='-o ', desc='file name to store read1 with on-target sequences')
    cmd.args['out2'] = Argument(prefix='-O ', level='optional', desc='file name to store read2 with on-target sequences')
    cmd.args['html'] = Argument(prefix='-h ', desc='html report file')
    cmd.args['json'] = Argument(prefix='-j ', desc='json report file')
    # 下面的outputs设置起初是为了能够生成wdl设置, 这里使用”{}“引用其他Argument对象作为输入
    cmd.outputs['out1'] = Output(value="{out1}", type='outfile')
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}")
    cmd.outputs['json'] = Output(value="{json}")
    return cmd


def cnvkit():
    cmd = Command()
    cmd.meta.name = 'cnvkit'
    cmd.meta.source = 'https://github.com/etal/cnvkit'
    cmd.meta.version = '0.9.9'
    cmd.meta.desc = 'detecting copy number variants and alterations genome-wide from high-throughput sequencing'
    cmd.runtime.image = 'etal/cnvkit:latest'
    cmd.runtime.tool = 'cnvkit.py batch'
    cmd.runtime.cpu = 2
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['seq_method'] = Argument(prefix='-m ', default='hybrid', range=['hybrid', 'amplicon', 'wgs'], desc='Sequencing assay type')
    cmd.args['segment_method'] = Argument(prefix='--segment-method ', default='cbs', range=['cbs', 'flasso', 'haar', 'none', 'hmm', 'hmm-tumor','hmm-germline'], desc='method used in segment step')
    cmd.args['drop_low_cov'] = Argument(prefix='--drop-low-coverage', type='bool', default=True, desc='Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor')
    cmd.args['processes'] = Argument(prefix='-p ', default=5, desc='Number of subprocesses used to running each of the BAM files in parallel')
    cmd.args['rscript_path'] = Argument(prefix='--rscript-path ', default='Rscript', desc='Use this option to specify a non-default R')
    cmd.args['normal'] = Argument(prefix='-n ', level='optional', desc='Normal samples (.bam) used to construct the pooled, paired, or flat reference.')
    cmd.args['genome'] = Argument(prefix='-f ', type='infile', desc='Reference genome, FASTA format (e.g. UCSC hg19.fa)')
    cmd.args['targets'] = Argument(prefix='-t ', type='infile', level='optional', desc='Target intervals (.bed or .list)')
    cmd.args['antitargets'] = Argument(prefix='-a ', type='infile', level='optional', desc='Antitarget intervals (.bed or .list)')
    cmd.args['annotate'] = Argument(prefix='--annotate ', type='infile', level='optional', desc='Use gene models from this file to assign names to the target regions. Format: UCSC refFlat.txt or ensFlat.txt file (preferred), or BED, interval list, GFF, or similar.')
    cmd.args['access'] = Argument(prefix='--access', type='infile', level='optional', desc='Regions of accessible sequence on chromosomes (.bed), as output by the <access> command.')
    cmd.args['output_reference'] = Argument(prefix='--output-reference ', level='optional', desc='Output filename/path for the new reference file being created')
    cmd.args['cluster'] = Argument(prefix='--cluster', type='bool', default=False, desc='Calculate and use cluster-specific summary stats in the reference pool to normalize samples')
    cmd.args['reference'] = Argument(prefix='-r ', type='infile', level='optional', desc='Copy number reference file (.cnn)')
    cmd.args['outdir'] = Argument(prefix='-d ', default='.', desc='output directory')
    cmd.args['scatter'] = Argument(prefix='--scatter', type='bool', default=True, desc='Create a whole-genome copy ratio profile as a PDF scatter plot.')
    cmd.args['diagram'] = Argument(prefix='--diagram', type='bool', default=True, desc='Create an ideogram of copy ratios on chromosomes as a PDF')
    cmd.args['tumor_bams'] = Argument(prefix='', type='infile', array=True, desc='Mapped sequence reads (.bam)')
    cmd.outputs['cnr_file'] = Output(value='*.cnr')
    return cmd


def stringtie():
    cmd = Command()
    cmd.meta.name = 'StringTie'
    cmd.meta.source = 'https://ccb.jhu.edu/software/stringtie/index.shtml'
    cmd.meta.version = '2.2.0'
    cmd.meta.desc = 'StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts.'
    cmd.runtime.image = 'nanozoo/stringtie:2.1.7--420b0db'
    cmd.runtime.tool_dir = '/opt/stringtie/stringtie-2.2.0.Linux_x86_64/'
    cmd.runtime.tool = 'stringtie'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 8 * 1024 ** 3
    cmd.args['out_gtf'] = Argument(prefix='-o ', desc='Sets the name of the output GTF file where StringTie will write the assembled transcripts. This can be specified as a full path, in which case directories will be created as needed.')
    cmd.args['L'] = Argument(prefix='-L', type='bool', default=False, desc='long reads processing mode; also enforces -s 1.5 -g 0 (default:false)')
    cmd.args['mix'] = Argument(prefix='--mix', type='bool', default=False, desc='mixed reads processing mode; both short and long read data alignments are expected (long read alignments must be given as the 2nd BAM/CRAM input file)')
    cmd.args['e'] = Argument(prefix='-e', type='bool', default=False, desc='mixed reads processing mode; both short and long read data alignments are expected (long read alignments must be given as the 2nd BAM/CRAM input file)')
    cmd.args['v'] = Argument(prefix='-v', type='bool', default=False, desc='Turns on verbose mode, printing bundle processing details.')
    cmd.args['p'] = Argument(prefix='-p ', type='int', default=4, desc='Specify the number of processing threads (CPUs) to use for transcript assembly.')
    cmd.args['gene_model'] = Argument(prefix='-G ', type='infile', level='optional', desc='Use a reference annotation file (in GTF or GFF3 format) to guide the assembly process. The output will include expressed reference transcripts as well as any novel transcripts that are assembled. This option is required by options -B, -b, -e, -C')
    cmd.args['fr-firststrand'] = Argument(prefix='--rf', type='bool', default=False, desc='Assumes a stranded library fr-firststrand.')
    cmd.args['fr-secondstrand'] = Argument(prefix='--fr', type='bool', default=False, desc='Assumes a stranded library fr-secondstrand.')
    cmd.args['ptf'] = Argument(prefix='--ptf ', type='infile', level='optional', desc='Loads a list of point-features from a text feature file <f_tab> to guide the transcriptome assembly. Accepted point features are transcription start sites (TSS) and polyadenylation sites (CPAS). There are four tab-delimited columns in the feature file. The first three define the location of the point feature on the cromosome (sequence name, coordinate and strand), and the last is the type of the feature (TSS or CPAS)')
    cmd.args['label'] = Argument(prefix='-l ', default='STRG', desc='Sets <label> as the prefix for the name of the output transcripts')
    cmd.args['f'] = Argument(prefix='-f ', type='float', default=0.01, desc='Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given locus. Lower abundance transcripts are often artifacts of incompletely spliced precursors of processed transcripts')
    cmd.args['m'] = Argument(prefix='-m ', type='int', default=200, desc='Sets the minimum length allowed for the predicted transcripts')
    cmd.args['A'] = Argument(prefix='-A ', level='optional', desc='Gene abundances will be reported (tab delimited format) in the output file with the given name.')
    cmd.args['C'] = Argument(prefix='-C ', level='optional', desc='StringTie outputs a file with the given name with all transcripts in the provided reference file that are fully covered by reads (requires -G)')
    cmd.args['a'] = Argument(prefix='-a ', type='int', default=10, desc="Junctions that don't have spliced reads that align across them with at least this amount of bases on both sides are filtered out")
    cmd.args['j'] = Argument(prefix='-j ', type='float', default=1.0, desc="There should be at least this many spliced reads that align across a junction (i.e. junction coverage). This number can be fractional, since some reads align in more than one place. A read that aligns in n places will contribute 1/n to the junction coverage.")
    cmd.args['t'] = Argument(prefix='-t', type='bool', default=False, desc="This parameter disables trimming at the ends of the assembled transcripts. By default StringTie adjusts the predicted transcript's start and/or stop coordinates based on sudden drops in coverage of the assembled transcript")
    cmd.args['c'] = Argument(prefix='-c ', type='int', default=1, desc="Sets the minimum read coverage allowed for the predicted transcripts. A transcript with a lower coverage than this value is not shown in the output")
    cmd.args['s'] = Argument(prefix='-s ', type='float', default=4.75, desc='Sets the minimum read coverage allowed for single-exon transcripts')
    cmd.args['conservative'] = Argument(prefix='--conservative', type='bool', default=False, desc='Assembles transcripts in a conservative mode. Same as -t -c 1.5 -f 0.05')
    cmd.args['g'] = Argument(prefix='-g ', type='int', default=50, desc='Minimum locus gap separation value. Reads that are mapped closer than this distance are merged together in the same processing bundle')
    cmd.args['B'] = Argument(prefix='-B', type='bool', default=False, desc='This switch enables the output of Ballgown input table files (*.ctab) containing coverage data for the reference transcripts given with the -G option.')
    cmd.args['b'] = Argument(prefix='-b ', level='optional', desc='Just like -B this option enables the output of *.ctab files for Ballgown, but these files will be created in the provided directory <path> instead of the directory specified by the -o option. Note: adding the -e option is recommended with the -B/-b options, unless novel transcripts are still wanted in the StringTie GTF output.')
    cmd.args['M'] = Argument(prefix='-M ', type='float', default=0.95, desc='Sets the maximum fraction of muliple-location-mapped reads that are allowed to be present at a given locus')
    cmd.args['x'] = Argument(prefix='-x ', level='optional', desc="Ignore all read alignments (and thus do not attempt to perform transcript assembly) on the specified reference sequences. Parameter <seqid_list> can be a single reference sequence name (e.g. -x chrM) or a comma-delimited list of sequence names (e.g. -x 'chrM,chrX,chrY'). This can speed up StringTie especially in the case of excluding the mitochondrial genome, whose genes may have very high coverage in some cases, even though they may be of no interest for a particular RNA-Seq analysis. The reference sequence names are case sensitive, they must match identically the names of chromosomes/contigs of the target genome against which the RNA-Seq reads were aligned in the first place.")
    cmd.args['u'] = Argument(prefix='-u', type='bool', default=False, desc='Turn off multi-mapping correction. In the default case this correction is enabled, and each read that is mapped in n places only contributes 1/n to the transcript coverage instead of 1')
    cmd.args['cram_ref'] = Argument(prefix='--cram-ref ', type='infile', level='optional', desc="for CRAM input files, the reference genome sequence can be provided as a multi-FASTA file the same chromosome sequences that were used when aligning the reads. This option is optional but recommended as StringTie can make use of some alignment/junction quality data (mismatches around the junctions) that can be more accurately assessed in the case of CRAM files when the reference genome sequence is also provided")
    cmd.args['merge'] = Argument(prefix='--merge', type='bool', default=False, desc='In the merge mode, StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts. This mode is used in the new differential analysis pipeline to generate a global, unified set of transcripts (isoforms) across multiple RNA-Seq samples.')
    cmd.args['bam'] = Argument(type='infile', array=True, desc='input bam(s) file path')
    cmd.outputs['out_gtf'] = Output(value='{out_gtf}')
    # result of gffcompare
    cmd.outputs['tmap'] = Output(value='*.tmap')
    cmd.outputs['refmap'] = Output(value='*.refmap')
    return cmd


def gffcompare():
    cmd = Command()
    cmd.meta.name = 'gffcompare'
    cmd.meta.source = 'https://github.com/gpertea/gffcompare'
    cmd.meta.version = '0.12.6'
    cmd.meta.desc = 'StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential transcripts.'
    cmd.runtime.image = 'nanozoo/stringtie:2.1.7--420b0db'
    cmd.runtime.tool_dir = '/opt/gffcompare/gffcompare-0.12.6.Linux_x86_64/'
    cmd.runtime.tool = 'gffcompare'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['gtfs'] = Argument(prefix='', type='infile', level='optional', array=True, desc='provide a list of (query) GTF files to process')
    cmd.args['outprefix'] = Argument(prefix='-o ', default='gffcmp', desc='All output files created by gffcompare will have this prefix (e.g. .loci, .tracking, etc.)')
    cmd.args['ref'] = Argument(prefix='-r ', type='infile', level='optional', desc='An optional “reference” annotation GFF file. Each sample is matched against this file, and sample isoforms are tagged as overlapping, matching, or novel where appropriate')
    cmd.args['R'] = Argument(prefix='-R', type='bool', default=False, desc='If -r was specified, this option causes gffcompare to ignore reference transcripts that are not overlapped by any transcript in one of input1.gt,...,inputN.gtf. Useful for ignoring annotated transcripts that are not present in your RNA-Seq samples and thus adjusting the "sensitivity" calculation in the accuracy report written in the file.')
    cmd.args['Q'] = Argument(prefix='-Q', type='bool', default=False, desc='If -r was specified, this option causes gffcompare to ignore input transcripts that are not overlapped by any transcript in the reference. Useful for adjusting the “precision” calculation in the accuracy report written in the file.')
    cmd.args['M'] = Argument(prefix='-M', type='bool', default=False, desc='discard (ignore) single-exon transfrags and reference transcripts (i.e. consider only multi-exon transcripts)')
    cmd.args['N'] = Argument(prefix='-N', type='bool', default=False, desc='discard (ignore) single-exon reference transcripts; single-exon transfrags are still considered, but they will never find an exact match')
    cmd.args['D'] = Argument(prefix='-D', type='bool', default=False, desc='discard "duplicate" (redundant) query transfrags (i.e. those with the same intron chain) within a single sample (and thus disable "annotation" mode)')
    cmd.args['genome'] = Argument(prefix='-s ', type='infile', level='optional', desc='path to genome sequences (optional); this will enable the "repeat" ('r') classcode assessment;repeats must be soft-masked (lower case) in the genomic sequence')
    cmd.args['e'] = Argument(prefix='-e ', default=100, desc='Maximum distance (range) allowed from free ends of terminal exons of reference transcripts when assessing exon accuracy. By default, this is 100')
    cmd.args['d'] = Argument(prefix='-d ', default=100, desc='Maximum distance (range) for grouping transcript start sites, by default 100')
    cmd.args['strict_match'] = Argument(prefix='--strict-match', type='bool', default=False, desc="transcript matching takes into account the -e range for terminal exons; code '=' is only assigned if transcript ends are; within that range, otherwise code '~' is assigned just for intron chain; match (or significant overlap in the case of single exon transcripts)")
    cmd.args['p'] = Argument(prefix='-p ', default='TCONS', desc='The name prefix to use for consensus/combined transcripts in the <outprefix>.combined.gtf file')
    cmd.args['C'] = Argument(prefix='-C', type='bool', default=False, desc='Discard the “contained” transfrags from the .combined.gtf output. By default, without this option, gffcompare writes in that file isoforms that were found to be fully contained/covered (with the same compatible intron structure) by other transfrags in the same locus, with the attribute “contained_in” showing the first container transfrag found.')
    cmd.args['A'] = Argument(prefix='-A', type='bool', default=False, desc='Like -C but will not discard intron-redundant transfrags if they start on a different 5-end exon (keep alternate transcript start sites)')
    cmd.args['X'] = Argument(prefix='-X', type='bool', default=False, desc="Like -C but also discard contained transfrags if transfrag ends stick out within the container's introns")
    cmd.args['K'] = Argument(prefix='-K', type='bool', default=False, desc='For -C/-A/-X, do NOT discard any redundant transfrag matching a reference')
    cmd.args['T'] = Argument(prefix='-T', type='bool', default=False, desc='Do not generate .tmap and .refmap files for each input file')
    cmd.outputs['combined_gtf'] = Output(value='{outprefix}.combined.gtf')
    cmd.outputs['annotated_gtf'] = Output(value='{outprefix}.annotated.gtf', desc='If a single GTF/GFF query file is provided as input and no specific options to remove "duplicated"/redundant transfrags were given (-D, -S, -C, -A, -X), GffCompare outputs a file called <outprefix>.annotated.gtf instead of <outprefix>.combined.gtf Its format is similar to the above, but preserves transcript IDs (so the -p option is ignored)')
    return cmd


def gffread():
    cmd = Command()
    cmd.meta.name = 'gffread'
    cmd.meta.source = 'https://github.com/gpertea/gffread'
    cmd.meta.version = '0.12.7'
    cmd.meta.desc = 'GFF/GTF utility providing format conversions, filtering, FASTA sequence extraction and more.'
    cmd.runtime.image = ''
    cmd.runtime.tool_dir = '/opt/gffread/gffread-0.12.7.Linux_x86_64/'
    cmd.runtime.tool = 'gffread'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['genome'] = Argument(prefix='-g ', type='infile', level='optional', desc='full path to a multi-fasta file with the genomic sequences for all input mappings OR a directory with single-fasta files')
    cmd.args['ids'] = Argument(prefix='--ids ', type='infile', level='optional', desc='discard records/transcripts if their IDs are not listed in <IDs.lst>')
    cmd.args['minlen'] = Argument(prefix='-l ', default=10, desc='discard transcripts shorter than <minlen> bases')
    cmd.args['r'] = Argument(prefix='-r ', type='infile', level='optional', desc='only show transcripts overlapping coordinate range <start>..<end>')
    cmd.args['R'] = Argument(prefix='-R ', type='infile', level='optional', desc='for -r option, discard all transcripts that are not fully contained within the given range')
    cmd.args['E'] = Argument(prefix='-E', type='bool', default=False, level='optional', desc='expose (warn about) duplicate transcript IDs and other potential problems with the given GFF/GTF records')
    cmd.args['T'] = Argument(prefix='-T', type='bool', default=False, level='optional', desc='main output will be GTF instead of GFF3')
    cmd.args['m'] = Argument(prefix='-m ', type='infile', level='optional', desc='<chr_replace> is a name mapping table for converting reference sequence names, having this 2-column format: <original_ref_ID> <new_ref_ID>')
    cmd.args['w'] = Argument(prefix='-w ', level='optional', desc='write a fasta file with spliced exons for each transcript')
    cmd.args['input_gff'] = Argument(prefix='', type='infile', desc="input gtf/gff file")
    cmd.outputs['transcript_fa'] = Output(value='{w}')
    return cmd


def TransDecoder_LongOrfs():
    cmd = Command()
    cmd.meta.name = 'TransDecoderLongOrfs'
    cmd.meta.source = 'https://github.com/TransDecoder/TransDecoder'
    cmd.meta.version = '5.5.0'
    cmd.meta.desc = 'Use TransDecoder predict the likely coding regions. TransDecoder identifies candidate coding regions within transcript sequences, such as those generated by de novo RNA-Seq transcript assembly using Trinity, or constructed based on RNA-Seq alignments BAM.'
    cmd.runtime.image = ''
    cmd.runtime.tool_dir = '/opt/TransDecoder/TransDecoder-TransDecoder-v5.5.0/'
    cmd.runtime.tool = 'TransDecoder.LongOrfs'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['t'] = Argument(prefix='-t ', type='infile', desc='transcript fasta file')
    cmd.args['m'] = Argument(prefix='-m ', type='int', default=64, desc='minimum protein length (default: 100)')
    cmd.args['gene_trans_map'] = Argument(prefix='--gene_trans_map ', type='infile', level='optional', desc='gene-to-transcript identifier mapping file (tab-delimited, gene_id<tab>trans_id<return> )')
    cmd.args['genetic_code'] = Argument(prefix='-G ', default='universal', desc='genetic code (default: universal; see PerlDoc; options: Euplotes, Tetrahymena, Candida, Acetabularia)')
    cmd.args['only_top_strand'] = Argument(prefix='-S', type='bool', default=False, desc='strand-specific (only analyzes top strand)')
    cmd.args['outdir'] = Argument(prefix='--output_dir ', default='LongOrfs', desc='path to intended output directory')
    cmd.outputs['outdir'] = Output(value='{outdir}', type='outdir')
    return cmd


def transdecoder_predict():
    cmd = Command()
    cmd.meta.name = 'TransDecoderPredict'
    cmd.meta.source = 'https://github.com/TransDecoder/TransDecoder'
    cmd.meta.version = '5.5.0'
    cmd.meta.desc = 'Use TransDecoder predict the likely coding regions. TransDecoder identifies candidate coding regions within transcript sequences, such as those generated by de novo RNA-Seq transcript assembly using Trinity, or constructed based on RNA-Seq alignments BAM.'
    cmd.runtime.image = ''
    cmd.runtime.tool_dir = '/opt/TransDecoder/TransDecoder-TransDecoder-v5.5.0/'
    cmd.runtime.tool = 'TransDecoder.Predict'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['t'] = Argument(prefix='-t ', type='infile', desc='transcript fasta file')
    cmd.args['retain_long_orfs_mode'] = Argument(prefix='--retain_long_orfs_mode ', default='dynamic', desc='In dynamic mode, sets range according to 1%FDR in random sequence of same GC content.')
    cmd.args['retain_long_orfs_length'] = Argument(prefix='--retain_long_orfs_length ', level='optional', default=1000000, desc="under 'strict' mode, retain all ORFs found that are equal or longer than these many nucleotides even if no other evidence marks it as coding (default: 1000000) so essentially turned off by default")
    cmd.args['retain_pfam_hits'] = Argument(prefix='--retain_pfam_hits ', type='infile', level='optional', desc='file name of domain table output file from running hmmscan to search Pfam')
    cmd.args['retain_blastp_hits'] = Argument(prefix='--retain_blastp_hits ', type='infile', level='optional', desc="blastp output in '-outfmt 6' format. ")
    cmd.args['single_best_only'] = Argument(prefix='--single_best_only', type='bool', default=False, desc='Retain only the single best orf per transcript (prioritized by homology then orf length)')
    cmd.args['output_dir'] = Argument(prefix='--output_dir ', type='indir', default='LongOrfs_dir', desc='output directory from the TransDecoder.LongOrfs step (default: basename( -t val ) + ".transdecoder_dir")')
    cmd.args['no_refine_starts'] = Argument(prefix='--no_refine_starts', type='bool', default=False, desc="start refinement identifies potential start codons for 5' partial ORFs using a PWM, process on by default")
    cmd.outputs['LongOrfs_dir'] = Output(value='{output_dir}')
    cmd.outputs['pep_file'] = Output(value='*.transdecoder.pep')
    return cmd


def diamond_makedb():
    cmd = Command()
    cmd.meta.name = 'diamond_makedb'
    cmd.meta.source = 'https://github.com/bbuchfink/diamond'
    cmd.meta.version = '2.0.14'
    cmd.meta.desc = 'DIAMOND is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data'
    cmd.runtime.image = ''
    cmd.runtime.tool_dir = '/opt/diamond/'
    cmd.runtime.tool = 'diamond makedb'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['threads'] = Argument(prefix='-p ', type='int', default=4, desc='number of CPU threads')
    cmd.args['db'] = Argument(prefix='-d ', desc='database file name prefix')
    cmd.args['in'] = Argument(prefix='--in ', desc='input reference file in FASTA format')
    cmd.args['taxonmap'] = Argument(prefix='--taxonmap ', level='optional', desc='protein accession to taxid mapping file')
    cmd.args['taxonnodes'] = Argument(prefix='--taxonnodes ', level='optional', desc='taxonomy nodes.dmp from NCBI')
    cmd.args['taxonnames'] = Argument(prefix='--taxonnames ', level='optional', desc='taxonomy names.dmp from NCBI')
    cmd.outputs['db'] = Output(value='{db}.dmnd')
    return cmd


def diamond_blastp():
    cmd = Command()
    cmd.meta.name = 'diamond_blastp'
    cmd.meta.source = 'https://github.com/bbuchfink/diamond'
    cmd.meta.version = '2.0.14'
    cmd.meta.desc = 'DIAMOND is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data'
    cmd.runtime.image = ''
    cmd.runtime.tool_dir = '/opt/diamond/'
    cmd.runtime.tool = 'diamond blastp'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['threads'] = Argument(prefix='-p ', type='int', default=4, desc='number of CPU threads')
    cmd.args['db'] = Argument(prefix='-d ', desc='database file name prefix')
    cmd.args['query'] = Argument(prefix='-q ', type='infile', desc='query peptide fasta file')
    cmd.args['unaligned'] = Argument(prefix='--un ', desc='file for unaligned queries')
    cmd.args['aligned'] = Argument(prefix='--al ', desc='file or aligned queries')
    cmd.args['out'] = Argument(prefix='-o ', desc='output of alignment result')
    cmd.args['max-target-seqs'] = Argument(prefix='-k ', default=25, desc='maximum number of target sequences to report alignments for (default=25)')
    cmd.args['compress'] = Argument(prefix='--compress ', default=0, range=[0,1], desc='compression for output files')
    cmd.args['evalue'] = Argument(prefix='-e ', default=0.001, desc='maximum e-value to report alignments')
    cmd.args['identity'] = Argument(prefix='--id ', default=50, desc='minimum identity% to report an alignment')
    cmd.args['query-cover'] = Argument(prefix='--query-cover ', default=50, desc='minimum query cover% to report an alignment')
    cmd.args['subject-cover'] = Argument(prefix='--subject-cover ', default=50, desc='minimum subject cover% to report an alignment')
    cmd.args['mode'] = Argument(prefix='--', default='sensitive', range=['fast', 'mid-sensitive', 'sensitive', 'more-sensitive', 'very-sensitive', 'ultra-sensitive'], desc='indicate to enable which searching mode to use')
    cmd.args['algo'] = Argument(prefix='--algo ', default='0', range=['0', '1', 'ctg'], desc='Seed search algorithm (0=double-indexed/1=query-indexed/ctg=contiguous-seed)')
    cmd.args['other_args'] = Argument(prefix='', level='optional', desc='other customised arguments you wish to use, such as "-other value"')
    cmd.outputs['aligned'] = Output(value='{aligned}')
    cmd.outputs['unaligned'] = Output(value='{unaligned}')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def mhcflurry_predict():
    cmd = Command()
    cmd.meta.name = 'mhcflurry'
    cmd.meta.source = 'https://github.com/openvax/mhcflurry'
    cmd.meta.version = '2.0.14'
    cmd.meta.desc = 'MHCflurry is an open source package for peptide/MHC I binding affinity prediction.'
    cmd.runtime.image = ''
    cmd.runtime.tool = 'mhcflurry-predict'
    cmd.runtime.cpu = 4
    cmd.runtime.memory = 6 * 1024 ** 3
    cmd.args['input_csv'] = Argument(prefix='', type='infile', desc='The input CSV file is expected to contain columns "allele", "peptide", and, optionally, "n_flank", and "c_flank".')
    cmd.args['out'] = Argument(prefix='--out ', desc='output csv')
    cmd.args['no-throw'] = Argument(prefix='--no-throw', type='bool', default=False, desc='Return NaNs for unsupported alleles or peptides instead of raising')
    cmd.args['always-include-best-allele'] = Argument(prefix='--always-include-best-allele', type='bool', default=False, desc='Always include the best_allele column even when it is identical to the allele column (i.e. all queries are monoallelic)')
    cmd.args['models'] = Argument(prefix='--models ', type='indir', desc='Directory containing models. Either a binding affinity predictor or a presentation predictor can be used.')
    cmd.args['affinity-only'] = Argument(prefix='--affinity-only', type='bool', default=False, desc='Affinity prediction only (no antigen processing or presentation)')
    cmd.args['no-flanking'] = Argument(prefix='--no-flanking', type='bool', default=False, desc='Do not use flanking sequence information even when available')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def mhcflurry_predict_scan():
    cmd = Command()
    cmd.meta.name = 'mhcflurry'
    cmd.meta.source = 'https://github.com/openvax/mhcflurry'
    cmd.meta.version = '2.0.14'
    cmd.meta.desc = 'MHCflurry is an open source package for peptide/MHC I binding affinity prediction.'
    cmd.runtime.image = ''
    cmd.runtime.tool = 'mhcflurry-predict-scan'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['input_file'] = Argument(prefix='', type='infile', level='optional', desc='Sequences to predict, csv or fasta')
    cmd.args['alleles'] = Argument(prefix='--alleles ', array=True, delimiter=',', desc='Alleles to predict, such as LA-A*02:01,HLA-A*03:01')
    cmd.args['sequences'] = Argument(prefix='--sequences ', level='optional', desc='Sequences to predict (exclusive with passing an input file)')
    cmd.args['sequence-id-column'] = Argument(prefix='--sequence-id-column ', default='sequence_id', desc="Input CSV column name for sequence IDs. Default: 'sequence_id'")
    cmd.args['sequence-column'] = Argument(prefix='--sequence-column ', default='sequence', desc="Input CSV column name for sequences. Default: 'sequence'")
    cmd.args['peptide-lengths'] = Argument(prefix='--peptide-lengths ', array=True, delimiter=',', default=[8, 9, 10, 11], desc='Peptide lengths to consider. Pass as START-END (e.g. 8-11) or a comma-separated list (8,9,10,11). When using START-END, the range is INCLUSIVE on both ends. Default: 8-11')
    cmd.args['results-all'] = Argument(prefix='--results-all', type='bool', default=False, desc='Return results for all peptides regardless of affinity, etc.')
    cmd.args['results-best'] = Argument(prefix='--results-best ', level='optional', range=['presentation_score', 'processing_score', 'affinity', 'affinity_percentile'], desc='Take the top result for each sequence according to the specified predicted quantity')
    cmd.args['results-filtered'] = Argument(prefix='--results-filtered ',  default='affinity', range=['presentation_score', 'processing_score', 'affinity', 'affinity_percentile'], desc='Filter results by the specified quantity')
    cmd.args['threshold-presentation-score'] = Argument(prefix='--threshold-presentation-score ', default=0.7, desc='Threshold if filtering by presentation score')
    cmd.args['threshold-processing-score'] = Argument(prefix='--threshold-processing-score ', default=0.5, desc='Threshold if filtering by processing score.')
    cmd.args['threshold-affinity'] = Argument(prefix='--threshold-affinity ', default=500, desc='Threshold if filtering by affinity.')
    cmd.args['threshold-affinity-percentile'] = Argument(prefix='--threshold-affinity-percentile ', default=2.0, desc='Threshold if filtering by affinity percentile.')
    cmd.args['out'] = Argument(prefix='--out ', desc='output csv')
    cmd.args['output-delimiter'] = Argument(prefix='--output-delimiter ', default=',')
    cmd.args['no-affinity-percentile'] = Argument(prefix='--no-affinity-percentile', type='bool', default=False, desc='Do not include affinity percentile rank')
    cmd.args['models'] = Argument(prefix='--models ', type='indir', desc='Directory containing presentation models, such as mhcflurry/4/2.0.0/models_class1_presentation/models')
    cmd.args['no-flanking'] = Argument(prefix='--no-flanking', type='bool', default=False, desc='Do not use flanking sequence information in predictions')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def RNAmining():
    cmd = Command()
    cmd.meta.name = 'RNAmining'
    cmd.meta.source = 'https://rnamining.integrativebioinformatics.me/about'
    cmd.meta.version = '1.0.4'
    cmd.meta.desc = 'This tool was implemented using XGBoost machine learning algorithm. Machine learning is a subfield of computer science that developed from the study of pattern recognition and computational learning theories in artificial intelligence. This tool operate through a model obtained from training data analyzes and produces an inferred function, which can be used for mapping new examples'
    cmd.runtime.image = ''
    cmd.runtime.tool = 'python3 /opt/RNAmining/rnamining.py'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['query'] = Argument(prefix='-f ', type='infile', desc='The filename with a sequence to predict')
    cmd.args['organism'] = Argument(prefix='-organism_name ', default='Homo_sapiens', range=['Escherichia_coli', 'Arabidopsis_thaliana', 'Drosophila_melanogaster', 'Homo_sapiens', 'Mus_musculus', 'Saccharomyces_cerevisiae'], desc='The name of the organism you want to predict/train')
    cmd.args['prediction_type'] = Argument(prefix='-prediction_type ', default='coding_prediction', range=['coding_prediction', 'ncRNA_functional_assignation'])
    cmd.args['outdir'] = Argument(prefix='-output_folder ', default='.', desc='output directory')
    cmd.outputs['outdir'] = Output(value='{outdir}', type='outdir')
    cmd.outputs['predictions'] = Output(value='{outdir}/predictions.txt')
    cmd.outputs['codings'] = Output(value='{outdir}/codings.txt')
    return cmd


def MixMHC2pred():
    cmd = Command()
    cmd.meta.name = 'MixMHC2pred'
    cmd.meta.source = 'https://github.com/GfellerLab/MixMHC2pred'
    cmd.meta.version = '1.2'
    cmd.meta.desc= 'MixMHC2pred is a predictor of HLA class II ligands and epitopes. MixMHC2pred is meant for scoring different peptides and prioritising the most likely HLA-II ligands and epitopes. As it is trained on naturally presented peptides, it does not output a predicted affinity value, simply a score.'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.3'
    cmd.runtime.tool = '/opt/MixMHC2pred/MixMHC2pred_unix'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['input'] = Argument(prefix='-i ', type='infile', desc='File listing all the peptides with one peptide per line. It can also be a fasta file (note that the peptide description given by the lines with ">" are not kept). Please note that even in fasta format, the input should consist in a list of peptides: MixMHC2pred is not cutting inputted proteins into shorter fragments that could be presented but it uses the input sequences as given in the file directly')
    cmd.args['output'] = Argument(prefix='-o ', default='mixMHC2pred.txt', desc='The name of the output file (including the directory). Peptides are kept in the same order than in the input file.')
    cmd.args['alleles'] = Argument(prefix='-a ', array=True, desc='List of HLA-II alleles to test. Use for example the nomenclature DRB1_03_01 for HLA-DRB1*03:01 and DPA1_01_03__DPB1_04_01 for HLA-DPA1*01:03-DPB1*04:01. The full list of alleles available and corresponding nomenclature is given in the file Alleles_list.txt.')
    cmd.outputs['output'] = Output(value='{output}')
    return cmd


def makeblastdb():
    cmd = Command()
    cmd.meta.name = 'makeblastdb'
    cmd.meta.source = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
    cmd.meta.version = '2.12.0+'
    cmd.meta.desc = 'BLAST finds regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance.'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.3'
    cmd.runtime.tool_dir = '/opt/ncbi-blast-2.12.0+/bin'
    cmd.runtime.tool = 'makeblastdb'
    cmd.runtime.cpu = 3
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['input_file'] = Argument(prefix='-in ', type='infile', desc='Input file/database name')
    cmd.args['input_type'] = Argument(prefix='-input_type ', range=['asn1_bin', 'asn1_txt', 'blastdb', 'fasta'], default='fasta', desc='Type of the data specified in input_file')
    cmd.args['dbtype'] = Argument(prefix='-dbtype ', range=['nucl', 'prot'], default='prot', desc='Molecule type of target db')
    cmd.args['out'] = Argument(prefix='-out ', default='reference', desc='Name of BLAST database to be created')
    cmd.outputs['out'] = Output('{out}')
    return cmd


def blastp():
    cmd = Command()
    cmd.meta.name = 'blastp'
    cmd.meta.source = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
    cmd.meta.version = '2.12.0+'
    cmd.meta.desc = 'BLAST finds regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance.'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.3'
    cmd.runtime.tool_dir = '/opt/ncbi-blast-2.12.0+/bin'
    cmd.runtime.tool = 'blastp'
    cmd.runtime.cpu = 4
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['query'] = Argument(prefix='-query ', type='infile', desc='input file')
    cmd.args['task'] = Argument(prefix='-task ', default='blastp-short', range=['blastp', 'blastp-fast', 'blastp-short'], desc='Task to execute')
    cmd.args['db'] = Argument(prefix='-db ', type='infile', desc='database file')
    cmd.args['out'] = Argument(prefix='-out ', desc='output file name')
    cmd.args['evalue'] = Argument(prefix='-evalue ', type='int', default=10, desc='Expectation value (E) threshold for saving hits')
    cmd.args['word_size'] = Argument(prefix='-word_size ', type='int', default=2, desc='Word size for wordfinder algorithm')
    cmd.args['gapopen'] = Argument(prefix='-gapopen ', type='int', default=5, desc='Cost to open a gap')
    cmd.args['comp_based_stats'] = Argument(prefix='-comp_based_stats ', default='2', desc='Use composition-based statistics.D or d: default (equivalent to 2 ); 0 or F or f: No composition-based statistics; 1: Composition-based statistics as in NAR 29:2994-3005, 2001; 2 or T or t : Composition-based score adjustment as in Bioinformatics; 21:902-911,; 2005, conditioned on sequence properties; 3: Composition-based score adjustment as in Bioinformatics 21:902-911,; 2005, unconditionally')
    cmd.args['outfmt'] = Argument(prefix='-outfmt ', type='int', default='7', desc='alignment view options')
    cmd.args['max_target_seqs'] = Argument(prefix='-max_target_seqs ', type='int', default=500, desc='Maximum number of aligned sequences to keep')
    cmd.args['ungapped'] = Argument(prefix='-ungapped', type='bool', default=False, desc='Perform ungapped alignment only?')
    cmd.args['num_threads'] = Argument(prefix='-num_threads ', type='int', default=4, desc='Number of threads (CPUs) to use in the BLAST search')
    cmd.args['max_hsps'] = Argument(prefix='-max_hsps ', type='int', level='optional', desc='Set maximum number of HSPs per subject sequence to save for each query')
    cmd.args['qcov_hsp_perc'] = Argument(prefix='-qcov_hsp_perc ', type='int', level='optional', desc='Percent query coverage per hsp')
    cmd.outputs['out'] = Output(value='{out}')
    return cmd


def netMHCIIPan():
    cmd = Command()
    cmd.meta.name = 'netMHCIIPan4'
    cmd.meta.source = 'https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.1'
    cmd.meta.version = '4.1'
    cmd.meta.desc = 'The NetMHCIIpan-4.1 server predicts peptide binding to any MHC II molecule of known sequence using Artificial Neural Networks (ANNs). It is trained on an extensive dataset of over 500.000 measurements of Binding Affinity (BA) and Eluted Ligand mass spectrometry (EL), covering the three human MHC class II isotypes HLA-DR, HLA-DQ, HLA-DP, as well as the mouse molecules (H-2). The introduction of EL data extends the number of MHC II molecules covered, since BA data covers 59 molecules and EL data covers 74. As mentioned, the network can predict for any MHC II of known sequence, which the user can specify as FASTA format. The network can predict for peptides of any length.'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.3'
    cmd.runtime.tool_dir = '/opt/netMHCIIpan-4.1/'
    cmd.runtime.tool = 'netMHCIIpan'
    cmd.runtime.cpu = 4
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.args['BA'] = Argument(prefix='-BA', type='bool', default=False, desc='Include BA predictions, default is EL only')
    cmd.args['context'] = Argument(prefix='-context', type='bool', default=False, desc='Predict with context encoding')
    cmd.args['tdir'] = Argument(prefix='-tdir ', default='.', desc='Temporary directory')
    cmd.args['alleles'] = Argument(prefix='-a ', array=True, delimiter=',', desc='HLA allele')
    cmd.args['inptype'] = Argument(prefix='-inptype ', range=['0', '1'], default='0', desc='Input type [0] FASTA [1] Peptide')
    cmd.args['rankS'] = Argument(prefix='-rankS ', default=1.0, desc='Threshold for strong binders (%Rank)')
    cmd.args['rankW'] = Argument(prefix='-rankW ', default=5.0, desc='Threshold for weak binders (%Rank)')
    cmd.args['length'] = Argument(prefix='-length ', type='int', array=True, delimiter=',', default=[15], desc='Peptide length (multiple length with ,). Used for FASTA input only')
    cmd.args['sort'] = Argument(prefix='-s', type='bool', default=True, desc='Sort output on descending affinity')
    cmd.args['unique'] = Argument(prefix='-u', type='bool', default=False, desc='Print unique binding core only')
    cmd.args['infile'] = Argument(prefix='-f ', type='infile', desc=' File with the input data')
    cmd.args['xls'] = Argument(prefix='-xls', type='bool', default=False, desc='save output to xlsfile')
    cmd.args['xlsfile'] = Argument(prefix='-xlsfile ', default='NetMHCIIpan_out.xls', desc='save output to xlsfile')
    cmd.args['stdout'] = Argument(prefix='> ', default='raw.output.txt', desc='output file name')
    cmd.outputs['out'] = Output(value='{xlsfile}')
    return cmd


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals())
