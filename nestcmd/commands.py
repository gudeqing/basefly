from .nestcmd import Argument, Output, Command, TmpVar, TopVar

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
    cmd.args['run-reference-proteome-similarity'] = Argument(prefix='--run-reference-proteome-similarity', default=False, desc="Blast peptides against the reference proteome.")
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
