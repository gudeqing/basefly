import json
import os
import math
import pandas as pd
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar
from utils.get_fastq_info import get_fastq_info
__author__ = 'gdq'

"""
选定N个marker
从multipcr技术出发，每个marker设计一对引物，靶向富集目标dna片段
NGS测序

下机数据分析：
fastq数据清洗，应该考虑剔除primer，如何剔除：
    1. 如果read以primer开头，可以用cutadapt软件做到
    2. 如果primer在read中间，就需要慎重，是不是因为目标序列太短，导致出现了接头序列，此时应该可以用fastp进行接头序列切除
      根据PRJNA508621这个项目的数据，primer确实都在中间。正向primer序列前面有26bp未知特殊序列，同一个样本中这26bp几乎（96%）都一样，
      反向primer同样出现在read2的中间，其前面存在20bp的未知特殊序列，且不同样本中包含的几乎都一样。
    结论:使用cutadaptor去除primer
    
fastq数据质控: 
    fastp, 可以中间加上--correction参数, enable base correction in overlapped
    
比对：
    1. 可以直接和marker的序列进行比对
    2. 可以和参考基因组进行比对
    3. 可能会存在一些read能够比对到marker序列，但是却在参考基因组比对环节比对到非marker的序列区域，该如何处理该read？
        应该进行比对质量比较？质量相当，则保留，和marker比对的质量明显低于参考基因组，则去除？
    结论: 只需要和参考基因组比对即可,不指望merge序列后能提高比对准确率带来的收益
    
基因型确定(变异检测)：
    1. （*）如果marker都是序列较短的，比如小于<100, 两端各加上25bp的primer， 那么单端150bp应该可以完整测到，可以像microhapulator一样，直接用read信息去定型，要求一条read必须包含一个marker的所有SNP位点
    2. 可以按照常规做法，进行snp的变异检测，比如可以用gatk-haplotype-caller进行检测，然后再写脚本分析marker对应的区域的定型信息
    3. 很明显，第一种策略最直接高效，而第二种可能更适用于随机打断+探针捕获的测序技术
    结论: 直接基于比对结果,编写mhcaller.py实现基因定型

嵌合率计算:
    1. 根据基因型定性和定量结果进行计算,最简单的方式是计算多个marker的平均值

讨论microhapulator的分析策略的可能改进方案：
    1. 考虑去除primer（可以参考 https://link.springer.com/article/10.1007/s00414-021-02509-y ）这个文章提到他们设计的引物可能覆盖到了目标区域，也就有可能为某些目标SNP位点带来人为的基因型
    2. 可以用考虑用fastp软件进行merge工作
    3. 考虑是否进行序列去重？可以在fastp阶段完成去重的工作, 但对multipcr+ngs技术,这里不大适合根据序列本身进行去重. 如果有umi可以根据umi去重
    4. 不能设置过高的通用的count过滤阈值，这一点microhapulator做的比较好。PCR循环次数太高，可能会导致某些marker扩增更快，而导致有些marker的测序深度比其他的高，但是我们的分析目的不在于比较不同marker的差异，而在于关心一个marker可能出现哪些及什么比例的基因型。
    5. 继续第四点讨论，microhapulator设置了一个动态的阈值，即某个基因型的最低支持read数量不能低于该marker的总reads的2%，但这明显成了灵敏度的一个上限了，即如果混样比例低于2:100是不允许的，如果我们要达到5/1000的灵敏度，必须设置更小的值。
    6. 要提高灵敏度，需要尽量降低测序错误,
    7. 进行嵌合率分析时，先确定贡献者数量; 
    8. ACR（allele count ratio）
    9. 并不是为某个病人定制panel，那么通常不能保证所有marker都是包含有效信息的，因为marker的基因型可能在受体和供体之间完全一致。那这些marker的测序数据能不能用于背景矫正呢？
    10. 在确定基因型时的问题：
        1. 我们要提前分别测定供体和受体的基因型吗？如果是，分析起来明显更直接容易，且可以利用这些信息去区分噪音得到最大的精确度
        2. 我们一定要知道供体和受体中至少一人的基因型吗？貌似是这样的，一次测序中，似乎很难区分供体和受体的身份，除非特殊情况如可以通过性别信息推断
        3. 分析流程中，应该设计：一定要给受体的基因型信息或者测序数据作为参考，一定要给定贡献者数量，即告知DNA来源于几个人，供体的基因型提供作为可选，但是如果有多个供体，而且还需要确定身份，则需要提供供体的基因信息或数据
嵌合率计算：
    1. 如何计算？
    Percent recipient chimerism was calculated according to the formula:
    for each informative SNP
        B/ (A + B) [Donor(AA) and Recipient(BB)] or
        2B/(A + B) [Donor(AA) and Recipient(AB)]
    else:
        and then the mean calculated
    2. 把错误的read重新分配计数，或许提高比例计算的准确度。
其他：
    1. 造血干细胞移植后的嵌合率分析，不仅要知道嵌合率，还得知道dna主要贡献源是患者还是供者，那么应该是需要提取正常患者的正常细胞作为对照的
    2. 是否要考虑2次移植
    3. 如果能找到已有的造血干细胞移植数据进行验证就好了:
        https://pubmed.ncbi.nlm.nih.gov/33582770/ (所适用的marker符合microhaplotype的定义）
        https://www.nature.com/articles/s41409-022-01615-8
        https://pubmed.ncbi.nlm.nih.gov/22884060/
        Chimerism Assay Using Single Nucleotide Polymorphisms Adjacent and in Linkage-Disequilibrium Enables Sensitive Disease Relapse Monitoring after Hematopoietic Stem-Cell Transplantation

"""


def Cutadapter():
    cmd = Command()
    cmd.meta.name = 'Cutadapter'
    cmd.meta.source = 'https://cutadapt.readthedocs.io/en/stable/guide.html'
    cmd.runtime.image = 'genomicpariscentre/cutadapt:latest'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'cutadapt'
    cmd.args['cores'] = Argument(prefix='-j ', default=cmd.runtime.cpu, desc='Number of CPU cores to use. Use 0 to auto-detect')
    # hard-cut arguments, 会在adapter trimming之前发生
    cmd.args['cut-r1'] = Argument(prefix='--cut ', level='optional', desc="Remove LEN bases from each read (or R1 if paired; use -U option for R2). If LENis positive, remove bases from the beginning. If LEN is negative, remove bases from the end. Can be used twice if LENs have different signs. Applied *before* adapter trimming.")
    cmd.args['cut-r2'] = Argument(prefix='-U ', level='optional', desc='Remove LENGTH bases from R2')
    # 两端低质量碱基去除，会在adapter trimming之前发生
    cmd.args['quality-cutoff'] = Argument(prefix='-q ', level='optional', array=True, delimiter=',', desc="Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.")
    # adapter trimming 相关参数
    cmd.args['adapter3_r1'] = Argument(prefix='-a file:', type='infile', level='optional', desc="Sequence of an adapter ligated to the 3' end (paired data: of the first read).The adapter and subsequent bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only found if it is a suffix of the read.")
    cmd.args['adapter3_r2'] = Argument(prefix='-A file:', type='infile', level='optional', desc="Sequence of an adapter ligated to the 3' end (paired data: of the second read).The adapter and subsequent bases are trimmed. If a '$' character is appended ('anchoring'), the adapter is only found if it is a suffix of the read.")
    cmd.args['adapter5_r1'] = Argument(prefix='-g file:', type='infile', level='optional', desc="Sequence of an adapter ligated to the 5' end (paired data: of the first read).; The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended ('anchoring'), the adapter is only; found if it is a prefix of the read.")
    cmd.args['adapter5_r2'] = Argument(prefix='-G file:', type='infile', level='optional', desc="Sequence of an adapter ligated to the 5' end (paired data: of the second read).; The adapter and any preceding bases are trimmed. Partial matches at the 5' end are allowed. If a '^' character is prepended ('anchoring'), the adapter is only; found if it is a prefix of the read.")
    cmd.args['error-rate'] = Argument(prefix='-e ', default=0.1, desc="Maximum allowed error rate (if 0 <= E < 1), or absolute number of errors for full-length adapter match (if E is an integer >= 1). Error rate = no. of errors divided by length of matching region. Default: 0.1 (10%)")
    cmd.args['no-indels'] = Argument(prefix='--no-indels', type='bool', default=False, desc='Allow only mismatches in alignments. Default: allow both mismatches and indels')
    cmd.args['count'] = Argument(prefix='--times ', default=1, desc='Remove up to COUNT adapters from each read')
    cmd.args['overlap'] = Argument(prefix='--overlap ', default=3, desc='Require MINLENGTH overlap between read and adapter for an adapter to be found.')
    cmd.args['action'] = Argument(prefix='--action ', default='trim', desc="What to do if a match was found. trim: trim adapter and up- or downstream sequence; retain: trim, but retain adapter; mask: replace with 'N' characters; lowercase: convert to lowercase; none: leave unchanged. Default: trim")
    cmd.args['pair-adapters'] = Argument(prefix='--pair-adapters', type='bool', default=False, desc="Treat adapters given with -a/-A etc. as pairs. Either both or none are remove from each read pair.")
    # 如果read的注释中包含length信息，如使用length=150标记的，可以通过指定下面的tag进行长度信息更新
    cmd.args['length-tag'] = Argument(prefix='--length-tag ', level='optional', desc="Search for TAG followed by a decimal number in the description field of the read. Replace the decimal number with the correct length of the trimmed read. For example, use --length-tag 'length=' to correct fields like 'length=123'.")
    # 修改read名称, 在提取UMI的时候应该有用
    cmd.args['strip-suffix'] = Argument(prefix='--strip-suffix ', level='optional', desc='Remove this suffix from read names if present. Can be given multiple times.')
    cmd.args['add-prefix'] = Argument(prefix='--prefix ', level='optional', desc="Add this prefix to read names. Use {name} to insert the name of the matching adapter.")
    cmd.args['add-suffix'] = Argument(prefix='--suffix ', level='optional', desc='Add this suffix to read names; can also include {name}')
    cmd.args['rename'] = Argument(prefix='--rename ', level='optional', desc="Rename reads using TEMPLATE containing variables such as {id}, {adapter_name}")
    # 去除plolyA和N
    cmd.args['trim-ploy-a'] = Argument(prefix='--poly-a', type='bool', default=False, desc='Trim poly-A tails')
    cmd.args['trim-n'] = Argument(prefix='--trim-n', type='bool', default=False, desc="Trim N's on ends of reads.")
    # 最后的过滤
    cmd.args['min-len'] = Argument(prefix='-m ', type='int', level='optional', desc='Discard reads shorter than LEN')
    cmd.args['max-n'] = Argument(prefix='--max-n ', level='optional', desc="Discard reads with more than COUNT 'N' bases. If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.")
    cmd.args['discard-untrimmed'] = Argument(prefix='--discard-untrimmed', type='bool', default=False, desc='Discard reads that do not contain an adapter.')
    # 输出
    cmd.args['out1'] = Argument(prefix='-o ', type='outstr', desc="Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input. Summary report is sent to standard output. Use '{name}' for demultiplexing (see docs). Default: write to standard output")
    cmd.args['out2'] = Argument(prefix='-p ', type='outstr', level='optional', desc='Write R2 to FILE')
    cmd.args['info-file'] = Argument(prefix='--info-file ', type='outstr', level='optional', desc='Write information about each read and its adapter matches into FILE.')
    # The options --untrimmed-output, --discard-trimmed and -discard-untrimmed are mutually exclusive.
    cmd.args['untrimmed-output'] = Argument(prefix='--untrimmed-output ', type='outstr', level='optional', desc=' Write reads that do not contain any adapter to FILE')
    cmd.args['untrimmed-paired-output'] = Argument(prefix='--untrimmed-paired-output ', type='outstr', level='optional', desc='Write second read in a pair to this FILE when no adapter was found.')
    # 输入
    cmd.args['read1'] = Argument(prefix='', type='infile', editable=False, desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='', type='infile', editable=False, level='optional', desc='read2 fastq file')
    cmd.outputs['out1'] = Output(value='{out1}')
    cmd.outputs['out2'] = Output(value='{out2}')
    cmd.outputs['info-file'] = Output(value='{info-file}')
    cmd.outputs['untrimmed-r1'] = Output(value='{untrimmed-output}')
    cmd.outputs['untrimmed-r2'] = Output(value='{untrimmed-paired-output}')
    return cmd


def Fastp():
    cmd = Command()
    cmd.meta.name = 'Fastp'
    cmd.meta.source = 'https://github.com/OpenGene/fastp'
    cmd.meta.version = '0.23.4'
    cmd.runtime.image = 'gudeqing/gatk4.3-bwa-fastp-gencore-mutscan:1.0'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'fastp'
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', editable=False, desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', editable=False, level='optional', desc='read2 fastq file')
    cmd.args['threads'] = Argument(prefix='-w ', default=cmd.runtime.cpu, desc='thread number')
    cmd.args['min_length'] = Argument(prefix='-l ', default=35, desc='reads shorter than length_required will be discarded')
    cmd.args['correction'] = Argument(prefix='--correction', type='bool', default=False, desc='enable base correction in overlapped regions')
    cmd.args['overlap_diff_percent_limit'] = Argument(prefix='--overlap_diff_percent_limit ', default=20, desc='The maximum number of mismatched bases to detect overlapped region of PE reads')
    cmd.args['dedup'] = Argument(prefix='--dedup ', type='bool', default=False, desc='enable deduplication to drop the duplicated reads/pairs')
    cmd.args['trim_front1'] = Argument(prefix='--trim_front1 ', level='optional', desc='trimming how many bases in front for read1')
    cmd.args['disable_quality_filtering'] = Argument(prefix='-Q', type='bool', default=False, desc='quality filtering is enabled by default. If this option is specified, quality filtering is disabled')
    cmd.args['disable_adapter_trimming'] = Argument(prefix='-A', type='bool', default=False, desc='adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled')
    cmd.args['out1'] = Argument(prefix='-o ', type='outstr', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', level='optional', type='outstr', desc='clean read2 output fastq file')
    cmd.args['html'] = Argument(prefix='-h ', type='outstr', desc='html report file')
    cmd.args['json'] = Argument(prefix='-j ', type='outstr', desc='json report file')
    cmd.outputs['out1'] = Output(value="{out1}")
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}")
    cmd.outputs['json'] = Output(value="{json}")
    return cmd


def BwaMem(sample, platform, bwamem2=False):
    cmd = Command()
    cmd.meta.name = 'BwaMem2' if bwamem2 else "BwaMem"
    cmd.runtime.image = 'gudeqing/gatk4.3-bwa-fastp-gencore-mutscan:1.0'
    if bwamem2:
        cmd.meta.version = '2.2.1_x64-linux'
        cmd.runtime.tool = 'bwa-mem2 mem -M -Y -v 3'
    else:
        cmd.meta.version = '0.7.17-r1188'
        cmd.runtime.tool = 'bwa mem -M -Y -v 3'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 8
    cmd.args['include_read_header'] = Argument(prefix='-C', type='bool', default=False, desc='Append FASTA/FASTQ comment to SAM output')
    cmd.args['readgroup'] = Argument(prefix='-R ', desc='read group info', value=f'"@RG\\tID:{sample}\\tSM:{sample}\\tPL:{platform}"')
    cmd.args['k'] = Argument(prefix='-K ', default=100000000, desc='process INT input bases in each batch regardless of nThreads (for reproducibility)')
    cmd.args['t'] = Argument(prefix='-t ', default=cmd.runtime.cpu, desc='number of threads to use in computation')
    cmd.args['ref'] = Argument(prefix='', type='infile', format='fasta', desc='reference fasta file')
    cmd.args['read1'] = Argument(prefix='', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='', level='optional', type='infile', desc='read2 fastq file')
    cmd.args['_fix'] = Argument(type='fix', value=f'| samtools view -O BAM --threads {cmd.runtime.cpu} - ')
    cmd.args['out'] = Argument(prefix='> ', value=f'{sample}.unmerged.bam', type='outstr', desc='output bam file')
    cmd.outputs['out'] = Output(value="{out}", format='bam')
    return cmd


def MergeSamFiles():
    cmd = Command()
    cmd.meta.name = 'MergeSamFiles'
    cmd.meta.version = '4.3.0.0'
    cmd.meta.desc = 'Merges multiple SAM and/or BAM files into a single file, and sorting and indexing are optional'
    cmd.runtime.image = 'gudeqing/gatk4.3-bwa-fastp-gencore-mutscan:1.0'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'gatk --java-options -Xmx16g MergeSamFiles'
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', multi_times=True, desc='SAM or BAM input file', format='bam')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', type='outstr', desc='SAM or BAM file to write merged result')
    cmd.args['CREATE_INDEX'] = Argument(prefix='--CREATE_INDEX ', default='true', range=['true', 'false'], desc='Whether to create a BAM index when writing a coordinate-sorted BAM file.')
    cmd.args['SORT_ORDER'] = Argument(prefix='--SORT_ORDER ', default='coordinate', range=['unsorted', 'queryname', 'coordinate', 'duplicate', 'unknown'], desc='Sort order of output file')
    cmd.outputs['out'] = Output(value='{OUTPUT}', format='bam')
    return cmd


def CollectSequencingArtifactMetrics():
    cmd = Command()
    cmd.meta.name = 'GetSeqErrorMetrics'
    cmd.meta.version = '4.3.0.0'
    cmd.meta.desc = 'Collect metrics to quantify single-base sequencing artifacts.'
    cmd.runtime.image = 'gudeqing/gatk4.3-bwa-fastp-gencore-mutscan:1.0'
    cmd.runtime.memory = 8 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'gatk --java-options -Xmx8g CollectSequencingArtifactMetrics'
    cmd.args['INPUT'] = Argument(prefix='--INPUT ', type='infile', desc='SAM or BAM input file', format='bam')
    cmd.args['OUTPUT'] = Argument(prefix='--OUTPUT ', type='outstr', desc='prefix of file to write the output to')
    cmd.args['ref'] = Argument(prefix='--REFERENCE_SEQUENCE ', type='infile', desc='Reference sequence file')
    cmd.args['db_snp'] = Argument(prefix='--DB_SNP ', type='infile', level='optional', desc='VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis')
    cmd.outputs['out'] = Output(value='{OUTPUT}.error_summary_metrics')
    cmd.outputs['bait_bias_detail_metrics'] = Output(value='{OUTPUT}.bait_bias_detail_metrics')
    cmd.outputs['bait_bias_summary_metrics'] = Output(value='{OUTPUT}.bait_bias_summary_metrics')
    cmd.outputs['pre_adapter_detail_metrics'] = Output(value='{OUTPUT}.pre_adapter_detail_metrics')
    cmd.outputs['pre_adapter_summary_metrics'] = Output(value='{OUTPUT}.pre_adapter_summary_metrics')
    return cmd


def TypingMicroHap():
    cmd = Command()
    cmd.meta.name = 'TypingMicroHap'
    cmd.meta.source = 'in-house'
    cmd.meta.desc = "Microhaplotyp typing tool"
    cmd.runtime.image = 'gudeqing/gatk4.3-bwa-fastp-gencore-mutscan:1.0'
    cmd.runtime.memory = 12 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'python'
    cmd.args['script'] = Argument(prefix='', type='infile', value=f'{script_path}/microhap/mhcaller.py', desc='script path')
    cmd.args['_sub_command'] = Argument(prefix='', type='fix', value='micro_hap_caller')
    cmd.args['bam'] = Argument(prefix='-bam_file ', type='infile', format='bam', desc='input bam file')
    cmd.args['genome'] = Argument(prefix='-genome_file ', type='infile', format='fasta', desc='input indexed genome fasta file')
    cmd.args['micro_haps'] = Argument(prefix='-micro_hap_file ', type='infile', format='txt', desc='micro-haplotype definition file')
    cmd.args['error_rate_file'] = Argument(prefix='-error_rate_file ', type='infile', level='optional', format='txt', desc='result from CollectSequencingArtifactMetrics')
    cmd.args['out_prefix'] = Argument(prefix='-out_prefix ', type='outstr', desc='output file prefix')
    cmd.outputs['out'] = Output(value='{out_prefix}.csv', format='csv', desc='output typing profile', report=True)
    cmd.outputs['out_json'] = Output(value='{out_prefix}.json', format='json', report=True)
    return cmd


def Chimerism():
    cmd = Command()
    cmd.meta.name = 'Chimerism'
    cmd.meta.source = 'in-house'
    cmd.meta.desc = "Chimerism calculation based on microhaplotype typing result"
    cmd.runtime.image = 'gudeqing/gatk4.3-bwa-fastp-gencore-mutscan:1.0'
    cmd.runtime.memory = 12 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'python'
    cmd.args['script'] = Argument(prefix='', type='infile', value=f'{script_path}/microhap/mhcaller.py', desc='script path')
    cmd.args['_sub_command'] = Argument(prefix='', type='fix', value='chemerism')
    cmd.args['donor_profile'] = Argument(prefix='-donor_profile ', type='infile', desc="donor's micro-haplotype typing result")
    cmd.args['recipient_profile'] = Argument(prefix='-recipient_profile ', type='infile', desc="recipient's micro-haplotype typing result")
    cmd.args['test_profile'] = Argument(prefix='-test_profile ', type='infile', desc="testing sample's micro-haplotype typing result")
    cmd.args['out'] = Argument(prefix='-out_json ', type='outstr', desc='output json file name')
    cmd.outputs['out'] = Output(value='{out}', format='json', report=True)
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'ctDNA'
    wf.meta.source = ""
    wf.meta.desc = \
    """
    1. 原始数据处理: primer切除和质控
    2. 比对
    3. 基因定型
    4. 嵌合率计算
    """
    wf.meta.version = "1.0"
    # 定义流程输入参数
    wf.init_argparser()
    # fastq 输入参数
    wf.add_argument('-fastq_info', nargs='+', required=True,
                    help='A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json]')
    wf.add_argument('-r1_name', default='(.*)_1.fastq.gz',
                    help="python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'")
    wf.add_argument('-r2_name', default='(.*)_2.fastq.gz',
                    help="python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'")
    wf.add_argument('-exclude_samples', default=tuple(), nargs='+', help='samples to exclude from analysis')
    # 其他输入文件
    wf.add_argument('-microhaps', help="micro-haplotype definition file for target region")
    wf.add_argument('-ref', default='/home/hxbio04/biosofts/MicroHapulator/microhapulator/data/hg38.fasta', help='reference fasta file')
    wf.add_argument('-db_snp', default='/home/hxbio04/dbs/SNPdb/dbsnp_146.hg38.vcf.gz', help='VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis')
    wf.add_argument('-donor_name', required=False, help='donor sample name, if donor profile not provided, we wil start from fastq file to get the profile')
    wf.add_argument('-recipient_name', required=False, help='recipient sample name, if recipient profile not provided, we wil start from fastq file to get the profile')
    wf.add_argument('-donor_profile', required=False, help='donor micro-hap typing result')
    wf.add_argument('-recipient_profile', required=False, help='donor micro-hap typing result')
    wf.add_argument('-forward_primer', required=True, help="read1 5'end primer")
    wf.add_argument('-reverse_primer', required=True, help="read2 5'end primer")
    # 收集参数
    wf.parse_args()
    top_vars = dict(
        ref=TopVar(value=wf.args.ref, type='infile'),
        microhaps=TopVar(value=wf.args.microhaps, type='infile'),
        donor_profile=TopVar(value=wf.args.donor_profile, type='infile'),
        recipient_profile=TopVar(value=wf.args.recipient_profile, type='infile'),
        donor_name=TopVar(value=wf.args.donor_name, type='str'),
        recipient_name=TopVar(value=wf.args.recipient_name, type='str'),
        forward_primer=TopVar(value=wf.args.forward_primer, type='infile'),
        reverse_primer=TopVar(value=wf.args.reverse_primer, type='infile'),
        db_snp=TopVar(value=wf.args.db_snp, type='infile'),
    )
    wf.add_topvars(top_vars)

    # 提取fastq信息
    fastq_info = get_fastq_info(fastq_info=wf.args.fastq_info, r1_name=wf.args.r1_name, r2_name=wf.args.r2_name)
    if len(fastq_info) <= 0:
        raise Exception('No fastq file found !')

    # ---fastq信息添加到topvar 方便流程转化---
    for sample, reads in fastq_info.items():
        if len(reads) == 2:
            r1s, r2s = reads
        else:
            r1s = reads[0]
            r2s = [None] * len(r1s)
        new_r1s = []
        new_r2s = []
        for ind, (r1, r2) in enumerate(zip(r1s, r2s), start=1):
            r1_name = sample + f'_R1_{ind}' if len(r1s) > 1 else sample + '_R1'
            r1_topvar = TopVar(value=r1, name=r1_name, type='infile')
            new_r1s.append(r1_topvar)
            r2_name = sample + f'_R2_{ind}' if len(r1s) > 1 else sample + '_R2'
            r2_topvar = TopVar(value=r2, name=r2_name, type='infile')
            new_r2s.append(r2_topvar)
            wf.topvars[r1_name] = r1_topvar
            wf.topvars[r2_name] = r2_topvar
        fastq_info[sample] = [new_r1s, new_r2s]
    # ----end of fastq info adding to topvars----

    # 一些验证
    donor_name = wf.topvars['donor_name'].value
    if donor_name and (donor_name not in fastq_info):
        raise Exception(f'{donor_name} is not found in samples')
    recipient_name = wf.topvars['recipient_name'].value
    if recipient_name and (recipient_name not in fastq_info):
        raise Exception(f'{recipient_name} is not found in samples')

    # 开始数据处理
    typing_tasks = []
    chimerism_tasks = []
    for sample, reads in fastq_info.items():
        # 跳过不需要分析的样本
        if sample in wf.args.exclude_samples:
            print(f'Skip {sample} for it is in excluded sample list')
            continue
        # fastq预处理和比对 考虑一个样本存在多对fastq的处理
        if sample in wf.args.exclude_samples:
            continue
        if len(reads) == 2:
            r1s, r2s = reads
        else:
            r1s = reads[0]
            r2s = [None]*len(r1s)

        # 如果一个样本有多组fastq，将分别处理，最后合并bam
        bwa_tasks = []
        for ind, (r1, r2) in enumerate(zip(r1s, r2s)):
            uniq_tag = f'{sample}-{ind}' if len(r1s) > 1 else sample
            # cutdapter去除primer
            cutadapter_task, args = wf.add_task(Cutadapter(), tag=uniq_tag)
            args['adapter5_r1'].value = wf.topvars['forward_primer']
            args['adapter5_r2'].value = wf.topvars['reverse_primer']
            args['pair-adapters'].value = True
            args['error-rate'].value = 0.25
            args['overlap'].value = 5
            # 丢掉不包含primer的
            args['discard-untrimmed'].value = True
            args['length-tag'].value = "length="
            args['read1'].value = r1
            args['read2'].value = r2
            args['info-file'].value = uniq_tag + '.primer.cutInfo.txt'
            args['out1'].value = uniq_tag + '.noPrimer.R1.fq.gz'
            args['out2'].value = uniq_tag + '.noPrimer.R2.fq.gz'

            # fastp QC
            fastp_task, args = wf.add_task(Fastp(), tag=uniq_tag)
            args['correction'].value = True
            args['disable_quality_filtering'].value = True
            args['disable_adapter_trimming'].value = True
            args['read1'].value = cutadapter_task.outputs['out1']
            args['out1'].value = f'{sample}.clean.R1.fq.gz'
            if r2 is not None:
                args['read2'].value = cutadapter_task.outputs['out2']
                args['out2'].value = f'{sample}.clean.R2.fq.gz'
            args['html'].value = f'{sample}.fastp.html'
            args['json'].value = f'{sample}.fastp.json'

            bwa_task, args = wf.add_task(BwaMem(uniq_tag, 'Illumina'), tag=uniq_tag)
            args['read1'].value = fastp_task.outputs['out1']
            args['read2'].value = fastp_task.outputs['out2']
            args['ref'].value = wf.topvars['ref']
            bwa_tasks.append(bwa_task)

        # 合并比对结果
        merge_bam_task, args = wf.add_task(MergeSamFiles(), tag=sample, depends=bwa_tasks)
        args['INPUT'].value = [x.outputs['out'] for x in bwa_tasks]
        args['SORT_ORDER'].value = 'coordinate'
        args['OUTPUT'].value = sample + '.sorted.bam'

        # get seq error metrics
        error_metrics_task, args = wf.add_task(CollectSequencingArtifactMetrics(), tag=sample)
        args['INPUT'].value = merge_bam_task.outputs['out']
        args['OUTPUT'].value = sample
        args['ref'].value = wf.topvars['ref']
        args['db_snp'].value = wf.topvars['db_snp']

        # typing
        typing_task, args = wf.add_task(TypingMicroHap(), tag=sample)
        args['genome'].value = wf.topvars['ref']
        args['bam'].value = merge_bam_task.outputs['out']
        args['micro_haps'].value = wf.topvars['microhaps']
        args['out_prefix'].value = sample + '.typing'
        args['error_rate_file'].value = error_metrics_task.outputs['out']
        typing_tasks.append(typing_task)

        if wf.topvars['donor_profile'].value:
            chimerism_task, args = wf.add_task(Chimerism(), tag=sample)
            chimerism_tasks.append(chimerism_task)
            args['donor_profile'].value = wf.topvars['donor_profile']
            args['recipient_profile'].value = wf.topvars['recipient_profile']
            args['test_profile'].value = typing_task.outputs['out']
            args['out'].value = sample + '.chimerism.json'

    # chimerism
    if donor_name and recipient_name:
        donor_typing_task = [x for x in typing_tasks if x.tag == donor_name][0]
        recipient_typing_task = [x for x in typing_tasks if x.tag == recipient_name][0]
        for each in typing_tasks:
            chimerism_task, args = wf.add_task(Chimerism(), tag=each.tag)
            chimerism_tasks.append(chimerism_task)
            args['donor_profile'].value = donor_typing_task.outputs['out']
            args['recipient_profile'].value = recipient_typing_task.outputs['out']
            args['test_profile'].value = each.outputs['out']
            args['out'].value = each.tag + '.chimerism.json'

    wf.run()
    if wf.success:
        result = dict()
        for task in chimerism_tasks:
            chimer_json = os.path.join(task.wkdir, task.outputs['out'].value)
            with open(chimer_json) as f:
                result[task.tag] = json.load(f)
        df = pd.DataFrame(result)
        exp_chimerism = []
        for each in df.columns:
            ratio_tag = each.rsplit('-', 1)[1]
            if 'v' in ratio_tag:
                x, y = ratio_tag.split('v')
                exp_chimerism.append(int(y)/(int(x) + int(y)))
            else:
                if ratio_tag in ["1", "3"]:
                    exp_chimerism.append(0)
                elif ratio_tag in ["2", "4"]:
                    exp_chimerism.append(1)
        df.loc['exp_chimerism'] = exp_chimerism
        df.to_csv(os.path.join(wf.wkdir, 'Report', "all.chimerism.csv"))


if __name__ == '__main__':
    pipeline()

