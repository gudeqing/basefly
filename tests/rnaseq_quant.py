import os
from basefly import Argument, Output, Command, Workflow, TopVar, get_fastq_info

__author__ = 'gdq'

"""
提示:
0. 写workflow时，参数赋值规范建议：args[X].value = TopVar[?] | task.Outputs[?] | TmpVar()
*. 如果不是为了写wdl流程，可以不使用TmpVar，直接赋值就ok
1. 一定要正确定义参数的类型, type is one of ['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix']
    其中‘fix'可以用于表示命令行中的固定字符串或固定参数, 如 “bwa xxx | samtools sort -" 中的‘| samtools sort -’ 可以用fix固定
2. 参数的添加顺序对于命令行的正确形成很重要，这里字典的有序性得到利用
3. 定义output时，value(或path）属性对应的值可以直接用{}引用cmd.args的key，

关于meta：
name: 定义命令行的名称，会参与具体task的name的形成，建议组成：[数字，字母，’-‘], 下划线会自动被替换为中划线’-‘
其他字段都是描述工具的开发作者(author)，链接(source)，版本号(version)，简介（desc)
"""

All_IN_ONE_IMAGE = "gudeqing/rnaseq_envs:1.4"


def fastp():
    cmd = Command()
    cmd.meta.name = 'Fastp'
    cmd.meta.source = 'https://github.com/OpenGene/fastp'
    cmd.meta.version = '0.23.2'
    cmd.meta.desc = 'A fast Fastq QC tool'
    cmd.runtime.image = 'gudeqing/rnaseq_envs:1.4'
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'fastp'
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', editable=False, desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', editable=False, level='optional', desc='read2 fastq file')
    cmd.args['threads'] = Argument(prefix='-w ', default=cmd.runtime.cpu, desc='thread number')
    cmd.args['min_length'] = Argument(prefix='-l ', default=35, desc='reads shorter than length_required will be discarded')
    cmd.args['correction'] = Argument(prefix='--correction', type='bool', default=True, desc='enable base correction in overlapped regions')
    cmd.args['overlap_diff_percent_limit'] = Argument(prefix='--overlap_diff_percent_limit ', default=10, desc='The maximum number of mismatched bases to detect overlapped region of PE reads')
    cmd.args['dedup'] = Argument(prefix='--dedup ', type='bool', default=False, desc='enable deduplication to drop the duplicated reads/pairs')
    cmd.args['trim_front1'] = Argument(prefix='--trim_front1 ', level='optional', desc='trimming how many bases in front for read1')
    cmd.args['out1'] = Argument(prefix='-o ', type='outstr', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', level='optional', type='outstr', desc='clean read2 output fastq file')
    cmd.args['html'] = Argument(prefix='-h ', type='outstr', desc='html report file')
    cmd.args['json'] = Argument(prefix='-j ', type='outstr', desc='json report file')
    cmd.outputs['out1'] = Output(value="{out1}")
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}", report=True)
    cmd.outputs['json'] = Output(value="{json}", report=True)
    return cmd


def salmon_index():
    # 定义一个command
    cmd = Command()
    cmd.meta.name = 'SalmonIndex'
    cmd.meta.source = 'https://github.com/COMBINE-lab/salmon/'
    cmd.meta.version = '1.5.2'
    cmd.meta.desc = 'build salmon index, transcript fasta from Gencode is recommended'
    cmd.runtime.image = All_IN_ONE_IMAGE
    cmd.runtime.memory = 6 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'salmon index'
    cmd.args['index'] = Argument(prefix='-i ', type='str', default='salmon_index', desc='salmon index directory name')
    cmd.args['transcripts'] = Argument(prefix='-t ', type='infile', desc='Transcript fasta file')
    cmd.args['gencode'] = Argument(prefix='--gencode', type='bool', default=True, desc='This flag will expect the input transcript; fasta to be in GENCODE format, and will split the transcript name at the first "|" character; These reduced names will be used in the output and when looking for these transcripts in a gene to transcript GTF.')
    cmd.args['threads'] = Argument(prefix='--threads ', default=cmd.runtime.cpu, desc='Number of threads to use during indexing')
    cmd.args['decoys'] = Argument(prefix='--decoys ', type='infile', level='optional', desc='Treat these sequences ids from the reference as the decoys that may have sequence homologous to some known transcript.')
    cmd.outputs['index_dir'] = Output(path='{index}', type='outdir')
    return cmd


def salmon_quant():
    # 定义一个command
    cmd = Command()
    cmd.meta.name = 'SalmonQuant'
    cmd.meta.source = 'https://github.com/COMBINE-lab/salmon/'
    cmd.meta.version = '1.5.2'
    cmd.meta.desc = 'transcript expression quantification'
    cmd.runtime.image = All_IN_ONE_IMAGE
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'salmon quant'
    cmd.args['libType'] = Argument(prefix='--libType ', default='A', desc='默认自动判定文库类型')
    cmd.args['indexDir'] = Argument(prefix='-i ', type='indir', desc='transcript fasta index directory')
    cmd.args['read1'] = Argument(prefix='-1 ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-2 ', type='infile', desc='read2 fastq file')
    cmd.args['outDir'] = Argument(prefix='-o ', type='outstr', default='quant', desc='output directory')
    cmd.args['gcBias'] = Argument(prefix='--gcBias ', type='bool', default=True, desc='perform gc Bias correction')
    cmd.outputs['transcript'] = Output(path="{outDir}" + "/quant.sf")
    cmd.outputs['outDir'] = Output(path="{outDir}", type='outdir')
    return cmd


def quant_merge():
    cmd = Command()
    cmd.meta.name = 'quantMerge'
    cmd.meta.source = 'https://github.com/COMBINE-lab/salmon/'
    cmd.meta.version = '1.5.2'
    cmd.meta.desc = 'Merge multiple quantification results into a single file'
    cmd.runtime.image = All_IN_ONE_IMAGE
    cmd.runtime.tool = 'salmon quantmerge'
    # 下面的quants参数对应的是目录，可以接收多个值作为输入
    cmd.args['quants'] = Argument(prefix="--quants ", array=True, type='indir', desc='salmon quant dir list')
    cmd.args['names'] = Argument(prefix='--names ', array=True, level='optional', desc='sample names')
    cmd.args['column'] = Argument(prefix='--column ', default='TPM', desc='indicate which column to merge, default: TPM')
    cmd.args['genes'] = Argument(prefix='--genes ', type='bool', default=False, desc='indicate to merge at gene level')
    cmd.args['out'] = Argument(prefix='--output ', type='outstr', default=f'merged.{cmd.args["column"].default}.txt', desc='output file name')
    cmd.outputs['result'] = Output(path="{out}", report=True)
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'RNASeqQuant'
    wf.meta.source = "basefly"
    wf.meta.version = "0.1"
    wf.meta.function = 'Transcript expression Quantification'
    wf.meta.desc = """
    Fast transcript quantification using fastp and salmon
    you are required to input:
    1. fastq information
    2. A transcriptome fasta file
    
    # 使用示例：
    1. 使用docker：
    python rnaseq_quant.py -fastq_info testdata/fastqs/ -transcripts testdata/transcripts/transcripts.fa -outdir Result --run --plot --docker
    
    2. 使用本地环境（在云平台可以根据docker镜像创建一个本地环境：docker run -v /home/hxbio04/basefly/:/home/hxbio04/basefly/ -w $PWD --rm -it gudeqing/rnaseq_envs:1.4)
    python rnaseq_quant.py -fastq_info testdata/fastqs/ -transcripts testdata/transcripts/transcripts.fa -outdir Result --run --plot
    
    3. 在云平台，假设挂载 （1）流程数据集（scripts_dataset）（2）待分析数据集（testdata_dataset） （3）存储参考数据库的数据集
    启动实例后： 
    python /enigma/datasets/*/rnaseq_quant.py -r1_name '(.*).R1.fastq.gz' -r2_name '(.*).R2.fastq.gz' -transcripts transcripts.fa -outdir Result --run --plot
    注意：
        1. 上面“/enigma/datasets/*/rnaseq_quant.py”中的”*“是通配符，在不确定流程数据集的具体路径的时候派上用场
        2. -outdir：指定输出目录，如果目录不存在会尝试创建，该目录是流程的工作目录，也是结果的输出目录
        3. 默认流程会去遍历整个/enigma/datasets/目录，找到能被参数r1_name和r2_name匹配到文件作为待分析数据文件
        4. -transcripts参数是指定分析时需要的参考数据库，如果不提供完整路径，默认会去遍历整个/enigma/datasets/目录，找到该参数指定的文件作为参考输入
        5. --plot表示运行时同时绘制流程运行状态图，此状态图路径为 outdir + "/state.svg"，建议能够让用户实时预览该文件以查看流程的运行情况
        6. outdir + "/logs/" 为流程运行日志文件目录，建议能够让用户实时查看该目录以了解流程的运行情况，方便debug
        7. outdir + "wf.*.running.*.log"路径为流程运行日志，建议能够让用户实时预览该文件以查看流程的运行情况，默认流程会print这个日志信息
    """

    # 定义流程输入参数
    # 初始化流程输入参数解析器
    wf.init_argparser()
    wf.add_argument('-fastq_info', nargs='+', default='/enigma/datasets/',
                    help='A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json]')
    wf.add_argument('-r1_name', default='(.*).R1.fastq.gz',
                    help="python regExp that describes the full name of read1 fastq file name. "
                         "It requires at least one pair small brackets, "
                         "and the string matched in the first pair brackets will be used as sample name. "
                         "Example: '(.*).R1.fq.gz'")
    wf.add_argument('-r2_name', default='(.*).R2.fastq.gz',
                    help="python regExp that describes the full name of read2 fastq file name. "
                         "It requires at least one pair small brackets,"
                         " and the string matched in the first pair brackets will be used as sample name. "
                         "Example: '(.*).R2.fq.gz'")
    wf.add_argument('-exclude_samples', default=tuple(), nargs='+', help='specify samples to be excluded from analysis')
    wf.add_argument('-transcripts', default='transcripts.fa', help="transcriptome fasta file, if only a file name (but not full path) provided, we will try to find it in '/enigma/datasets/'")
    # 开始收集参数
    wf.parse_args()

    # 进行具体的流程编排任务
    # 提取fastq信息
    fastq_info = get_fastq_info(fastq_info=wf.args.fastq_info, r1_name=wf.args.r1_name, r2_name=wf.args.r2_name)
    if len(fastq_info) <= 0:
        raise Exception('No fastq file found !')

    # check输入文件transcripts
    if not os.path.exists(wf.args.transcripts):
        # /enigma/datasets/搜索文件
        for root, dirs, files in os.walk('/enigma/datasets/', followlinks=True):
            if wf.args.transcripts in files:
                wf.args.transcripts = os.path.join(root, wf.args.transcripts)
                break

    # 考虑将流程的输入文件信息放入到TopVar, 可以方便将来把流程转换为其他格式的流程，当然，你也可以选择不这么做
    wf.topvars = dict(
        transcripts=TopVar(value=wf.args.transcripts, type='infile'),
    )

    # index
    index_task, args = wf.add_task(salmon_index())
    args['transcripts'].value = wf.topvars['transcripts']

    # fastp and quant
    merge_depends = []
    for sample, (r1, r2) in fastq_info.items():
        # 向流程中添加task
        fastp_task, args = wf.add_task(fastp(), tag=sample, parent_wkdir='fastp')
        args['read1'].value = os.path.abspath(r1[0])  # 如果read1有多个文件，也仅仅使用第一个作为输入，暂不设计处理多个的情况
        args['read2'].value = os.path.abspath(r2[0])
        args['out1'].value = f'{sample}.clean.R1.fq'
        args['out2'].value = f'{sample}.clean.R2.fq'
        args['html'].value = f'{sample}.fastp.html'
        args['json'].value = f'{sample}.fastp.json'

        # 添加第二个task
        task, args = wf.add_task(salmon_quant(), tag=sample, parent_wkdir='salmon')
        task.depends = [fastp_task.task_id]
        args['read1'].value = fastp_task.outputs["out1"]
        args['read2'].value = fastp_task.outputs["out2"]
        args['indexDir'].value = index_task.outputs['index_dir']
        args['outDir'].value = sample
        # 收集多个目标task，作为其他的输入
        merge_depends.append(task.task_id)

    # merge transcript TPM, 下面添加的任务名称会是cmd.meta.name + '-' + str(tag)
    task, args = wf.add_task(quant_merge(), tag='TPM')
    # 由于同一个command被使用多次，这里要记得重新命名, 这样生成的wdl task才不会漏掉
    task.cmd.meta.name = 'MergeTranscriptTPM'
    task.depends = merge_depends
    # 针对当前任务，column 和 genes 参数需要固定，不允许修改,用 'fix' 表示
    args['column'].value = 'TPM'
    args['column'].type = 'fix'
    args['genes'].value = False
    args['genes'].type = 'fix'
    args['quants'].value = [wf.tasks[task_id].outputs['outDir'] for task_id in task.depends]

    # merge transcript Count
    task, args = wf.add_task(quant_merge(), tag='Count')
    task.cmd.meta.name = 'MergeTranscriptCount'
    task.depends = merge_depends
    # 针对当前任务，column 和 genes 参数需要固定，不允许修改
    args['column'].value = 'NumReads'
    args['column'].type = 'fix'
    args['genes'].value = False
    args['genes'].type = 'fix'
    args['quants'].value = [wf.tasks[task_id].outputs['outDir'] for task_id in task.depends]
    # 由于生成wdl的task时，优先使用参数的默认值，所以下面直接修改默认值
    args['out'].default = f'merged.{args["column"].value}.txt'

    # 流程编排结束，运行流程
    wf.run()


if __name__ == '__main__':
    pipeline()
