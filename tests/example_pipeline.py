import sys; sys.path.append('..')
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, TmpVar, get_fastq_info

"""
注意:
pycharm里设置code completion 允许“suggest variable and parameter name”, 可以极大方便流程编写
0. 写workflow时，参数赋值规范建议：args[X].value = TopVar[?] | task.Outputs[?] | TmpVar()
*. 如果不是为了写wdl流程，可以不使用TmpVar，直接赋值就ok
1. 一定要正确定义参数的类型, type is one of ['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix']
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
"""


def fastp():
    cmd = Command()
    cmd.meta.name = 'fastp'
    cmd.runtime.image = 'gudeqing/fastp:0.21.0'
    cmd.runtime.tool = 'fastp'
    # 可以直接用访问属性的方式添加参数，这个得益于使用Munch对象而不是原生字典
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', desc='read2 fastq file')
    # 当然，可以直接用字典的方式添加参数
    cmd.args['out1'] = Argument(prefix='-o ', type='str', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', type='str', desc='clean read2 output fastq file')
    # 下面的outputs设置起初是为了能够生成wdl设置,
    cmd.outputs['out1'] = Output(path="{out1}")  # 这里使用”{}“引用其他Argument对象作为输入
    cmd.outputs['out2'] = Output(path="{out2}")
    return cmd


def salmon():
    # 定义一个command
    cmd = Command()
    cmd.meta.name = 'salmon'
    cmd.meta.desc = 'transcript expression quantification'
    cmd.runtime.image = "combinelab/salmon:latest"
    cmd.runtime.memory = 2*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'salmon quant'
    cmd.args['libType'] = Argument(prefix='--libType ', default='A')
    cmd.args['indexDir'] = Argument(prefix='-i ', type='indir', desc='transcript fasta index directory')
    cmd.args['read1'] = Argument(prefix='-1 ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-2 ', type='infile', desc='read2 fastq file')
    cmd.args['outDir'] = Argument(prefix='-o ', type='str', default='quant', desc='output directory')
    cmd.args['gcBias'] = Argument(prefix='--gcBias ', type='bool', default=True, desc='perform gc Bias correction')
    cmd.outputs['transcript'] = Output(path="{outDir}" + "/quant.sf")
    cmd.outputs['outDir'] = Output(path="{outDir}", type='outdir')
    return cmd


def quant_merge():
    cmd = Command()
    cmd.meta.name = 'quantMerge'
    cmd.meta.desc = 'Merge multiple quantification results into a single file'
    cmd.runtime.image = "combinelab/salmon:latest"
    cmd.runtime.tool = 'salmon quantmerge'
    # 下面的quants参数对应的是目录，所以type='indir'
    cmd.args['quants'] = Argument(prefix="--quants ", array=True, type='indir', desc='salmon quant dir list')
    cmd.args['names'] = Argument(prefix='--names ', array=True, level='optional')
    cmd.args['column'] = Argument(prefix='--column ', default='TPM')
    cmd.args['genes'] = Argument(prefix='--genes ', type='bool', default=False)
    cmd.args['out'] = Argument(prefix='--output ', default=f'merged.{cmd.args["column"].default}.txt')
    cmd.outputs['result'] = Output(path="{out}")
    return cmd


def pipeline():
    fastq_info = get_fastq_info(fastq_dirs=('testdata/',), r1_name='(.*).R1.fastq', r2_name='(.*).R2.fastq')
    top_vars = dict(
        index_dir=TopVar(value='testdata/index/', type='indir'),
    )

    wf = Workflow(top_vars=top_vars)
    wf.meta.name = 'PipelineExample'
    wf.meta.desc = 'This is a simple pipeline for fast gene/transcript quantification. '
    wf.meta.desc += 'workflow = [fastq -> Fastp -> Salmon]'

    merge_depends = []
    for sample, (r1, r2) in fastq_info.items():
        r1 = r1[0]
        r2 = r2[0]
        # 向流程中添加task
        task, args = wf.add_task(fastp(), name='fastp-'+sample)
        # 可随意带入任何信息，如样本信息
        task.sample = sample
        # 给task分组信息，同一批次循环中的task属于同一组，这对于wdl的scatter转换非常重要
        task.group = 'batch1'
        task_id = task.task_id
        args['read1'].value = TmpVar(name='read1', value=r1, type='infile')
        args['read2'].value = TmpVar(name='read2', value=r2, type='infile')
        args['out1'].value = TmpVar(name='~{sample}.clean.R1.fq', value=f'{sample}.clean.R1.fq', type='str')
        args['out2'].value = TmpVar(name='~{sample}.clean.R2.fq', value=f'{sample}.clean.R2.fq', type='str')

        depend_task = task
        task, args = wf.add_task(salmon(), name='salmon-'+sample)
        task.depends = [task_id]
        task.sample = sample
        task.group = 'batch1'
        args['read1'].value = depend_task.outputs["out1"]
        args['read2'].value = depend_task.outputs["out2"]
        args['indexDir'].value = top_vars['index_dir']
        args['outDir'].value = TmpVar(name="sample", value=sample, type='str')
        # 上面的sample只是一个普通的字符串，所以必须添加wdl属性用以辅助wdl生成
        merge_depends.append(task.task_id)

    # merge transcript TPM
    task, args = wf.add_task(quant_merge())
    # 由于同一个command被使用多次，这里要记得重新命名，不然生成的wdl中的task可能命名重复
    task.cmd.meta.name = 'MergeTranscriptTPM'
    task.depends = merge_depends
    # 针对当前任务，column 和 genes 参数需要固定，不允许修改,用 'fix' 表示，便于wdl的正确生成。
    args['column'].value = 'TPM'
    args['column'].type = 'fix'
    args['genes'].value = False
    args['genes'].type = 'fix'
    args['quants'].value = [wf.tasks[task_id].outputs['outDir'] for task_id in task.depends]

    # merge transcript Count
    task, args = wf.add_task(quant_merge())
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

    for task_id, task in wf.tasks.items():
        pass
        # print(task.task_id)
        # print(task.outputs)
        # print(task.cmd.format_cmd(wf.tasks))
        # for line in task.argo_template(wf.tasks):
        #     print(line)

    wf.to_wdl(f'{wf.meta.name}.wdl')
    # wf.to_argo_worflow(f'{wf.meta.name}.yaml')
    wf.to_nestcmd(outdir='look')


if __name__ == '__main__':
    pipeline()
