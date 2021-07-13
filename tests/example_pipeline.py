from nestcmd.nestcmd import Argument, Output, Command, Workflow
"""
pycharm里设置code completion 允许“suggest variable and parameter name”, 可以极大方便流程编写
要想生成正确的wdl:
1. 一定要正确定义参数的类型，比如是否为”infile“ or "indir",
 type is one of ['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix']
 其中‘fix'可以用于添加管道符等固定字符串
2. 参数的添加顺序对于命令行的正确形成很重要，这里，字典是有序性得到利用
"""


def fastp():
    cmd = Command()
    cmd.meta.name = 'fastp'
    cmd.runtime.image = 'gudeqing/fastp:0.21.0'
    cmd.runtime.tool = 'fastp'
    # 可以直接用访问属性的方式添加参数，这个得益于使用Munch对象而不是原生字典
    cmd.args.read1 = Argument(prefix='-i ', type='infile', desc='read1 fastq file')
    cmd.args.read2 = Argument(prefix='-I ', type='infile', desc='read2 fastq file')
    # 当然，可以直接用字典的方式添加参数
    cmd.args['out1'] = Argument(prefix='-o ', type='str', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', type='str', desc='clean read2 output fastq file')
    # 下面的outputs设置纯粹是为了能够生成wdl设置
    cmd.outputs['out1'] = Output(path="~{out1}")
    cmd.outputs['out2'] = Output(path="~{out2}")
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
    cmd.args.libType = Argument(prefix='--libType ', default='A')
    cmd.args['indexDir'] = Argument(prefix='-i ', type='indir', desc='transcript fasta index directory')
    cmd.args['read1'] = Argument(prefix='-1 ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-2 ', type='infile', desc='read2 fastq file')
    cmd.args['outDir'] = Argument(prefix='-o ', type='str', default='quant', desc='output directory')
    cmd.args['gcBias'] = Argument(prefix='--gcBias ', type='bool', default=True, desc='perform gc Bias correction')
    cmd.outputs['transcript'] = Output(path="~{outDir}" + "/quant.sf", locate='quant')
    cmd.outputs['outDir'] = Output(path="~{outDir}", type='Directory')
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
    cmd.outputs['result'] = Output(path="~{out}", locate='quant')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'PipelineExample'
    wf.meta.desc = 'This is a simple pipeline for fast gene/transcript quantification. '
    wf.meta.desc += 'workflow = [fastq -> Fastp -> Salmon]'
    index_dir = 'testdata/index/'

    # init
    def init_func():
        samples = ['s1', 's2']
        read1s = ['testdata/s1.R1.fastq', 'testdata/s2.R1.fastq']
        read2s = ['testdata/s1.R2.fastq', 'testdata/s2.R2.fastq']
        return zip(samples, read1s, read2s)

    merge_depends = []
    for sample, r1, r2 in init_func():
        # task_name = 'quant_' + sample
        task, args = wf.add_task(fastp())
        # 带入样本信息
        task.sample = sample
        # 给task分组信息，用于wdl转换时的判读
        task.group = 'batch1'
        task_id = task.task_id
        args['read1'].value = r1
        # 给参数添加wdl属性,目的是为了正确生成wdl流程，如果无需生成wdl流程，则无需添加wdl属性。
        # 下面代码中的each表示wdl中对array进行循环时的临时变量, each结构：{sample:{r1:r2}}
        args['read1'].wdl = "each.right.left"
        args['read2'].value = r2
        args['read2'].wdl = "each.right.right"
        args['out1'].value = f'{sample}.clean.R1.fq'
        args['out1'].wdl = "~{each.left}.clean.R1.fq"
        args['out2'].value = f'{sample}.clean.R2.fq'
        args['out2'].wdl = "~{each.left}.clean.R2.fq"
        # task的outputs对象不会用于wdl的生成,所以无需添加wdl属性
        task.outputs['out1'] = Output(path=f'{sample}.clean.R1.fq')
        task.outputs['out2'] = Output(path=f'{sample}.clean.R2.fq')
        # print(task.cmd.format_cmd())

        depend_task = task
        task, args = wf.add_task(salmon())
        task.depends = [task_id]
        task.sample = sample
        task.group = 'batch1'
        args['read1'].value = depend_task.outputs["out1"].path
        args['read1'].wdl = f"{depend_task.cmd.meta.name}.out1"
        args['read2'].value = depend_task.outputs["out2"].path
        args['read2'].wdl = f"{depend_task.cmd.meta.name}.out2"
        args['indexDir'].value = index_dir
        args['outDir'].value = sample
        args['outDir'].wdl = "each.left"
        task.outputs['outDir'] = Output(path=sample, type='Directory')
        task.outputs['transcript'] = Output(path=sample + '/' + 'quant.sf')
        task.outputs['gene'] = Output(path=sample + '/' + 'quant_gene.sf')
        merge_depends.append(task.task_id)

    # merge transcript TPM
    task, args = wf.add_task(quant_merge())
    task.cmd.meta.name = 'MergeTranscriptTPM'
    task.depends = merge_depends
    # 针对当前任务，column 和 genes 参数需要固定，不允许修改,用 'fix' 表示，便于wdl的正确生成。
    args['column'].value = 'TPM'
    args['column'].type = 'fix'
    args['genes'].value = False
    args['genes'].type = 'fix'
    args['quants'].value = [wf.tasks[task_id].outputs['outDir'].path for task_id in task.depends]
    args['quants'].wdl = f'{wf.tasks[merge_depends[0]].cmd.meta.name}.outDir'
    task.outputs['result'] = Output(path=f'merged.{args["column"].value}.txt', type='File')

    # merge transcript Count
    task, args = wf.add_task(quant_merge())
    task.cmd.meta.name = 'MergeTranscriptCount'
    task.depends = merge_depends
    # 针对当前任务，column 和 genes 参数需要固定，不允许修改
    args['column'].value = 'NumReads'
    args['column'].type = 'fix'
    args['genes'].value = False
    args['genes'].type = 'fix'
    args['quants'].value = [wf.tasks[task_id].outputs['outDir'].path for task_id in task.depends]
    args['quants'].wdl = f'{wf.tasks[merge_depends[0]].cmd.meta.name}.outDir'
    args['out'].default = f'merged.{args["column"].value}.txt'
    task.outputs['result'] = Output(path=f'merged.{args["column"].value}.txt', type='File')

    for _, task in wf.tasks.items():
        # print(task.task_id)
        print(task.outputs)
        print(task.cmd.format_cmd())

    wf.to_wdl(f'{wf.meta.name}.wdl')


if __name__ == '__main__':
    pipeline()
