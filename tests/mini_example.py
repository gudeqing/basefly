import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar


# 定义如下函数时，最好不要带参数, 如果带参数，在给参数赋值时，建议使用value而不是default
def stepA():
    cmd = Command()
    cmd.meta.name = 'Module-A'
    cmd.meta.desc = 'This module is ude to collect information'
    cmd.runtime.image = 'docker/whalesay'
    cmd.runtime.tool = 'echo'
    cmd.args['content'] = Argument(type='str', default="'hello cow boy'", desc='the content of input')
    cmd.args['outfile'] = Argument(prefix='> ', type='outstr', default='outfile.txt', desc='output file name')
    cmd.outputs['outfile'] = Output(type='outfile', value='{outfile}')
    return cmd


def stepB(out_prefix='outfile2'):
    cmd = Command()
    cmd.meta.name = 'Module-B'
    cmd.meta.desc = 'This tool is used to get ouput2'
    cmd.meta.source = 'in-house development'
    cmd.runtime.image = 'docker/whalesay'
    cmd.runtime.tool = 'cat'
    cmd.args['infile'] = Argument(prefix='', type='infile', desc='input file')
    cmd.args['outfile'] = Argument(prefix='> ', type='outstr', value=f'{out_prefix}.txt', desc='output file')
    cmd.outputs['outfile'] = Output(type='outfile', value='~{outfile}')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'MiniPipeline'
    wf.meta.desc = 'This is a mini pipeline example'
    wf.meta.source = 'in-house development'
    wf.init_argparser()
    wf.add_argument('-input', help='input information for the mini workflow')
    wf.parse_args()

    wf.topvars = {
        'content': TopVar(value=wf.args.input, type='str', name='content'),
    }

    # 王流程中添加task，tag用于形成task的具体名称：name= cmd.meta.name +'-' + tag
    task_a, args = wf.add_task(stepA(), tag='firstStep')
    # 参数的赋值必须调用'.value'属性赋值，赋值对象可以是TopVar或Output或其他
    args['content'].value = wf.topvars['content']
    # 增加更多的task
    for idx in [1, 2]:
        task_b, args = wf.add_task(stepB(), name=f'B-{idx}', parent_wkdir='B-tasks', depends=[task_a])
        args['infile'].value = task_a.outputs['outfile']
    # 运行流程
    wf.run()
    # 下面可以补充其他代码执行总结或清理工作
    if wf.success:
        print('do something else')


if __name__ == '__main__':
    pipeline()



