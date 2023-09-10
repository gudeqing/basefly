import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar


def stepA():
    cmd = Command()
    cmd.meta.name = 'Module-A'
    cmd.meta.desc = 'This module is ude to collect information'
    cmd.runtime.image = 'docker/whalesay'
    cmd.runtime.tool = 'echo'
    cmd.args['content'] = Argument(value="'hello cow boy'")
    cmd.args['outfile'] = Argument(prefix='> ', value='outfile.txt')
    cmd.outputs['outfile'] = Output(type='outfile', value='~{outfile}')
    return cmd


def stepB():
    cmd = Command()
    cmd.meta.name = 'Module-B'
    cmd.meta.desc = 'This tool is used to get ouput2'
    cmd.meta.source = 'in-house development'
    cmd.runtime.image = 'docker/whalesay'
    cmd.runtime.tool = 'cat'
    cmd.args['infile'] = Argument(prefix='', type='infile')
    cmd.args['outfile'] = Argument(prefix='> ', value='outfile2.txt')
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

    wf.topvars = dict(
        content=TopVar(value=wf.args.input, type='str')
    )

    task_a, args = wf.add_task(stepA(), tag='firstStep')
    args['content'].value = wf.topvars['content']

    for idx in [1, 2]:
        task_b, args = wf.add_task(stepB(), name=f'B-{idx}', parent_wkdir='B-tasks', depends=[task_a])
        args['infile'].value = task_a.outputs['outfile']

    wf.run()


if __name__ == '__main__':
    pipeline()



