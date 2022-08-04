import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar


def stepA():
    cmd = Command()
    cmd.meta.name = 'A'
    cmd.runtime.image = 'docker/whalesay'
    cmd.runtime.tool = 'echo'
    cmd.args['content'] = Argument(value="'hello cow boy'")
    cmd.args['outfile'] = Argument(prefix='> ', value='outfile.txt')
    cmd.outputs['outfile'] = Output(type='outfile', value='~{outfile}')
    return cmd


def stepB():
    cmd = Command()
    cmd.meta.name = 'B'
    cmd.runtime.image = 'docker/whalesay'
    cmd.runtime.tool = 'cat'
    cmd.args['infile'] = Argument(prefix='', type='infile')
    cmd.args['outfile'] = Argument(prefix='> ', value='outfile2.txt')
    cmd.outputs['outfile'] = Output(type='outfile', value='~{outfile}')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'mini_pipeline'
    wf.meta.desc = 'This is a simple pipeline'
    wf.init_argparser()
    wf.add_argument('-input', help='input information')
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



