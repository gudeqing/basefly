import sys; sys.path.append('..')
from nestcmd.nestcmd import Argument, Output, Command, Workflow, TopVar


def stepA():
    cmd = Command()
    cmd.meta.name = 'A'
    cmd.runtime.image = 'docker/whalesay'
    cmd.runtime.tool = 'cowsay'
    cmd.args['content'] = Argument(value="'hello cow boy'")
    cmd.args['outfile'] = Argument(prefix='> ', value='outfile.txt')
    cmd.outputs['outfile'] = Output(type='File', value='~{outfile}')
    return cmd


def stepB():
    cmd = Command()
    cmd.meta.name = 'B'
    cmd.runtime.image = 'docker/whalesay'
    cmd.runtime.tool = 'cat'
    cmd.args['infile'] = Argument(prefix='', type='infile')
    cmd.args['outfile'] = Argument(prefix='> ', value='outfile2.txt')
    cmd.outputs['outfile'] = Output(type='File', value='~{outfile}')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'mini_pipeline'
    wf.meta.desc = 'This is a simple pipeline'

    task_a, args = wf.add_task(stepA(), name='A')
    task_b, args = wf.add_task(stepB(), name='B', depends=[task_a.task_id])
    args['infile'].value = task_a.outputs['outfile']

    wf.to_argo_worflow('mini_wf.yaml')
    wf.to_nestcmd(outdir='look', run=True)


if __name__ == '__main__':
    pipeline()



