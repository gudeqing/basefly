import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar, ToWdlTask
from utils.get_fastq_info import get_fastq_info
__author__ = 'gdq'


def vcf2bed():
    cmd = Command()
    cmd.meta.name = 'vcf2bed'
    cmd.runtime.image = ''
    cmd.runtime.tool = ''
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', level='optional', desc='read2 fastq file')
    cmd.args['threads'] = Argument(prefix='-w ', default=7, desc='thread number')
    cmd.args['other_args'] = Argument(prefix='', default='', desc="other arguments you want to use, such as '-x val'")
    cmd.args['out1'] = Argument(prefix='-o ', type='str', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', level='optional', type='str', desc='clean read2 output fastq file')
    cmd.args['html'] = Argument(prefix='-h ', type='str', desc='html report file')
    cmd.args['json'] = Argument(prefix='-j ', type='str', desc='json report file')
    # 定义输出
    cmd.outputs['out1'] = Output(value="{out1}")  # 这里使用”{}“引用其他Argument对象作为输入
    cmd.outputs['out2'] = Output(value="{out2}")
    cmd.outputs['html'] = Output(value="{html}")
    cmd.outputs['json'] = Output(value="{json}")
    return cmd
