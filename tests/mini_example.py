import os
script_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import sys; sys.path.append(script_path)
from basefly.basefly import Argument, Output, Command, Workflow, TopVar

# 参考代码的注释信息，了解如何写basefly流程


# 定义如下函数时，最好不要带参数, 如果带参数，在给参数赋值时，建议使用value而不是default
def stepA():
    # 初始化tool
    cmd = Command()
    cmd.meta.name = 'ModuleA'  # 工具名称，将被用作task名称的前缀，不能有重复，最好不用特殊字符串，不建议包含‘-’或‘_'
    cmd.meta.desc = 'This module is used to collect information'
    cmd.meta.version = 'version_number'  # 建议填写软件版本号
    cmd.meta.author = 'who_create_the_tool'

    # 定义运行环境
    cmd.runtime.image = 'docker/whalesay'  # 如果不使用docker，则可以省略
    cmd.runtime.tool_dir = ''   # 通常填写工具所在目录，可以省略
    cmd.runtime.tool = 'echo'  # 通常是命令行中所使用的tool的名称，如"gatk" 或 "bwa mem"

    # 定义参数信息，注意顺序是很重要的
    # 参数描述的越详细越好，有利于一键生成有意义的流程说明文档
    cmd.args['content'] = Argument(type='str', default="'hello cow boy'", desc='the content of input')
    cmd.args['outfile'] = Argument(prefix='> ', type='outstr', default='outfile.txt', desc='output file name')

    # 定义输出，可以用’{}‘或’~{}‘引用前面定义的参数，后续解析命令行时将自动从cmd.args抓取信息进行实际信息替换
    # 你不需要定义工作目录作为输出，默认会给每一个task定义一个名为“_wkdir_”的输出对象，每个task也会默认有wkdir属性
    cmd.outputs['outfile'] = Output(type='outfile', value='{outfile}')

    # 记得返回cmd
    return cmd


def stepB(out_prefix='outfile2'):
    # 注意：这个函数，带入了参数，认为是不规范的行为，但是被允许，如果你不打算将流程转化为其他格式流程，不用考虑潜在影响
    cmd = Command()
    cmd.meta.name = 'ModuleB'
    cmd.meta.desc = 'This tool is used to get ouput2'
    cmd.meta.source = 'in-house development'
    cmd.runtime.image = 'docker/whalesay'
    cmd.runtime.tool = 'cat'
    cmd.args['infile'] = Argument(prefix='', type='infile', desc='input file')
    # 注意：默认前缀为空，前缀可以包含一个“{}”引用一次参数的值或多个“{v}"多次引用参数值，如你可以定以 prefix = '@{v}:@{v}' 以处理参数前缀和参数值交融的情景
    # 注意：前缀末尾可以包含空格符，因为空格也算是前缀的一部分
    cmd.args['outfile'] = Argument(prefix='> ', type='outstr', value=f'{out_prefix}.txt', desc='output file')
    cmd.outputs['outfile'] = Output(type='outfile', value='~{outfile}')
    return cmd


def pipeline():
    wf = Workflow()
    wf.meta.name = 'MiniPipeline'
    wf.meta.desc = 'This is a mini pipeline example'
    wf.meta.source = 'in-house development'
    wf.meta.version = '1.1'
    wf.init_argparser()  # 初始化流程参数，使得流程具备add_argument功能
    # 添加当前流程特有的输入参数
    # 这里本质调用的是python的”argparse“包，因此可以参见argparse的说明文档
    wf.add_argument('-input', help='input information for the mini workflow')
    # 收集参数
    wf.parse_args()

    # 进一步整理流程最开始的输入信息,放到wf.topvars字典, 这样方便流程的转化，也方便后续输入参数检查和引用，强烈建议如此处理输入
    wf.topvars = {
        'content': TopVar(value=wf.args.input, type='str', name='content'),
    }

    # 往流程中添加第一个task，如下，tag用于形成task的具体名称：name= cmd.meta.name +'-' + tag
    task_a, args = wf.add_task(stepA(), tag='firstStep')
    # 参数的赋值必须调用'.value'属性赋值，赋值对象可以是TopVar或Output或其他
    args['content'].value = wf.topvars['content']

    # 增加更多的task
    for idx in [1, 2]:
        # 默认depends是空的列表，你可以显示添加depends信息，如果不输入，basefly将自动根据输入信息添加depends信息
        # parent_wkdir参数可以帮你把类似的任务归类到同一个目录进行执行，让结果目录看起来更简洁
        task_b, args = wf.add_task(stepB(), name=f'B-{idx}', parent_wkdir='B-tasks', depends=[task_a])

        # 参数赋值时，一定要记得用‘.value'属性，否则会出错
        args['infile'].value = task_a.outputs['outfile']

    # 运行流程, 必须添加下面这一行，否则流程无效
    wf.run()

    # 下面可以补充其他代码执行整理工作
    if wf.success:
        print('do something else')


if __name__ == '__main__':
    pipeline()
