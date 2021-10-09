import os
import re
import json
from uuid import uuid4
from dataclasses import dataclass, field
# from typing import Any, List, Dict, Literal
from typing import Any, List, Dict
from typing_extensions import Literal
# shadow default dict
# from munch import Munch as dict

__author__ = 'gdq'

"""
设计思路
1. 定义argument,runtime,outputs,meta
    argument: 全面描述参数，最终可以根据参数列表自动组装成命令行
    runtime：运行环境描述，如软件，docker image，memory，cpu限制等信息
    outputs: 描述程序输出
    meta: 描述程序/command的信息，如名称，作者，版本等信息
2. 由上述对象进一步定义Command对象, 即Command = {Argument, RunTime, Output, Meta}
3. 定义Task对象：
    1. 给Command的参数具体赋值
    2. 加上depend信息
    3. 继承Command的outputs
4. 定义workflow对象： 
    1. 由task对象可构成workflow对象，workflow = dict(task1=task_body, task2=task_body, ...)
    2. 当然也可以给workflow对象添加outputs对象, task的outputs继承自Command的outputs
    
5. 定义方法: 依据Task对象生成具体的cmd信息。
6. 定义方法: 将Command/Task对象转换成wdl脚本
7. 定义方法将workflow转换成wdl流程等、argo_workflow、nestcmd（我自定义的流程）

注意：
1. python3.6以后的 dict有序性对于cmd的正确解析非常重要，因为定义Argument的顺序非常重要，很多命令行要求参数有序。
3. Argument支持“多值且可重复”参数，如下，但是没有办法转化成wdl的，因为wdl不支持这么复杂的参数
    Argument(prefix='--x ', array=True, multi_times=True, default=[[1, 2], ['c', 'y']], delimiter=',')
"""


@dataclass()
class Argument:
    name: str = '?'
    value: Any = None
    # prefix 可以是如 ’-i '或 'i=', 对于前者, 空格一定要带上
    prefix: str = ''
    # type is one of ['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix']
    # fix 类型表示该参数并不是真正的参数，其为固定的字符串. 例如其可以用来表示管道符如‘| samtools sort’
    type: Literal['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix'] = 'str'
    level: Literal['required', 'optional'] = 'required'
    # for bool type, default is one of ['false', true']
    default: Any = None
    range: Any = None
    array: bool = False
    # delimiter：当一个参数可以输入多个值，这几个值的分隔符由此指定
    delimiter: str = ' '
    # 指示一个参数是否可以多次使用
    multi_times: bool = False
    # format 是描述输入文件的格式用的
    format: str = None
    # order字段是为了参数排序设计
    order: int = 0
    desc: str = 'This is description of the argument.'

    def __post_init__(self):
        # type 类型检查和矫正
        if type(self.default) == int and self.type == 'str':
            # print(f'{self.name}: type of default value is not agree with type specified!')
            self.type = 'int'
        elif type(self.default) == float and self.type == 'str':
            # print(f'{self.name}: type of default value is not agree with type specified!')
            self.type = 'float'
        elif type(self.default) == bool:
            self.type = 'bool'

        if self.type == 'bool':
            # 对于布尔参数，其一定为必要参数类型,可选范围为True或False
            self.level = 'required'
            self.range = {True, False}
            if not self.default:
                # 对于没有默认值的bool参数，强行赋值为false，即该参数默认不参与命令行的形成
                self.default = False
        elif self.type == 'fix':
            if not self.value:
                self.value = self.default
            else:
                self.default = self.value
        # 如果某个参数可以接收多个值，那么default，value的类型强制设为列表
        if self.array:
            if not self.multi_times:
                self.__annotations__['default'] = List[Any]
                self.__annotations__['value'] = List[Any]
            else:
                self.__annotations__['default'] = List[List[Any]]
                self.__annotations__['value'] = List[List[Any]]
        # 对于有集合候选的参数，先对默认值检查一番
        if type(self.range) == set:
            if self.default not in self.range:
                raise Exception(f'default value is not in {self.range}')


@dataclass()
class RunTime:
    image: str = None
    # 软件所在目录
    tool_dir: str = ''
    # 软件名称
    tool: str = ''
    # 执行环境设置, 以下为命令运行所需的最少计算资源
    memory: int = 1000
    cpu: int = 2
    # 执行环境设置，以为命令行所需的最大计算资源
    max_memory: int = 0
    max_cpu: int = 0
    # 运行时间上线
    timeout: int = 3600*24

@dataclass()
class Output:
    path: str = None  # will be deprecated, replace with value
    value: Any = None
    # out_id: str = field(default_factory=uuid4)
    type: Literal['str', 'int', 'float', 'bool', 'outfile', 'outdir'] = 'outfile'
    # 设计report参数，用于指定该参数最终是否作为流程的outputs
    report: bool = True
    # 由command形成task之后，可以带入task_id
    task_id: str = None
    name: str = None

    def __post_init__(self):
        if self.value is None:
            # value属性是后来修改增加的，本来用来替代path，为了向前兼容，保留path属性
            self.value = self.path


@dataclass()
class Meta:
    name: str = None
    desc: str = 'This is description of the tool/workflow.'
    author: str = 'unknown'
    source: str = 'source URL for the tool'
    version: str = 'unknown'


@dataclass()
class Command:
    meta: Meta = field(default_factory=Meta)
    runtime: RunTime = field(default_factory=RunTime)
    args: Dict[str, Argument] = field(default_factory=dict)
    # 下面支持用wdl的定义’~‘符号, 当前脚本要求所有命令行的输出皆为文件
    outputs: Dict[str, Output] = field(default_factory=dict)

    def format_cmd(self, wf_tasks=None):
        # 该函数应该在参数被具体赋值后才能调用
        if not self.args:
            raise Exception(f'Command {self.meta.name} has no args !')
        cmd = ''
        if self.runtime.tool_dir and not self.runtime.tool_dir.endswith('/'):
            self.runtime.tool_dir += '/'
        cmd += self.runtime.tool_dir + self.runtime.tool

        for arg_name, arg in self.args.items():
            # 对于极端情况:如果不小心定义了"有默认值的非必须参数"，仅当明确赋值时才参与命令行的形成
            if arg.level == 'required':
                arg_value = arg.value or arg.default
            else:
                arg_value = arg.value

            if arg_value is None:
                if arg.level == 'required':
                    raise Exception(f'No value found for {arg_name}!')
                else:
                    # 对于非必须参数，且没有赋值的参数，直接跳过，不参与命令行的形成
                    continue

            # 当参数值为output or topVar类型时，需如下特殊处理得到实际的参数值
            if type(arg_value) == Output:
                value_dict = dict()
                for k, v in wf_tasks[arg_value.task_id].cmd.args.items():
                    if type(v.value) != TmpVar:
                        value_dict[k] = v.value or v.default
                    else:
                        value_dict[k] = v.value.value
                arg_value = arg_value.value.replace('~', '').format(**value_dict)
            elif type(arg_value) == TopVar or type(arg_value) == TmpVar:
                arg_value = arg_value.value
            elif type(arg_value) == list:
                arg_value = arg_value.copy()
                for ind, each in enumerate(arg_value):
                    if type(each) == Output:
                        value_dict = dict()
                        for k, v in wf_tasks[each.task_id].cmd.args.items():
                            if type(v.value) != TmpVar:
                                value_dict[k] = v.value or v.default
                            else:
                                value_dict[k] = v.value.value
                        # 当前设计中可以用’~{target}‘引用其他参数的形成outputs的值,这里需要替换回来
                        # '~{}'也是wdl的语法, 其实起初这个功能是为wdl设置的
                        arg_value[ind] = each.value.replace('~', '').format(**value_dict)
                    elif type(each) == TopVar or type(each) == TmpVar:
                        arg_value[ind] = each.value

            # 对于可以接收多个值的参数
            if arg.array:
                if not arg.multi_times:
                    arg_value = arg.delimiter.join([str(x) for x in arg_value])
                else:
                    arg_value = [arg.delimiter.join([str(x) for x in v]) for v in arg_value]

            # 处理bool值参数
            if arg.type == "bool" or type(arg.value) == bool:
                if arg_value:
                    cmd += ' ' + arg.prefix
                else:
                    # 如果bool类型且value为false，则该参数不参与命令行形成
                    continue
            else:
                if not arg.multi_times:
                    cmd += ' ' + arg.prefix + str(arg_value)
                else:
                    cmd += ' ' + arg.prefix + (' ' + arg.prefix).join(arg_value)
        return cmd

    def format_wdl_task(self, outfile=None, wdl_version='development'):
        # 形成wdl格式的task非常复杂，因此单独提供ToWdlTask处理，该函数只有当你仅需要某个task的wdl版本时需要
        if not outfile:
            outfile = f'{self.meta.name}.wdl'
        ToWdlTask(self, outfile, wdl_version)


@dataclass()
class Task:
    cmd: Command
    name: str = None
    task_id: str = field(default_factory=uuid4)
    depends: List[str] = field(default_factory=list)

    def __post_init__(self):
        # task name
        if self.name is None:
            self.name = self.cmd.meta.name + '-' + str(self.task_id)

        # 为每一个output带入
        for key in self.cmd.outputs.keys():
            self.cmd.outputs[key].task_id = self.task_id
            self.cmd.outputs[key].name = key
        # 继承 cmd 的outputs
        self.outputs = self.cmd.outputs

    def argo_template(self, wf_tasks):
        # 仅输入文件需要作处理，其他参数如数字或字符串已经在command中硬编码好了
        # 输入文件需要处理，是因为需要申明有这样一个文件存在
        template_name = self.name or self.task_id
        lines = [f'- name: {template_name.replace("_", "-")}']
        artifacts = [k for k, v in self.cmd.args.items() if (v.type == 'infile' or v.type == 'indir')]
        if artifacts:
            lines += [' ' * 2 + 'inputs:']
            lines += [' ' * 4 + 'artifacts:']
            for each in artifacts:
                if type(self.cmd.args[each].value) != list:
                    lines += [' ' * 6 + f'- name: {each}']
                    lines += [' ' * 8 + f'path: {self.cmd.args[each].value.value}']
                else:
                    for i, v in enumerate(self.cmd.args[each].value):
                        lines += [' ' * 6 + f'- name: {each}_{i}']
                        lines += [' ' * 8 + f'path: {v.value}']

        # container info
        lines += [' '*2 + 'container:']
        lines += [' '*4 + f'image: {self.cmd.runtime.image}']
        lines += [' '*4 + 'command: [sh, -c]']
        lines += [' '*4 + f'args: ["{self.cmd.format_cmd(wf_tasks=wf_tasks)}"]']

        # outputs
        if self.cmd.outputs:
            value_dict = {k: v.value or v.default for k, v in self.cmd.args.items()}
            lines += [' '*2 + 'outputs:']
            lines += [' '*4 + 'artifacts:']
            for k, v in self.cmd.outputs.items():
                lines += [' ' * 6 + f'- name: {k}']
                lines += [' ' * 8 + f'path: {v.value.format(**value_dict)}']
        return [' '*2+x for x in lines]


@dataclass()
class TopVar:
    """
    该对象用于描述输入流程的起始输入对象, 这样设计的目的是为了方便流程的转化
    """
    value: Any
    name: str = 'notNamed'
    type: Literal['str', 'int', 'float', 'bool', 'infile', 'indir'] = 'infile'

    def __post_init__(self):
        if self.type in ['infile', 'indir']:
            # 对输入文件的路径进行绝对化
            self.value = os.path.abspath(self.value)


@dataclass()
class TmpVar:
    """
    该对象用于描述流程中如循环时的临时变量, 纯粹是为wdl的scatter模式设计, name属性将作为wdl文件中的传递值。
    例如：如循环中某个变量为A，则把输出文件名定义为out_file = ’~{A}.txt‘.
    那么，这时脚本中可以写: outfile = TmpVar(name="~{A}.txt", type='str') , 这里name中“~{}”是wdl传递变量值的语法。
    """
    value: Any
    name: str = 'notNamed'
    type: Literal['str', 'int', 'float', 'bool', 'infile', 'indir'] = 'str'

    def __post_init__(self):
        if self.type in ['infile', 'indir']:
            # 对输入文件的路径进行绝对化
            self.value = os.path.abspath(self.value)
        if '~{' in self.name:
            self.name = f'"{self.name}"'


@dataclass()
class Workflow:
    meta: Meta = field(default_factory=Meta)
    tasks: Dict[str, Task] = field(default_factory=dict)
    outputs: Dict[str, Output] = field(default_factory=dict)
    top_vars: Dict[str, TopVar] = field(default_factory=dict)

    def __post_init__(self):
        for k, v in self.top_vars.items():
            # 将key作为var的名称
            v.name = k

    def add_task(self, cmd: Command, depends: list = None, name: str = None):
        task = Task(cmd=cmd, depends=depends, name=name)
        self.tasks[task.task_id] = task
        return task, task.cmd.args

    def to_wdl(self, outfile):
        ToWdlWorkflow(self).write_wdl(outfile)

    def to_argo_worflow(self, outfile):
        # 有待完善
        lines = ['apiVersion: argoproj.io/v1alpha1']
        lines += ['kind: Workflow']
        lines += ['metadata:']
        lines += [' ' * 2 + f'generateName: {self.meta.name.lower()}']
        lines += ['spec:']
        # entry point
        lines += [' '*2 + 'entrypoint: main']
        artifacts = [k for k, v in self.top_vars.items() if v.type in ['infile', 'indir']]
        if artifacts:
            lines += [' ' * 2 + 'arguments:']
            lines += [' ' * 4 + 'artifacts:']
            for each in artifacts:
                if type(self.top_vars[each].value) != list:
                    lines += [' ' * 4 + f'- name: {each}']
                    lines += [' ' * 6 + f'path: {self.top_vars[each].value}']
                else:
                    for i, v in enumerate(self.top_vars[each].value):
                        lines += [' ' * 4 + f'- name: {each}_{i}']
                        lines += [' ' * 6 + f'path: {v}']

        # DAG templates
        lines += ['']
        lines += [' '*2 + 'templates:']
        lines += [' '*2 + '- name: main']
        lines += [' '*4 + 'dag:']
        lines += [' '*6 + 'tasks:']
        for task_id, task in self.tasks.items():
            task.name = task.name.replace('_', '-')
            lines += [' '*6 + f'- name: {task.name}']
            if task.depends:
                lines += [' '*8 + 'dependencies: ' + str([self.tasks[x].name for x in task.depends]).replace("'", '')]
            lines += [' '*8 + f'template: {task.name}']
            args = task.cmd.args
            # 通过TopVar或者Output传递的文件参数
            artifacts = [k for k, v in args.items() if (v.type == 'infile' or v.type == 'indir')]
            if artifacts:
                lines += [' '*8 + f'arguments:']
                lines += [' '*10 + f'artifacts:']
                for each in artifacts:
                    if type(args[each].value) != list:
                        lines += [' ' * 10 + f'- name: {each}']
                        if type(args[each].value) in [TopVar, TmpVar]:
                            lines += [' ' * 12 + f'from: "{{{{workflow.artifacts.{args[each].value.name}}}}}"']
                        elif type(args[each].value) == Output:
                            depend_task = self.tasks[args[each].value.task_id].name
                            depend_name = args[each].value.name
                            lines += [' ' * 12 + f'from: "{{{{tasks.{depend_task}.outputs.artifacts.{depend_name}}}}}"']
                        else:
                            # 只能假设来自topVar了
                            lines += [' ' * 12 + f'from: "{{{{workflow.artifacts.{each}}}}}"']
                    else:
                        for i, v in enumerate(args[each].value):
                            lines += [' ' * 10 + f'- name: {each}_{i}']
                            if type(v) in [TopVar, TmpVar]:
                                lines += [' ' * 12 + f'from: "{{{{workflow.artifacts.{v.name}}}}}"']
                            elif type(v) == Output:
                                depend_task = self.tasks[v.task_id].name
                                depend_name = v.name
                                lines += [' ' * 12 + f'from: "{{{{tasks.{depend_task}.outputs.artifacts.{depend_name}}}}}"']
                            else:
                                # 只能假设来自topVar了
                                lines += [' ' * 12 + f'from: "{{{{workflow.artifacts.{each}_{i}}}}}"']

        for task_id, task in self.tasks.items():
            lines += ['']
            lines += task.argo_template(self.tasks)

        # write argo workflow
        with open(outfile, 'w') as f:
            for line in lines:
                f.write(line+'\n')

    def to_nestcmd(self, outdir, run=False, no_docker=False, threads=3, retry=1, time_wait_resource=15,
                   no_monitor_resource=False, no_check_resource=False):
        """
        生成nestcmd格式的workflow并且允许直接本地执行
        """
        import configparser
        wf = configparser.ConfigParser()
        wf.optionxform = str
        outdir = os.path.abspath(outdir)
        wf['mode'] = dict(
            outdir=outdir,
            threads=threads,
            retry=retry,
            monitor_resource=not no_monitor_resource,
            monitor_time_step=2,
            check_resource_before_run=not no_check_resource,
        )

        for task_id, task in self.tasks.items():
            cmd_wkdir = os.path.join(outdir, task.name)
            mount_vols = {cmd_wkdir}
            for k, v in task.cmd.args.items():
                if type(v.value) in [list, tuple]:
                    values = v.value
                else:
                    values = [v.value]
                for value in values:
                    if type(value) == Output and value.type in ['outfile', 'outdir']:
                        if not value.value.startswith('${{mode:'):
                            value.value = os.path.join("${{mode:outdir}}", self.tasks[value.task_id].name, value.value)
                        mount_vols.add(os.path.join(outdir, self.tasks[value.task_id].name))
                    elif (type(value) == TopVar or type(value) == TmpVar) and value.type in ['infile', 'indir']:
                        if value.type == 'infile':
                            file_dir = os.path.dirname(value.value)
                            mount_vols.add(os.path.abspath(file_dir))
                        elif value.type == 'indir':
                            mount_vols.add(os.path.abspath(value.value))
                    elif v.type in ['infile', 'indir'] and type(value) == str:
                        # 直接从参数的类型来判定输入文件，这样，即使没有定义TopVar或者TmpVar等对象，也可以顺利生成nestcmd
                        if v.type == 'infile':
                            mount_vols.add(os.path.abspath(os.path.dirname(value)))
                        elif v.type == 'indir':
                            mount_vols.add(os.path.abspath(value))

            wf[task.name] = dict(
                depend=','.join(self.tasks[x].name for x in task.depends) if task.depends else '',
                cmd=task.cmd.format_cmd(self.tasks),
                mem=task.cmd.runtime.memory,
                cpu=task.cmd.runtime.cpu,
                max_mem=task.cmd.runtime.max_memory,
                max_cpu=task.cmd.runtime.max_cpu,
                timeout=task.cmd.runtime.timeout,
                image='' if no_docker else (task.cmd.runtime.image or ''),
                wkdir=cmd_wkdir,
                mount_vols=';'.join(mount_vols)
            )
        os.makedirs(outdir, exist_ok=True)
        outfile = os.path.join(outdir, f'{self.meta.name}.ini')
        with open(outfile, 'w') as configfile:
            wf.write(configfile)
        if run:
            from .runner import RunCommands
            if os.path.exists(os.path.join(outdir, 'cmd_state.txt')):
                RunCommands(outfile, timeout=time_wait_resource).continue_run()
            else:
                RunCommands(outfile, timeout=time_wait_resource).parallel_run()


class ToWdlTask(object):
    type_conv_dict = {
        'str': 'String',
        'int': 'Int',
        'float': 'Float',
        'bool': 'Boolean',
        'infile': 'File',
        'outfile': 'File',
        'indir': 'Directory',
        'outdir': 'Directory',
    }

    def __init__(self, command: Command, outfile=None, wdl_version='development'):
        self.cmd = command
        self.task_name = None
        self.wdl_version = wdl_version
        self.wdl = self.format_wdl_task(outfile, wdl_version)

    def get_task_meta(self):
        return self.cmd.meta.__dict__.copy()

    def get_runtime(self):
        runtime = self.cmd.runtime.__dict__.copy()
        if self.cmd.runtime.memory is None:
            runtime.pop('memory')
        if self.cmd.runtime.cpu is None:
            runtime.pop('cpu')
        runtime.pop('tool')
        runtime.pop('tool_dir')
        if self.wdl_version in ["development", "1.0"]:
            runtime['docker'] = runtime.pop('image')
        return runtime

    def get_parameter_meta(self):
        arg_meta = dict()
        for arg_name, detail in self.cmd.args.items():
            if detail.type == 'fix':
                continue
            detail = detail.__dict__.copy()
            keys = ['prefix', 'type', 'level', 'default', 'range', 'array', 'desc']
            arg_meta[arg_name] = {k: v for k, v in detail.items() if k in keys}
        return arg_meta

    def get_inputs_and_cmd(self):
        inputs = []
        cmd = []
        for arg_name, detail in self.cmd.args.items():
            arg_info = ''
            detail = detail.__dict__
            # define type of arg
            if detail['type'] == 'fix':
                if detail['value'] is False:
                    pass
                elif detail['value'] is True:
                    cmd += [detail['prefix']]
                else:
                    cmd += [detail['prefix'] + (detail['value'] or detail['default'])]
                continue
            elif detail['type'] == 'bool':
                arg_info = 'Boolean'
            elif detail['type'] == 'int':
                arg_info = 'Int'
            elif detail['type'] == 'float':
                arg_info = 'Float'
            elif detail['type'] == 'str':
                arg_info = 'String'
            elif detail['type'] == 'infile':
                arg_info = 'File'
            elif detail['type'] == 'indir':
                arg_info = 'Directory'
            else:
                raise Exception(f'unexpected type {arg_info}')

            if detail['array'] or detail['multi_times']:
                arg_info = f'Array[{arg_info}]'

            if detail['level'] == 'optional' and detail['type'] != 'bool':
                arg_info += '?'

            # add arg name
            arg_info += ' ' + arg_name

            # get default value
            if detail['type'] == "bool":
                if detail['default']:
                    arg_info += ' = true'
                else:
                    arg_info += ' = false'
            elif detail['type'] == "str":
                if detail['default'] is not None:
                    arg_info += ' = "{}"'.format(detail['default'])
            else:
                if detail['default'] is not None:
                    if detail['type'] != 'indir':
                        arg_info += ' = {}'.format(detail['default'])
                    else:
                        # 目前cromwell并不支持给Directory参数赋默认值
                        pass

            inputs.append(arg_info)

            # format cmd
            if detail['multi_times']:
                if detail['array']:
                    print('多值且可重复使用参数，wdl如何表示我也不知道！')
                    cmd += ['多值且可重复使用参数，wdl如何表示我也不知道！']
                else:
                    if not detail['prefix']:
                        raise Exception(f'可重复使用参数必须要有前缀:{arg_info}')
                    if detail['type'] == 'bool':
                        raise Exception(f'可重复使用参数对应的值不能是bool类型:{arg_info}')
                    cmd += ['~{sep=" " ' + 'prefix(' + f'"{detail["prefix"]}", ' + arg_name + ')}']
            else:
                if not detail['array']:
                    if detail['prefix'] == '':
                        cmd += ['~{' + arg_name + '}']
                    else:
                        if detail['type'] != 'bool':
                            cmd += ['~{' + f'"{detail["prefix"]}"' + ' + ' + arg_name + '}']
                        else:
                            cmd += ['~{' + f'if {arg_name} then "{detail["prefix"]} " else ""' + '}']
                else:
                    delimiter = detail['delimiter']
                    if detail['prefix'] == '':
                        cmd += ['~{sep=' + f'"{delimiter}" ' + arg_name + '}']
                    else:
                        if detail['type'] != 'bool':
                            prefix = '~{' + f'if defined({arg_name}) then "{detail["prefix"]} " else ""' + '}'
                            cmd += [prefix + '~{sep=' + f'"{delimiter}" ' + arg_name + '}']
                        else:
                            raise Exception('多值参数对应的值不能是bool类型！')
        return inputs, cmd

    def get_outputs(self):
        outputs = []
        for name, v in self.cmd.outputs.items():
            if '~' not in v.path:
                value = v.path.replace('{', '~{')
            else:
                value = v.path
            outputs += [self.type_conv_dict[v.type] + ' ' + name + ' = ' + f'"{value}"']
        return outputs

    def format_wdl_task(self, outfile=None, wdl_version='development'):
        """
        version = 1.0

        task cmd_name {
            input {
                xxx
            }
            command <<<
                xxx
            .>>>
            output {
                xxx
            }
            runtime {
                xxx
            }
            meta {
                xxx
            }
            parameter_meta {
                xxx
            }
        }
        :return: *.wdl file
        """
        data = self.cmd
        self.task_name = data.meta.name
        lines = f'task {data.meta.name}' + '{\n'
        # input
        lines += ' ' * 4 + 'input {\n'
        inputs, cmd = self.get_inputs_and_cmd()
        for line in inputs:
            lines += ' ' * 8 + line + '\n'
        lines += ' ' * 8 + '# for runtime\n'
        for k, v in self.get_runtime().items():
            if k in ['cpu', 'time_minutes']:
                typ = 'Int'
            else:
                typ = 'String'
                v = f'"{v}"'
            lines += ' ' * 8 + f'{typ} {k} = {v}' + '\n'
        lines += ' ' * 4 + '}\n\n'

        # command
        lines += ' ' * 4 + 'command <<<\n'
        lines += ' ' * 8 + 'set -e \n'
        base_cmd = ''
        tool_dir = data.runtime.tool_dir
        if tool_dir and not tool_dir.endswith('/'):
            tool_dir += '/'
        base_cmd += tool_dir + data.runtime.tool
        lines += ' ' * 8 + base_cmd + ' \\\n'
        for line in cmd:
            lines += ' ' * 8 + line + ' \\\n'
        lines = lines[:-2] + '\n'
        lines += ' ' * 4 + '>>>\n\n'

        # output
        lines += ' ' * 4 + 'output {\n'
        for line in self.get_outputs():
            lines += ' ' * 8 + line + '\n'
        lines += ' ' * 4 + '}\n\n'

        # runtime
        lines += ' ' * 4 + 'runtime {\n'
        for k, v in self.get_runtime().items():
            # lines += ' ' * 8 + f'{k}: "{v}"' + '\n'
            lines += ' ' * 8 + f'{k}: {k}' + '\n'
        lines += ' ' * 4 + '}\n\n'

        # task meta
        lines += ' ' * 4 + 'meta {\n'
        for k, v in self.get_task_meta().items():
            lines += ' ' * 8 + f'{k}: "{v}"' + '\n'
        lines += ' ' * 4 + '}\n\n'

        # parameter_meta
        lines += ' ' * 4 + 'parameter_meta {\n'
        for k, v in self.get_parameter_meta().items():
            v_info = '{' + ', '.join(f'{k}: "{v}"' for k, v in v.items()) + '}'
            lines += ' ' * 8 + f'{k}: {v_info}' + '\n'
        lines += ' ' * 4 + '}\n\n'

        # write wdl file
        lines += '}\n'
        if outfile:
            with open(outfile, 'w', encoding='utf-8') as f:
                f.write(f'version {wdl_version}\n\n')
                f.write(lines)
        else:
            return lines


class ToWdlWorkflow(object):
    """
    workflow的group信息的依据是：是否可以在同一个循环中并发
    如何知道task输入的依赖信息：要求给参数加一个wdl属性，对应使用wdl语法, 或者根据output对象推断，
    由于WDL语法并非十分结构化，难以正确形成，生成的流程需人工核查调整
    """
    type_conv_dict = {
        'str': 'String',
        'int': 'Int',
        'float': 'Float',
        'bool': 'Boolean',
        'infile': 'File',
        'outfile': 'File',
        'indir': 'Directory',
        'outdir': 'Directory',
    }

    def __init__(self, wf: Workflow):
        self.wf = wf

    def group_task(self):
        group = dict()
        for tid, task in self.wf.tasks.items():
            if hasattr(task, 'group'):
                group.setdefault(task.group, []).append(tid)
            else:
                group[tid] = [tid]
        return group

    def get_group_cmd_lst(self, task_ids):
        tasks = [self.wf.tasks[x] for x in task_ids]
        cmd_lst = []
        cmd_names = []
        for t in tasks:
            if t.cmd.meta.name not in cmd_names:
                cmd_names.append(t.cmd.meta.name)
                cmd_lst.append(t.cmd)
        return cmd_lst

    def format_call_cmds(self, cmds: List[Command], scatter=False):
        cmd_used_times = dict()
        lines = ''
        if scatter:
            lines += ' '*4 + 'scatter (each in keys(getFastqInfo.fastq_info)) { \n'
            lines += ' '*8 + "String sample = each\n"
            lines += ' '*8 + "File read1 = getFastqInfo.fastq_info[each][0][0]\n"
            lines += ' '*8 + "File read2 = getFastqInfo.fastq_info[each][1][0]\n"
            space_increase = 2
        else:
            lines = ''
            space_increase = 1
        for cmd in cmds:
            task_name = cmd.meta.name
            cmd_used_times.setdefault(task_name, 0)
            cmd_used_times[task_name] += 1
            if cmd_used_times[task_name] > 1:
                task_name = task_name + str(cmd_used_times[task_name])
            lines += ' '*4*space_increase + f'call {task_name} ' + '{\n'
            lines += ' '*4*(space_increase+1) + 'input: \n'
            for arg_name, detail in cmd.args.items():
                if hasattr(detail, 'wdl'):
                    if '~' in detail.wdl:
                        lines += ' '*4*(space_increase+1) + arg_name + ' = "' + detail.wdl + '"' + ',\n'
                    else:
                        lines += ' '*4*(space_increase+1) + arg_name + ' = ' + detail.wdl + ',\n'
                # 如果参数的value为output对象，则需要转换成wdl形式
                elif type(detail.value) == Output:
                    task_name = self.wf.tasks[detail.value.task_id].cmd.meta.name
                    detail.wdl = task_name + '.' + detail.value.name
                    lines += ' '*4*(space_increase+1) + arg_name + ' = ' + detail.wdl + ',\n'
                elif type(detail.value) == TopVar or type(detail.value) == TmpVar:
                    lines += ' '*4*(space_increase+1) + arg_name + ' = ' + detail.value.name + ',\n'
                elif type(detail.value) == list:
                    # print(detail.value)
                    wdl_str_set = set()
                    for each in detail.value:
                        if type(each) == Output:
                            task_name = self.wf.tasks[each.task_id].cmd.meta.name
                            wdl_str = task_name + '.' + each.name
                            wdl_str_set.add(wdl_str)
                        elif type(each) == TopVar or type(each) == TmpVar:
                            wdl_str_set.add(each.name)
                    # 如果wdl_str_set长度=1，说明这个list里面存储的是并发结果
                    if len(wdl_str_set) == 1:
                        wdl_str_lst = wdl_str_set.pop()
                    else:
                        wdl_str_lst = list(wdl_str_set)
                    lines += ' ' * 4 * (space_increase + 1) + arg_name + ' = ' + str(wdl_str_lst) + ',\n'

            lines = lines[:-2] + '\n'
            lines += ' '*4*space_increase + '}\n\n'
        if scatter:
            lines += ' '*4 + '}\n\n'
        return lines

    def write_wdl(self, outfile):
        wdl = 'version development\n\n'
        wdl += 'workflow pipeline {\n'

        # 这一部分是针对fastq数据特殊设计的
        wdl += " "*4 + "input {\n"
        if self.wf.top_vars:
            for k, v in self.wf.top_vars.items():
                var_value = v.value
                if v.type == 'infile':
                    var_type = 'File'
                    var_value = f'"{v.value}"'
                elif v.type == 'int':
                    var_type = 'Int'
                elif v.type == 'float':
                    var_type = 'Float'
                elif v.type == 'bool':
                    var_type = 'Boolean'
                    var_value = str(v.value).lower()
                elif v.type == 'indir':
                    var_type = 'Directory'
                    var_value = f'"{v.value}"'
                elif v.type == 'str':
                    var_type = 'String'
                    var_value = f'"{v.value}"'
                else:
                    var_type = 'String'
                    var_value = f'"{v.value}"'
                if var_type == 'Directory':
                    # print('warn: wdl目前不支持给Directory变量赋默认值')
                    wdl += ' ' * 4 * 2 + f'{var_type} {k}\n'
                else:
                    wdl += ' '*4*2 + f'{var_type} {k} = {var_value}\n'

        wdl += ' '*4 + '}\n\n'
        wdl += ' '*4*1 + 'call getFastqInfo{}\n'

        all_cmds = []
        scattered = []
        for grp, tids in self.group_task().items():
            cmd_lst = self.get_group_cmd_lst(tids)
            wdl += self.format_call_cmds(cmd_lst, scatter=len(tids) > 1)
            repeat = set(x.meta.name for x in cmd_lst) & set(x.meta.name for x in all_cmds)
            if repeat:
                raise Exception(f'cmd {repeat} is duplicated, '
                                f'did you forget to give a group name for cmd that called in the same loop? '
                                f'Or, you should rename cmd for different tasks !')
            all_cmds += cmd_lst
            scattered += [True]*len(cmd_lst) if len(tids) > 1 else [False]*len(cmd_lst)

        # add workflow meta
        wdl += ' '*4 + 'meta {\n'
        for k, v in self.wf.meta.__dict__.items():
            wdl += ' ' * 4*2 + f'{k}: "{v}"' + '\n'
        wdl += ' ' * 4 + '}\n\n'

        # add output section
        output_lst = []
        for cmd, is_scattered in zip(all_cmds, scattered):
            name = cmd.meta.name
            if is_scattered:
                output_lst += [f'Array[{self.type_conv_dict[v.type]}] {(name+"_"+k).replace(".", "_")} = {name}.{k}' for k, v in cmd.outputs.items()]
            else:
                output_lst += [f'{self.type_conv_dict[v.type]} {(name+"_"+k).replace(".", "_")} = {name}.{k}' for k, v in cmd.outputs.items()]
        wdl += ' ' * 4 + 'output{\n'
        for line in output_lst:
            wdl += ' ' * 4*2 + line + '\n'
        wdl += ' ' * 4 + '}\n\n'

        # end of workflow
        wdl += '}\n\n'

        # add get_fastq_info
        wdl += '\n' + get_fastq_info_wdl_str() + '\n'

        # format_tasks
        for cmd in all_cmds:
            wdl += ToWdlTask(cmd).wdl
            wdl += '\n'

        # write wdl
        with open(outfile, 'w') as f:
            f.write(wdl)


def get_fastq_info_wdl_str():
    return """
task getFastqInfo{
    input {
        Array[Directory]? fastq_dirs
        Array[File]? fastq_files
        String r1_name = '(.*).read1.fastq.gz'
        String r2_name = '(.*).read2.fastq.gz'
        String docker = 'gudeqing/getfastqinfo:1.0'
    }

    command <<<
        set -e
        python /get_fastq_info.py \\
            ~{if defined(fastq_dirs) then "-fastq_dirs " else ""}~{sep=" " fastq_dirs} \\
            ~{if defined(fastq_files) then "-fastq_files " else ""}~{sep=" " fastq_files} \\
            -r1_name '~{r1_name}' \\
            -r2_name '~{r2_name}' \\
            -out fastq.info.json
    >>>

    output {
        Map[String, Array[Array[File]]] fastq_info = read_json("fastq.info.json")
        File fastq_info_json = "fastq.info.json"
    }

    runtime {
        docker: docker
    }

    parameter_meta {
        fastq_dirs: {desc: "directory list, target fastq files should be in these directories. All target files in 'fastq_files' or 'fastq_dirs' will be used", level: "optional", type: "indir", range: "", default: ""}
        fastq_files: {desc: "target fastq file list. 'fastq_files' or 'fastq_dirs' must be provided.", level: "optional", type: "infile", range: "", default: ""}
        r1_name: {desc: "python regExp that describes the full name of read1 fastq file name. It requires at least one pair of small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'", level: "required", type: "str", range: "", default: ""}
        r2_name: {desc: "python regExp that describes the full name of read2 fastq file name. It requires at least one pair of small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'", level: "required", type: "str", range: "", default: ""}
    }
}
    """


# ----tools------
def get_fastq_info(fastq_files: tuple = None, fastq_dirs: tuple = None, out='fastq.info.json',
                   r1_name: str = "(.*).R1.fq.gz", r2_name: str = "(.*).R2.fq.gz",
                   link_data=False, add_s_to_numeric_name=False, middle2underscore=False):
    """
    :param fastq_files: target fastq file list. 'fastq_files' or 'fastq_dirs' must be provided.
    :param fastq_dirs: directory list, target fastq files should be in these directories. All target files in 'fastq_files' or 'fastq_dirs' will be used.
    :param r1_name: python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'
    :param r2_name: python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'
    :param link_data: bool to indicate if to make soft links for fastq files
    :param out: output file that contains three columns: [sample_name, read1_abs_path, read2_abs_path]
    :param add_s_to_numeric_name: bool value to indicate if to add a 'S' letter at the head of the sample name that startswith numeric string.
    :param middle2underscore: bool value to indicate if to transform '-' letter to '_' letter for a sample name.
    :return: result_dict： {sample: [[r1, r1'], [r2, r2']], ...}
    """
    if not (fastq_dirs or fastq_files):
        raise Exception("At least one of 'fastq_files' or 'fastq_dirs' must be provided.")

    result_dict = dict()

    if fastq_files:
        for each in fastq_files:
            name = os.path.basename(each)
            directory = os.path.dirname(each)
            is_read1 = True
            match = re.fullmatch(r1_name, name)
            if not match:
                match = re.fullmatch(r2_name, name)
                is_read1 = False
            if match:
                # first matched group is sample name
                sample = match.groups()[0]
                result_dict.setdefault(sample, [[], []])
                if is_read1:
                    if each not in result_dict[sample][0]:
                        result_dict[sample][0].append(each)
                    else:
                        print(f'warn: duplicated path found for {each}, and we will only keep the first one!')
                else:
                    if each not in result_dict[sample][1]:
                        result_dict[sample][1].append(each)
                    else:
                        print(f'warn: duplicated path found for {each}, and we will only keep the first one!')

    if fastq_dirs:
        for path in fastq_dirs:
            path = os.path.abspath(path)
            for root, dirs, files in os.walk(path):
                for each in files:
                    is_read1 = True
                    match = re.fullmatch(r1_name, each)
                    if not match:
                        match = re.fullmatch(r2_name, each)
                        is_read1 = False
                    if match:
                        # first matched group is sample name
                        sample = match.groups()[0]
                        result_dict.setdefault(sample, [[], []])
                        file_path = os.path.join(root, each)
                        if is_read1:
                            if file_path not in result_dict[sample][0]:
                                result_dict[sample][0].append(file_path)
                            else:
                                print(f'warn: duplicated path found for {file_path}, and we will only keep the first one!')
                        else:
                            if file_path not in result_dict[sample][1]:
                                result_dict[sample][1].append(file_path)
                            else:
                                print(f'warn: duplicated path found for {file_path}, and we will only keep the first one!')

    new_result = dict()
    if link_data:
        os.mkdir('rawdata')
        os.chdir('rawdata')
    for sample, lst in result_dict.items():
        read1 = sorted(lst[0])
        read2 = sorted(lst[1])
        if middle2underscore:
            sample = sample.replace('-', '_')
        if add_s_to_numeric_name:
            if sample.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                sample = 'S' + sample
        new_result[sample] = [read1, read2]
        if link_data:
            # make link
            os.mkdir(sample)
            for each in read1:
                os.symlink(each, os.path.join(sample, os.path.basename(each)))
            for each in read2:
                os.symlink(each, os.path.join(sample, os.path.basename(each)))

    if out.endswith('.json'):
        with open(out, 'w') as f:
            json.dump(new_result, f, indent=2)
    else:
        with open(out, 'w') as f:
            for k, v in new_result.items():
                read1 = ';'.join(v[0])
                read2 = ';'.join(v[1])
                f.write(f'{k}\t{read1}\t{read2}\n')

    return new_result
