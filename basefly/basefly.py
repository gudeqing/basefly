import glob
import os
import shutil
import json
import sys
import re
import argparse
import textwrap
from functools import partial
from uuid import uuid4, UUID
from dataclasses import dataclass, field
try:
    from typing import Any, List, Dict, Literal
except:
    from typing import Any, List, Dict
    from typing_extensions import Literal
# shadow default dict
# from munch import Munch as dict
from .runner import run_wf
from .runner import RunCommands

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
7. 定义方法将workflow转换成wdl流程等（放弃）、argo_workflow（待完善）、自定义的流程(不妨命名为basefly)

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
    # type is one of ['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix', 'outstr']
    # fix 类型表示该参数并不是真正的参数，其为固定的字符串. 例如其可以用来表示管道符如‘| samtools sort’，方便拼接多个命令
    # outstr 类型表示该参数是用于形成输出文件名或路径的参数，引入该类型的目的是为了输出的wf.args.json更准确
    type: Literal['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix', 'outstr'] = 'str'
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
    # editable字段用于表示该参数是否适合作为流程参数进行修改，如果不可以，则不会显示在参数的配置文件中wf.args.json
    editable: bool = True

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
            self.level = 'optional'
            self.editable = False
            if self.value is None:
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
        if type(self.range) in {set, list}:
            if self.default not in self.range:
                raise Exception(f'default value is not in {self.range}')

        if self.type == 'outstr':
            self.editable = False


class RunTime:
    image: str = None
    # 软件所在目录
    tool_dir: str = ''
    # 软件名称
    tool: str = ''
    # 执行环境设置, 以下为命令运行所需的最少计算资源, memory的单位为字节
    memory: int = 1024
    cpu: int = 2
    # 执行环境设置，以为命令行所需的最大计算资源
    max_memory: int = 0
    max_cpu: int = 0
    # 运行时间上限默认7天
    timeout: int = 3600*24*7
    # docker cmd prefix
    docker_cmd_prefix2: str = 'docker run --rm --privileged --user `id -u`:`id -g` -i --entrypoint /bin/bash'
    docker_cmd_prefix: str = 'docker run --rm --privileged -i --entrypoint /bin/bash'


@dataclass()
class Output:
    path: str = None  # will be deprecated, replace with value
    value: Any = None
    # out_id: str = field(default_factory=uuid4)
    type: Literal['str', 'int', 'float', 'bool', 'outfile', 'outdir'] = 'outfile'
    # 设计report参数，用于指定该参数最终是否作为流程的outputs
    report: bool = False
    # 由command形成task之后，可以带入task_id
    task_id: str = None
    name: str = None
    desc: str = None
    format: str = None

    def __post_init__(self):
        if self.value is None:
            # value属性是后来修改增加的，本来用来替代path，为了向前兼容，保留path属性
            self.value = self.path
        else:
            self.path = self.value


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
    
    def __post_init__(self):
        # other_args是本设计自留的特殊参数,可以用来传递用户从来没有定义但是软件本身确实包含的参数
        self.args['other_args'] = Argument(prefix='', default='', level='optional', desc='This argument is designed to provide any arguments that are not wrapped in Command')

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
                arg_value = arg.value if arg.value is not None else arg.default
            else:
                arg_value = arg.value

            if arg_value is None:
                if arg.level == 'required':
                    raise Exception(f'No value found for required argument {arg_name} in {self.meta.name}')
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
                try:
                    arg_value = arg_value.value.replace('~', '').format(**value_dict)
                except KeyError as e:
                    print(e, f'=> failed to format value of {arg_name} with value {arg_value.value}')
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
                        # 当前设计中可以用’~{other_arg_name}‘或者‘{other_arg_name}’引用其他参数形成outputs的值
                        arg_value[ind] = each.value.replace('~', '').format(**value_dict)
                    elif type(each) == TopVar or type(each) == TmpVar:
                        arg_value[ind] = each.value
                arg_value = [x for x in arg_value if x is not None]
                if len(arg_value) == 0:
                    continue

            # 对于可以接收多个值的参数
            if arg.array:
                if not arg.multi_times:
                    arg_value = arg.delimiter.join([str(x) for x in arg_value if x is not None])
                else:
                    # 这里可以处理多值且可重复使用的参数，形如："--i 1 3 -i 2 5", 因此做如下处理
                    if all([type(v) == list or type(v) == tuple for v in arg_value]):
                        arg_value = [arg.delimiter.join([str(x) for x in v if x is not None]) for v in arg_value]
                    else:
                        raise Exception('如果一个参数可以接受多个值，而且可以重复使用，必须传入嵌套列表作为值！')
                if len(arg_value) == 0:
                    continue

            if not arg.array and arg.multi_times:
                arg_value = [str(x) for x in arg_value]

            # 处理bool值参数
            if arg.type == "bool" or type(arg.value) == bool:
                if arg_value:
                    cmd += ' ' + arg.prefix
                else:
                    # 如果bool类型且value为false，则该参数不参与命令行形成
                    continue
            else:
                if arg_value is not None:
                    if not arg.multi_times:
                        if '{}' not in arg.prefix:
                            cmd += ' ' + arg.prefix + str(arg_value)
                        else:
                            cmd += ' ' + arg.prefix.format(arg_value)
                    else:
                        if arg_value:
                            cmd += ' ' + arg.prefix + (' ' + arg.prefix).join(arg_value)
        return cmd

    def format_wdl_task(self, outfile=None, wdl_version='development'):
        # 形成wdl格式的task非常复杂，因此单独提供ToWdlTask处理，该函数只有当你仅需要某个task的wdl版本时需要
        if not outfile:
            outfile = f'{self.meta.name}.wdl'
        ToWdlTask(self, outfile, wdl_version)

    def run_now(self, wkdir, wf_tasks=None, docker=False):
        import subprocess
        print(f'Running {self.meta.name} ...')
        os.makedirs(wkdir, exist_ok=True)
        if self.runtime.image and docker:
            # get mount volumes
            mount_vols = {wkdir}
            for k, v in self.args.items():
                if type(v.value) in [list, tuple]:
                    values = v.value
                else:
                    values = [v.value]
                for value in values:
                    if (type(value) == TopVar or type(value) == TmpVar) and value.type in ['infile', 'indir'] and value.value is not None:
                        if value.type == 'infile':
                            file_dir = os.path.dirname(value.value)
                            mount_vols.add(os.path.abspath(file_dir))
                        elif value.type == 'indir':
                            mount_vols.add(os.path.abspath(value.value))
                    elif (type(value) == str) and (v.type in ['infile', 'indir']):
                        if v.type == 'infile':
                            mount_vols.add(os.path.abspath(os.path.dirname(value)))
                        elif v.type == 'indir':
                            mount_vols.add(os.path.abspath(value))

            with open(os.path.join(wkdir, 'cmd.sh'), 'w') as f:
                # f.write('set -o pipefail\n')
                f.write(self.format_cmd(wf_tasks) + '\n')
                f.write(f'chown -R {os.getuid()}:{os.getgid()} {wkdir}' + '\n')
            docker_cmd = self.runtime.docker_cmd_prefix
            for each in mount_vols:
                docker_cmd += f' -v {each}:{each} '
            docker_cmd += f'-w {wkdir} {self.runtime.image} cmd.sh'
            subprocess.check_call(docker_cmd, cwd=wkdir, shell=True)
        else:
            subprocess.check_call(self.format_cmd(wf_tasks), cwd=wkdir, shell=True, executable='/bin/bash')
        # format output
        value_dict = dict()
        for k, v in self.args.items():
            if type(v.value) == TopVar:
                v.value = v.value.value
            value_dict[k] = v.value or v.default
        invalid_outs = []
        for name, out in self.outputs.items():
            out.value = out.value.replace('~', '').format(**value_dict)
            out.value = os.path.join(wkdir, out.value)
            if out.type in ['outfile', 'outdir'] and not os.path.exists(out.value):
                invalid_outs.append(name)
        for each in invalid_outs:
            # print('drop invalid outputs of', self.meta.name, self.outputs[each].value)
            self.outputs.pop(each)


@dataclass()
class Task:
    cmd: Command
    name: str = None
    tag: str = None
    parent_wkdir: str = ""
    task_id: str = field(default_factory=uuid4)
    depends: List[str] = field(default_factory=list)
    wkdir: str = ''

    def __post_init__(self):
        # task name
        if self.name is None:
            if self.tag:
                self.name = self.cmd.meta.name + '-' + str(self.tag)
            else:
                # self.name = self.cmd.meta.name + '-' + str(self.task_id)
                self.name = self.cmd.meta.name

        # 为每一个output带入
        for key in self.cmd.outputs.keys():
            self.cmd.outputs[key].task_id = self.task_id
            self.cmd.outputs[key].name = key
        # 继承 cmd 的outputs
        self.outputs = self.cmd.outputs
        for ind, each in enumerate(self.depends):
            if hasattr(each, 'task_id'):
                self.depends[ind] = each.task_id
        self.depends = [x for x in self.depends if x is not None]

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
            if self.value is not None:
                if os.path.exists(self.value):
                    self.value = os.path.abspath(self.value)
                else:
                    raise FileExistsError(self.value)


@dataclass()
class TmpVar:
    """
    该对象用于描述流程中如循环时的临时变量, 纯粹是为wdl的scatter模式设计, name属性将作为wdl文件中的传递值。
    例如：如循环中某个变量为A，则把输出文件名定义为out_file = ’~{A}.txt‘.
    那么，这时脚本中可以写: outfile = TmpVar(name="~{A}.txt", type='str') , 这里name中“~{}”是wdl传递变量值的语法。
    后来发现想直接生成完整的wdl流程几乎不可能，该函数将来可能会被删除，不建议使用
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
    topvars: Dict[str, TopVar] = field(default_factory=dict)
    argparser = None
    args = None
    success = False
    add_argument = None
    wkdir: str = None
    task_order: int = 0

    def __post_init__(self):
        for k, v in self.topvars.items():
            # 将key作为var的名称
            v.name = k

    def run(self):
        """
        生成类config格式的workflow(->wf.ini)并且允许直接本地执行
        """
        self.finalize()
        import configparser
        wf = configparser.ConfigParser()
        wf.optionxform = str
        parameters = self.argparser.parse_args()
        outdir = os.path.abspath(parameters.outdir)
        self.wkdir = outdir
        self.parameters = parameters

        if parameters.skip:
            self.skip_steps(parameters.skip, skip_depend=not parameters.no_skip_depend)

        if parameters.to_cwl:
            self.to_cwl_workflow(outdir)

        if parameters.rerun_steps:
            rerun_steps = self.get_all_depends(parameters.rerun_steps)
            print('rerun the following steps:', rerun_steps)
        else:
            rerun_steps = tuple()

        if parameters.update_args:
            self.update_args(parameters.update_args)

        if parameters.assume_success_steps:
            all_task_names = set(x.name for x in self.tasks.values())
            all_cmd_names = set(x.cmd.meta.name for x in self.tasks.values())
            assume_success_tasks = set()
            for each in parameters.assume_success_steps:
                matched = []
                if each in all_task_names:
                    matched.append(assume_success_tasks)
                elif each.endswith('*'):
                    matched += [x for x in all_task_names if x.startswith(each[:-1])]
                elif each in all_cmd_names:
                    matched += [x for x in self.tasks.values() if x.cmd.meta.name == each]
                if not matched:
                    raise Exception(f'{each} matches no task, you may check the task name by "--list_task"')
                assume_success_tasks.update(matched)
        else:
            assume_success_tasks = set()

        wf['mode'] = dict(
            outdir=outdir,
            threads=parameters.threads,
            retry=parameters.retry,
            monitor_resource=parameters.monitor_resource,
            monitor_time_step=3,
            check_resource_before_run=not parameters.no_check_resource_before_run,
        )

        for task_id, task in self.tasks.items():
            cmd_wkdir = os.path.join("${mode:outdir}", task.parent_wkdir, task.name)
            task.wkdir = os.path.join(self.wkdir, task.parent_wkdir, task.name)
            mount_vols = {cmd_wkdir}

            # format output
            value_dict = dict()
            for k, v in task.cmd.args.items():
                value_dict[k] = v.value or v.default
            for _, out in task.outputs.items():
                out.value = out.value.replace('~', '').format(**value_dict)

            if task.name in assume_success_tasks:
                task.depends = None
                wf[task.name] = dict(
                    depend='',
                    cmd="using previous result",
                    mem=task.cmd.runtime.memory,
                    cpu=task.cmd.runtime.cpu,
                    max_mem=task.cmd.runtime.max_memory,
                    max_cpu=task.cmd.runtime.max_cpu,
                    timeout=task.cmd.runtime.timeout,
                    image='' if not parameters.docker else (task.cmd.runtime.image or ''),
                    wkdir=cmd_wkdir,
                    mount_vols=';'.join(mount_vols)
                )
                continue

            # 参数赋值处理
            for k, v in task.cmd.args.items():
                if type(v) != Argument:
                    raise Exception(f'Wrong assignment of value to {k} with {v}')
                if v.level == 'required' and v.value is None:
                    v.value = v.default
                if type(v.value) in [list, tuple]:
                    values = v.value
                else:
                    values = [v.value]
                for value in values:
                    if type(value) == Output and value.type in ['outfile', 'outdir']:
                        if not value.value.startswith('${{mode:'):
                            try:
                                value.value = os.path.join("${{mode:outdir}}", self.tasks[value.task_id].parent_wkdir, self.tasks[value.task_id].name, value.value)
                            except Exception as e:
                                print(e, value, task.name)
                                print(f'you may skip this step {task.name} by --assume_success_step')
                        # mount_vols.add(os.path.join(outdir, self.tasks[value.task_id].name))
                        mount_vols.add(os.path.join("${mode:outdir}", self.tasks[value.task_id].parent_wkdir, self.tasks[value.task_id].name))

                    elif (type(value) == TopVar or type(value) == TmpVar) and value.type in ['infile', 'indir'] and value.value is not None:
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
                image='' if not parameters.docker else (task.cmd.runtime.image or ''),
                docker_cmd_prefix=task.cmd.runtime.docker_cmd_prefix,
                wkdir=cmd_wkdir,
                mount_vols=';'.join(mount_vols)
            )

        if parameters.list_cmd:
            self.list_cmd()
        elif parameters.show_cmd:
            self.show_cmd(parameters.show_cmd)
        elif parameters.list_task:
            self.list_task()
        elif parameters.dry_run:
            os.makedirs(outdir, exist_ok=True)
            self.generate_docs(os.path.join(outdir, f'{self.meta.name}.ReadMe.md'))
            self.dump_args(out=os.path.join(outdir, 'wf.static.args.json'))
            with open(os.path.join(outdir, "wf.run.cmd.txt"), 'w') as f:
                args = []
                for each in sys.argv:
                    if {'(', '{', ';'} & set(each):
                        args.append('"' + each + '"')
                    else:
                        args.append(each)
                f.write('python ' + ' '.join(args) + '\n')
                # print(sys.argv)
                f.write('>>>Argument Detail\n')
                f.write('{}\n'.format(dict(parameters.__dict__.items())))
            outfile = os.path.join(outdir, f'{self.meta.name}.ini')
            with open(outfile, 'w') as configfile:
                wf.write(configfile)
            # 仅仅为了生成流程图
            RunCommands(outfile, draw_state_graph=True)
        elif parameters.run:
            os.makedirs(outdir, exist_ok=True)
            self.dump_args(out=os.path.join(outdir, 'wf.static.args.json'))
            with open(os.path.join(outdir, "wf.run.cmd.txt"), 'w') as f:
                args = []
                for each in sys.argv:
                    if {'(', '{', ';'} & set(each):
                        args.append('"' + each + '"')
                    else:
                        args.append(each)
                f.write('python ' + ' '.join(args) + '\n')
                # print(sys.argv)
                f.write('>>>Argument Detail\n')
                f.write('{}\n'.format(dict(parameters.__dict__.items())))
            outfile = os.path.join(outdir, f'{self.meta.name}.ini')
            with open(outfile, 'w') as configfile:
                wf.write(configfile)

            wf = run_wf(
                outfile,
                timeout=parameters.wait_resource_time,
                assume_success_steps=tuple(assume_success_tasks),
                plot=parameters.plot,
                rerun_steps=rerun_steps
            )
            self.success = wf.failed == 0

            # 通过软连接汇总输出目录
            outputs_dir = os.path.join(outdir, 'Outputs')
            os.makedirs(outputs_dir, exist_ok=True)
            shutil.copyfile(outfile, os.path.join(outputs_dir, f'{self.meta.name}.ini'))
            shutil.copyfile(os.path.join(outdir, 'wf.static.args.json'), os.path.join(outputs_dir, 'wf.static.args.json'))
            shutil.copyfile(os.path.join(outdir, 'wf.run.cmd.txt'), os.path.join(outputs_dir, 'wf.run.cmd.txt'))
            shutil.copyfile(os.path.join(outdir, 'state.svg'), os.path.join(outputs_dir, 'state.svg'))
            for name, out in self.outputs.items():
                if '${{' in out.value:
                    src_dir = out.value.replace('${{mode:outdir}}', outdir)
                else:
                    if out.task_id in self.tasks:
                        src_dir = os.path.join(self.tasks[out.task_id].wkdir, out.value)
                        if src_dir.endswith('/.'):
                            src_dir = src_dir[:-1]
                    else:
                        continue
                scan_dir = glob.glob(src_dir)
                if scan_dir:
                    src_dir = scan_dir[0]
                if os.path.exists(src_dir):
                    targets = [src_dir]
                    parent_dir = self.tasks[out.task_id].wkdir
                    if os.path.isfile(src_dir):
                        tmp_lst = os.listdir(parent_dir)
                        if 'cmd.sh' in tmp_lst:
                            targets.append(os.path.join(parent_dir, 'cmd.sh'))
                        if 'docker.cmd.sh' in tmp_lst:
                            targets.append(os.path.join(parent_dir, 'docker.cmd.sh'))
                    for src_dir in targets:
                        # print('Found expected output: ', src_dir)
                        final_out_dir = os.path.join(outdir, 'Outputs', self.tasks[out.task_id].name)
                        os.makedirs(final_out_dir, exist_ok=True)
                        dst_path = os.path.join(final_out_dir, os.path.basename(src_dir))
                        if dst_path.endswith('/.'):
                            dst_path = dst_path[:-1]
                        dst_path = dst_path.rstrip('/')
                        # 删除已经存在的结果
                        if os.path.exists(dst_path):
                            if os.path.isfile(dst_path):
                                os.remove(dst_path)
                            else:
                                if not os.path.islink(dst_path):
                                    shutil.rmtree(dst_path)
                                else:
                                    os.remove(dst_path)
                        # 如果输出结果是文件，则创建硬链接，否则软连接
                        if os.path.isfile(src_dir):
                            os.link(src_dir, dst_path)
                        else:
                            os.symlink(src_dir.rstrip('/'), dst_path)
                else:
                    print('Failed to found expected output: ', src_dir)
        else:
            print('No actions, you may provide one action parameter: --run, --dry_run, --list_cmd, --list_task, -show_cmd')

    def add_topvars(self, var_dict):
        for k, v in var_dict.items():
            # 将key作为var的名称
            v.name = k
            if v.type in ['infile', 'indir'] and v.value:
                v.value = os.path.abspath(v.value)
        self.topvars.update(var_dict)

    def add_task(self, cmd: Command, depends: list = [], parent_wkdir: str ='', name: str = None, tag: str = None):
        self.task_order += 1
        task = Task(cmd=cmd, depends=depends.copy(), name=name, tag=tag, parent_wkdir=parent_wkdir)
        existed_names = {x.name for x in self.tasks.values()}
        if task.name in existed_names:
            raise Exception(f'{task.name} duplicated, please rename it')
        # 更新工作目录
        task.wkdir = os.path.join(self.wkdir, task.parent_wkdir, task.name)
        self.tasks[task.task_id] = task
        return task, task.cmd.args

    def get_task_by_name(self, matcher):
        for task_id, task in self.tasks.items():
            if task.name.startswith(matcher):
                return task
        raise Exception(f'cannot found task {matcher}')

    def dump_args(self, out='arguments.json'):
        """
        该函数的目的是为了输出流程中每个步骤的参数json模板文件，方便客户后续修改参数
        如果一个command被调用多次，则仅记录第一次被调用的参数
        """
        if self.parameters.dry_run:
            print('如果后续需要输入本次生成的wf.static.args.json，强烈建议仅保留被修改的参数，同时注意一个command会被重复使用的情况')
        cmd_names = set()
        arg_value_dict = dict()
        for tid, task in self.tasks.items():
            cmd_name = task.cmd.meta.name
            if cmd_name not in cmd_names:
                tmp_dict = arg_value_dict.setdefault(cmd_name, dict())
                for arg_name, arg in task.cmd.args.items():
                    if arg.editable:
                        # 对于必须参数，又没有默认值的，那么其一定是在流程中动态赋值的参数
                        if arg.level == 'required' and (arg.default is None):
                            continue
                        # 对于赋值为列表的参数, 且赋值特殊，也不能输出
                        if type(arg.value) == list:
                            if type(arg.value[0]) in {TopVar, TmpVar, Output}:
                                continue
                        # 对赋值对象为TopVar的参数，其一定是流程最开始的输入参数，也强制归属为动态赋值参数
                        if type(arg.value) in {TopVar, TmpVar, Output}:
                            # tmp_dict[arg_name] = arg.value.value
                            continue
                        else:
                            tmp_dict[arg_name] = arg.value or arg.default
                        # 对于非必须参数，如果没有赋值，则将赋值设置为None
                        if arg.level == 'optional' and arg.value is None:
                            tmp_dict[arg_name] = None

            cmd_names.add(cmd_name)
        # print(arg_value_dict)
        with open(out, 'w') as f:
            json.dump(arg_value_dict, f, indent=2)

    def update_args(self, arg_json_file):
        with open(arg_json_file) as f:
            cfg = json.load(f)
        for tid, task in self.tasks.items():
            for arg_name, arg in task.cmd.args.items():
                cmd_name = task.cmd.meta.name
                if (cmd_name in cfg) and (arg_name in cfg[cmd_name]):
                    arg.value = cfg[cmd_name][arg_name]

    def skip_steps(self, steps, skip_depend=True):
        """
        :param steps: list containing cmd.meta.name or task.name
        :param skip_depend: if also skip steps that depend on the steps
        :return:
        """
        task_copy = self.tasks.copy()
        skip_tasks = [tid for tid, x in self.tasks.items() if x.name in steps or x.cmd.meta.name in steps]
        # pop task
        _ = [self.tasks.pop(x) for x in skip_tasks]
        if skip_depend:
            # find out depending task
            total_deduced_skips = list()
            while True:
                to_be_skip = list()
                for tid in self.tasks.keys():
                    # print(self.tasks[tid].depends, self.tasks[tid].name)
                    if set(self.tasks[tid].depends) - set(self.tasks.keys()):
                        to_be_skip.append(tid)
                        total_deduced_skips.append(tid)
                _ = [self.tasks.pop(x) for x in to_be_skip]
                # judge if all dependencies are available
                if all(len(set(self.tasks[x].depends) - set(self.tasks.keys())) == 0 for x in self.tasks):
                    break
            skip_tasks += total_deduced_skips
        print(f'total {len(skip_tasks)} tasks will be skipped', [task_copy[tid].name for tid in skip_tasks])

    def get_all_depends(self, steps):
        """
        给定一个或多个step的名称或者前缀(要求字符串末尾加*标志），找到下游所有可能依赖他们的steps
        :param steps:
        :return: 列表，包含起始步骤和所有依赖他们的步骤
        """
        all_tasks = self.tasks.copy()
        step_names = [x.name for _, x in all_tasks.items()]
        all_matched_steps = []
        for each in steps:
            match_steps = []
            if each in step_names:
                match_steps.append(each)
            else:
                if each.endswith('*'):
                    match_steps += [x for x in step_names if x.startswith(each[:-1])]
            if match_steps:
                all_matched_steps += match_steps
            else:
                raise Exception(f'{each} matches no task, you may check the task name by "--list_task"')
        init_tasks = [tid for tid, x in self.tasks.items() if x.name in all_matched_steps]
        followed_tasks = list()
        # 从所有的tasks中删掉起始步骤
        _ = [all_tasks.pop(x) for x in init_tasks]
        # 通过循环找到依赖的步骤
        while True:
            tmp_followed = list()
            for tid in all_tasks.keys():
                # print(self.tasks[tid].depends, self.tasks[tid].name)
                # 检查当前task的依赖步骤是否全部在更新后的tasks中
                if set(all_tasks[tid].depends) - set(all_tasks.keys()):
                    # 有些依赖没有找到，说明当前步骤依赖曾经删除的步骤
                    tmp_followed.append(tid)
                    followed_tasks.append(tid)
            if tmp_followed:
                # 从all task中删除 tmp_followed
                _ = [all_tasks.pop(x) for x in tmp_followed]
            else:
                # 说明所有task的依赖都存在，可以停止循环
                break
        target_tids = set(init_tasks + followed_tasks)
        target_names = [self.tasks[tid].name for tid in target_tids]
        return tuple(target_names)

    def list_cmd(self):
        print(sorted({tsk.cmd.meta.name for _, tsk in self.tasks.items()}))

    def list_task(self):
        print(sorted([tsk.name for _, tsk in self.tasks.items()]))

    def show_cmd(self, cmd_name):
        for _, tsk in self.tasks.items():
            if tsk.cmd.meta.name == cmd_name:
                print(tsk.cmd.format_cmd(self.tasks))
                break

    def init_argparser(self):
        if len(sys.argv) <= 1:
            exit('Please provide at least one argument, use -h for help')
        parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description=self.meta.desc
        )
        wf_args = parser.add_argument_group('Arguments for controlling running mode')
        wf_args.add_argument('-outdir', metavar='workdir', default=os.path.join(os.getcwd(), 'Result'), help='结果目录')
        wf_args.add_argument('--run', default=False, action='store_true', help="默认不运行流程, 仅生成流程; 如果outdir目录已经存在cmd_state.txt文件，则自动续跑")
        wf_args.add_argument('--docker', default=False, action='store_true', help="指示是否使用docker镜像创造运行环境")
        wf_args.add_argument('--plot', default=False, action='store_true', help="该参数用于通知是否对分析流程进行实时绘图，生成的流程图可以显示任务之间的依赖关系,也可以显示每个任务的参数信息，同时还能根据每个节点的颜色判断任务的运行状态以及任务完成需要的时间'")
        wf_args.add_argument('-threads', metavar='max-workers', default=3, type=int, help="允许的最大并行的任务数量, 默认3")
        wf_args.add_argument('-update_args', metavar='update-args', required=False, help="输入参数配置文件, json格式. 流程每次运行时都会输出wf.args.json, 由于其中包含一些样本信息参数，不可直接使用，但可以根据此模板修改")
        wf_args.add_argument('-skip', metavar=('step1', 'task3'), default=list(), nargs='+', help='指定要跳过的步骤或具体task,空格分隔,默认程序会自动跳过依赖他们的步骤, 使用--list_cmd or --list_task可查看候选')
        wf_args.add_argument('--no_skip_depend', default=False, action='store_true', help="当使用skip参数时, 如果同时指定该参数，则不会自动跳过依赖的步骤")
        wf_args.add_argument('-rerun_steps', metavar=('task3', 'task_prefix'), default=list(), nargs='+', help="指定需要重跑的步骤，不论其是否已经成功完成，空格分隔, 这样做的可能原因可能是: 重新设置了命令参数. 使用--list_task可查看候选,可使用task的前缀，并且以'*'结尾，将自动匹配符合前缀条件的所有task")
        wf_args.add_argument('-assume_success_steps', metavar=('task_name', 'task_cmd_name'), default=list(), nargs='+', help="假定哪些步骤已经成功运行，不论其是否真的已经成功完成，空格分隔, 这样做的可能原因: 利用之前已经成功运行的结果(需要把之前的运行结果放到当前结果目录). 使用--list_task可查看候选, 可使用task的前缀，并且以'*'结尾，将自动匹配符合前缀条件的所有task,也可以使用cmd.meata.name")
        wf_args.add_argument('-retry', metavar='max-retry', default=1, type=int, help='某步骤运行失败后再尝试运行的次数, 默认1次. 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini')
        wf_args.add_argument('--list_cmd', default=False, action="store_true", help="仅仅显示当前流程包含的主步骤, 不会显示指定跳过的步骤")
        wf_args.add_argument('-show_cmd', metavar='cmd-query', help="提供一个cmd名称,输出该cmd的样例")
        wf_args.add_argument('--list_task', default=False, action="store_true", help="仅仅显示当前流程包含的详细步骤, 且已经排除指定跳过的步骤")
        wf_args.add_argument('--monitor_resource', default=False, action='store_true', help='是否监控每一步运行时的资源消耗, 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini')
        wf_args.add_argument('-wait_resource_time', metavar='wait-time', default=900, type=int, help="等待资源的时间上限, 默认每次等待时间为900秒, 等待时间超过这个时间且资源不足时判定任务失败")
        wf_args.add_argument('--no_check_resource_before_run', default=False, action='store_true', help="指示运行某步骤前检测指定的资源是否足够, 如不足, 则该步骤失败; 如果设置该参数, 则运行前不检查资源. 如需对某一步设置不同的值,可运行前修改pipeline.ini. 如需更改指定的资源, 可在运行流程前修改pipeline.ini")
        wf_args.add_argument('--dry_run', default=False, action='store_true', help='不运行流程，仅仅输出markdown格式的流程说明文档和流程配置文件')
        wf_args.add_argument('--to_cwl', default=False, action='store_true', help='输出cwl格式的流程草稿')
        self.argparser = parser
        # for user defined arguments
        # self.add_argument = partial(self.argparser.add_argument, required=True)
        self.add_argument = self.argparser.add_argument

    def parse_args(self):
        if self.argparser is None:
            self.init_argparser()
        self.args = self.argparser.parse_args()
        self.wkdir = self.args.outdir

    def finalize(self):
        for task_id, task in self.tasks.items():
            # depends统一转为task_id
            for ind, each in enumerate(task.depends):
                if hasattr(each, 'task_id'):
                    task.depends[ind] = each.task_id
                elif type(each) != UUID:
                    raise Exception(f'valid "depends" for task {task.name} should be UUID object or Task Object but not "{each}"')

            # 根据输入数据矫正忘了添加的依赖
            for k, v in task.cmd.args.items():
                if type(v.value) == Output:
                    if v.value.task_id not in task.depends:
                        # print(f'You may forgot to add dependency {self.tasks[v.value.task_id].name} for {task.name}. We will fixed it.')
                        task.depends.append(v.value.task_id)
                elif type(v.value) == list:
                    for each in v.value:
                        if type(each) == Output:
                            if each.task_id not in task.depends:
                                task.depends.append(each.task_id)
                else:
                    pass
            # 为了生成cwl，下面的代码移动到run中
            # # format output
            # value_dict = dict()
            # for k, v in task.cmd.args.items():
            #     value_dict[k] = v.value or v.default
            # for _, out in task.outputs.items():
            #     out.value = out.value.replace('~', '').format(**value_dict)

            # 添加workflow的output
            for _name, out in task.outputs.items():
                if out.report:
                    self.outputs[task.name+'.'+_name] = out

    def generate_docs(self, out):
        """
        文档结构：
        1.	引言
        1.1.	目的
        1.2.	项目背景
        1.3.	术语定义
        1.4.	参考资料
        2.	流程概述
        2.1.	流程名称
        2.2.	流程说明
        2.3.	流程主要功能
        3.	流程总体结构
        3.1.	模块功能汇总表
        3.2.	流程总体逻辑结构图
        3.3.	流程运行控制：
                    输入参数
                    输出文件
        3.4.	流程运行环境
        4.	模块
        4.1.	模块说明
        4.2.	主要功能
        4.3.	模块输入参数
        4.4.	模块输出
        """
        tool_names = []
        tools = []
        for task_id, task in self.tasks.items():
            if task.cmd.meta.name not in tool_names:
                tool_names.append(task.cmd.meta.name)
                tools.append(task.cmd)
        contents = []
        contents += ['## 系统概述']
        contents += ['### 系统名称']
        contents += [f'* {self.meta.name}']
        contents += ['### 系统说明']
        contents += [f'* 系统版本号：{self.meta.version}']
        contents += [f'* 适用范围：{self.meta.source}']
        contents += [f'* 参考来源：{self.meta.source}']
        contents += ['### 系统主要功能']
        contents += [f'{self.meta.desc}']
        contents += ['## 系统总体结构']
        contents += ['### 模块列表']
        contents += ["| name | desc | source | vesion |"]
        contents += ["| :--- | :---: | :---: | :---: |"]
        for cmd in tools:
            contents += [f"|{cmd.meta.name}|{cmd.meta.desc}|{cmd.meta.source}|{cmd.meta.version}|"]
        contents += ['### 系统总体结构图']
        contents += ["本系统为数据分析流程，流程结构图可能随参数有所变化，一下仅为一个典型示例图"]
        contents += [f'![分析流程的有向无循环(DAG)图](./state.svg "{self.meta.name}")']
        contents += ['## 系统运行环境']
        contents += ['该系统主要内容是分析流程，需要在Linux操作系统下运行，需安装python > 3.7的环境，环境中需要安装graphviz和xcmds等包. '
                     '一般情况下，该环境完全由docker容器提供，包括流程运行每个环节所需要的软件或工具，具体所需软件工具可参考前面章节的”模块列表“.']
        contents += ['## 系统运行控制']
        contents += ['* 系统跟随容器一起运行，理论上没有时间限制，主要取决于容器运行周期，容器运行周期又取决于流程本身的具体任务执行情况，'
                     '分析任务完成后, 容器自动停止。针对系统所包含的分析流程而言，有如下参数控制流程的运行，以方便用户能够根据具体情况分析数据']
        for arg in self.argparser.__dict__['_actions']:
            contents += [f'+ **{arg.dest}**: {arg.help}']

        for tool in tools:
            contents += [f'## 模块 {tool.meta.name}']
            contents += [f'### 模块说明']
            contents += [f'* 简介: {tool.meta.desc}']
            contents += [f'* 参考：{tool.meta.source}']
            contents += [f'### 模块运行环境']
            contents += [f'* 镜像：{tool.runtime.image}']
            contents += [f'* CPU: {tool.runtime.cpu}']
            contents += [f'* Memory: {tool.runtime.memory}']
            contents += [f'### 模块输入文件参数']
            for name, arg in tool.args.items():
                if not name.startswith('_') and (arg.type in ['infile', 'indir']):
                    contents += [f'#### {name}']
                    contents += [f'+ 类型: {arg.type}']
                    contents += [f'+ 默认值: {arg.default}']
                    contents += [f'+ 描述: {arg.desc}']
            contents += [f'### 模块普通参数']
            for name, arg in tool.args.items():
                if not name.startswith('_') and  (arg.type not in ['infile', 'indir']):
                    contents += [f'#### {name}']
                    contents += [f'+ 类型: {arg.type}']
                    contents += [f'+ 默认值: {arg.default}']
                    contents += [f'+ 描述: {arg.desc}']
            contents += [f'### 模块输出']
            for key, output in tool.outputs.items():
                contents += [f'#### {key}']
                contents += [f'+ *type*: {output.type}']
                contents += [f'+ *value example*: {output.value}']
                contents += [f'+ *description*: {output.desc}']

        with open(out, 'w') as f:
            for each in contents:
                if each.startswith('#'):
                    f.write('\n')
                f.write(each + '\n')

    def to_wdl_tasks(self, outfile=None):
        # 输出每一个Command的WDL版本
        with open(outfile or f'{self.meta.name}.tasks.wdl', 'w') as f:
            seen = set()
            for _, task in self.tasks.items():
                if task.cmd.meta.name not in seen:
                    seen.add(task.cmd.meta.name)
                    wdl_str = ToWdlTask(task.cmd).wdl
                    f.write(wdl_str + '\n')

    def to_wdl_workflow(self, outfile=None):
        print('目前不能保证能生产完全正确的WDL流程, 甚至可能出错，不再维护该功能. 建议使用to_wdl_task后自己编写workflow')
        ToWdlWorkflow(self).write_wdl(outfile)

    def to_argo_workflow(self, outfile=None):
        print('目前生成argo workflow 不一定正确')
        outfile = outfile or f'{self.meta.name}.argo.yaml'
        lines = ['apiVersion: argoproj.io/v1alpha1']
        lines += ['kind: Workflow']
        lines += ['metadata:']
        lines += [' ' * 2 + f'generateName: {self.meta.name.lower()}']
        lines += ['spec:']
        # entry point
        lines += [' '*2 + 'entrypoint: main']
        artifacts = [k for k, v in self.topvars.items() if v.type in ['infile', 'indir']]
        if artifacts:
            lines += [' ' * 2 + 'arguments:']
            lines += [' ' * 4 + 'artifacts:']
            for each in artifacts:
                if type(self.topvars[each].value) != list:
                    lines += [' ' * 4 + f'- name: {each}']
                    lines += [' ' * 6 + f'path: {self.topvars[each].value}']
                else:
                    for i, v in enumerate(self.topvars[each].value):
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

    def to_cwl_tool(self, cmd: Command):
        convert_type: dict[str, str] = {
            'str': 'string',
            'outstr': 'string',
            'infile': 'File',
            'outfile': 'File',
            'indir': 'Directory',
            'outdir': 'Directory',
            'bool': 'boolean',
            'int': 'int',
            'float': 'float',
        }
        contents = '#!/usr/bin/env cwl-runner\n'
        contents += 'cwlVersion: v1.0\n'
        contents += 'class: CommandLineTool\n'
        contents += f'label: "{cmd.meta.name}"\n'
        contents += f'doc: |\n'
        for each in textwrap.wrap(cmd.meta.desc):
            contents += ' '*4 + each + '\n'
        contents += '\n'

        contents += 'hints:\n'
        contents += ' '*2 + 'SoftwareRequirement:\n'
        contents += ' '*4 + 'packages:\n'
        contents += ' '*6 + f'{cmd.meta.name}:\n'
        contents += ' '*8 + f'specs: ["{cmd.meta.source}"]\n'
        contents += ' '*8 + f'version: ["{cmd.meta.version}"]\n'
        contents += ' '*2 + 'DockerRequirement:\n'
        contents += ' ' * 4 + f"dockerPull: {cmd.runtime.image}\n"
        contents += ' '*2 + "ResourceRequirement:\n"
        contents += ' '*4 + f"coresMin: {cmd.runtime.cpu}\n"
        # RAM的单位为M
        contents += ' '*4 + f"ramMin: {int(cmd.runtime.memory/1024)}\n"

        contents += '\n'
        base_cmd = cmd.runtime.tool_dir + cmd.runtime.tool
        base_cmd = base_cmd.split()
        if base_cmd:
            contents += f'baseCommand: {base_cmd}\n'

        contents += 'inputs:\n'
        pos = 0
        for arg_name, arg in cmd.args.items():
            if arg.type == 'fix':
                continue
            if arg_name in ('out', 'in', 'inputs', 'outputs'):
                # arg名称和cwl保留字段冲突，加x以区别
                arg_name = arg_name + 'x'
                arg.name = arg_name
            pos += 1
            contents += ' '*2 + f'{arg_name}:\n'
            arg_type = convert_type[arg.type]
            if arg.array and (not arg.multi_times):
                arg_type += '[]'
            if arg.level == 'optional':
               arg_type += '?'

            if arg.multi_times:
                if arg.array:
                    print('! 对于多值且可以重复使用的参数，可能生成结果存在问题')
                contents += ' '*4 + 'type:\n'
                contents += ' '*6 + 'type: array\n'

                # 当输入的文件是bam文件或vcf.gz文件时，需要加secondaryFiles
                if arg.value is not None and arg.type == 'infile':
                    if arg.format == 'bam':
                        contents += ' ' * 6 + 'secondaryFiles:\n'
                        contents += ' ' * 8 + '- .bai\n'
                    elif arg.format == 'vcf.gz':
                        contents += ' ' * 6 + 'secondaryFiles:\n'
                        contents += ' ' * 8 + '- .tbi\n'

                contents += ' '*6 + f'items: {arg_type}\n'
                contents += ' '*6 + 'inputBinding:\n'
                contents += ' '*8 + f'prefix: {arg.prefix.strip()}\n'
                if arg.prefix.endswith(' '):
                    contents += ' ' * 8 + 'separate: true\n'
                else:
                    contents += ' ' * 8 + 'separate: false\n'
                contents += ' '*4 + 'inputBinding:\n'
                contents += ' '*6 + f'position {pos}\n'
            else:
                contents += ' '*4 + f'type: {arg_type}\n'

                # 当输入的文件是bam文件或vcf.gz文件时，需要加secondaryFiles
                if arg.value is not None and arg.type == 'infile':
                    if arg.format == 'bam':
                        contents += ' ' * 4 + 'secondaryFiles:\n'
                        contents += ' ' * 6 + '- .bai\n'
                    elif arg.format == 'vcf.gz':
                        contents += ' ' * 4 + 'secondaryFiles:\n'
                        contents += ' ' * 6 + '- .tbi\n'

                contents += ' '*4 + 'inputBinding:\n'
                contents += ' '*6 + f'position: {pos}\n'
                contents += ' '*6 + f'prefix: {arg.prefix.strip()}\n'
                if '[' in arg_type:
                    contents += ' '*6 + f'itemSeparator: "{arg.delimiter}"\n'
                if arg.prefix.endswith(' '):
                    contents += ' '*6 + 'separate: true\n'
                else:
                    contents += ' '*6 + 'separate: false\n'
        contents += '\n'

        contents += 'outputs:\n'
        for out_name, out_obj in cmd.outputs.items():
            if out_name in ('out', 'in', 'inputs', 'outputs'):
                # arg名称和cwl保留字段冲突，加x以区别
                out_name = out_name + 'x'
            contents += ' '*2 + f'{out_name}:\n'
            contents += ' '*4 + f'type: {convert_type[out_obj.type]}\n'
            contents += ' '*4 + 'outputBinding:\n'
            out_expr = out_obj.value
            # 目前仅考支持输出表达式包含一个参数名称变量
            matched_arg_name = re.findall('{(.*?)}', out_expr)
            if matched_arg_name:
                matched_arg_name = matched_arg_name[0]
                if matched_arg_name in ('out', 'in', 'inputs', 'outputs'):
                    out_expr = out_expr.replace('{'+matched_arg_name+'}', '$(inputs.'+matched_arg_name+'x)')
                else:
                    out_expr = out_expr.replace('{'+matched_arg_name+'}', '$(inputs.'+matched_arg_name+')')
            contents += ' '*6 + f'glob: {out_expr}'

        return contents

    def to_cwl_workflow(self, outdir):
        # 生成流程草稿
        os.makedirs(outdir, exist_ok=True)
        convert_type: dict[str, str] = {
            'str': 'string',
            'outstr': 'string',
            'infile': 'File',
            'outfile': 'File',
            'indir': 'Directory',
            'outdir': 'Directory',
            'bool': 'boolean',
            'int': 'int',
            'float': 'float',
        }
        contents = '#!/usr/bin/env cwl-runner\n'
        contents += 'cwlVersion: v1.0\n'
        contents += 'class: Workflow\n'

        contents += 'inputs:\n'
        for name, top_var in self.topvars.items():
            contents += ' '*2 + f'{name}: {convert_type[top_var.type]}\n'
        contents += '\n'

        contents += 'outputs:\n'
        for name, out_obj in self.outputs.items():
            contents += ' ' * 2 + f'{name}:\n'
            contents += ' ' * 4 + f'type: {convert_type[out_obj.type]}\n'
            contents += ' ' * 4 + f'outputSource: {name.replace(".", "/", 1)}\n'
        contents += '\n'

        contents += 'steps:\n'
        for task_id, task in self.tasks.items():
            contents += ' '*2 + f'{task.name}\n'
            contents += ' '*4 + f'run: {task.cmd.meta.name}.tool.cwl\n'
            contents += ' '*4 + 'in:\n'
            # 参数赋值
            for arg_name, arg in task.cmd.args.items():
                if type(arg.value) == TopVar:
                    contents += ' ' * 6 + f'{arg_name}: {arg.value.name}\n'
                if type(arg.value) in {Output}:
                    contents += ' '*6 + f'{arg_name}:\n'
                    depend_task_name = self.tasks[arg.value.task_id].name
                    contents += ' '*8 + f'source: {depend_task_name}/{arg.value.name}\n'
            # 定义输出
            out_names = []
            for out_name, out_obj in task.cmd.outputs.items():
                out_names.append(out_name)
            contents += ' '*4 + f'out: [{",".join(out_names)}]\n'
            contents += '\n'

        # 输出流程
        outfile = os.path.join(outdir, f'{self.meta.name}.wf.cwl')
        with open(outfile, 'w') as f:
            f.write(contents)
        # 输出tools
        added_tools = []
        for task_id, task in self.tasks.items():
            cmd = task.cmd
            name = cmd.meta.name
            if name not in added_tools:
                added_tools.append(name)
                outfile = os.path.join(outdir, f'{name}.tool.cwl')
                with open(outfile, 'w') as f:
                    f.write(self.to_cwl_tool(cmd))
        # 新建2个空文件
        for each in ['LICENSE', 'R.md']:
            outfile = os.path.join(outdir, each)
            with open(outfile, 'w') as f:
                f.write('')


class ToWdlTask(object):
    type_conv_dict = {
        'str': 'String',
        'outstr': 'String',
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
            if '~' not in v.path and '{' in v.path:
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
    不再打算继续优化该函数
    """
    type_conv_dict = {
        'str': 'String',
        'outstr': 'String',
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
        if self.wf.topvars:
            for k, v in self.wf.topvars.items():
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
        wdl += '\n' + self.get_fastq_info_wdl_str() + '\n'

        # format_tasks
        for cmd in all_cmds:
            wdl += ToWdlTask(cmd).wdl
            wdl += '\n'

        # write wdl
        with open(outfile, 'w') as f:
            f.write(wdl)

    @staticmethod
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

