"""
设计思路
1. 定义argument,runtime,outputs,meta
    argument: 全面描述参数，最终可以根据参数列表自动组装成命令行
    runtime：运行环境描述，如软件，docker image，memory，cpu限制等信息
    outputs: 描述程序输出
    meta: 描述程序/command的信息，如名称，作者，版本等信息
2. 由上述对象进一步定义Command对象, 即Command = {Argument, RunTime, Output, Meta}
3. 定义Task对象：
    1. 给Command的参数赋值
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
