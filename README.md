# 简介
> Basefly主要目的是用python定义一种简便但严格的方式描述分析流程，和WDL、CWL的目的类似，旨在方便流程的编写。

## Basefly流程开发示例:
1. [mini_example](./tests/mini_example.py)
2. [gene expression quantification: fastp+salmon](./tests/rnaseq_quant.py)

## 流程运行
### basefly流程自带一个runner，支持多个任务的并发运行和可视化。缺点：流程中所有任务的命令一定都是在计算开始之前编排好的，这意味着每一步的输出文件路径和数量是要预先知道的，可以理解为流程是预先编译好的，而不是动态执行的。

### basefly自带runner的帮助文档可以参考如下bismark流程
```
$ python bismark.py -h
Warn: cannot import pygraphviz and state graph will not be drawn!
usage: bismark.py [-h] [--run] [--docker] [--plot] [-update_args update-args] [-threads max-workers] [-outdir workdir] [-skip step1 [task3 ...]] [--no_skip_depend]
                  [-rerun_steps task3 [task_prefix ...]] [-assume_success_steps task3 [task_prefix ...]] [-retry max-retry] [--list_cmd] [-show_cmd cmd-query] [--list_task]
                  [--monitor_resource] [-wait_resource_time wait-time] [--no_check_resource_before_run] -fastq_info FASTQ_INFO [FASTQ_INFO ...] [-r1_name R1_NAME] [-r2_name R2_NAME]
                  [-exclude_samples EXCLUDE_SAMPLES [EXCLUDE_SAMPLES ...]] [-genome_folder GENOME_FOLDER] [--dedup]

Bismark is a program to map bisulfite treated sequencing reads to a genome of interest and perform methylation calls in a single step. The output can be easily imported into a genome
viewer, such as SeqMonk, and enables a researcher to analyse the methylation levels of their samples straight away. It's main features are: * Bisulfite mapping and methylation calling in
one single step * Supports single-end and paired-end read alignments * Supports ungapped, gapped or spliced alignments * Alignment seed length, number of mismatches etc. are adjustable *
Output discriminates between cytosine methylation in CpG, CHG and CHH context

optional arguments:
  -h, --help            show this help message and exit
  -fastq_info FASTQ_INFO [FASTQ_INFO ...]
                        A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json] (default: None)
  -r1_name R1_NAME      python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair
                        brackets will be used as sample name. Example: '(.*).R1.fq.gz' (default: (.*).R1.fastq)
  -r2_name R2_NAME      python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair
                        brackets will be used as sample name. Example: '(.*).R2.fq.gz' (default: (.*).R2.fastq)
  -exclude_samples EXCLUDE_SAMPLES [EXCLUDE_SAMPLES ...]
                        samples to exclude from analysis (default: ())
  -genome_folder GENOME_FOLDER
                        The path to the folder containing the genome to be bisulfite converted. (default: None)
  --dedup               This step is recommended for whole-genome bisulfite samples, but should not be used for reduced representation libraries such as RRBS, amplicon or target
                        enrichment libraries. (default: False)

Arguments for controlling running mode:
  --run                 运行流程，默认不运行, 仅生成流程，如果outdir目录已经存在cmd_state.txt文件，则自动续跑 (default: False)
  --docker              default won't use docker even if docker image is provided. (default: False)
  --plot                generate directed acyclic graph for whole workflow timely (default: False)
  -update_args update-args
                        输入参数配置文件, 其包含流程所有软件需要的参数，该配置文件设置的参数值将是流程中最后实际用的值 (default: None)
  -threads max-workers  允许的最大并行的cmd数目, 默认5 (default: 5)
  -outdir workdir       分析目录或结果目录 (default: /home/danny.gu/PycharmProjects/basefly/bismark/Result)
  -skip step1 [task3 ...]
                        指定要跳过的步骤或具体task,空格分隔,默认程序会自动跳过依赖他们的步骤, 使用--list_cmd or --list_task可查看候选 (default: [])
  --no_skip_depend      当使用skip参数时, 如果同时指定该参数，则不会自动跳过依赖的步骤 (default: False)
  -rerun_steps task3 [task_prefix ...]
                        指定需要重跑的步骤，不论其是否已经成功完成，空格分隔, 这样做的可能原因可以是: 你重新设置了参数. 使用--list_task可查看候选，也可以使用task的前缀指定属于同一个步骤的task (default: [])
  -assume_success_steps task3 [task_prefix ...]
                        假定哪些步骤已经成功运行，不论其是否真的已经成功完成，空格分隔, 这样做的可能原因: 利用之前已经成功运行的结果. 使用--list_task可查看候选，也可以使用task的前缀指定属于同一个步骤的task (default: [])
  -retry max-retry      某步骤运行失败后再尝试运行的次数, 默认1次. 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini (default: 1)
  --list_cmd            仅仅显示当前流程包含的主步骤, 且已经排除指定跳过的步骤 (default: False)
  -show_cmd cmd-query   提供一个cmd名称,输出该cmd的样例 (default: None)
  --list_task           仅仅显示当前流程包含的详细步骤, 且已经排除指定跳过的步骤 (default: False)
  --monitor_resource    是否监控每一步运行时的资源消耗, 如需对某一步设置不同的值, 可在运行流程前修改pipeline.ini (default: False)
  -wait_resource_time wait-time
                        等待资源的时间上限, 默认每次等待时间为15秒, 等待时间超过这个时间且资源不足时判定任务失败 (default: 60)
  --no_check_resource_before_run
                        指示运行某步骤前检测指定的资源是否足够, 如不足, 则该步骤失败; 如果设置该参数, 则运行前不检查资源. 如需对某一步设置不同的值,可运行前修改pipeline.ini. 如需更改指定的资源, 可在运行流程前修改pipeline.ini (default: False)

```
----
# 设计思路如下
## 定义Argument对象，由如下属性描述
1. name: str = '?'
    > 描述参数名称
2. value: Any = None
    > 描述参数取值
3. prefix: str = ''
   > 参数前缀，如‘-i ’， ‘--i’
4. type: Literal['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix'] = 'str'
   > 描述参数类型, type is one of ['str', 'int', 'float', 'bool', 'infile', 'indir', 'fix']
   >> + str: 字符串类型
   >> + int: 整数类型
   >> + float: 整数类型
   >> + bool: bool类型，表示该参数无值，如果设置为True，则命令行中将出现该参数的前缀，否则不出现
   >> + infile：表示参数对应的值是输入文件路径
   >> + indir：表示参数对应的只是输入文件目录
   >> + fix： 表示该参数并不是真正的参数，其为固定的字符串，我们需要给value直接赋值. 例如其可以用来表示管道符如‘| samtools sort’，如此可以让一个command可以串联执行多个程序
5. level: Literal['required', 'optional'] = 'required'
   > 描述参数是否为必需参数
6. default: Any = None
   > 描述参数的默认值, 如果指定level="optional", 则默认值会被忽略
7. range: Any = None
   > 描述参数的取值范围，用列表表示，如range['a', 'b', 'c']
8. array: bool = False
   > 描述参数是否可以赋值多个值，如果True，给value赋值时应该使用列表， 进一步，如果此时还有array=True, 则赋值是嵌套列表
9. delimiter: str = ' '
    > 当一个参数可以输入多个值（即array=True)，形成命令行时，这几个值的分隔符由此指定
10. multi_times: bool = False
    > 指示一个参数是否可以重复多次使用，如果True，给value赋值时应该使用列表, 进一步，如果此时还有array=True, 则赋值是嵌套列表
11. format: str = None
    > 指示输入文件的格式，用后缀表示，仅仅起提示作用
12. order: int = 0
    > order字段是为了参数排序设计，暂无用处，可以用于将来的说明文档
13. desc: str = 'This is description of the argument.'
    > 参数的说明

## 定义Runtime对象，其由如下属性
1. image: str = None
   > docker镜像地址
2. tool_dir: str = ''
   > 命令行的起始字符串，如软件所在目录
3. tool: str = ''
   > 软件名称，命令的行的拼接顺序是 tool_dir + tool + arguments
4. memory: int = 1024
   > 执行环境设置, 指定最少计算资源, memory的单位为字节
5. cpu: int = 2
   > 执行环境设置, 指定最少cpu数量
6. max_memory: int = 0
   > 执行环境设置, 指定最大计算资源, memory的单位为字节
7. max_cpu: int = 0
   > 执行环境设置, 指定最大cpu数量
8. timeout: int = 3600*24*7
   > 最长运行时间，如果运行时间超过该值，则将强制终止该命令

### 定义Output对象，其由如下属性描述
1. value: Any = None
   > 输出值，可以是文件路径或目录路径
2. type: Literal['str', 'int', 'float', 'bool', 'outfile', 'outdir'] = 'outfile'
   > "outfile"表示是输出文件，"outdir"表示是输出目录
3. report: bool = True
   > 设计report参数，用于指定该参数最终是否作为流程的outputs, 暂无用处
4. task_id: str = None
   > 由command形成task之后，可以带入task_id
5. name: str = None
   > 名称
6. desc: str = None
   > 描述信息

## 定义Meta对象，由如下属性描述
> 该对象用于存储Command或者workflow的描述信息
1. name: str = None
   > 描述command或workflow的名称
2. desc: str = 'This is description of the tool/workflow.'
   > command的功能概述
3. author: str = 'unknown'
   > 工具的作者
4. source: str = 'source URL for the tool'
   > 工具的来源，如github地址或工具主页
5. version: str = 'unknown'
   > 记录版本信息

## 定义Command对象，其由如下对象组成
>  Command = {Argument, RunTime, Output, Meta}
* Argument: 全面描述参数，最终可以根据参数列表自动组装成命令行
* Runtime：运行环境描述，如软件，docker image，memory，cpu限制等信息
* Outputs: 描述程序的输出
* Meta: 描述程序的信息，如名称，作者，版本等信息

## 定义Task对象，主要由配备具体参数的Command实例组成:
> Task = {Command, task_id, depends, name}
>> * task_id是自动生成的uuid4对象
>> * depends是列表，记录其依赖的Task或task_id
>> * name是task的名称，最好不要有重复

## 定义TopVar对象
> 用于记录流程的输入的文件或参数信息，这个设计的目的是为了方便流程的转化
1. value: Any
2. name: str = 'notNamed'
3. type: Literal['str', 'int', 'float', 'bool', 'infile', 'indir'] = 'infile'

## 定义workflow对象，主要由Task组成的字典构成
> workflow = {Meta, TopVar, dict(task_id=Task, task_id=Task, ...)}

----

