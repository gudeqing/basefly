## 设计思路
1. 定义argument, runtime, outputs, meta
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
7. 定义方法将workflow转换为第三方格式流程: (1) WDL（暂停开发） (2) argo_workflow（待测试）
