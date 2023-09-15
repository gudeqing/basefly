# coding=utf-8
__author__ = 'gudeqing'

import re
import time
import os
import shutil
import configparser
import psutil
import queue
import logging
from subprocess import PIPE
import threading
from threading import Timer, Lock
import weakref
import atexit
import signal

try:
    import pygraphviz as pgv
except Exception as e:
    print('Warn: cannot import pygraphviz and state graph will not be drawn!')
    pgv = None

PROCESS_local = weakref.WeakKeyDictionary()
TIMERS = weakref.WeakKeyDictionary()

# 终止流程时，如何kill容器，可以考虑运行docker时给容器一个特别的名称，然后在流程退出前执行docker kill命令进行清理
Container_Marker = 'bfly'+ str(time.time())

@atexit.register
def _kill_processes_when_exit():
    # register函数位于atexit模块，用于在程序退出时运行，进行必要的清理等
    # print("....Ending....")
    living_processes = list(PROCESS_local.items())
    while living_processes:
        for proc, cmd_name in living_processes:
            if psutil.pid_exists(proc.pid):
                print('Stop running task {} with pid={}:'.format(cmd_name, proc.pid))
                subprocs = list(proc.children(recursive=True))[::-1]
                for ind, subproc in enumerate(subprocs, start=1):
                    if psutil.pid_exists(subproc.pid):
                        print(f'...terminate {ind}th children process with pid={subproc.pid}')
                        subproc.kill()
                proc.kill()
            PROCESS_local.pop(proc)
        living_processes = list(PROCESS_local.items())

    # 清理docker容器
    try:
        proc = psutil.Popen('docker ps --format "{{.ID}}:{{.Names}}"', shell=True, stderr=PIPE, stdout=PIPE)
        stdout, stderr = proc.communicate()
        for line in stdout.split():
            cid, name = line.decode().split(':', 1)
            if name.endswith(Container_Marker):
                print('kill running container:', line)
                os.system(f"docker kill {cid}")
    finally:
        print('finish cleaning up, but you are recommended to check it manually')


def shutdown(signum, frame):
    print('\nYou are Stopping the BaseFly workflow, and we will help to kill the derived processes!')
    exit(0)


# kill signal will be captured
signal.signal(signal.SIGTERM, shutdown)
signal.signal(signal.SIGINT, shutdown)


def set_logger(name='workflow.log', logger_id='x'):
    logger = logging.getLogger(logger_id)
    logger.propagate = False
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(name, mode='w+')
    fh.setLevel(logging.INFO)
    sh = logging.StreamHandler()
    sh.setLevel(logging.WARNING)
    # fmt = '%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
    fmt = '%(asctime)s: %(message)s'
    format_str = logging.Formatter(fmt)  # 设置日志格式
    fh.setFormatter(format_str)
    logger.addHandler(sh)
    logger.addHandler(fh)
    return logger


class Command(object):
    def __init__(self, cmd, name, timeout=3600*24*7, outdir=os.getcwd(), image=None, docker_cmd_prefix=None, mount_vols=None, wkdir=None,
                 monitor_resource=True, monitor_time_step=2, logger=None, **kwargs):
        self.name = name
        self.cmd = cmd
        self.image = image
        self.docker_cmd_prefix = docker_cmd_prefix or 'docker run --rm --privileged -i --entrypoint /bin/bash'
        self.mount_vols = mount_vols
        self.proc = None
        self.stdout = None
        self.stderr = None
        self.timeout = int(timeout)
        self.used_time = 0
        self.max_used_mem = 0
        self.max_used_cpu = 0
        self.threads_num = 0
        self.monitor = monitor_resource
        self.monitor_time_step = int(monitor_time_step)
        self.outdir = outdir
        self.wkdir = wkdir
        self.times = 0
        if not logger:
            self.logger = set_logger(name=os.path.join(self.outdir, 'command.log'))
        else:
            self.logger = logger

    def _monitor_resource(self):
        while psutil.pid_exists(self.proc.pid) and self.proc.is_running():
            try:
                self.max_used_mem = self.proc.memory_full_info().vms
                self.max_used_cpu = self.proc.cpu_percent(interval=0.5)
                for subproc in self.proc.children(recursive=True):
                    if psutil.pid_exists(subproc.pid):
                        # 获取进程占用的memory信息
                        memory = subproc.memory_full_info().vms
                        # 获取cpu信息
                        used_cpu = subproc.cpu_percent(interval=1)
                        # print(self.name, subproc.pid, memory, used_cpu, subproc.num_threads())
                        if memory > self.max_used_mem:
                            self.max_used_mem = memory
                        if used_cpu > self.max_used_cpu:
                            self.max_used_cpu = used_cpu
                # 获取主进程的线程数量
                self.threads_num = self.proc.num_threads()
                time.sleep(self.monitor_time_step)
            except Exception as e:
                # print('Failed to capture cpu/mem info for: ', e)
                break
            finally:
                # print(self.proc.pid, self.max_used_cpu, self.max_used_mem, self.threads_num)
                pass

    def run(self):
        if not self.wkdir:
            cmd_wkdir = os.path.join(self.outdir, self.name)
        else:
            cmd_wkdir = self.wkdir
        if os.path.exists(cmd_wkdir):
            try:
                print(f'Remove or Rename existed workdir {cmd_wkdir} to prevent potential error')
                try:
                    shutil.rmtree(cmd_wkdir)
                except:
                    os.rename(cmd_wkdir, cmd_wkdir+'_old')
            except Exception as e:
                print(f'Failed to remove/rename {cmd_wkdir}: {e}')

        os.makedirs(cmd_wkdir, exist_ok=True)
        with open(os.path.join(cmd_wkdir, 'cmd.sh'), 'w') as f:
            f.write('set -o pipefail\n')
            f.write(self.cmd + '\n')

        if self.image:
            docker_cmd = self.docker_cmd_prefix
            docker_cmd += f' --name {self.name}-{Container_Marker}'
            for each in self.mount_vols.split(';'):
                docker_cmd += f' -v {each}:{each} '
            docker_cmd += f'-w {cmd_wkdir} {self.image} cmd.sh'
            with open(os.path.join(cmd_wkdir, 'docker.cmd.sh'), 'w') as f:
                f.write((docker_cmd+'\n'))

        start_time = time.time()
        self.logger.warning("Step: {}".format(self.name))

        # submit task
        if self.image:
            self.logger.info("> {}".format(docker_cmd))
            self.proc = psutil.Popen(docker_cmd, shell=True, stderr=PIPE, stdout=PIPE, cwd=cmd_wkdir)
        else:
            self.logger.info("> {}".format(self.cmd))
            self.proc = psutil.Popen(self.cmd, shell=True, stderr=PIPE, stdout=PIPE, cwd=cmd_wkdir, executable='/bin/bash')
            # add ‘exec’ will cause cmd to inherit the shell process, instead of having the shell launch a child process
            # self.proc = psutil.Popen("exec " + self.cmd, shell=True, stderr=PIPE, stdout=PIPE, cwd=cmd_wkdir)

        # tracing process
        PROCESS_local[self.proc] = self.name

        if self.monitor:
            thread = threading.Thread(target=self._monitor_resource, daemon=True)
            thread.start()

        # timer = Timer(self.timeout, self.proc.kill)
        # TIMERS[timer] = self.name
        # try:
        #     timer.start()
        #     self.stdout, self.stderr = self.proc.communicate()
        #     if self.monitor:
        #         thread.join()
        # finally:
        #     timer.cancel()
        # 使用自带参数timeout计时
        try:
            self.stdout, self.stderr = self.proc.communicate(timeout=self.timeout)
        except psutil.TimeoutExpired:
            self.proc.kill()
            self.stdout, self.stderr = self.proc.communicate()

        self._write_log()

        if self.image:
            # 无论是否运行成功，修改结果文件权限
            docker_cmd = "docker run --rm --privileged -i"
            # 由于有些镜像会导致出错如”/bin/chown: cannot execute binary file“，所以采用固定镜像来执行
            docker_cmd += f' -v {cmd_wkdir}:{cmd_wkdir} bash:latest '
            docker_cmd += f'chown -R {os.getuid()}:{os.getgid()} {cmd_wkdir}'
            os.system(docker_cmd)

        end_time = time.time()
        self.used_time = round(end_time - start_time, 2)

    def _write_log(self):
        log_dir = os.path.join(self.outdir, 'logs')
        if not os.path.exists(log_dir):
            try:
                os.mkdir(log_dir)
            except FileExistsError:
                pass

        prefix = os.path.join(self.outdir, 'logs', self.name+'.'+str(self.proc.pid))
        if self.stderr:
            with open(prefix+'.stderr.txt', 'wb') as f:
                f.write(self.stderr)
        if self.stdout:
            with open(prefix+'.stdout.txt', 'wb') as f:
                f.write(self.stdout)
        if self.monitor:
            with open(prefix+'.resource.txt', 'w') as f:
                f.write('max_cpu (cpu_percent*0.01): {:.2f}\n'.format(self.max_used_cpu*0.01))
                f.write('max_mem (Virtual Memory Size; M): {:.2f}\n'.format(self.max_used_mem/1024**2))
                f.write('thread_num (num_threads): {}\n'.format(self.threads_num))


class CommandNetwork(object):
    def __init__(self, cmd_config):
        # self.parser = configparser.ConfigParser()
        self.parser = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
        self.parser.read(cmd_config, encoding='utf-8')
        self.pool_size = self.parser.getint('mode', 'threads')
        self.outdir = os.path.abspath(self.parser.get('mode', 'outdir'))

    def names(self):
        sections = self.parser.sections()
        # mode section is not cmd
        sections.pop(sections.index('mode'))
        return sections

    def orphans(self):
        independent_cmds = list()
        for name in self.names():
            if 'depend' not in self.parser[name]:
                independent_cmds.append(name)
            else:
                depend = self.parser[name]['depend'].strip()
                if not depend:
                    independent_cmds.append(name)
                else:
                    for each in depend.split(','):
                        if each not in self.names():
                            raise Exception(f'Step "{each}" is not in your pipeline! A spelling mistake?')
        return independent_cmds

    def get_dependency(self, name):
        if 'depend' not in self.parser[name]:
            return []
        else:
            depend = self.parser[name]['depend'].strip()
            if not depend:
                return []
            else:
                return [x.strip() for x in depend.split(',')]

    def get_cmd_description_dict(self, name):
        tmp_dict = dict(self.parser[name])
        tmp_dict['name'] = name
        if 'cpu' not in tmp_dict:
            tmp_dict['cpu'] = 0
        if 'mem' not in tmp_dict:
            tmp_dict['mem'] = 0
        if 'depend' not in tmp_dict:
            tmp_dict['depend'] = None
        if 'retry' not in tmp_dict:
            tmp_dict['retry'] = self.parser.getint('mode', 'retry')
        else:
            tmp_dict['retry'] = self.parser.getint(name, 'retry')
        if 'monitor_resource' not in tmp_dict:
            tmp_dict['monitor_resource'] = self.parser.getboolean('mode', 'monitor_resource')
        else:
            tmp_dict['monitor_resource'] = self.parser.getboolean(name, 'monitor_resource')
        if 'timeout' not in tmp_dict:
            tmp_dict['timeout'] = 3600*24*1
        else:
            tmp_dict['timeout'] = self.parser.getint(name, 'timeout')
        if 'monitor_time_step' not in tmp_dict:
            tmp_dict['monitor_time_step'] = self.parser.getint('mode', 'monitor_time_step')
        else:
            tmp_dict['monitor_time_step'] = self.parser.getint(name, 'monitor_time_step')
        if 'check_resource_before_run' not in tmp_dict:
            tmp_dict['check_resource_before_run'] = self.parser.getboolean('mode', 'check_resource_before_run')
        else:
            tmp_dict['check_resource_before_run'] = self.parser.getboolean(name, 'check_resource_before_run')
        if 'wkdir' not in tmp_dict:
            tmp_dict['wkdir'] = os.path.join(self.outdir, name)
        else:
            tmp_dict['wkdir'] = self.parser.get(name, 'wkdir')
        return tmp_dict


class CheckResource(object):
    @staticmethod
    def available_mem():
        # memory
        return psutil.virtual_memory().available*0.98

    @staticmethod
    def available_cpu():
        total = psutil.cpu_count()
        return int(total - total*psutil.cpu_percent()*0.01)

    def is_enough(self, cpu, mem, timeout=10):
        start_time = time.time()
        enough_num = 0
        while True:
            if float(cpu) <= self.available_cpu() \
                    and float(mem) <= self.available_mem():
                enough_num += 1
                if enough_num >= 3:
                    return True
                if enough_num >= 1 and timeout <= 10:
                    return True
            if time.time() - start_time >= timeout:
                return False
            time.sleep(3)


class StateGraph(object):
    def __init__(self, state):
        if type(state) != dict:
            state_dict = dict()
            with open(state) as f:
                header = f.readline().strip('\n').split('\t')
                for line in f:
                    lst = line.strip('\n').split('\t')
                    state_dict[lst[0]] = dict(zip(header[1:], lst[1:]))
            state = state_dict
        self.state = state
        self.graph = pgv.AGraph(directed=True, rankdir='LR')
        self.used_colors = dict()
        self.color_dict = dict(
            success='#7FFF00',
            failed='#FFD700',
            running='#9F79EE',
            queueing='#87CEFF',
            killed='red',
            outdoor='#A8A8A8',
        )

    def _add_nodes(self):
        for node, cmd_info in self.state.items():
            status = cmd_info['state']
            node_detail = node.split('_', 1)
            if status in self.color_dict:
                color = self.color_dict[status]
            else:
                color = '#A8A8A8'
            self.used_colors[status] = color
            used_time = cmd_info['used_time']
            if isinstance(used_time, str):
                if used_time == 'unknown':
                    pass
                else:
                    try:
                        float(used_time)
                        node_detail.append(used_time+'s')
                    except ValueError:
                        node_detail.append(used_time)
            elif float(used_time) <= 0:
                pass
            else:
                node_detail.append(str(used_time) + 's')
            self.graph.add_node(
                node,
                # 谷歌浏览器可以正常显示tooltip
                # tooltip=cmd_info['cmd'].replace(' ', '\n').replace('\\', ''),
                tooltip=cmd_info['cmd'].replace('\\', ''),
                shape="box",
                style="rounded, filled",
                fillcolor=color,
                color="mediumseagreen",
                label='\n'.join(node_detail)
            )

    def _add_edges(self):
        for target in self.state:
            sources = self.state[target]['depend'].strip()
            if sources:
                sources = sources.split(',')
                edges = zip(sources, [target]*len(sources))
                if self.state[target]['state'] == 'success':
                    color = 'green'
                elif self.state[target]['state'] == 'running':
                    color = '#836FFF'
                else:
                    color = '#4D4D4D'
                self.graph.add_edges_from(edges, color=color)
            else:
                self.graph.add_edge('Input', target, color='green')

    def _add_legend(self):
        subgraph = self.graph.add_subgraph(name='cluster_sub', label='Color Legend', rank='max')
        subgraph.graph_attr['color'] = 'lightgrey'
        subgraph.graph_attr['style'] = 'filled'
        subgraph.graph_attr['ratio'] = 'compress'
        for node, color in self.used_colors.items():
            subgraph.add_node(
                node,
                shape="note",
                style="filled",
                fillcolor=color,
                color="mediumseagreen",
            )
        nodes = list(self.used_colors.keys())
        for ind in range(len(nodes)):
            if ind <= len(nodes) - 2:
                subgraph.add_edge(nodes[ind], nodes[ind+1], style='invis')

    def draw(self, img_file='state.svg'):
        self._add_nodes()
        self._add_edges()
        self._add_legend()
        img_fmt = os.path.splitext(img_file)[1][1:]
        self.graph.draw(path=img_file, format=img_fmt, prog='dot')


class RunCommands(CommandNetwork):
    __LOCK__ = Lock()

    def __init__(self, cmd_config, timeout=10, logger=None, draw_state_graph=True):
        super().__init__(cmd_config)
        self.end = False
        self.ever_queued = set()
        self.queue = self.__init_queue()
        self.state = self.__init_state()
        self.task_number = len(self.state)
        self.success = 0
        self.failed = 0
        self.is_continue = False
        # wait resource time limit
        self.timeout = timeout
        if not logger:
            os.makedirs(self.outdir, exist_ok=True)
            files = os.listdir(self.outdir)
            order = 1
            for each in files:
                if re.fullmatch(r'wf.\d.running.*.log', each):
                    if os.path.getsize(os.path.join(self.outdir, each)) > 1:
                        order += 1
            self.logger = set_logger(name=os.path.join(self.outdir, f'wf.{order}.running.{time.time()}.log'))
        else:
            self.logger = logger
        # draw state graph
        self.draw_state_graph = draw_state_graph if pgv else False
        self._draw_state()

    def __init_queue(self):
        cmd_pool = queue.Queue()
        for each in self.orphans():
            cmd_pool.put(each)
            self.ever_queued.add(each)
        return cmd_pool

    def __init_state(self):
        state_dict = dict()
        for name in self.names():
            state_dict[name] = dict()
            fields = ['state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
            for each in fields:
                state_dict[name][each] = 'unknown'
            state_dict[name]['cmd'] = self.parser[name]['cmd']
            state_dict[name]['depend'] = ','.join(self.get_dependency(name))
            state_dict[name]['times'] = 0
        return state_dict

    def _update_queue(self):
        success = set(x for x in self.state if self.state[x]['state'] == 'success')
        failed = set(x for x in self.state if self.state[x]['state'] == 'failed')
        waiting = set(self.names()) - self.ever_queued
        if not waiting:
            self.queue.put(None)
            return
        for each in waiting:
            dependency = set(self.get_dependency(each))
            if dependency & failed:
                self.ever_queued.add(each)
                self.state[each]['state'] = 'failed'
                self.state[each]['used_time'] = 'FailedDependencies'
                self.logger.warning(each + ' cannot be started for some failed dependencies!')
            if not (dependency - success):
                self.ever_queued.add(each)
                self.queue.put(each, block=True)
        # 按照self.names排序，也即原始流程中的顺序排序
        current_queue = list(self.queue.queue)
        reorder_queue = queue.Queue()
        for each in self.names():
            if each in current_queue:
                reorder_queue.put(each, block=True)
        self.queue = reorder_queue

    def _update_state(self, cmd=None, killed=False):
        if cmd is not None:
            cmd_state = self.state[cmd.name]
            if cmd.proc is None:
                cmd_state['state'] = 'failed'
                cmd_state['used_time'] = 'NotEnoughResource'
                cmd_state['times'] += 1
                self.logger.warning(cmd.name + ' cannot be started for not enough resource!')
            else:
                cmd_state['state'] = 'success' if cmd.proc.returncode == 0 else 'failed'
                cmd_state['used_time'] = cmd.used_time
                cmd_state['mem'] = cmd.max_used_mem
                cmd_state['cpu'] = cmd.max_used_cpu
                cmd_state['pid'] = cmd.proc.pid
                cmd_state['times'] += 1
        success = set(x for x in self.state if self.state[x]['state'] == 'success')
        self.success = len(success)
        failed = set(x for x in self.state if self.state[x]['state'] == 'failed')
        self.failed = len(failed)
        running_or_queueing = self.ever_queued - success - failed
        waiting = set(self.names()) - self.ever_queued
        tmp_dict = {y: x for x, y in PROCESS_local.items()}
        for each in running_or_queueing:
            try:
                if each in tmp_dict:
                    self.state[each]['pid'] = tmp_dict[each].pid
                    if tmp_dict[each].is_running():
                        if killed:
                            self.state[each]['state'] = 'killed'
                        else:
                            self.state[each]['state'] = 'running'
                else:
                    self.state[each]['state'] = 'queueing'
            except Exception as e:
                pass
        for each in waiting:
            self.state[each]['state'] = 'outdoor'

    def _write_state(self):
        outfile = os.path.join(self.outdir, 'cmd_state.txt')
        if time.localtime().tm_min % 5 == 0:
            back_file = os.path.join(self.outdir, 'bak.cmd_state.txt')
            if os.path.exists(outfile):
                os.rename(outfile, back_file)
        with open(outfile, 'w') as f:
            fields = ['name', 'state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
            f.write('\t'.join(fields)+'\n')
            for name in self.state:
                content = '\t'.join([str(self.state[name][x]) for x in fields[1:]])
                f.write(name+'\t'+content+'\n')

    def _draw_state(self):
        if self.draw_state_graph:
            outfile = os.path.join(self.outdir, 'state.svg')
            back_file = os.path.join(self.outdir, 'bak.state.svg')
            if os.path.exists(outfile):
                os.rename(outfile, back_file)
            StateGraph(self.state).draw(outfile)

    def _update_status_when_exit(self):
        # print('final update status')
        self._update_state(killed=True)
        self._write_state()
        self._draw_state()

    def single_run(self):
        while True:
            if self.queue.empty():
                time.sleep(3)
                with self.__LOCK__:
                    self._update_queue()
                    self._write_state()
                    self._draw_state()
                continue
            name = self.queue.get(block=True)
            if name is None:
                self.queue.put(None)
                self.end = True
                break
            tmp_dict = self.get_cmd_description_dict(name)
            if 'outdir' in tmp_dict:
                tmp_dict.pop('outdir')
            if 'logger' in tmp_dict:
                tmp_dict.pop('logger')

            if self.state[name]['times'] <= int(tmp_dict['retry']):
                cmd = Command(**tmp_dict, outdir=self.outdir, logger=self.logger)
                enough = True
                success = False
                if tmp_dict['check_resource_before_run']:
                    if not CheckResource().is_enough(tmp_dict['cpu'], tmp_dict['mem'], self.timeout):
                        enough = False
                        self.logger.warning(f"No resource to start task {name} at {self.state[name]['times']}'th time!")
                        if self.state[name]['times'] < int(tmp_dict['retry']):
                            # 把任务在再放回任务队列
                            self.queue.put(name, block=True)
                if enough:
                    if self.state[name]['times'] > 0:
                        self.logger.info('It is {}th time to retry {}'.format(self.state[name]['times'], cmd.name))
                    self.state[cmd.name]['state'] = 'running'
                    with self.__LOCK__:
                        self._draw_state()
                    cmd.run()
                    if cmd.proc.returncode == 0:
                        success = True
                    else:
                        if self.state[name]['times'] < int(tmp_dict['retry']):
                            # 如果不是最后一次执行任务，失败的任务再放回任务队列
                            self.queue.put(name, block=True)
                # 最后一次尝试运行，且执行失败，则判定任务执行失败
                if (not success) and self.state[name]['times'] == int(tmp_dict['retry']):
                    self.logger.warning(f'Failed to complete task {name}!')
                    if cmd.stderr:
                        self.logger.warning(cmd.stderr)
                # 只有本次执行任务成功，或者最是后一次尝试执行任务时，才去更新状态
                if success or self.state[name]['times'] == int(tmp_dict['retry']):
                    if not enough:
                        self.logger.warning('Local resource is Not enough for {}!'.format(cmd.name))
                    with self.__LOCK__:
                        self._update_state(cmd)
                        self._update_queue()
                        self._write_state()
                        self._draw_state()
                else:
                    # 中间失败时递增执行次数
                    with self.__LOCK__:
                        self.state[name]['times'] += 1

    def parallel_run(self, assume_success_steps=tuple()):
        atexit.register(self._update_status_when_exit)
        if assume_success_steps:
            for each in assume_success_steps:
                self.ever_queued.add(each)
                self.state[each]['state'] = 'success'
                self.state[each]['used_time'] = 'unknown'
            self.queue = queue.Queue()
            self._update_queue()
        start_time = time.time()
        pool_size = self.parser.getint('mode', 'threads')
        threads = list()
        for _ in range(pool_size):
            thread = threading.Thread(target=self.single_run, daemon=True)
            threads.append(thread)
            thread.start()
            # 每隔3秒提交一个任务，一定程度可以避免同时提交多个消耗资源比较大的任务。
            time.sleep(3)

        # update state
        time.sleep(2)
        with self.__LOCK__:
            self._update_state()
            self._write_state()
            self._draw_state()
        # join threads
        _ = [x.join() for x in threads]
        percent = f'{self.success/self.task_number:.2%}'
        failed = self.task_number - self.success
        self.logger.warning("Total time used for current running: {}s".format(time.time() - start_time))
        total_used_time = sum(float(self.state[x]['used_time']) for x in self.state if str(self.state[x]['used_time']).replace('.', '').isnumeric())
        self.logger.warning(f'Finished {percent}: Success={self.success}, Failed={failed}, Total={self.task_number}')
        self.logger.warning("Total Time Used for Success Tasks Without Parallel Running: {:.2f} minutes".format(total_used_time/60))
        return self.success, len(self.state)

    def continue_run(self, rerun_steps=tuple(), assume_success_steps=tuple()):
        self.ever_queued = set()
        # 使用已有状态信息更新状态
        existed_state_file = os.path.join(self.outdir, 'cmd_state.txt')
        if not os.path.exists(existed_state_file):
            raise Exception('We found no cmd_state.txt file in {}!'.format(self.outdir))
        with open(existed_state_file, 'r') as f:
            _ = f.readline()
            header = ['name', 'state', 'used_time', 'mem', 'cpu', 'pid', 'depend', 'cmd']
            for line in f:
                line_lst = line.strip().split('\t')
                task_info = dict(zip(header, line_lst))
                task_name = task_info['name']

                if task_name not in self.state:
                    self.logger.warning(line_lst[0] + ' was skipped for it is no longer in the modified pipeline.ini')
                    continue

                if task_name in rerun_steps:
                    # 这里不会用旧信息更新任务信息，因此重跑的任务使用的是最新的命令
                    task_info['state'] = 'unknown'
                    self.state[task_name]['state'] = 'unknown'

                if task_name in assume_success_steps:
                    task_info['state'] = 'success'

                if task_info['state'] == 'success':
                    self.ever_queued.add(task_name)
                    if task_info['cmd'] != self.state[task_name]['cmd']:
                        # 命令行的更改可能只是字符串层面的更改，而没有实质的任务性质更改
                        print(f'we noticed that command of {task_name} changed, but we will not rerun it !')
                        print('Old cmd:', task_info['cmd'])
                        print('New cmd:', self.state[task_name]['cmd'])
                    # 对于被已经判定成功的task，使用旧的信息进行更新
                    task_info.pop("name")
                    self.state[task_name].update(task_info)

        failed = set(self.names()) - self.ever_queued
        if failed:
            self.is_continue = True
            self.logger.warning('Continue to run {} tasks: {}'.format(len(failed), failed))
            self.queue = queue.Queue()
            self._update_queue()
            self._draw_state()
            self.parallel_run()
        else:
            self.logger.warning('Nothing to continue run since all steps are in success status')


def draw_state(cmd_state, out='state.svg'):
    StateGraph(cmd_state).draw(img_file=out)


def run_wf(wf, plot=False, timeout=300, rerun_steps:tuple=None, assume_success_steps:tuple=None):
    """
    :param wf: pipeline configuration file
    :param plot: if set, running state will be visualized if pygraphviz installed
    :param timeout: time to wait for enough resource to initiate a task, default 300
    :param rerun_steps: tell which finished tasks need to be rerun. By default the runner will start from failed steps.
    :return:
    """
    workflow = RunCommands(wf, timeout=timeout, draw_state_graph=plot)
    assume_success_tasks = []
    for each in assume_success_steps:
        if each in workflow.state.keys():
            assume_success_tasks.append(each)
        else:
            if each.endswith("*"):
                assume_success_tasks += [x for x in workflow.state.keys() if x.startswith(each[:-1])]
            else:
                print(f'{each} is not a valid task name and we will skip this one !')

    rerun_tasks = []
    for each in rerun_steps:
        if each in workflow.state.keys():
            rerun_tasks.append(each)
        else:
            if each.endswith("*"):
                rerun_tasks += [x for x in workflow.state.keys() if x.startswith(each[:-1])]
            else:
                print(f'{each} is not a valid task name and we will skip this one !')

    state = os.path.join(workflow.outdir, 'cmd_state.txt')
    state_bak = os.path.join(workflow.outdir, 'cmd_state_history.txt')
    if os.path.exists(state):
        shutil.copyfile(state, state_bak)
        workflow.continue_run(rerun_steps=tuple(rerun_tasks), assume_success_steps=tuple(assume_success_tasks))
    else:
        workflow.parallel_run(assume_success_steps=tuple(assume_success_tasks))
    return workflow


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['run_wf', 'draw_state'])

