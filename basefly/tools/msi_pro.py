from basefly import Command, Argument, Output


def msi_pro():
    cmd = Command()
    cmd.meta.name = 'MSIPro'
    cmd.meta.source = 'https://github.com/xjtu-omics/msisensor-pro'
    cmd.meta.function = 'MSI calculation for tumor only sample'
    cmd.meta.version = '1.2.0'
    cmd.meta.desc = """
    This module evaluate MSI using tumor only sample
    Example:
    1. msisensor-pro pro -d /path/to/reference.list -i 0.1 -t /path/to/case1_tumor_sorted.bam -o /path/to/case1_output
    2. msisensor-pro pro -d /path/to/reference.list_baseline -t /path/to/case1_tumor_sorted.bam -o /path/to/case1_output
    """
    cmd.runtime.image = "pengjia1110/msisensor-pro:latest"
    cmd.runtime.tool = "msisensor-pro pro"
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.args['sites'] = Argument(prefix='-d ', type='infile', desc='homopolymer and microsatellite file')
    cmd.args['tumor'] = Argument(prefix='-t ', type='infile', desc='tumor bam file with index')
    cmd.args['out_prefix'] = Argument(prefix='-o ', type='outstr', desc='output prefix')
    cmd.args['bed'] = Argument(prefix='-e ', type='infile', level='optional', desc='bed file, optional')
    cmd.args['threshold'] = Argument(prefix='-i ', type='float', level='optional', default=0.1, desc='minimal threshold for instable sites detection (just for tumor only data), default=0.1')
    cmd.args['coverage'] = Argument(prefix='-c ', type='int', default=15, desc='coverage threshold for msi analysis, WXS: 20; WGS: 15')
    cmd.args['region'] = Argument(prefix='-r', level='optional', desc='choose one region, format example: 1:10000000-20000000')
    cmd.args['min_homo_size'] = Argument(prefix='-p ', type='int', default=8, desc='minimal homopolymer size for distribution analysis, default=8')
    cmd.args['max_homo_size'] = Argument(prefix='-m ', type='int', default=50, desc='maximal homopolymer size for distribution analysis, default=50')
    cmd.args['min_ms_size'] = Argument(prefix='-s ', type='int', default=5, desc='minimal microsatellite size for distribution analysis, default=5')
    cmd.args['max_ms_size'] = Argument(prefix='-w ', type='int', default=40, desc='maximal microsatellite size for distribution analysis, default=40')
    cmd.args['span_size'] = Argument(prefix='-u ', type='int', default=500, desc='span size around window for extracting reads, default=500')
    cmd.args['threads'] = Argument(prefix='-b ', type='int', default=4, desc='threads number for parallel computing, default=2')
    cmd.args['only_ms'] = Argument(prefix='-y ', type='int', level='optional', default=0, desc='output microsatellite only, 0: no; 1: yes, default=0')
    cmd.args['only_homo'] = Argument(prefix='-x ', type='int', level='optional', default=0, desc='output homopolymer only, 0: no; 1: yes, default=0')
    cmd.args['out_zero_cov'] = Argument(prefix='-0 ', type='int', level='optional', default=0, desc='output site have no read coverage, 1: no; 0: yes, default=0')
    cmd.outputs['out_all'] = Output(value='{out_prefix}_all', type='outfile', desc='all MSI sites')
    cmd.outputs['out_unstable'] = Output(value='{out_prefix}_unstable', type='outfile', desc='unstable MSI sites')
    cmd.outputs['out_distribution'] = Output(value='{out_prefix}_dis', type='outfile', desc='MS distribution')
    cmd.outputs['out'] = Output(value='{out_prefix}', type='outfile', desc='MSI statistic summary')
    return cmd


if __name__ == '__main__':
    msi_pro().run_on_terminal()
