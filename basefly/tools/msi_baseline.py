from basefly import Command, Argument, Output


def msi_baseline():
    cmd = Command()
    cmd.meta.name = 'MSIScan'
    cmd.meta.source = 'https://github.com/xjtu-omics/msisensor-pro'
    cmd.meta.function = 'create msi baseline'
    cmd.meta.version = '1.2.0'
    cmd.meta.desc = """
    This module build baseline for MSI detection with pro module using only tumor sequencing data. 
    To achieve it, you need sequencing data from normal samples(-i)
    Example:
    msisensor-pro baseline -d /path/to/reference.list -i /path/to/configure.txt -o /path/to/baseline/directory 
    """
    cmd.runtime.image = "pengjia1110/msisensor-pro:latest"
    cmd.runtime.tool = "msisensor-pro baseline"
    cmd.runtime.memory = 5 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.args['sites'] = Argument(prefix='-d ', type='infile', desc='homopolymer and microsatellite file')
    cmd.args['bam_files'] = Argument(prefix='-i ', type='infile', desc='configure files for building baseline (text file), first column is case name, second column is bam file path')
    cmd.args['outdir'] = Argument(prefix='-o ', type='outstr', desc='output directory')
    cmd.args['min_coverage'] = Argument(prefix='-c ', type='int', default=15, desc='coverage threshold for msi analysis, WXS: 20; WGS: 15')
    cmd.args['min_fraction'] = Argument(prefix='-l ', type='float', default=0.5, desc='a site with a ratio of deteced in all samples less than this parameter will be removed in following analysis')
    cmd.args['threads'] = Argument(prefix='-b ', type='int', default=1, desc='threads number for parallel computing')
    cmd.args['only_ms'] = Argument(prefix='-y ', type='int', level='optional', default=0, desc='output microsatellite only, 0: no; 1: yes, default=0')
    cmd.args['only_homo'] = Argument(prefix='-x ', type='int', level='optional', default=0, desc='output homopolymer only, 0: no; 1: yes, default=0')
    cmd.outputs['out'] = Output(value='{outdir}', type='outdir', desc='output directory')
    return cmd


if __name__ == '__main__':
    msi_baseline().run_on_terminal()



