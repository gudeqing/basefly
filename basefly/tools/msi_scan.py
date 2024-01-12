from basefly import Argument, Output, Command, Workflow, TopVar


def msi_scan():
    cmd = Command()
    cmd.meta.name = 'MSIScan'
    cmd.meta.source = 'https://github.com/xjtu-omics/msisensor-pro'
    cmd.meta.function = 'scan microsatellites information'
    cmd.meta.version = '1.2.0'
    cmd.meta.desc = """
    This module scan the reference genome to get microsatellites information.
    example: msisensor-pro scan -d /path/to/reference.fa -o /path/to/reference.list
    """
    cmd.runtime.image = "pengjia1110/msisensor-pro:latest"
    cmd.runtime.tool = "msisensor-pro scan"
    cmd.runtime.memory = 5 * 1024 **3
    cmd.runtime.cpu = 4
    cmd.args['fasta'] = Argument(prefix='-d ', type='infile', desc='reference genome sequences file, *.fasta or *.fa format')
    cmd.args['out'] = Argument(prefix='-o ', type='outstr', desc='output homopolymers and microsatellites file')
    cmd.args['min_homo_size'] = Argument(prefix='-l ', type='int', default=8, desc='minimal homopolymer(repeat unit length = 1) size, default=8')
    cmd.args['max_homo_size'] = Argument(prefix='-m ', type='int', default=50, desc='maximal homopolymer size, default=50')
    cmd.args['context_length'] = Argument(prefix='-c ', type='int', default=5, desc='context length, default=5')
    cmd.args['max_len'] = Argument(prefix='-s ', type='int', default=5, desc='maximal length of microsatellite, default=5')
    cmd.args['min_repeat'] = Argument(prefix='-r ', type='int', default=5, desc='minimal repeat times of microsatellite(repeat unit length>=2), default=5')
    cmd.args['only_homo'] = Argument(prefix='-p ', type='int', level='optional', default=0, desc='output homopolymer only, 0: no; 1: yes, default=0')
    cmd.outputs['out'] = Output(value='{out}', desc='output homopolymers and microsatellites file')
    return cmd


if __name__ == '__main__':
    Workflow().to_cwl_tool(cmd=msi_scan(), write_out=True)
    msi_scan().run_on_terminal()


