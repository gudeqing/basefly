from basefly import Argument, Output, Command, Workflow, TopVar


def fastp():
    cmd = Command()
    cmd.meta.name = 'Fastp'
    cmd.meta.source = 'https://github.com/OpenGene/fastp'
    cmd.meta.version = '0.23.4'
    cmd.meta.function = 'fastq QC, adapter trimming'
    cmd.meta.desc = """
    fastp is a tool used in bioinformatics for the quality control and preprocessing of raw sequence data. 
    fastp is known for its speed and efficiency, and it can process data in parallel, making it suitable for large datasets.
    fastp provides several key functions:
    * It can filter out low-quality reads, which are sequences that have a high probability of containing errors. This is done based on quality scores that are assigned to each base in a read.
    * It can trim adapter sequences, which are artificial sequences added during the preparation of sequencing libraries and are not part of the actual sample's genome.
    * It provides comprehensive quality control reports, including information on sequence quality, GC content, sequence length distribution, and more.
    """
    cmd.runtime.image = 'gudeqing/gatk4.3-bwa-fastp-gencore-mutscan:1.0'
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.runtime.tool = 'fastp'
    cmd.args['read1'] = Argument(prefix='-i ', type='infile', editable=False, desc='read1 fastq file')
    cmd.args['read2'] = Argument(prefix='-I ', type='infile', editable=False, level='optional', desc='read2 fastq file')
    cmd.args['threads'] = Argument(prefix='-w ', default=cmd.runtime.cpu, desc='thread number')
    cmd.args['min_length'] = Argument(prefix='-l ', default=15, desc='reads shorter than length_required will be discarded')
    cmd.args['correction'] = Argument(prefix='--correction', type='bool', default=False, desc='enable base correction in overlapped regions')
    cmd.args['overlap_diff_percent_limit'] = Argument(prefix='--overlap_diff_percent_limit ', default=10, desc='The maximum number of mismatched bases to detect overlapped region of PE reads')
    cmd.args['dedup'] = Argument(prefix='--dedup ', type='bool', default=False, desc='enable deduplication to drop the duplicated reads/pairs')
    cmd.args['trim_front1'] = Argument(prefix='--trim_front1 ', level='optional', desc='trimming how many bases in front for read1')
    cmd.args['trim_tail1'] = Argument(prefix='--trim_tail1 ', level='optional', desc='trimming how many bases in tail for read1')
    cmd.args['fix_mgi_id'] = Argument(prefix='--fix_mgi_id', type='bool', default=False, desc='the MGI FASTQ ID format is not compatible with many BAM operation tools, enable this option to fix it.')
    cmd.args['enable_umi'] = Argument(prefix='--umi', type='bool', default=False, desc='enable unique molecular identifier (UMI) preprocessing')
    cmd.args['umi_loc'] = Argument(prefix='--umi_loc ', level='optional', desc='specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])')
    cmd.args['umi_len'] = Argument(prefix='--umi_len ', type='int', level='optional', desc='if the UMI is in read1/read2, its length should be provided (int [=0])')
    cmd.args['umi_prefix'] = Argument(prefix='--umi_prefix ', level='optional', desc='if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG).')
    cmd.args['umi_skip'] = Argument(prefix='--umi_skip ', level='optional', desc='if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0')
    cmd.args['out1'] = Argument(prefix='-o ', level='optional', type='outstr', desc='clean read1 output fastq file')
    cmd.args['out2'] = Argument(prefix='-O ', level='optional', type='outstr', desc='clean read2 output fastq file')
    cmd.args['html'] = Argument(prefix='--html ', type='outstr', desc='html report file')
    cmd.args['json'] = Argument(prefix='--json ', type='outstr', desc='json report file')
    cmd.outputs['out1'] = Output(value="{out1}", desc='output clean read1 file')
    cmd.outputs['out2'] = Output(value="{out2}", desc='output clean read2 file')
    cmd.outputs['html'] = Output(value="{html}", desc='output html report file')
    cmd.outputs['json'] = Output(value="{json}", desc='output json report file')
    return cmd


if __name__ == '__main__':
    # Workflow().to_cwl_tool(cmd=hisat_genotype(), write_out=True)
    fastp().run_on_terminal()

