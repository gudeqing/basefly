from basefly import Command, Argument, Output


def VardictPaired():
    cmd = Command()
    cmd.meta.name = 'VardictPaired'
    cmd.meta.source = 'https://github.com/AstraZeneca-NGS/VarDictJava'
    cmd.meta.version = '1.8.3'
    cmd.meta.function = 'somatic variant calling in paired mode using vardict'
    cmd.meta.desc = "VarDictJava is a variant discovery program written in Java and Perl."
    cmd.runtime.image = 'hydragenetics/vardict:1.8.3'
    cmd.runtime.memory = 16 * 1024 ** 3
    cmd.runtime.cpu = 8
    cmd.runtime.tool = 'java -Xmx16g -jar /usr/local/lib/VarDict-1.8.3.jar'
    cmd.args['sample'] = Argument(prefix='-N ', desc='sample name')
    cmd.args['bam'] = Argument(prefix='-b ', type='infile', array=True, delimiter='\\|', format='bam', desc='The indexed BAM files, tumor|normal')
    cmd.args['genome'] = Argument(prefix='-G ', type='infile', format='fasta', desc='The reference fasta. Should be indexed (.fai).')
    cmd.args['threads'] = Argument(prefix='-th ', default=cmd.runtime.cpu, desc='Threads count.')
    cmd.args['min-freq'] = Argument(prefix='-f ', default="0.0001", desc='The threshold for allele frequency')
    cmd.args['chromosome'] = Argument(prefix='-c ', default=1, desc='The column of chromosome')
    cmd.args['region_start'] = Argument(prefix='-S ', default=2, desc='The column of region start')
    cmd.args['region_end'] = Argument(prefix='-E ', default=3, desc='The column of region end')
    cmd.args['gene'] = Argument(prefix='-g ', default=4, desc='The column of gene name')
    cmd.args['fisher'] = Argument(prefix='--fisher', type='bool', default=True, desc='Fisher exact test.')
    cmd.args['mfreq'] = Argument(prefix='-mfreq ', default=0.25, desc="The variant frequency threshold to determine variant as good in case of monomer MSI. Default: 0.25")
    cmd.args['nmfreq'] = Argument(prefix='-nmfreq ', default=0.1, desc='The variant frequency threshold to determine variant as good in case of non-monomer MSI')
    cmd.args['read_position_filter'] = Argument(prefix='-P ', default=5, desc="If the mean variants position is less that specified, it's considered false")
    cmd.args['min-reads'] = Argument(prefix='-r ', default=2, desc='The minimum number of variant reads')
    cmd.args['nosv'] = Argument(prefix='--nosv', type='bool', default=True,)
    cmd.args['count_overlap_by_first_read'] = Argument(prefix='-UN', type='bool', default=False, desc='Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using first read only')
    cmd.args['count_overlap_by_forward_read'] = Argument(prefix='-u', type='bool', default=True, desc='Indicate unique mode, which when mate pairs overlap, the overlapping part will be counted only once using foward read only')
    cmd.args['bed'] = Argument(prefix='', type='infile', desc='region or bed file')
    cmd.args['_fix'] = Argument(type='fix', value='| var2vcf_paired.pl -A -p 5 -q 22.5 -d 5 -v 2 -f 0.00001 ', desc='pipe to another script')
    cmd.args['names'] = Argument(prefix='-N ', array=True, delimiter='\\|', desc='The sample name(s).  If only one name is given, the matched will be simply names as "name-match".')
    cmd.args['output'] = Argument(prefix='> ', type='outstr', desc='output vcf name')
    cmd.outputs['out'] = Output(value='{output}', format='bam')
    return cmd


if __name__ == '__main__':
    VardictPaired().run_on_terminal()
