from basefly import Command, Argument, Output

"""
https://www.sciencedirect.com/science/article/pii/S2001037021003688
MSI classification using the various available algorithms has been accomplished in different types of cancer using WGS, WES, RNAseq, or targeted sequencing data [81].
The performance of some algorithms has been inter-compared in some studies. 
For example, the accuracy of MSIsensor, MANTIS and mSINGS using their default setting were evaluated in COAD/READ, UCEC and STAD cohorts downloaded from Cancer Genomics Hub.
MSIsensor and MANTIS displayed equally good sensitivity, which was clearly better than mSINGS in all cohorts [69]. 
In another study from Jia et al. MSIsensor and MANTIS outperformed mSINGS but has similar accuracy as MSIsensor-pro using 1532 TCGA normal-tumor paired whole-exome sequencing data from 3 cancer types [71]. 
So far, no guidance is available regarding the choice of algorithm for MSI detection. 
The decision depends largely on type of sequencing data, availability of paired normal sample, type of cancer as well as type of DNA specimen. 
MSI detection has been achieved on NGS data obtained from fresh-frozen tissue [72], [82], [83], FFPE tissue [72], [82], [83], [84] and liquid biopsy samples [17], [29], [63], [65], [85]. 
Most algorithms were designed based on NGS data obtained from tissue samples which usually contain high tumor purity.
MSIsensor, MANTIS, mSINGS, MSIsensor-pro and MSI NGS calling algorithms show reduced power for calling MSI at tumor purities < 10% [71],
thereby restricting their application in liquid biopsy where mutant DNA can be masked by presence of excessive wild-type alleles. 
In contrast, some of the algorithms aim to detect samples bearing low mutation levels and can push the limit of MSI detection down to 0.05% tumor purity (Table 1).
"""


def msi_paired():
    cmd = Command()
    cmd.meta.name = 'MSIPaired'
    cmd.meta.source = 'https://github.com/xjtu-omics/msisensor-pro'
    cmd.meta.function = 'MSI calculation for tumor normal paired samples'
    cmd.meta.version = '1.2.0'
    cmd.meta.desc = """
    This module evaluate MSI using the difference between normal and tumor length distribution of microsatellites
    Example:
    msisensor-pro msi -d /path/to/reference.list -n /path/to/case1_normal_sorted.bam -t /path/to/case1_tumor_sorted.bam -o /path/to/case1_output
    """
    cmd.runtime.image = "pengjia1110/msisensor-pro:latest"
    cmd.runtime.tool = "msisensor-pro msi"
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.cpu = 4
    cmd.args['sites'] = Argument(prefix='-d ', type='infile', desc='homopolymer and microsatellite file')
    cmd.args['normal'] = Argument(prefix='-n ', type='infile', desc='normal bam file with index')
    cmd.args['tumor'] = Argument(prefix='-t ', type='infile', desc='tumor bam file with index')
    cmd.args['out_prefix'] = Argument(prefix='-o ', type='outstr', desc='output prefix')
    cmd.args['bed'] = Argument(prefix='-e ', type='infile', level='optional', desc='bed file, optional')
    cmd.args['fdr'] = Argument(prefix='-f ', type='float', default=0.05, desc='FDR threshold for somatic sites detection, default=0.05')
    cmd.args['coverage'] = Argument(prefix='-c ', type='int', default=15, desc='coverage threshold for msi analysis, WXS: 20; WGS: 15')
    cmd.args['normalize'] = Argument(prefix='-z ', type='int', default=0, desc='coverage normalization for paired tumor and normal data, 0: no; 1: yes, default=0')
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
    cmd.outputs['out_somatic'] = Output(value='{out_prefix}_somatic', type='outfile', desc='somatic MSI sites')
    cmd.outputs['out_germline'] = Output(value='{out_prefix}_germline', type='outfile', desc='germline MSI sites')
    cmd.outputs['out_distribution'] = Output(value='{out_prefix}_dis', type='outfile', desc='MS distribution')
    cmd.outputs['out'] = Output(value='{out_prefix}', type='outfile', desc='MSI statistic summary')
    return cmd


if __name__ == '__main__':
    msi_paired().run_on_terminal(to_cwl=True)

