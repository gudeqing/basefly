from basefly import Argument, Output, Command, Workflow, TopVar

"""
# Download the data package
cd HLA-LA/graphs
wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz

# Index the graph
Finally, pre-compute the graph index structure - this can take a few hours and might take up to 40G of memory:
../bin/HLA-LA --action prepareGraph --PRG_graph_dir ../graphs/PRG_MHC_GRCh38_withIMGT

Interpreting the output from HLA*PRG:LA
For each sample with ID $mySampleID, the main output file is ../working/$mySampleID/hla/R1_bestguess_G.txt. Columns:

Locus: Locus (HLA-A, HLA-B, etc..).
Chromosome: 1 or 2. Calls are not phased across loci!
Allele: Called allele at G group resolution (exons 2/3 for class I genes, exon 2 for class II genes).
Q1: Quality score (~ probability), between 0 and 1.
Q2: Ignore.
AverageCoverage: Locus average coverage.
CoverageFirstDecile: First decile coverage.
MinimumCoverage: Minimum column coverage for locus.
proportionkMersCovered: Proportion k-mers belonging to the called allele observed in locus input data.
LocusAvgColumnError: Average alignment error for this locus.
NColumns_UnaccountedAllele_fGT0.2: Number of columns with evidence for the presence of alleles not accounted for by the called alleles.
perfectG: Perfect translation from HLA*PRG:LA call to G grop resolution? (1 = Yes).
"""


def HLA_LA():
    cmd = Command()
    cmd.meta.name = 'HLALA'
    cmd.meta.version = '1.0.3'
    cmd.meta.source = 'https://hxbiogit.asia/DiltheyLab/HLA-LA'
    cmd.meta.function = "HLA typing"
    cmd.meta.desc = """
    HLA*LA carries out HLA typing based on a population reference graph and employs a new linear projection method to align reads to the graph
    """
    cmd.runtime.image = 'taniguti/hla-la:latest'
    cmd.runtime.image = 'skchronicles/genome-seek_hla:v0.1.0'
    cmd.runtime.cpu = 4
    cmd.runtime.memory = 10 * 1024 ** 3
    cmd.runtime.tool = 'HLA-LA.pl'
    cmd.args['bam'] = Argument(prefix='--BAM ', type='infile', desc='input bam file')
    cmd.args['graph'] = Argument(prefix='--customGraphDir ', type='indir', default='/home/hxbio04/dbs/HLA-LA/graphs/PRG_MHC_GRCh38_withIMGT/', desc='graph indexed directory')
    cmd.args['sample_id'] = Argument(prefix='--sampleID ', type='outstr', desc='sample unique id')
    cmd.args['threads'] = Argument(prefix='--maxThreads ', type='int', default=7, desc='threads number')
    cmd.args['working_dir'] = Argument(prefix='--workingDir ', default='.', desc='working directory')
    cmd.outputs['out_G'] = Output(value='{sample_id}/hla/R1_bestguess_G.txt')
    return cmd


if __name__ == '__main__':
    Workflow().to_cwl_tool(cmd=HLA_LA(), write_out=True)
    HLA_LA().run_on_terminal()
