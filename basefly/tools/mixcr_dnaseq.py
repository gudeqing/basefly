from basefly import Argument, Output, Command, Workflow, TopVar

"""
https://www.takarabio.com/about/bioview-blog/research-news/tcr-seq-methods-strengths-weaknesses-and-rankings
Most TCR-seq methods can be broadly categorized based on the starting material. 
DNA-based methods are thought to be better suited for quantification of individual TCR clones, with a single template per cell, but this can become expensive.
RNA-based methods are more sensitive and provide quantitative estimates of TCR abundance and gene expression, 
as well as the ability to easily incorporate unique molecular identifiers (UMIs), reducing PCR bias and enhancing accurate identification of variants and rare mutations.

"""


def mixcr_dnaseq():
    cmd = Command()
    cmd.meta.name = 'mixcr'
    cmd.meta.desc = """
    该软件需要license！
    MiXCR is a universal software for fast and accurate analysis of raw T- or B- cell receptor repertoire sequencing data. 
    https://docs.milaboratories.com/mixcr/reference/mixcr-analyze/#generic-non-targeted-shotgun-data-rna-seq
    """
    cmd.meta.source = "https://github.com/milaboratory/mixcr"
    cmd.runtime.image = '?'
    cmd.runtime.memory = 5*1024**3
    cmd.runtime.cpu = 2
    cmd.runtime.tool = 'mixcr analyze shotgun'
    cmd.args['species'] = Argument(prefix='--species ', default='hsa', desc="Possible values: hsa (or HomoSapiens), mmu (or MusMusculus), rat, spalax, alpaca, lamaGlama, mulatta (Macaca Mulatta), fascicularis (Macaca Fascicularis) or any species from IMGT ® library")
    cmd.args['material'] = Argument(prefix='--starting-material ', default='dna', desc='Starting material. Possible values: rna, dna')
    cmd.args['receptor'] = Argument(prefix='--receptor-type ', level='optional', default='tcr', desc='')
    cmd.args['contig-assembly'] = Argument(prefix='contig-assembly', type='bool', default=True)
    cmd.args['threads'] = Argument(prefix='-t ', default=4, desc='threads number')
    cmd.args['only-productive'] = Argument(prefix='--only-productive', type='bool', default=True, desc='Filter out-of-frame sequences and clonotypes with stop-codons in clonal sequence export')
    cmd.args['impute-germline-on-export'] = Argument(prefix='--impute-germline-on-export', type='bool', default=True, desc='Impute germline on export')
    cmd.args['report'] = Argument(prefix='--report ', desc='')
    cmd.args['read1'] = Argument(prefix='', type='infile', desc='input read1 file')
    cmd.args['read2'] = Argument(prefix='', type='infile', level='optional')
    cmd.args['out_prefix'] = Argument(prefix='', desc='output prefix')
    cmd.outputs['report'] = Output(value='{report}')
    cmd.outputs['TRAD'] = Output(value='{out_prefix}.clonotypes.TRAD.tsv')
    cmd.outputs['TRB'] = Output(value='{out_prefix}.clonotypes.TRB.tsv')
    cmd.outputs['IGH'] = Output(value='{out_prefix}.clonotypes.IGH.tsv')
    cmd.outputs['IGK'] = Output(value='{out_prefix}.clonotypes.IGk.tsv')
    cmd.outputs['IGL'] = Output(value='{out_prefix}.clonotypes.IGL.tsv')


