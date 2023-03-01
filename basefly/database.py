class database():
    """
    定义一个结果输出数据库
    预制如下表：

    基因表达量表
    转录本表达量表
    融合基因表
    GATK: 小snp/indel表
    mixcr: TCR clone结果表
    quantiseq: 结果表
    HLA定型结果表
    先抗原结果表
    """
    def __init__(self, name=None, password=None):
        pass

    def sample_info(self):
        meta_fields = [
            'Sample_id',
            'Sample_name',
            'Patient_id',
            'Age',
            'Gender',
            'Disease',
            'DiseaseSubtype',
        ]
        pass

    def gene_info(self):
        pass

    def make_gene_exp_table(self, feature_type='gene', exp_type='TPM'):
        pass

    def make_starfusion_table(self):
        pass

    def make_somatic_snp_indel_table(self):
        pass

    def make_germline_snp_indel_table(self):
        pass

    def make_tcr_clone_table(self):
        pass

    def make_hla_type_table(self):
        pass

    def make_quantiseq_table(self):
        pass

    def make_mhcflurry_table(self):
        pass

    def make_netmhciipan_table(self):
        pass