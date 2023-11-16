import pymongo
from pysam import VariantFile, VariantHeader


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


Variant_DB = {
    'sample_id',
    'sample_name',
    'patient_id',
    'workflow',
    'time',
    'cancer'
}


def parse_vep_csq(row_record, csq_header):
    # csq_header = fr.header.info['CSQ'].description.split('Format: ')[1]
    csq_dict_list = []
    for each in row_record.info['CSQ']:
        csq_dict = dict(zip(csq_header.split('|'), each.split('|')))
        for key, value in csq_dict.items():
            if type(value) == str and '%3D' in value:
                csq_dict[key] = value.replace('%3D', '=')
        csq_dict_list.append(csq_dict)
    return csq_dict_list


def parse_info_field(row_record, csq_header=None):
    info_dict = dict(row_record.info)
    if 'CSQ' in info_dict:
        csq_dict_lst = parse_vep_csq(info_dict['CSQ'], csq_header=csq_header)
        info_dict['CSQ'] = csq_dict_lst
    return info_dict


def create_variant_db(vcfs:list, sources:list, canonical_refseq=None, host='0.0.0.0', port=27017, db_name='small_variants', collection_name='hgvs'):
    client = pymongo.MongoClient(host=host, port=port)
    db = client[db_name]
    coll = db[collection_name]

    if len(vcfs) > 1 and len(sources) == 1:
        sources = sources*len(vcfs)
    elif len(vcfs) == len(sources):
        pass
    else:
        raise Exception("vcfs and sources are not properly matched")

    if canonical_refseq is not None:
        g2t = dict([x.strip().split()[:2] for x in open(canonical_refseq) if len(x.strip().split()) > 1])
    else:
        g2t = dict()

    for vcf, source in zip(vcfs, sources):
        with VariantFile(vcf) as f:
            main_keys = ['source', 'chrom', 'start', 'var_id', 'gene', 'hgvs', 'ref', 'alt', 'canonical']
            hgvs_keys = ['gene', 'transcript', 'exon', 'chgvs', 'phgvs']
            for record in f:
                if record.info['AAChange_refGene'][0] != '.':
                    chrom = record.chrom
                    start = record.pos
                    ref = record.ref
                    alt = record.alts[0]
                    var_id = record.id
                    canonical = 'None'
                    hgvs_values = []
                    for hgvs in record.info['AAChange_refGene']:
                        detail = hgvs.split(':')
                        gene = detail[0]
                        if len(detail) == 5:
                            g, t, e, c, p = detail
                            hgvs_values = detail
                        elif len(detail) == 4:
                            g, t, e, c = detail
                            hgvs_values = [g, t, e, c, 'None']
                        else:
                            print(f'{detail} is out of expectation, the stored info maybe Wrong!')
                            hgvs_values = detail

                        if gene in g2t:
                            ct = g2t[gene]
                            if ct.startswith(detail[1]):
                                canonical = hgvs

                        hgvs_dict = dict(zip(hgvs_keys, hgvs_values))
                        main_values = [source, chrom, start, var_id, gene, hgvs_dict, ref, alt, canonical]
                        coll.insert_one(dict(zip(main_keys, main_values)))
                else:
                    print('skip line with empty "AAChange_refGene"')


def query_variant_db(query, out='query_result.vcf', host='0.0.0.0', port=27017, db_name='variants', coll_name='hgvs'):
    """
    支持的header有 ['gene', 'transcript', 'exon', 'chgvs', 'phgvs']
    """
    import pymongo
    client = pymongo.MongoClient(host=host, port=port)
    db = client[db_name]
    coll = db[coll_name]

    with open(query) as f, open(out, 'w') as fw:
        for line in vcf_header():
            fw.write(line+'\n')

        header = f.readline().strip().split()
        header = ['hgvs.'+x if x in ('transcript', 'exon', 'chgvs', 'phgvs') else x for x in header]
        new_rows = []
        for line in f:
            condition = dict(zip(header, line.strip().split()))
            docs = list(coll.find(condition))
            if not docs:
                print('failed query:', line.strip())
                continue
            # sources = [doc['source'] for x in docs]
            # mutation_set = {(doc['chrom'], doc['start'], doc['ref'], doc['alt']) for x in docs}
            filtered_docs = [doc for doc in docs if doc["canonical"] != "None"]
            if filtered_docs:
                docs = filtered_docs
            # if len(docs) > 1:
                # print(f'{len(docs)} records for {line.strip()}')
            for doc in docs:
                lst = [doc[x] for x in ['chrom', 'start', 'var_id', 'ref', 'alt']]
                lst.append('.')
                lst.append("PASS")
                query = f'Query={condition}'
                if doc["canonical"] == 'None':
                    # lst.append(f'{query};Source={doc["source"]};AAChange_refGene={doc["hgvs"]}')
                    hgvs = ':'.join(doc['hgvs'].values())
                    lst.append(f'AAChange_refGene={doc["hgvs"]}')
                else:
                    # lst.append(f'{query};Source={doc["source"]};AAChange_refGene={doc["canonical"]}')
                    lst.append(f'AAChange_refGene={doc["canonical"]}')
                new_line = '\t'.join((str(x) for x in lst))+'\n'
                if new_line not in new_rows:
                    fw.write(new_line)
                    new_rows.append(new_line)


def del_db(host='0.0.0.0', port=27017, db_name='variants', coll_name='hgvs', condition=None):
    import pymongo
    client = pymongo.MongoClient(host=host, port=port)
    if db_name not in client.list_database_names():
        exit(f'databse {db_name} is not found')
    db = client[db_name]

    colls = db.list_collection_names()
    print('Find collections:', colls)
    if coll_name not in colls:
        exit(f'{coll_name} is not found')
    coll = db[coll_name]

    if condition is None:
        print(f'You are dropping the whole collection {coll_name}')
        password = input('Password:')
        if password != 'yes':
            exit('Password is Invalid')
        else:
            coll.drop()
    else:
        query = dict(x.split('=') for x in condition.split(','))
        print(f'You are deleting docs by querying {query}')
        password = input('Password:')
        if password != 'yes':
            exit('Password is Invalid')
        else:
            if len(list(coll.find(query))) > 0:
                coll.delete_many(query)
            else:
                print('Nothing matched and no deletion')
