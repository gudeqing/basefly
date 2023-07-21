import pylab as p
import pysam


def add_contig_header(vcf, out, ref_dict=None, header_txt=None):
    """
    如果不提供header信息，将自动添加hg19的contig信息
    :param vcf:
    :param ref_dict:
    :param header_txt:
    :return:
    """
    if ref_dict:
        contig_info = []
        with open(ref_dict) as f:
            _ = f.readline()
            for line in f:
                lst = line.strip().split()
                name = lst[1].split(':')[1]
                length = lst[2].split(':')[1]
                md5 = lst[3].split(':')[1]
                ref_file = lst[4].split(':')[2]
                contig_info.append([name, length])
        contig_info = [f'##reference={ref_file}'] + contig_info
    elif header_txt:
        with open(header_txt) as f:
            contig_info = [x.strip() for x in f.readlines()]
    else:
        contig_info = [
            f'##reference=hg19',
            "##contig=<ID=chr1,length=249250621>",
            "##contig=<ID=chr2,length=243199373>",
            "##contig=<ID=chr3,length=198022430>",
            "##contig=<ID=chr4,length=191154276>",
            "##contig=<ID=chr5,length=180915260>",
            "##contig=<ID=chr6,length=171115067>",
            "##contig=<ID=chr7,length=159138663>",
            "##contig=<ID=chr8,length=146364022>",
            "##contig=<ID=chr9,length=141213431>",
            "##contig=<ID=chr10,length=135534747>",
            "##contig=<ID=chr11,length=135006516>",
            "##contig=<ID=chr12,length=133851895>",
            "##contig=<ID=chr13,length=115169878>",
            "##contig=<ID=chr14,length=107349540>",
            "##contig=<ID=chr15,length=102531392>",
            "##contig=<ID=chr16,length=90354753>",
            "##contig=<ID=chr17,length=81195210>",
            "##contig=<ID=chr18,length=78077248>",
            "##contig=<ID=chr19,length=59128983>",
            "##contig=<ID=chr20,length=63025520>",
            "##contig=<ID=chr21,length=48129895>",
            "##contig=<ID=chr22,length=51304566>",
            "##contig=<ID=chrMT,length=16569>",
            "##contig=<ID=chrX,length=155270560>",
            "##contig=<ID=chrY,length=59373566>",
        ]

    # 读取vcf
    vcf = pysam.VariantFile(vcf)
    header = vcf.header.copy()
    # 检查是否存在contig信息
    if header.contigs.__len__() >= 1:
        raise Exception('contig info is not empty!')
    # 读取一个记录判断染色体序号是否为数字
    record = next(vcf)
    if str(record.contig).isnumeric() and not ref_dict:
        contig_info = [x.replace('=chr', '=') for x in contig_info]
    # update header
    for line in contig_info:
        if type(line) == list:
            header.contigs.add(line[0], length=line[1],)
        else:
            header.add_line(line)
    # 输出新的vcf
    vcf_out = pysam.VariantFile(out, "w", header=header)
    vcf_out.write(record)
    for r in vcf:
        vcf_out.write(r)
    vcf_out.close()
    return out


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', type=Path, help='path to vcf')
    parser.add_argument('-ref_dict', type=Path, required=False, help='path to reference dict file. Recommended!')
    parser.add_argument('-header_txt', type=Path, required=False, help='path to a file contains vcf contig header lines')
    parser.add_argument('-out', help='output file name')
    args = parser.parse_args()
    add_contig_header(vcf=args.vcf, out=args.out, ref_dict=args.ref_dict, header_txt=args.header_txt)
