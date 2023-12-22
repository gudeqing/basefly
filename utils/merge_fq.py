import os


def merge_fq(fq_lst, outdir, group='tumor'):
    # 准备输出目录
    r1_outdir = os.path.join(outdir, group + '_R1')
    r2_outdir = os.path.join(outdir, group + '_R2')
    os.makedirs(r1_outdir)
    os.makedirs(r2_outdir)
    # 使用第一对fastq的原名作为合并后文件的名称
    fq_r1_name = os.path.basename(fq_lst[0])
    fq_r2_name = os.path.basename(fq_lst[1])

    r1 = []
    r2 = []
    if len(fq_lst) == 2:
        print(f'Just one pair fastq for group {group}, no need to merge')
        os.system(f'cp {fq_lst[0]} {r1_outdir}')
        os.system(f'cp {fq_lst[1]} {r2_outdir}')
    if len(fq_lst) % 2 != 0:
        raise Exception('fastq number is not even, please confirm the fastq input info')
    for idx, fq in enumerate(fq_lst):
        # 偶数为r1, 奇数为r2
        if idx % 2 == 0:
            r1.append(fq)
        else:
            r2.append(fq)
    os.system(f"cat {' '.join(r1)} > {r1_outdir}/{fq_r1_name}")
    os.system(f"cat {' '.join(r2)} > {r2_outdir}/{fq_r2_name}")


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('outdir', help='fastq outdir')
    parser.add_argument('fastqs', help='fastq info, sep by comma and white space for intra-group and between group')
    args = parser.parse_args()
    fastq_group_lst = args.fastqs.split()
    normal = None
    if len(fastq_group_lst) == 1:
        tumor = fastq_group_lst[0]
    elif len(fastq_group_lst) == 2:
        tumor, normal = fastq_group_lst
    else:
        raise Exception('too many group of fastq files !')
    merge_fq(tumor.split(','), args.outdir, group='tumor')
    if normal:
        merge_fq(normal.split(','), args.outdir, group='normal')


