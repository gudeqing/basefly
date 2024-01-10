# coding=utf-8
import os
import re
import json
__author__ = 'gudeqing'


def get_fastq_info(fastq_info:tuple, pair_info=None, out='fastq.info.json',
                   r1_name="(.*).R1.fq.gz", r2_name="(.*).R2.fq.gz",
                   link_data=False, add_s_to_numeric_name=False, middle2underscore=False):
    """
    :param fastq_info: a list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json]
    :param pair_info: 'pair info file that contains two columns without any header: [tumor_name, normal_name]
    :param r1_name: python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'
    :param r2_name: python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'
    :param link_data: bool to indicate if to make soft links for fastq files
    :param out: output file that contains three columns: [sample_name, read1_abs_path, read2_abs_path]
    :param add_s_to_numeric_name: bool value to indicate if to add a 'S' letter at the head of the sample name that startswith numeric string.
    :param middle2underscore: bool value to indicate if to transform '-' letter to '_' letter for a sample name.
    :return: result_dict: {sample: [[r1, r1'], [r2, r2']], ...}
    """
    result_dict = dict()
    fastq_dirs = []
    fastq_files = []
    for each in fastq_info:
        if os.path.isdir(each):
            fastq_dirs.append(os.path.abspath(each))
        elif os.path.isfile(each):
            if each.endswith(('.fq', 'fq.gz', 'fastq', 'fastq.gz')):
                fastq_files.append(os.path.abspath(each))
            elif each.endswith('.json'):
                with open(each) as f:
                    result_dict.update(json.load(f))
            else:
                with open(each) as f:
                    for line in f:
                        lst = line.strip().split('\t')
                        tmp = result_dict.setdefault(lst[0], list())
                        tmp.append(lst[1].split(';'))
                        if len(lst) >= 3:
                            tmp.append(lst[2].split(';'))

    if r1_name == r2_name:
        raise Exception('read1 filename == read2 filename ?!')

    if fastq_files:
        for each in fastq_files:
            name = os.path.basename(each)
            directory = os.path.dirname(each)
            is_read1 = True
            match = re.fullmatch(r1_name, name)
            if not match:
                match = re.fullmatch(r2_name, name)
                is_read1 = False
            if match:
                # first matched group is sample name
                sample = match.groups()[0]
                result_dict.setdefault(sample, [[], []])
                if is_read1:
                    if each not in result_dict[sample][0]:
                        result_dict[sample][0].append(each)
                    else:
                        print(f'warn: duplicated path found for {each}, and we will only keep the first one!')
                else:
                    if each not in result_dict[sample][1]:
                        result_dict[sample][1].append(each)
                    else:
                        print(f'warn: duplicated path found for {each}, and we will only keep the first one!')

    if fastq_dirs:
        for path in fastq_dirs:
            for root, dirs, files in os.walk(path):
                for each in files:
                    is_read1 = True
                    match = re.fullmatch(r1_name, each)
                    if not match:
                        match = re.fullmatch(r2_name, each)
                        is_read1 = False
                    if match:
                        # first matched group is sample name
                        sample = match.groups()[0]
                        result_dict.setdefault(sample, [[], []])
                        file_path = os.path.join(root, each)
                        if is_read1:
                            if file_path not in result_dict[sample][0]:
                                result_dict[sample][0].append(file_path)
                            else:
                                print(f'warn: duplicated path found for {file_path}, and we will only keep the first one!')
                        else:
                            if file_path not in result_dict[sample][1]:
                                result_dict[sample][1].append(file_path)
                            else:
                                print(f'warn: duplicated path found for {file_path}, and we will only keep the first one!')

    new_result = dict()
    if link_data:
        os.mkdir('rawdata')
        os.chdir('rawdata')
    for sample, lst in result_dict.items():
        read1 = sorted(lst[0])
        read2 = sorted(lst[1])
        if middle2underscore:
            sample = sample.replace('-', '_')
        if add_s_to_numeric_name:
            if sample.startswith(('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')):
                sample = 'S' + sample
        new_result[sample] = [read1, read2]
        if link_data:
            # make link
            os.mkdir(sample)
            for each in read1:
                os.symlink(each, os.path.join(sample, os.path.basename(each)))
            for each in read2:
                os.symlink(each, os.path.join(sample, os.path.basename(each)))

    if pair_info:
        with open(pair_info) as fr, open(out+'.pair', 'w') as f:
            for line in fr:
                tumor, normal = line.strip().split()
                if tumor in new_result and normal in new_result:
                    tr1 = ';'.join(new_result[tumor][0])
                    tr2 = ';'.join(new_result[tumor][1])
                    nr1 = ';'.join(new_result[normal][0])
                    nr2 = ';'.join(new_result[normal][1])
                    lst = [tumor, tr1, tr2, normal, nr1, nr2]
                    f.write('\t'.join(lst) + '\n')
                else:
                    print(f'{tumor} or {normal} fastq is not found !')

    if out.endswith('.json'):
        with open(out, 'w') as f:
            json.dump(new_result, f, indent=2)
    else:
        with open(out, 'w') as f:
            for k, v in new_result.items():
                read1 = ';'.join(v[0])
                read2 = ';'.join(v[1])
                f.write('{k}\t{read1}\t{read2}\n'.format(k=k, read1=read1, read2=read2))

    return new_result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-fastq_info', required=True, nargs='+', help="A list with elements from [fastq file, fastq parent dir, fastq_info.txt, fastq_info.json].")
    parser.add_argument('-r1_name', required=True, help="python regExp that describes the full name of read1 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R1.fq.gz'")
    parser.add_argument('-r2_name', required=True, help="python regExp that describes the full name of read2 fastq file name. It requires at least one pair small brackets, and the string matched in the first pair brackets will be used as sample name. Example: '(.*).R2.fq.gz'")
    parser.add_argument('-out', required=False, default='fastq.info.json', help='output file that contains three columns: [sample_name, read1_abs_path, read2_abs_path]')
    parser.add_argument('-pair_info', required=False, help='pair info file that contains two columns without any header: [tumor, normal]')
    parser.add_argument('--add_s_to_numeric_name', action='store_true', default=False, help="if to add a 'S' letter at the head of the sample name that startswith numeric string.")
    parser.add_argument('--middle2underscore', action='store_true', default=False, help="if to transform '-' letter to '_' letter for a sample name")
    parser.add_argument('--link', action='store_true', default=False, help="if to transform '-' letter to '_' letter for a sample name")
    args = parser.parse_args()
    get_fastq_info(
        fastq_info=args.fastq_info,
        r1_name=args.r1_name, r2_name=args.r2_name, pair_info=args.pair_info,
        link_data=args.link, out=args.out,
        add_s_to_numeric_name=args.add_s_to_numeric_name,
        middle2underscore=args.middle2underscore
    )
