#!/usr/bin/env python

import sys
import os
from pysam import AlignmentFile, VariantFile
import tempfile
import csv
from subprocess import Popen, PIPE


def generate_region_list(hash, wkdir='.'):
    fh = tempfile.NamedTemporaryFile('w', dir=wkdir, delete=False)
    writer = csv.writer(fh, delimiter='\t')
    for chr_, positions in hash.items():
        for pos in sorted(positions.keys()):
            writer.writerow([chr_, pos, pos])
    fh.close()
    return fh.name


def filter_sites_in_hash(region_list, bam_file, ref_fasta, sample, output_dir, insertion_centric, map_qual, base_qual):
    bam_readcount_cmd = ['bam-readcount', '-f', ref_fasta, '-l', region_list, '-w', '0', '-b', str(base_qual), '-q', str(map_qual)]
    if insertion_centric:
        bam_readcount_cmd.append('-i')
        output_file = os.path.join(output_dir, sample + '.bam_readcount.indel.tsv')
    else:
        output_file = os.path.join(output_dir, sample + '.bam_readcount.snv.tsv')
    bam_readcount_cmd.append(bam_file)
    execution = Popen(bam_readcount_cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = execution.communicate()
    if execution.returncode == 0:
        with open(output_file, 'wb') as output_fh:
            output_fh.write(stdout)
    else:
        sys.exit(stderr)


def run(vcf_file, sample, ref_fasta, bam_file, output_dir='.', min_base_qual=20, min_mapping_qual=0, not_only_pass=False):
    rc_for_indel = {}
    rc_for_snp = {}

    has_chr_prefix = False
    with AlignmentFile(bam_file) as bam:
        one_contig = bam.get_reference_name(0)
        if one_contig.startswith('chr'):
            has_chr_prefix = True

    with VariantFile(vcf_file) as vcf:
        for variant in vcf:
            if list(variant.filter)[0].lower() != 'pass' and (not not_only_pass):
                continue
            ref = variant.ref
            contig = variant.contig
            if contig[:3] != 'chr' and has_chr_prefix:
                contig = 'chr' + contig
            pos = variant.pos
            for var in variant.alts:
                if len(ref) > 1 or len(var) > 1:
                    #it's an indel or mnp
                    if len(ref) == len(var) or (len(ref) > 1 and len(var) > 1):
                        sys.stderr.write("Complex variant or MNP will be skipped: %s\t%s\t%s\t%s\n" % (contig, pos, ref, var))
                        continue
                    elif len(ref) > len(var):
                        #it's a deletion
                        pos += 1
                    else:
                        # it's an insertion
                        pass
                    rc_for_indel.setdefault(contig, dict()).setdefault(pos, list())
                    rc_for_indel[contig][pos].append((ref, var))
                else:
                    #it's a SNP
                    rc_for_snp.setdefault(contig, dict()).setdefault(pos, list())
                    rc_for_snp[contig][pos].append((ref, var))

    if len(rc_for_snp.keys()) > 0:
        region_file = generate_region_list(rc_for_snp, wkdir=output_dir)
        filter_sites_in_hash(region_file, bam_file, ref_fasta, sample, output_dir, False, min_mapping_qual, min_base_qual)
    else:
        output_file = os.path.join(output_dir, sample + '.bam_readcount.snv.tsv')
        open(output_file, 'w').close()

    if len(rc_for_indel.keys()) > 0:
        region_file = generate_region_list(rc_for_indel, wkdir=output_dir)
        filter_sites_in_hash(region_file, bam_file, ref_fasta, sample, output_dir, True, min_mapping_qual, min_base_qual)
    else:
        output_file = os.path.join(output_dir, sample + '.bam_readcount.indel.tsv')
        open(output_file, 'w').close()


if __name__ == '__main__':
    from xcmds import xcmds
    xcmds.xcmds(locals(), include=['run'])
