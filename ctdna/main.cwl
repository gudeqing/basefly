cwlVersion: v1.2
class: Workflow
requirements:
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  project_label:
    type: string?
  fq:
    type:
      type: array
      items:
        type: array
        items: File
    default:
      [
        [
          { "class": File, "path": /home/hxbio04/data/uCaler_AML_Panel/data/uCaler_AML_Panel_v1_Mix_Illumina/tumor-down0.1_T_R1.fastq.gz },
          { "class": File, "path": /home/hxbio04/data/uCaler_AML_Panel/data/uCaler_AML_Panel_v1_Mix_Illumina/tumor-down0.1_T_R2.fastq.gz }
        ],
#        [
#          { "class": File, "path": /home/hxbio04/data/uCaler_AML_Panel/data/uCaler_AML_Panel_v1_Promega_male_G1471_Illumina/LQCL230117-9_N_R1.fastq.gz },
#          { "class": File, "path": /home/hxbio04/data/uCaler_AML_Panel/data/uCaler_AML_Panel_v1_Promega_male_G1471_Illumina/LQCL230117-9_N_R2.fastq.gz }
#        ],
      ]
  bed:
    type: File
    default:
      class: File
      path: /home/hxbio04/data/uCaler_AML_Panel/uCaler_AML_Panel_v1.0.hg19.sorted.bed

  ref:
    type: File
    default:
      class: File
      path: /home/hxbio04/dbs/hs37d5/hs37d5.fa
    secondaryFiles:
      - .fai?
      - ^.dict?
      - ".0123?"
      - .ann?
      - .bwt.2bit.64?
      - .pac?
      - .amb?
      - .bwt?
      - .sa?
  ref_dict:
    type: File
    default:
      class: File
      path: /home/hxbio04/dbs/hs37d5/hs37d5.dict
  StatSeqError_tumor_script:
    type: File
    default:
      class: File
      path: /home/hxbio04/basefly/utils/stat_3bases_error.py
  AddVcfContig_tumor_script:
    type: File
    default:
      class: File
      path: /home/hxbio04/basefly/utils/add_vcf_contig.py
  VcfFilter_vardict_tumor_script:
    type: File
    default:
      class: File
      path: /home/hxbio04/basefly/utils/vcf_filter.py

steps:
  step1:
    run: ctDNApaired.cwl
    in:
      fq: fq
      bed: bed
      ref: ref
      ref_dict: ref_dict
      StatSeqError_tumor_script: StatSeqError_tumor_script
      AddVcfContig_tumor_script: AddVcfContig_tumor_script
      VcfFilter_vardict_tumor_script: VcfFilter_vardict_tumor_script
    when: $(inputs.fq.length == 2)
    out: [
      Bamdst_preUMI_tumor_outdir,
      Bamdst_final_tumor_outdir,
      ABRA2_tumor_out,
      ABRA2_tumor_out_bai,
      StatSeqError_tumor_stat_per_site,
      StatSeqError_tumor_context_error_rate,
      Bamdst_preUMI_normal_outdir,
      Bamdst_final_normal_outdir,
      ABRA2_normal_out,
      ABRA2_normal_out_bai,
      StatSeqError_normal_stat_per_site,
      StatSeqError_normal_context_error_rate,
      VcfFilter_vardict_tumor_final_vcf,
      VcfFilter_vardict_tumor_final_xls,
      VcfFilter_vardict_tumor_discarded_vcf,
      VcfFilter_vardict_tumor_log,
      VcfFilter_vardict_tumor_snv_metrics,
      MultiQC_outdir
    ]

  step2:
    run: ctDNAsingle.cwl
    in:
      fq: fq
      bed: bed
      ref: ref
      ref_dict: ref_dict
      StatSeqError_tumor_script: StatSeqError_tumor_script
      AddVcfContig_tumor_script: AddVcfContig_tumor_script
      VcfFilter_vardict_tumor_script: VcfFilter_vardict_tumor_script
    when: $(inputs.fq.length == 1)
    out: [
      ABRA2_tumor_out,
      ABRA2_tumor_out_bai,
      Bamdst_preUMI_tumor_outdir,
      Bamdst_final_tumor_outdir,
      StatSeqError_tumor_stat_per_site,
      StatSeqError_tumor_context_error_rate,
      VcfFilter_vardict_tumor_final_vcf,
      VcfFilter_vardict_tumor_final_xls,
      VcfFilter_vardict_tumor_discarded_vcf,
      VcfFilter_vardict_tumor_log,
      VcfFilter_vardict_tumor_snv_metrics,
      MultiQC_outdir
    ]

  step3:
    run: vep.cwl
    in:
      vcf_file1: step1/VcfFilter_vardict_tumor_final_vcf
      vcf_file2: step2/VcfFilter_vardict_tumor_final_vcf
      # datebase: datebase
      # dbNSFP: dbNSFP
      # custom_cosmic: custom_cosmic
      # custom_gnomad: custom_gnomad
      # SpliceAI_raw_snv: SpliceAI_raw_snv
      # SpliceAI_raw_indel: SpliceAI_raw_indel
      # custom_bed: custom_bed
      # dbscSNV: dbscSNV
      # custom_clinvar: custom_clinvar
    out: [vep_file_output]

  step4:
    run: snv-anno.cwl
    in:
      vcf_file: step3/vep_file_output
      input_project_label: project_label
    out: [out_annovar]

outputs:
  out_snv:
    type: File
    outputSource: step4/out_annovar

  Bamdst_preUMI_tumor_outdir:
    type: Directory
    outputSource:
      - step1/Bamdst_preUMI_tumor_outdir
      - step2/Bamdst_preUMI_tumor_outdir
    pickValue: first_non_null

  ABRA2_tumor_out:
    type: File
    outputSource:
      - step1/ABRA2_tumor_out
      - step2/ABRA2_tumor_out
    pickValue: first_non_null

  ABRA2_tumor_out_bai:
    type: File
    outputSource:
      - step1/ABRA2_tumor_out_bai
      - step2/ABRA2_tumor_out_bai
    pickValue: first_non_null

  Bamdst_final_tumor_outdir:
    type: Directory
    outputSource:
      - step1/Bamdst_final_tumor_outdir
      - step2/Bamdst_final_tumor_outdir
    pickValue: first_non_null

  StatSeqError_tumor_stat_per_site:
    type: File
    outputSource:
      - step1/StatSeqError_tumor_stat_per_site
      - step2/StatSeqError_tumor_stat_per_site
    pickValue: first_non_null

  StatSeqError_tumor_context_error_rate:
    type: File
    outputSource:
      - step1/StatSeqError_tumor_context_error_rate
      - step2/StatSeqError_tumor_context_error_rate
    pickValue: first_non_null

  Bamdst_preUMI_normal_outdir:
    type: Directory
    outputSource:
      - step1/Bamdst_preUMI_normal_outdir

  ABRA2_normal_out:
    type: File
    outputSource:
      - step1/ABRA2_normal_out

  ABRA2_normal_out_bai:
    type: File
    outputSource:
      - step1/ABRA2_normal_out_bai

  Bamdst_final_normal_outdir:
    type: Directory
    outputSource:
      - step1/Bamdst_final_normal_outdir

  StatSeqError_normal_stat_per_site:
    type: File
    outputSource:
      - step1/StatSeqError_normal_stat_per_site

  StatSeqError_normal_context_error_rate:
    type: File
    outputSource:
      - step1/StatSeqError_normal_context_error_rate

  VcfFilter_vardict_tumor_final_vcf:
    type: File
    outputSource:
      - step1/VcfFilter_vardict_tumor_final_vcf
      - step2/VcfFilter_vardict_tumor_final_vcf
    pickValue: first_non_null

  VcfFilter_vardict_tumor_final_xls:
    type: File
    outputSource:
      - step1/VcfFilter_vardict_tumor_final_xls
      - step2/VcfFilter_vardict_tumor_final_xls
    pickValue: first_non_null

  VcfFilter_vardict_tumor_discarded_vcf:
    type: File
    outputSource:
      - step1/VcfFilter_vardict_tumor_discarded_vcf
      - step2/VcfFilter_vardict_tumor_discarded_vcf
    pickValue: first_non_null

  VcfFilter_vardict_tumor_log:
    type: File
    outputSource:
      - step1/VcfFilter_vardict_tumor_log
      - step2/VcfFilter_vardict_tumor_log
    pickValue: first_non_null

  VcfFilter_vardict_tumor_snv_metrics:
    type: File
    outputSource:
      - step1/VcfFilter_vardict_tumor_snv_metrics
      - step2/VcfFilter_vardict_tumor_snv_metrics
    pickValue: first_non_null

  MultiQC_outdir:
    type: Directory
    outputSource:
      - step1/MultiQC_outdir
      - step2/MultiQC_outdir
    pickValue: first_non_null
