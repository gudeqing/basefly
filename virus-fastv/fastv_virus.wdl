version development

workflow fastvirus{
    call virus_fastv{}
}


task virus_fastv{
    input {
        Directory fastq_dir
        Directory tool_dir
        String r1_name = '(.*)_R1.fastq.gz'
        String r2_name = '(.*)_R2.fastq.gz'
        File kmer = '/data_analysis/software/fastv-0.8.1/data/SARS-CoV-2.kmer.fa'
        File genome = '/data_analysis/software/fastv-0.8.1/data/SARS-CoV-2.genomes.fa'
    }

    command <<<
        set -e
        export LC_ALL="en_US.utf8"
        python3 \
        ~{tool_dir}/virus-fastv/virus_fastv.py \
        ~{"-k " + kmer} \
        ~{"-g " + genome} \
        ~{"-fastq_info " + fastq_dir} \
        ~{"-r1_name " + "'" + r1_name + "'"} \
        ~{"-r2_name " + "'" + r2_name + "'"} \
        ~{"-tool_dir " + tool_dir} \
        --run

    >>>

    output {
        Directory result_dir = "Result"
    }

    runtime {
        docker: "registry-xdp-v3-yifang.xdp.basebit.me/basebitai/zheng_fastv:002"
    }

    meta {
        name: "pipeline"
        desc: "this is a pipeline wrapper using wdl"
        version: "002"
    }
}