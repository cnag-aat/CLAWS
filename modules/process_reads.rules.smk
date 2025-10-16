from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule trim_galore:
  input:
    read1 = "illumina.1.fastq.gz",
    read2 = "illumina.2.fastq.gz",
  output:
    trim1 = "illumina.trimed.1.fastq.gz",
    trim2 = "illumina.trimed.2.fastq.gz",
  params:
    outdir = "illumina_trim",
    opts = "--gzip -q 20 --paired --retain_unpaired",
  threads: 4
  conda:
    "../envs/trim_galore0.6.10.yaml"
  shell:
    "mkdir -p {params.outdir};"
    "cd {params.outdir}; "
    "trim_galore -j {threads} {params.opts} {input.read1} {input.read2} ;"
    "b=`basename {input.read1} .1.fastq.gz`;"
    "ln -s $b.1_val_1.fq.gz {output.trim1};"
    "ln -s $b.2_val_2.fq.gz {output.trim2};"

rule concat_reads:
  input:
    fastqs = ["reads1.fastq.gz", "reads2.fastq.gz"]
  output:
    final_fastq = "all_reads.fastq.gz"
  log:
    "{dir}logs/" + str(date) + ".j%j.concat.{ext}.out",
    "{dir}logs/" + str(date) + ".j%j.concat.{ext}.err"
  benchmark:
    "{dir}logs/" + str(date) + ".concat.benchmark.{ext}.txt"
  conda:
    '../envs/ass_base.yaml'
  threads: 2
  shell:
    "zcat {input.fastqs} | pigz -p {threads} -c  > {output.final_fastq};"

rule build_meryl_db:
  input:
    fastq = "all_reads.fastq.gz"
  output:
     out_dir = directory("meryl_db.meryl")
  params:
    kmer = 21
  conda:
    "../envs/merqury1.3.yaml"
  threads: 4
  shell:
    "meryl k={params.kmer} count output {output.out_dir} {input.fastq};"

rule concat_meryl:
  input:
    input_run = "meryl_db"
  output:
    meryl_all = directory("meryl_db_all"),
    histogram = "meryl.hist"
  params:
    kmer = 21,
  conda:
    "../envs/merqury1.3.yaml"
  threads: 4
  shell:
    "meryl union-sum output {output.meryl_all} {input.input_run};"
    "meryl histogram {output.meryl_all} > {output.histogram};"

rule genomescope2:
  input:
    histogram = "meryl.hist"
  output:
    outdir = os.getcwd(),
    summary = "summary.txt",
    log_plot = "log_plot.png",
    transformed_linear = "transformed_linear_plot.png"
  params:
    ploidy = 2,
    kmer = 21,
    opts = ""
  conda:
    "../envs/genomescope2.0.yaml"
  threads: 1
  shell:
    "genomescope2 -i {input.histogram} -o {output.outdir} -p {params.ploidy} -k {params.kmer} {params.opts};"

rule smudgeplot:
  input:
    histogram = "meryl.hist",
    meryl = "meryl_db",
  output:
    plot = "smudgeplot_smudgeplot.png"
  params:
    dir = "smudgeplot_k21"
  conda:
    "../envs/merqury1.3.yaml"
  threads: 2
  shell:
    "mkdir -p {params.dir}; cd {params.dir}; "
    "smudgeplot.py cutoff {input.histogram} U > upper.cutoff;"
    "smudgeplot.py cutoff {input.histogram} L > lower.cutoff;"
    "up=$(cat upper.cutoff);"
    "low=$(cat lower.cutoff);"
    "meryl print less-than $up greater-than $low {input.meryl} > meryl.U$up.L$low.dump;"
    "smudgeplot.py hetkmers -o meryl.U$up.L$low < meryl.U$up.L$low.dump;"
    "smudgeplot.py plot meryl.U$up.L$low''_coverages.tsv;"

rule long_ranger:
  input: 
    mkfastq_dir = "input_dir/"
  output:
    fastq_out = "preprocess_10X_linkedreads/10X.barcoded.fastq.gz",
    sum_out = "preprocess_10X_linkedreads/10X.summary.csv"
  params:
    path = "/scratch/project/devel/aateam/src/10X/longranger-2.2.2",
    outdir = "preprocess_10X_linkedreads/",
    sample = "raw_10X",
    rmcmd = ""
  threads: 8
  shell:
      "export PATH={params.path}:$PATH;"
      "mkdir -p {params.outdir};"
      "cd {params.outdir};"
      "longranger basic --id={params.sample} --sample={params.sample} --fastqs={input.mkfastq_dir} --localcores={threads};"
      "cp {params.sample}/outs/barcoded.fastq.gz {output.fastq_out};"
      "cp {params.sample}/outs/summary.csv {output.sum_out};"
      "{params.rmcmd}"

rule filtlong:
  input:
    reads = "ont_reads.fastq.gz"
  output:
    outreads = "ont_reads.filtlong.fastq.gz"
  params:
    minlen = 1000,
    min_mean_q = 80,
    opts = ""
  threads: 8
  log:
    "{dir}logs/" + str(date) + ".j%j.filtlong.out",
    "{dir}logs/" + str(date) + ".j%j.filtlong.err"
  benchmark:
    "{dir}logs/" + str(date) + ".filtlong.benchmark.txt"
  conda:
    '../envs/filtlong0.2.1.yaml'
  shell:
    "filtlong --version;"
    "filtlong --min_length {params.minlen} --min_mean_q {params.min_mean_q} {params.opts} {input.reads} | pigz -p {threads} -c > {output.outreads};"

rule nanoplot:
  input:
    fastq = "reads.ont.fastq.gz"
  output:
    stats = "nanostats_out/NanoStats.txt",
    yield_len = "nanostats_out/Yield_By_Length.png",
    read_len = "nanostats_out/WeightedHistogramReadlength.png"
  params:
    outdir = "nanostats_out/",
  log:
    "{dir}logs/" + str(date) + ".NanoStats.out",
    "{dir}logs/" + str(date) + ".NanoStats.err"
  benchmark:
    "{dir}logs/" + str(date) + ".NanoStats.benchmark.txt"
  threads: 4
  conda:
    '../envs/nanoplot1.40.0'
  shell:
    "mkdir -p {params.outdir};" 
    "cd {params.outdir}; "
    "NanoPlot -t {threads} --plots dot --fastq {input.fastq} -o .; "

rule Kraken2:
  input:
    read = "reads.fastq.gz",
    database = "/scratch/groups/assembly/shared/databases/kraken2/bacteria_db/",
    kmers = "/scratch/groups/assembly/shared/databases/kraken2/bacteria_db/database100mers.kmer_distrib"
  output:
    report = "kraken2.report",
    readsout = "kraken2.seqs.out",
    abundance =  "bracken_abundance.txt",
  params:
    additional = "",
    prefix = "assembly"
  log:
    "logs/" + str(date) + ".kraken.out",
    "logs/" + str(date) + ".kraken.err",
  benchmark:
    "logs/" + str(date) + ".kraken.benchmark.txt",
  threads: 4
  conda:
    '../envs/kraken2.1.2.yaml'
  shell:
     "kraken2 --threads {threads} --db {input.database}  --use-names --report {output.report} {params.additional} {input.read} > {output.readsout}; "
     "est_abundance.py -i {output.report} -k {input.kmers} -l S -t 10 -o {output.abundance}; "