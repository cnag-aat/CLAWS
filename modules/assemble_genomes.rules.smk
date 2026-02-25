from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

if not os.path.exists("logs"):
  os.makedirs("logs")

rule flye:
  input:
    reads = os.getcwd() + "/ontreads.fastq.gz",
  output:
    assembly = os.getcwd() + "/flye/flye.assembly.fasta",
    gfa = os.getcwd() + "flye/assembly_graph.gfa",
    gfa_plot = os.getcwd() + "flye/assembly_graph.gfa.png"
  params:
    outdir = os.getcwd() + "/flye/",
    readtype = "nano-raw",
    pol_iterations = 2,
    other_flye_opts = ""  #include here genome size in pipeline
  threads: 24
  log:
    "logs/" + str(date) + ".flye.out",
    "logs/" + str(date) + ".flye.err",
  benchmark:
    "logs/" + str(date) + ".flye.benchmark.txt"
  conda:
    '../envs/flye2.9.5.yaml'
  envmodules: 
    'Mesa/21.1.7-GCCcore-11.2.0'
  shell:
    "mkdir -p {params.outdir}out;"
    "cd {params.outdir};"
    "echo 'Running command: flye --{params.readtype} {input.reads} -o {params.outdir}out -t {threads} -i {params.pol_iterations} {params.other_flye_opts}';"
    "flye --{params.readtype} {input.reads} -o {params.outdir}out -t {threads} -i {params.pol_iterations} {params.other_flye_opts};"
    "ln -s {params.outdir}out/assembly.fasta {output.assembly};"
    "ln -s {params.outdir}out/assembly_graph.gfa {output.gfa};"
    "Bandage image {output.gfa} {output.gfa_plot};"

rule nextdenovo:
  input:
    reads = os.getcwd() + "/ontreads.fastq.gz"
  output:
    assembly = os.getcwd() + "/nextDenovo/nextdenovo.assembly.fasta"
  params:
    outdir = os.getcwd() + "/nextDenovo/",
    config = "nextdenovo.cfg"
  threads: 24
  log:
    "logs/" + str(date) + ".nextdenovo.out",
    "logs/" + str(date) + ".nextdenovo.err",
  benchmark:
    "logs/" + str(date) + ".nextdenovo.benchmark.txt"
  envmodules:
    "NextDenovo/2.5.0"
  shell:
    "mkdir -p {params.outdir};"
    "ls {input.reads} > {params.outdir}long_reads.fofn;"
    "cp {params.config} {params.outdir};"
    "nextDenovo {params.config};"
    "ln -s {params.outdir}03.ctg_graph/nd.asm.fasta {output.assembly};"

rule hifiasm:
  input:
    reads = os.getcwd() + "/ontreads.fastq.gz",
    h1 = "hic_reads.1.fastq.gz",
    h2 = "hic_reads.2.fastq.gz"
  output:
      fastas = [os.getcwd() + "/hifiasm/hfsm.asm.bp.hap1.fasta", os.getcwd() + "/hifiasm/hfsm.asm.bp.hap2.fasta"],
      gfas = [os.getcwd() + "/hifiasm/hfsm.asm.bp.hap1.gfa", os.getcwd() + "/hifiasm/hfsm.asm.bp.hap2.gfa"],
      gfa_plots = [os.getcwd() + "hifiasm/hfsm.asm.bp.hap1.gfa.gfa.png", os.getcwd() + "hifiasm/hfsm.asm.bp.hap2.gfa.png"]
  params:
    outdir = os.getcwd() + "/hifiasm/",
    out_base = "hfsm.asm",
    phase = "",
    teloseq = "TTAGGG",
    other_hfsm_opts = " --ont ",
  threads: 24
  log:
    "logs/" + str(date) + ".hifiasm.out",
    "logs/" + str(date) + ".hifiasm.err",
  benchmark:
    "logs/" + str(date) + ".hifiasm.benchmark.txt"
  envmodules: 
    'Mesa/21.1.7-GCCcore-11.2.0'
  conda:
    '../envs/hifiasm0.24.0-r702.yaml'
  shell:
    "mkdir -p {params.outdir};"
    "cd {params.outdir};"
    "echo 'Running command: hifiasm {params.other_hfsm_opts} {params.teloseq} {params.phase}  -o {params.out_base} -t {threads} {input.reads}';"
    "hifiasm {params.other_hfsm_opts}  {params.teloseq} {params.phase} -o {params.out_base} -t {threads} {input.reads};"
    "for i in {output.fastas}; do b=`basename $i \".fa\"`; ln -s $b.p_ctg.gfa $b.gfa; awk \'/^S/{{print \">\"$2\"\\n\"$3}}\' $b.gfa > $b.fa; Bandage image $b.gfa $b.gfa.png;done"

