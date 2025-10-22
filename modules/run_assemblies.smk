from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module assembly_workflow:
  snakefile: "../modules/assemble_genomes.rules.smk"

working_dir = config["Outputs"]["base_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]

keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

assembled = []

##0. Define path for files and variables

flye_dir = config["Outputs"]["flye_dir"]
nextdenovo_dir = config["Outputs"]["nextdenovo_dir"]
hifiasm_dir = config["Outputs"]["hifiasm_dir"]
flye_assembly = config["Outputs"]["flye_out"]
nextdenovo_assembly = config["Outputs"]["nextdenovo_out"]
hifiasm_assemblies = config["Outputs"]["hifiasm_out"]

if re.search("nano", config["Parameters"]["lr_type"]):
  reads2assemble = ONT_filtered
elif re.search("pacbio", config["Parameters"]["lr_type"]):
  reads2assemble = hifi_reads
  if re.search(".bam", hifi_reads):
      hifi_base = os.path.basename(hifi_reads).replace(".bam", "")
      reads2assemble = config["Outputs"]["preprocess_lr"] + hifi_base + "_bam.fastq"

##1. Run assemblers
if config["Parameters"]["run_flye"] == True:
  assembled.append(flye_assembly)
  flye_base = os.path.splitext(os.path.basename(flye_assembly))[0]
  fl_gfa = flye_dir + flye_base + ".gfa"
  use rule flye from assembly_workflow with:
    input:
      reads = reads2assemble,
    output:
      assembly = flye_assembly,
      gfa_plot = report(flye_dir + "assembly_graph.gfa.png",
            caption="../report/flye.rst",
            category = "Evaluate assemblies",
            subcategory = flye_dir),
      gfa = fl_gfa
    params:
      outdir = flye_dir,
      readtype = config["Parameters"]["lr_type"],
      pol_iterations = config["Flye"]["Flye polishing iterations"],
      other_flye_opts = config["Flye"]["options"]
    log:
      flye_dir + "logs/" + str(date) + ".j%j.flye.out",
      flye_dir + "logs/" + str(date) + ".j%j.flye.err"
    benchmark:
      flye_dir + "logs/" + str(date) + ".flye.benchmark.txt",
    conda:
      '../envs/flye2.9.5.yaml'
    threads: config["Flye"]["Flye cores"]

if config["Parameters"]["run_hifiasm"] == True:
  gfas = []
  gfa_plots = []
  for i in hifiasm_assemblies:
    assembled.append(i)
    gfas.append(hifiasm_dir + os.path.splitext(os.path.basename(i))[0] + ".gfa")
    gfa_plots.append(hifiasm_dir + os.path.splitext(os.path.basename(i))[0] + ".gfa.png")

  use rule hifiasm from assembly_workflow with:
    input:
      reads = reads2assemble,
      h1 = hic_pe1 if config["Hifiasm"]["phase"] == True
           else [],
      h2 = hic_pe2 if config["Hifiasm"]["phase"] == True
           else []
    output:
      fastas = hifiasm_assemblies,
      gfas = gfas,
      gfa_plots = gfa_plots
    params:
      outdir = hifiasm_dir,
      out_base = config["Parameters"]["base_name"] + ".hfsm.asm",
      teloseq = " --telo-m " + config["Parameters"]["telo_repeat"] if config["Parameters"]["telo_repeat"]
                else "",
      phase = " --h1 " + hic_pe1 + " --h2 " + hic_pe2 if config["Hifiasm"]["phase"] == True
           else "",
      other_hfsm_opts = config["Hifiasm"]["options"] ,
    log:
      hifiasm_dir + "logs/" + str(date) + ".j%j.hifiasm.out",
      hifiasm_dir + "logs/" + str(date) + ".j%j.hifiasm.err"
    benchmark:
      hifiasm_dir + "logs/" + str(date) + ".hifiasm.benchmark.txt",
    conda:
      '../envs/hifiasm0.24.0-r702.yaml'
    threads: config["Hifiasm"]["Hifiasm cores"]

if config["Parameters"]["run_nextdenovo"] == True:
  assembled.append(nextdenovo_assembly)
  use rule nextdenovo from assembly_workflow with:
    input:
      reads = reads2assemble
    output:
      assembly = nextdenovo_assembly
    params:
      outdir = nextdenovo_dir,
      config = config["Parameters"]["ndconfFile"]
    log:
      nextdenovo_dir + "logs/" + str(date) + ".j%j.nextdenovo.out",
      nextdenovo_dir + "logs/" + str(date) + ".j%j.nextdenovo.err"
    benchmark:
      nextdenovo_dir + "logs/" + str(date) + ".j%j.nextdenovo.benchmark.txt"
    envmodules:
      "NextDenovo/2.5.0"
    threads: config["Nextdenovo"]["Nextdenovo cores"]