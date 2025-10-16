from datetime import datetime
import os
import re
import subprocess

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module preprocess_workflow:
  snakefile: "../modules/process_reads.rules.smk"

working_dir = config["Outputs"]["base_dir"]
scripts_dir = config["Inputs"]["scripts_dir"]

keepfiles = config["Parameters"]["keep_intermediate"]
base = config["Parameters"]["base_name"]

#0. Obtain inputs

fastqs = {}
ont_reads = ""
ont_fastqs = {}
kraken_ins = {}

if config["Inputs"]["ONT_reads"]:
  ont_reads = config["Inputs"]["ONT_reads"]      
ONT_filtered = config["Inputs"]["ONT_filtered"] 
ont_fastqs["filtered_ont"] = ONT_filtered

if config["Inputs"]["ONT_dir"] and config["Wildcards"]["ONT_wildcards"]:
  ont_dir = config["Inputs"]["ONT_dir"]
  ont_list = config["Wildcards"]["ONT_wildcards"].split(',')
  if ont_reads == "":
    ont_reads = config["Outputs"]["filtlong_dir"] + "reads.ont.fastq.gz"  
  extensions = ["fastq.gz"]
  for i in extensions:
    fastqs["ont."+i] = []
    for file in ont_list:
      fastqs["ont."+i].append(ont_dir + file + "." + i)

if ont_reads != "":
  ont_fastqs["raw_ont"] = ont_reads
  if not os.path.exists(ONT_filtered) and not os.path.exists(os.path.dirname(ONT_filtered) + "/logs"):
    os.makedirs(os.path.dirname(ONT_filtered) + "/logs")

pe1_reads = config["Inputs"]["ILLUMINA_pe1"]
pe2_reads = config["Inputs"]["ILLUMINA_pe2"]

if config["Inputs"]["illumina_dir"] and config["Wildcards"]["illumina_wildcards"]:
  illumina_dir = config["Inputs"]["illumina_dir"]
  illumina_list = config["Wildcards"]["illumina_wildcards"].split(',')

  extensions = ["1.fastq.gz", "2.fastq.gz"]
  illumina_processed = config["Inputs"]["processed_illumina"]
  ill_out_dir = os.path.dirname(os.path.dirname(illumina_processed)) + "/"

  if pe1_reads == None:
    pe1_reads = ill_out_dir + "reads.illumina.1.fastq.gz"
    pe2_reads = ill_out_dir + "reads.illumina.2.fastq.gz"

  if not os.path.exists(ill_out_dir + "logs/"):
    os.makedirs(ill_out_dir + "logs/")
  
  for i in extensions:
    fastqs["illumina." + i] = []
    for file in illumina_list:
      fastqs["illumina." + i].append(illumina_processed + file + ".trimmed." + i)

r10X_reads = config["Inputs"]["ILLUMINA_10X"]
longranger_inputs = {}
longranger_dir = {}
if config["Inputs"]["processed_10X"] != None or len(config["Inputs"]["raw_10X"])>0:
  r10X_list = config["Wildcards"]["10X_wildcards"].split(',')
  r10X_dir = config["Inputs"]["processed_10X"]

  if len(config["Inputs"]["raw_10X"]) >0:
    for d in config["Inputs"]["raw_10X"]:
      l = config["Inputs"]["raw_10X"][d].split(',')
      for i in l:
        #print (d, i)
        longranger_inputs[i] = i
        longranger_dir[i] = d
  extensions = ["barcoded.fastq.gz"]
  r10X_reads = r10X_dir + "reads.illumina10X.barcoded.fastq.gz"
  if not os.path.exists(r10X_dir + "logs/"):
    os.makedirs(r10X_dir + "logs/")
  for i in extensions:
    fastqs["illumina10X." + i] = []
    for file in r10X_list:
      fastqs["illumina10X." + i].append(r10X_dir + file + ".lr." + i)

if config["Parameters"]["run_kraken2"] == True:
  dbname = os.path.basename(config["Kraken2"]["database"])

  kraken_ins[os.path.dirname(ONT_filtered) + "/Kraken/filtered_ont/"+dbname+"/filtlong_"+dbname] = ONT_filtered
  logs = os.path.dirname(ONT_filtered)+"/Kraken/filtered_ont/"+dbname + "/filtlong_"+dbname+"_logs"
  if not os.path.exists(logs):
      os.makedirs(logs)
  if ont_reads != "":
    kraken_ins[os.path.dirname(ONT_filtered) + "/Kraken/raw_ont/"+dbname+"/raw_ont_"+dbname] = ont_reads
    logs = os.path.dirname(ont_reads)+"/Kraken/raw_ont/"+dbname + "/raw_ont_"+dbname+"_logs"
    if not os.path.exists(logs):
      os.makedirs(logs)
  if pe1_reads != None:
    kraken_ins[ill_out_dir + "Kraken/"+dbname+"/illumina_"+dbname] = []
    kraken_ins[ill_out_dir + "Kraken/"+dbname+"/illumina_"+dbname].append(pe1_reads)
    kraken_ins[ill_out_dir + "Kraken/"+dbname+"/illumina_"+dbname].append(pe2_reads)
    logs= ill_out_dir + "Kraken/"+dbname+"/illumina_"+dbname+"_logs"
    if not os.path.exists(logs):
      os.makedirs(logs)
  
  if r10X_reads != None:
    kraken_ins[r10X_dir + "Kraken/"+dbname+"/10X_"+dbname] = []
    kraken_ins[r10X_dir + "Kraken/"+dbname+"/10X_"+dbname].append(pe1_reads)
    kraken_ins[r10X_dir + "Kraken/"+dbname+"/10X_"+dbname].append(pe2_reads)
    logs= r10X_dir + "Kraken/"+dbname+"/10X_"+dbname+"_logs"
    if not os.path.exists(logs):
      os.makedirs(logs)

genomescope_dir = ""
if config["Finalize"]["Merqury db"]:
  genomescope_dir = os.path.dirname(config["Finalize"]["Merqury db"]) + "/genomescope2_k" + str(config["Finalize"]["Meryl K"])
  meryl_loc = os.path.dirname(config["Finalize"]["Merqury db"]) + "/tmp_meryl/"
  if not os.path.exists(config["Finalize"]["Merqury db"]) and not os.path.exists(meryl_loc):
    os.makedirs(meryl_loc)
  reads_loc = {}
  meryl_dbs = []

  if config["Inputs"]["processed_illumina"] != None:
    illumina_processed = config["Inputs"]["processed_illumina"]
    illumina_list = config["Wildcards"]["illumina_wildcards"].split(',')
    extensions = [".trimmed.1.fastq.gz", ".trimmed.2.fastq.gz"]
    for i in illumina_list:
      for e in extensions:
        reads_loc[i + e] = illumina_processed + i + e
        if "illumina" in config["Finalize"]["Meryl reads"]:
            meryl_dbs.append(i + e)
  elif config["Inputs"]["processed_10X"] != None:
    extensions = [".lr.barcoded.fastq.gz"]
    for i in r10X_list:
      for e in extensions:
        reads_loc[i + e] = r10X_dir + i + e
        if "illumina" in config["Finalize"]["Meryl reads"]:
            meryl_dbs.append(i + e)
  else:
    if pe1_reads != None:
      reads_loc[os.path.basename(pe1_reads)] = pe1_reads
      if "illumina" in config["Finalize"]["Meryl reads"]:
          meryl_dbs.append(os.path.basename(pe1_reads))
    if pe2_reads != None:
      reads_loc[os.path.basename(pe2_reads)] = pe2_reads
      if "illumina" in config["Finalize"]["Meryl reads"]:
          meryl_dbs.append(os.path.basename(pe2_reads))
    if r10X_reads != None:
      reads_loc[os.path.basename(r10X_reads)] = r10X_reads
      if "illumina" in config["Finalize"]["Meryl reads"]:
          meryl_dbs.append(os.path.basename(r10X_reads))

  if ONT_filtered:
    if "ont" in config["Finalize"]["Meryl reads"]:
        meryl_dbs.append(os.path.basename(ONT_filtered))
    reads_loc[os.path.basename(ONT_filtered)] = ONT_filtered

extra_filtlong_opts = config["Filtlong"]["options"]
if extra_filtlong_opts == None:
  extra_filtlong_opts = ""

extra_genomescope_opts = config["Parameters"]["genomescope_additional_options"]
if extra_genomescope_opts == None:
  extra_genomescope_opts = ""

if config['Inputs']['HiC_dir'] and config["Wildcards"]["HiC_wildcards"]:
  hic_readsd = config["Inputs"]["HiC_dir"]
  hic_list = config["Wildcards"]["HiC_wildcards"].split(',')

  if len(hic_list) > 1:
    hic_concatd = config['Outputs']['concat_HiC_dir']  
    if not os.path.exists(hic_concatd + "logs/"):
      os.makedirs(hic_concatd + "logs/")
  
    extensions = ["1.fastq.gz", "2.fastq.gz"]
    for i in extensions:
      fastqs["hic."+i] = []
      for file in hic_list:
        fastqs["hic."+i].append(hic_readsd + file + "." + i)  
    hic_pe1 = hic_concatd + "reads.hic.1.fastq.gz"
    hic_pe2 = hic_concatd + "reads.hic.2.fastq.gz"
  else:
    hic_pe1 = hic_readsd + hic_list[0] + ".1.fastq.gz"
    hic_pe2 = hic_readsd + hic_list[0] + ".2.fastq.gz"

#1. Preprocess reads

if config["Inputs"]["processed_illumina"] != None and config["Inputs"]["illumina_dir"] != None:
  use rule trim_galore from preprocess_workflow with:
    input:
      read1 = config["Inputs"]["illumina_dir"] + "{file}.1.fastq.gz",
      read2 = config["Inputs"]["illumina_dir"] + "{file}.2.fastq.gz"
    output:
      trim1 = illumina_processed + "{file}.trimmed.1.fastq.gz",
      trim2 = illumina_processed + "{file}.trimmed.2.fastq.gz",
      report1 = report(illumina_processed + "{file}.1.fastq.gz_trimming_report.txt",
                caption = "../report/trimgalore1.rst",
                category = "Process reads",
                subcategory = "Illumina"),
      report2 = report(illumina_processed + "{file}.2.fastq.gz_trimming_report.txt",
                caption = "../report/trimgalore2.rst",
                category = "Process reads",
                subcategory = "Illumina")
    params:
      outdir = illumina_processed,
      opts = config["Trim_Galore"]["options"]
    log: 
      os.path.dirname(os.path.dirname(illumina_processed)) + "/logs/" + str(date) + ".j%j.{file}.trim_galore.out",
      os.path.dirname(os.path.dirname(illumina_processed)) + "/logs/" + str(date) + ".j%j.{file}.trim_galore.err",
    benchmark:
      os.path.dirname(os.path.dirname(illumina_processed)) + "/logs/" + str(date) + ".{file}.trim_galore.benchmark.txt"
    conda:
      "../envs/trim_galore0.6.10.yaml"
    threads: config["Trim_Galore"]["Trim_Illumina_cores"]

if len(fastqs) > 0:
  use rule concat_reads from preprocess_workflow with:
    input:
      fastqs = lambda wildcards: fastqs[wildcards.ext]
    output:
      final_fastq = "{dir}reads.{ext}"
    log:
      "{dir}logs/" + str(date) + ".j%j.concat.{ext}.out",
      "{dir}logs/" + str(date) + ".j%j.concat.{ext}.err"
    benchmark:
      "{dir}logs/" + str(date) + ".concat.benchmark.{ext}.txt"
    conda:
      '../envs/ass_base.yaml'
    threads: config["Parameters"]["concat_cores"]  

if ont_reads != "":
  use rule filtlong from preprocess_workflow with:
    input:
      reads = ont_reads
    output:
      outreads = ONT_filtered
    params:
      minlen = config["Filtlong"]["Filtlong minlen"],
      min_mean_q = config["Filtlong"]["Filtlong min_mean_q"],
      opts = extra_filtlong_opts
    log:
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".j%j.filtlong.out",
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".j%j.filtlong.err"
    benchmark:
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".filtlong.benchmark.txt"
    conda:
      '../envs/filtlong0.2.1.yaml'
    threads: config["Parameters"]["concat_cores"]  

if len(longranger_inputs) > 0:

  use rule long_ranger from preprocess_workflow with: 
    input: 
      mkfastq_dir =  lambda wildcards: longranger_dir[wildcards.bname]
    output:
      fastq_out = r10X_dir + "{bname}.lr.barcoded.fastq.gz",
      sum_out = r10X_dir + "{bname}.lr.barcoded.summary.csv"
    params:
      path = config["Parameters"]["longranger_path"],
      outdir = r10X_dir,
      sample = lambda wildcards: longranger_inputs[wildcards.bname],
      rmcmd =  lambda wildcards: "echo 'Removing longranger run dir: " + \
               r10X_dir + longranger_inputs[wildcards.bname] + "'; rm -r " +\
               r10X_dir + longranger_inputs[wildcards.bname] if keepfiles == False else "" 
    log:
      r10X_dir + "logs/" + str(date) + ".j%j.longranger.{bname}.out",
      r10X_dir + "logs/" + str(date) + ".j%j.longranger.{bname}.err"
    benchmark:
      r10X_dir + "logs/" + str(date) + ".longranger.{bname}.benchmark.txt",
    threads: config["Parameters"]["longranger_cores"] 

if config["Finalize"]["Merqury db"]:
  if not os.path.exists(config["Finalize"]["Merqury db"]):
    use rule build_meryl_db from preprocess_workflow with:
      input:
        fastq = lambda wildcards: reads_loc[wildcards.db]
      output:
        out_dir = directory(meryl_loc + "{db}.meryl")  
      params:
        kmer = config["Finalize"]["Meryl K"],
      log:
        logs_dir + str(date) + ".j%j.build_meryl.{db}.out",
        logs_dir + str(date) + ".j%j.build_meryl.{db}.err" 
      benchmark:
        logs_dir + str(date) + ".build_meryl.{db}.benchmark.txt"
      conda:
        "../envs/merqury1.3.yaml"
      threads:
        config["Finalize"]["Meryl threads"]

    use rule concat_meryl from preprocess_workflow with:
      input:
        input_run = lambda wildcards: expand(rules.build_meryl_db.output.out_dir, db=meryl_dbs)
      output:
        meryl_all = directory(config["Finalize"]["Merqury db"]),
        histogram = os.path.dirname(config["Finalize"]["Merqury db"]) + "/meryl.hist"
      params:
        kmer = config["Finalize"]["Meryl K"],
      log:
        logs_dir + str(date) + ".j%j.concat_meryl.out",
        logs_dir + str(date) + ".j%j.concat_meryl.err" 
      benchmark:
        logs_dir + str(date) + ".concat_meryl.benchmark.txt" 
      threads:
        config["Finalize"]["Meryl threads"]

  if config["Parameters"]["run_smudgeplot"]:
    use rule smudgeplot from preprocess_workflow with:
      input:
        histogram = os.path.dirname(config["Finalize"]["Merqury db"]) + "/meryl.hist",
        meryl = config["Finalize"]["Merqury db"],
      output:
        plot = smudgeplot_dir + "/smudgeplot_smudgeplot.png"
      params:
        dir = smudgeplot_dir
      conda:
        "../envs/merqury1.3.yaml"
      log:
        logs_dir + str(date) + ".j%j.smudgeplot.out",
        logs_dir + str(date) + ".j%j.smudgeplot.err" 
      benchmark:
        logs_dir + str(date) + ".smudgeplot.benchmark.txt"  

  use rule genomescope2 from preprocess_workflow with:
    input:
      histogram = os.path.dirname(config["Finalize"]["Merqury db"]) + "/meryl.hist"
    output:
      outdir = directory(genomescope_dir),
      summary =  report(genomescope_dir + "/summary.txt",
                 caption = "../report/genomescope.rst",
                 category = "Process reads",
                 subcategory = "Illumina"),
      log_plot =  report (genomescope_dir + "/log_plot.png",
                 caption = "../report/genomescope.rst",
                 category = "Process reads",
                 subcategory = "Illumina"),
      transformed_linear = report (genomescope_dir + "/transformed_linear_plot.png",
                 caption = "../report/genomescope.rst",
                 category = "Process reads",
                 subcategory = "Illumina")
    params:
      ploidy = config["Parameters"]["ploidy"],
      kmer =  config["Finalize"]["Meryl K"],
      opts = extra_genomescope_opts
    conda:
      "../envs/genomescope2.0.yaml"
    log:
      logs_dir + str(date) + ".j%j.genomescope.out",
      logs_dir + str(date) + ".j%j.genomescope.err" 
    benchmark:
      logs_dir + str(date) + ".genomescope.benchmark.txt" 

#2. EVALUATE INPUT READS
if ont_reads != "":
  use rule nanoplot from preprocess_workflow with:
    input:
      fastq = lambda wildcards: ont_fastqs[wildcards.prefix]
    output:
      stats = report(config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/NanoStats.txt",
            caption="../report/nanostats.rst",
            category = "Process reads",
            subcategory = "{prefix}"),
      yield_len = report(config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/Yield_By_Length.png",
                caption="../report/nanostats.rst",
                category = "Process reads",
                subcategory = "{prefix}"),
      read_len = report(config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/WeightedHistogramReadlength.png",
               caption= "../report/nanostats.rst",
               category = "Process reads",
               subcategory = "{prefix}"),
    params:
      outdir = config["Outputs"]["filtlong_dir"] + "nanostats/{prefix}/"
    log:
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".j%j.NanoStats.{prefix}.out",
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".j%j.NanoStats.{prefix}.err"
    benchmark:
      config["Outputs"]["filtlong_dir"] + "logs/" + str(date) + ".NanoStats.{prefix}.benchmark.txt"
    conda:
      '../envs/nanoplot1.40.0.yaml'
    threads: config["Parameters"]["concat_cores"]

if config["Parameters"]["run_kraken2"] == True:
  use rule Kraken2 from preprocess_workflow with: 
    input:
      read = lambda wildcards: kraken_ins[wildcards.base],
      database = config["Kraken2"]["database"],
      kmers = config["Kraken2"]["kmer_dist"]
    output:
      report = report("{base}.kraken2.report.txt",
               caption="../report/kraken.rst",
               category = "Process reads",
               subcategory = "Kraken reports"),
      readsout = "{base}.kraken2.seqs.out",
      abundance = "{base}.bracken_abundance.txt"
    params:
      additional = lambda wildcards: config["Kraken2"]["additional_opts"] + " --paired " \
                   if re.search("illumina", wildcards.base) \ 
                   else config["Kraken2"]["additional_opts"],
      prefix = lambda wildcards: os.path.basename(wildcards.base)
    log:
      "{base}_logs/" + str(date) + ".j%j.kraken.out",
      "{base}_logs/" + str(date) + ".j%j.kraken.err",
    benchmark:
       "{base}_logs/" + str(date) + ".j%j.kraken.benchmark.txt",
    conda:
      '../envs/kraken2.1.2.yaml'
    threads: config["Kraken2"]["threads"]