from datetime import datetime
import os
import re
import subprocess

module purging_workflow:
  snakefile: "../modules/purging.rules.smk"
module scaffolding_workflow:
  snakefile: "../modules/scaffolding_rules.smk"
module hic_workflow:
    snakefile: "../modules/hic.rules.smk"

##0. Define path for files and variables
input_assemblies = {}
postpolish = []
scripts_dir = config["Inputs"]["scripts_dir"]

for i in config["Inputs"]["Assemblies for postpolishing"]:
  postpolish.append(i)

if config["Parameters"]["run_purgedups"] == True:
  for i in config["Inputs"]["Assemblies for postpolishing"]:
    if config["Hifiasm"]["external_purging"] or not re.search("hfsm", i):
      input_assemblies[working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/"] = i
      base_postpolish = os.path.splitext(os.path.basename(i))[0]
      if not os.path.exists(working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/logs"):
        os.makedirs(working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/logs")
      mappings_dir = working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/mappings"
      if not os.path.exists(mappings_dir):
        os.makedirs(mappings_dir)    
      postpolish.append( working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base_postpolish + ".pgd.fa")
      minimap2[base_postpolish] = i

tigmint_assemblies = {}
if config["Parameters"]["run_tigmint"] == True:
  for i in config["Inputs"]["Assemblies for postpolishing"]: 
    if not re.search("hfsm", i):
      step = config["Inputs"]["Assemblies for postpolishing"][i]
      base = os.path.splitext(os.path.basename(i))[0]
      assembly = i
      if config["Parameters"]["run_purgedups"] == True:
        pstep = step.split('_')[0]
        nstep = pstep.replace('s','')
        cstep = float(nstep) + 1
        step = "s0" + str(cstep) + "_" + step.split('_')[0].replace('s','p')
        assembly =  working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".pgd.fa"

      tigmint_assemblies[working_dir + step + "_Scaffolding_tigmint_ARKS_LINKS/"] = assembly
      base = os.path.splitext(os.path.basename(assembly))[0]
      postpolish.append( working_dir + step + "_Scaffolding_tigmint_ARKS_LINKS/" + base + ".10X.scffs.fa")
      if not os.path.exists(working_dir +step + "_Scaffolding_tigmint_ARKS_LINKS/logs"):
        os.makedirs(working_dir + step + "_Scaffolding_tigmint_ARKS_LINKS/logs")

hic_assemblies = {}
asslength = {}
pretext_files = []
pretext_in = []
tpf_files = []
yahs_mq= config['HiC']["yahs_mq"]
if config['Inputs']['HiC_dir'] and config["Wildcards"]["HiC_wildcards"]:
  hic_out_dir = {}
  if config['HiC']['deepseq'] == False:
    i = config['HiC']['assembly_qc']
    name = os.path.splitext(os.path.basename(i))[0]
    hic_assemblies[name] = i
    hic_out_dir[name] = config['Outputs']['hic_qc_dir']
  elif config["Parameters"]["run_yahs"]:
    if config["Parameters"]["run_tigmint"] == False and config["Parameters"]["run_purgedups"] == False:
      for i in config["Inputs"]["Assemblies for postpolishing"]:
        name = os.path.splitext(os.path.basename(i))[0]
        hic_assemblies[name] = i 
    elif config["Parameters"]["run_tigmint"] == True:
      for i in postpolish:
        if "10X" in i:
          name = os.path.splitext(os.path.basename(i))[0]
          hic_assemblies[name] = i
    elif config["Parameters"]["run_purgedups"] == True:
      for i in postpolish:
        if re.search("pgd", i):
          name = os.path.splitext(os.path.basename(i))[0]
          hic_assemblies[name] = i 
        elif config["Hifiasm"]["external_purging"] == False:
          if re.search("hfsm", i):
            name = os.path.splitext(os.path.basename(i))[0]
            hic_assemblies[name] = os.path.dirname(i) + "/" + name + ".fa"

    for i in config["Inputs"]["Assemblies for postpolishing"]:
      step = config["Inputs"]["Assemblies for postpolishing"][i]
      base = os.path.splitext(os.path.basename(i))[0]
      assembly = i
      pstep = step.split('_')[0]
      nstep = pstep.replace('s','')
      if re.search("hfsm", i):
        if config["Parameters"]["run_purgedups"] == True and config["Hifiasm"]["external_purging"] == True:
          cstep = float(nstep) + 1
          step = "s0" + str(cstep) + "_" + step.split('_')[0].replace('s','p')
          assembly =  working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".pgd.fa"
      else:
        if config["Parameters"]["run_tigmint"] == True and config["Parameters"]["run_purgedups"] == True:
          cstep = float(nstep) + 2
          step = "s0" + str(cstep) + "_" + step.split('_')[0].replace('s','p')
          assembly =  working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".pgd.10X.scffs.fa"
        elif config["Parameters"]["run_purgedups"] == True:
          cstep = float(nstep) + 1
          step = "s0" + str(cstep) + "_" + step.split('_')[0].replace('s','p')
          assembly =  working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".pgd.fa"
        elif config["Parameters"]["run_tigmint"] == True:
          cstep = float(nstep) + 1
          step = "s0" + str(cstep) + "_" + step.split('_')[0].replace('s','p')
          assembly =  working_dir + config["Inputs"]["Assemblies for postpolishing"][i] + "_run_purgedups/" + base + ".10X.scffs.fa"
      name = os.path.splitext(os.path.basename(assembly))[0]
      hic_out_dir[name] = working_dir + step + "_HiC_scaffolding/"
      postpolish.append(hic_out_dir[name] + name + ".yhs_mq" + str(yahs_mq) + ".fa")
      tpf_files.append(hic_out_dir[name] + name + ".yhs_mq" + str(yahs_mq) + ".fa.tpf")  

  
  for i in config["Inputs"]["Curated Assemblies"]:
    assembly = i     
    name = os.path.splitext(os.path.basename(assembly))[0]
    curated_assemblies[name] = config["Inputs"]["Curated Assemblies"][i]
    hic_out_dir[name] = config["Inputs"]["Curated Assemblies"][i] + "/"
    tpf_files.append(hic_out_dir[name] + name + ".fa.tpf")
  
  hic_bams = {}
  pretext_lrmap = {}
  for i in hic_out_dir:
    name = i
    hic_bams[name] = hic_out_dir[i] + "mappings/" + name + ".full_hic.bam"
    pretext_lrmap[name] = hic_out_dir[i] + "mappings/" + name + "_minimap2.bam"
    if not os.path.exists(hic_out_dir[i] + "/logs"):
      os.makedirs(hic_out_dir[i] + "/logs")
      
    if config['Parameters']['run_yahs']:  
    # if config['Parameters']['run_yahs'] and i not in curated_assemblies:
      name = i + ".yhs_mq" + str(yahs_mq)
    # elif i in curated_assemblies:
      # name = i
    hic_assemblies[name] = hic_out_dir[i] + name + ".fa"
    hic_bams[name] = hic_out_dir[i] + "mappings/" + name + ".full_hic.bam"
    pretext_lrmap[name] = hic_out_dir[i] + "mappings/" + name + "_minimap2.bam"
  for i in hic_assemblies:
    asslength[i] = os.path.dirname(hic_assemblies[i]) + "/"+ i + ".genome"
    pretext_in.append(os.path.dirname(hic_assemblies[i]) + "/" + i + ".fa")
    minimap2[i] = os.path.dirname(hic_assemblies[i]) + "/" + i + ".fa"

    for mq in config['HiC']['MQ']:
      if i in hic_out_dir:
        pretext_files.append(hic_out_dir[i] + "in_pretext/" + i + "_mq" + str(mq) + ".extensions.pretext")
        if not os.path.exists(hic_out_dir[i] + "in_pretext/logs"):
          os.makedirs(hic_out_dir[i] + "in_pretext/logs")
      else:
        pretext_files.append(os.path.dirname(hic_assemblies[i]) + "/out_pretext/" + i + "_mq" + str(mq) + ".extensions.pretext")
        if not os.path.exists(os.path.dirname(hic_assemblies[i]) + "/out_pretext/logs"):
          os.makedirs(os.path.dirname(hic_assemblies[i]) + "/out_pretext/logs")

#1- Run purgedups 
if config["Parameters"]["run_purgedups"] == True or config["Hifiasm"]["external_purging"]:
  calcuts = ""
  if config["Purge_dups"]["calcuts_options"] != None:
    calcuts = config["Purge_dups"]["calcuts_options"]

  use rule purge_dups from purging_workflow with:
    input:
      assembly_in = lambda wildcards: input_assemblies[wildcards.dir + "/"],
      reads = reads2assemble,
      mapping = "{dir}/mappings/{base_in}_minimap2.allreads.paf.gz"
    output:
      assembly_out = "{dir}/{base_in}.pgd.fa",
      plot = report ("{dir}/{base_in}.PB.cov.png",
             caption="../report/purgedups.rst",
             category = "Purgedups")
    params:
      scripts_dir = scripts_dir,
      base = "{base_in}",
      dir = "{dir}",
      calcuts_opts = calcuts
    log:
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.purge_dups.out",
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.purge_dups.err"
    benchmark:
      "{dir}/logs/" + str(date) + ".{base_in}.purge_dups.benchmark.txt"
    conda:
      "../envs/purge_dups1.2.6.yaml"
    threads: config["Purge_dups"]["purgedups_cores"]

#2- Run tigmint
if config["Parameters"]["run_tigmint"] == True:
  tigmint_opts = ""
  if config["scaffolding_10X"]["tigmint_options"] != None:
    tigmint_opts = config["scaffolding_10X"]["tigmint_options"]

  if config["Inputs"]["processed_10X"] != None or len(config["Inputs"]["raw_10X"])>0:
    r10X_dir = config["Inputs"]["processed_10X"]
    r10X_reads = r10X_dir + "reads.illumina10X.barcoded.fastq.gz"

  use rule scaffolding_10X from scaffolding_workflow with:
    input:
      assembly = lambda wildcards: tigmint_assemblies[wildcards.dir + "/"],
      reads = r10X_reads
    output:
      scaffolded = "{dir}/{base_in}.10X.scffs.fa"
    params:
      base = "{base_in}",
      dir = "{dir}",
      opts = tigmint_opts,
      scripts = scripts_dir,
      rmcmd =  lambda wildcards: "echo 'Removing bam file'; "  \
            # "rm $base_ass.$base_reads.sortbx.bam;" \
             "echo 'Removing tmp dir'; " \
             "rm -r $TMPDIR;"
             if keepfiles == False else "" 
    log:
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.10X_scaffolding.out",
      "{dir}/logs/" + str(date) + ".j%j.{base_in}.10X_scaffolding.err",
    benchmark:
      "{dir}/logs/" + str(date) + ".{base_in}.10X_scaffolding.benchmark.txt",
    threads: config["scaffolding_10X"]["tigmint_cores"]

#3- Hi-C scaffolding
if config['Inputs']['HiC_dir'] and config["Wildcards"]["HiC_wildcards"]:
  use rule assembly_prepare from hic_workflow with:
    input:
      lambda wildcards: hic_assemblies[wildcards.name]
    output:
      glen = "{dir}/{name}.genome",
      bwa = "{dir}/{name}.fa.bwt.2bit.64",
      faidx = "{dir}/{name}.fa.fai"   
    params:
      scripts = scripts_dir,
      workdir = lambda wildcards: os.path.dirname(hic_assemblies[wildcards.name]),
    log:
      "{dir}/logs/" + str(date) + ".j%j.rule_assembly.prepare.{name}.out",
      "{dir}/logs/"  + str(date) + ".j%j.rule_assembly.prepare.{name}.err"  
    benchmark:
      "{dir}/logs/" + str(date) + ".rule_assembly.prepare.benchmark.{name}.txt"
    threads: 2

  if config['Parameters']['run_yahs']:
    use rule run_yahs from hic_workflow with:
      input:
        mappedsort = "{directory}/in_pretext/mappings/{name}.CM.mq" + str(config['HiC']["yahs_mq"]) + ".bam",
        sla = lambda wildcards: hic_assemblies[wildcards.name],
        index = lambda wildcards: hic_assemblies[wildcards.name] + ".fai"
      output:
        outyahs = "{directory}/{name}.yhs_mq" + str(yahs_mq)+ ".fa",
        agp = "{directory}/{name}.yhs_mq" + str(yahs_mq) + ".agp"  
      params:
        name = "{name}",
        yahsdir = "{directory}/run_yahs/",
        mq = config['HiC']["yahs_mq"],
        contig_ec = " --no-contig-ec" if config['HiC']["yahs_contig_ec"] == False
                    else "",
        yahsopts = config["HiC"]['yahs_opts'],
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_yahs.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_yahs.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_yahs.{name}.benchmark.txt"
      threads:  config['HiC']['yahs_cores']

if config['HiC']['get_pretext']:
  use rule get_tpf from hic_workflow with:
    input:
      fasta = "{directory}/{name}.fa"
    output:
      tpf = "{directory}/{name}.fa.tpf",
    params:
      scripts_dir = scripts_dir,
      dir = lambda wildcards: wildcards.directory
    log:
      "{directory}/logs/" + str(date) + ".j%j.get_tpf.{name}.out",
      "{directory}/logs/" + str(date) + ".j%j.get_tpf.{name}.err"
    benchmark:
      "{directory}/logs/" + str(date) + ".get_tpf.{name}.benchmark.txt"
    threads:  2

  use rule generate_pretext from hic_workflow with:
    input:
      mapbam = "{directory}/mappings/{name}.CM.mq{mq}.bam"
    output:
      pret = "{directory}/{name}_mq{mq}.pretext", 
      hr_pret = "{directory}/{name}_mq{mq}.HR.pretext", 
      ptd = "{directory}/snapshots/PTdone{mq}.{name}.txt"
    wildcard_constraints:
      mq="\d+",
    params:
      scripts_dir = scripts_dir,
      outd = "{directory}",
      name = "{name}",
      mq = "{mq}",
      sort = config['HiC']['sort_pretext']
    log:
      "{directory}/logs/" + str(date) + ".j%j.rule_pretext.{mq}.{name}.out",
      "{directory}/logs/" + str(date) + ".j%j.rule_pretext.{mq}.{name}.err"
    benchmark:
      "{directory}/logs/" + str(date) + ".rule_pretext.{mq}.{name}.txt"
    threads: 2

  use rule add_extensions_pretext from hic_workflow with:
    input:
      tel = lambda wildcards: telo_bgs[wildcards.name],
      gaps = "{directory}/{name}.gaps.bg",
      ontcov = lambda wildcards: "{directory}/{name}.LRcoverage.bg" if len(minimap2) > 0 else "",
      pret = "{directory}/{name}_mq{mq}.pretext",  
      hr_pret = "{directory}/{name}_mq{mq}.HR.pretext",
    output:
      pretext = "{directory}/{name}_mq{mq}.extensions.pretext", 
      hr_pretext = "{directory}/{name}_mq{mq}.extensions.HR.pretext", 
    wildcard_constraints:
      mq="\d+",
    params:
      scripts_dir = scripts_dir,
      outd = "{directory}",
      mq = "{mq}",
    log:
      "{directory}/logs/" + str(date) + ".j%j.rule_add_ext.{mq}.{name}.out",
      "{directory}/logs/" +  str(date) + ".j%j.rule_add_ext.{mq}.{name}.err"
    benchmark:
      "{directory}/logs/" + str(date) + ".rule_add_ext.{mq}.{name}.benchmark.txt"
    threads: 2
