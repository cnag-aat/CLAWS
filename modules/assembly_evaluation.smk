from datetime import datetime
import os
import re
import subprocess

module eval_workflow:
  snakefile: "../modules/evaluate_assemblies.rules.smk"

##0. Define path for files and variables

inputs = config["Inputs"]["Assemblies for polishing"]

assemblies = []
if config["Finalize"]["final Evaluations"] == True:
  for i in inputs:
    assemblies.append(i)
  for i in polished:
    assemblies.append(i) 
  for i in postpolish:
    assemblies.append(i)
  for i in assembled:
    assemblies.append(i)
  for i in pretext_in:
    b=os.path.splitext(i)[0]
    if not b + ".fasta" in assemblies:
      assemblies.append(i) 

in_files = {}
gfa_files = {}
evals_dir = {}
BuscoSummaries = []
StatsFiles = []
telo_bgs = {}
haps = {}
hap1_files = {}
hap2_files = {}
diploid_fasta = ""
telomeres = []
MerqurySummaries = []
MerquryQV = []
MerquryDups = []
Merqury_in = []
lrtype = config["Parameters"]["lr_type"]

if config["Finalize"]["BUSCO lineage"]:
  buscodb = os.path.basename(config["Finalize"]["BUSCO lineage"])

if config["Parameters"]["run_hifiasm"] == True:
  for file in gfas:
    ass_base = os.path.splitext(os.path.basename(file))[0]
    basedirname = os.path.basename(os.path.dirname(file))
    evalassdir =  eval_dir + basedirname + "/"
    gfa_files[evalassdir + ass_base] = file

if config["Parameters"]["run_flye"] == True:
  ass_base = os.path.splitext(os.path.basename(fl_gfa))[0]
  basedirname = os.path.basename(os.path.dirname(fl_gfa))
  evalassdir =  eval_dir + basedirname + "/"
 #gfa_files[evalassdir + ass_base] = flye_assembly

for file in assemblies:
  ass_base = os.path.splitext(os.path.basename(file))[0]
  basedirname = os.path.basename(os.path.dirname(file))
  if basedirname == "hypo" or basedirname == "nextpolish":
    basedirname = os.path.basename(os.path.dirname(os.path.dirname(file)))
  evalassdir =  eval_dir + basedirname + "/"
  if not os.path.exists(evalassdir + "logs/"):
    os.makedirs(evalassdir + "logs/")
  if not os.path.exists(evalassdir + "stats"):
    os.makedirs(evalassdir + "stats")
  buscodir = evalassdir + "busco/"
  if not os.path.exists(buscodir) and config["Finalize"]["BUSCO lineage"]:
    os.makedirs(buscodir)

  in_files[evalassdir + ass_base] = file
  evals_dir[evalassdir + ass_base] = evalassdir
  StatsFiles.append(evalassdir + "stats/" + ass_base + ".stats.txt")
  if config["Parameters"]["telo_repeat"]:
    telo_bgs[ass_base] = evalassdir + "telomeres/" + ass_base + "." + config["Parameters"]["telo_repeat"] + "_telo.bg"
    telomeres.append(evalassdir + "telomeres/" + ass_base + "." + config["Parameters"]["telo_repeat"] + "_telo.bg")
  if config["Finalize"]["BUSCO lineage"]:
    BuscoSummaries.append(buscodir + ass_base + "." + buscodb + ".short_summary.txt")

  if config["Finalize"]["Merqury db"]:
    if "hap" in file:
      hap_base = ass_base.split(".hap")
      fullbase = hap_base[0] 
      if not re.search("pgd", file) and not re.search("yhs", file):
        ev_dir = evalassdir 
        merqdir = evalassdir + "merqury/" + hap_base[0] + ".haps"
        if evalassdir + hap_base[0] + ".haps" in in_files and not file in in_files[evalassdir + hap_base[0] + ".haps"]:
          in_files[evalassdir + hap_base[0] + ".haps"].append(file)
          haps[merqdir + "/" + hap_base[0] + ".haps"].append(file)
        elif not evalassdir + hap_base[0] + ".haps" in in_files:
          in_files[evalassdir + hap_base[0] + ".haps"] = [file]
          haps[merqdir + "/" + hap_base[0] + ".haps"] = [file]
      else:
        hap_end = hap_base[1].split(".")
        for i in range(1, len(hap_end)):
          fullbase+= "." + hap_end[i]
        step_dir = os.path.basename(evalassdir.rstrip("/"))
        if re.search("^s([0-9]+)\.([0-9]+)_p([0-9]+)", step_dir):
          hap_step = step_dir.split(".")
          tmp_step = hap_step[1].split("_")
          if len(hap_end) == 3:
            tmp2 = hap_step[2].split("_")[1:]
            merqdir = os.path.dirname(evalassdir.rstrip("/")) + "/" + hap_step[0] + "_" + tmp_step[1] + "." + "_".join(tmp2) + "_haps/merqury/" + fullbase + ".haps"
            ev_dir = os.path.dirname(evalassdir.rstrip("/")) + "/" + hap_step[0] + "_" + tmp_step[1] + "." + "_".join(tmp2) + "_haps/"
          else:
            merqdir = os.path.dirname(evalassdir.rstrip("/")) + "/" + hap_step[0] + "_" + tmp_step[1] + "." + hap_step[2] + "_haps/merqury/" + fullbase + ".haps"
            ev_dir = os.path.dirname(evalassdir.rstrip("/")) + "/" + hap_step[0] + "_" + tmp_step[1] + "." + hap_step[2] + "_haps/"
        else:
           merqdir = evalassdir + "merqury/" +  fullbase + ".haps"
           ev_dir = evalassdir.rstrip("/")
        if not os.path.exists(ev_dir + "/logs/"):
          os.makedirs(ev_dir + "/logs/")
        if ev_dir +fullbase + ".haps" in in_files and not file in in_files[ev_dir + fullbase + ".haps"]:
          in_files[ev_dir + fullbase + ".haps"].append(file)
          haps[merqdir + "/" + fullbase + ".haps"].append(file)
        elif not ev_dir + fullbase + ".haps" in in_files:
          in_files[ev_dir + fullbase + ".haps" ] = [file]
          haps[merqdir + "/" + fullbase + ".haps"] = [file]
      if "yhs" in fullbase or ass_base in curated_assemblies:
        if "hap1" in file:
          hap1_files[ev_dir + "diploid/" + fullbase] = file
        elif "hap2" in file:
          hap2_files[ev_dir + "diploid/" + fullbase] = file
        if ev_dir + "diploid/" + fullbase in hap1_files and ev_dir + "diploid/" + fullbase in hap2_files:
          if not os.path.exists(ev_dir + "diploid/logs/"):
            os.makedirs(ev_dir + "diploid/logs")
          diploid_fasta = ev_dir + "diploid/" + fullbase + ".diploid.fa"
          if config["Parameters"]["telo_repeat"]:
            telo_bgs[fullbase + ".diploid"] = ev_dir + "diploid/telomeres/" + fullbase + ".diploid." + config["Parameters"]["telo_repeat"] + "_telo.bg"
            telomeres.append(ev_dir + "diploid/telomeres/" + fullbase + ".diploid." + config["Parameters"]["telo_repeat"] + "_telo.bg")
          in_files[ev_dir + "diploid/" + fullbase + ".diploid" ] = diploid_fasta
          hic_assemblies[fullbase + ".diploid"] = diploid_fasta
          pretext_lrmap[fullbase + ".diploid"] = ev_dir + "diploid/mappings/" + fullbase + ".diploid_minimap2.bam"
          hic_bams[fullbase + ".diploid"] = ev_dir + "diploid/mappings/" + fullbase + ".diploid.full_hic.bam"
          asslength[fullbase + ".diploid"] = ev_dir + "diploid/"+ fullbase + ".diploid.genome"
          pretext_in.append(diploid_fasta)
          minimap2[fullbase + ".diploid"] = diploid_fasta
          tpf_files.append(diploid_fasta + ".tpf")

          # for mq in config['HiC']['MQ']:
          pretext_files.append(ev_dir + "diploid/in_pretext/" + fullbase + ".diploid_mq0.extensions.pretext")
          if not os.path.exists(ev_dir + "diploid/in_pretext/logs"):
            os.makedirs(ev_dir + "diploid/in_pretext/logs")
        if ass_base in curated_assemblies:
          if not os.path.exists(merqdir) and config["Finalize"]["Merqury db"]:
            os.makedirs(merqdir)
          MerqurySummaries.append(merqdir + "/" + fullbase + ".haps.completeness.stats")
          MerquryQV.append(merqdir + "/" + fullbase + ".haps.qv")
          MerquryDups.append(merqdir + "/" + fullbase + ".haps.false_duplications.txt")
    else:
      merqdir = evalassdir + "merqury/" + ass_base 
      if not os.path.exists(merqdir) and config["Finalize"]["Merqury db"]:
        os.makedirs(merqdir)
      MerqurySummaries.append(merqdir + "/" + ass_base + ".completeness.stats")
      MerquryQV.append(merqdir + "/" + ass_base + ".qv")
      MerquryDups.append(merqdir + "/" + ass_base + ".false_duplications.txt")

for base in haps:
  ploi = 0
  for i in haps[base]:
    if "hap1" in i:
      ploi+= 1
    elif "hap2" in i:
      ploi+=1
  if ploi == 2:
    MerqurySummaries.append(base + ".completeness.stats")
    MerquryQV.append(base + ".qv")
    MerquryDups.append(base +  ".false_duplications.txt")
    if not os.path.exists(os.path.dirname(base)):
      os.makedirs(os.path.dirname(base))
  else:
    dir = os.path.dirname(os.path.dirname(base))
    name = os.path.splitext(os.path.basename(file))[0]
    if not os.path.exists(dir + "/" + name):
      os.makedirs(dir + "/" + name)
    MerqurySummaries.append(dir + "/" + name + "/" + name + ".completeness.stats")
    MerquryQV.append(dir + "/" + name + "/" + name + ".qv")
    MerquryDups.append(dir + "/" + name + "/" + name +  ".false_duplications.txt")

#1- Perform alignments
if len(bwa) >0:
  use rule align_illumina from eval_workflow with:
    input:
      genome = lambda wildcards: bwa[wildcards.name],
      reads = [pe1_reads, pe2_reads] if r10X_reads == None else r10X_reads
    output:
      mapping = "{directory}/mappings/{name}_bwa.bam",
      stats = report("{directory}/{name}_bwa.stats.txt",
              caption = "../report/bwa.rst",
              category = "Read Mapping",
              subcategory = "{name}_bwa")
    log:
      "{directory}/logs/" + str(date) + ".j%j.bwa.{name}.out",
      "{directory}/logs/" + str(date) + ".j%j.bwa.{name}.err",
    benchmark:
      "{directory}/logs/" + str(date) + ".bwa.{name}.benchmark.txt",
    conda:
      "../envs/bwa-mem2.2.1.yaml"
    threads: config["Parameters"]["BWA_cores"]

if len(minimap2) > 0:
  use rule align_lr from eval_workflow with:
    input:
      genome = lambda wildcards: minimap2[wildcards.name],
      reads = lambda wildcards: ont_reads if wildcards.ext == "minimap2.allreads.paf.gz" and lrtype == "nano-raw"
                                else reads2assemble,
    output:
      mapping = "{directory}/mappings/{name}_{ext}"
    params:
      align_opts = lambda wildcards:"ax" if wildcards.ext == "minimap2.bam" else "x",
      split = "" if gsize < 4000 else " --split-prefix {name}_tmp ",
      type = lambda wildcards:"map-ont" if lrtype == "nano-raw" else "map-hifi",
      # tmp = "{directory}/mappings/{name}_{ext}.tmp",
      compress_cmd = lambda wildcards : "samtools view -Sb - | " \
                     "samtools sort -@ " + str(config["Parameters"]["minimap2_cores"]) +" -o " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext +";" +\
                     "samtools index -c " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext  \
                     if wildcards.ext == "minimap2.bam" else \
                     "gzip -c > " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext         
    wildcard_constraints:
      ext = "minimap2.(.+)"
    log:
      "{directory}/logs/" + str(date) + ".j%j.{ext}.{name}.out",
      "{directory}/logs/" + str(date) + ".j%j.{ext}.{name}.err",
    benchmark:
      "{directory}/logs/" + str(date) + ".{ext}.{name}.benchmark.txt",
    conda:
      "../envs/minimap2.24.yaml"
    threads: config["Parameters"]["minimap2_cores"]

if len(hic_assemblies) > 0:

  use rule align_hic_chromap from eval_workflow with:
    input:
        ass = lambda wildcards: hic_assemblies[wildcards.name],
        read1 = hic_pe1,
        read2 = hic_pe2
    output:
        mapped = "{directory}/mappings/{name}.CM.mq0.sorted.bam",
        unmapped = "{directory}/mappings/{name}.unmapped_hic.bam"
    params:
        name = '{name}',
        outd = '{directory}/mappings/',
        options = config['HiC']['align_opts'],
        mem = "4G",
        tmpd = "{directory}/mappings/{name}.chromap_tmp_mq0"
    log:
        "{directory}/logs/" + str(date) + ".j%j.rule_chromap.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_chromap.{name}.err"
    benchmark:
        "{directory}/logs/" + str(date) + ".rule_chromap.benchmark.{name}.txt"
    threads: config["Parameters"]["BWA_cores"]
    
  if not config['HiC']['deepseq']: 
    use rule read_screening from eval_workflow with:
      input:
        umapped = "{directory}/mappings/{name}.unmapped_hic.bam",
      output:
        ufasta = "{directory}/blast/unmapped_hic.{name}.{readsblast}.fasta",
        subufasta = "{directory}/blast/unmapped_hic.{name}.{readsblast}_reads.fasta",
        blastoutbscore = "{directory}/blast/unmapped_hic.{name}.{readsblast}_unmapped_reads_vs_nt_25cul1_1e25.megablast.sorted_by_bitscore.out",
        blastoutbscorethits = "{directory}/blast/unmapped_hic.{name}.{readsblast}_unmapped_reads_vs_nt_25cul1_1e25.megablast.sorted_by_bitscore.tophits",
        blastoutorganisms = "{directory}/blast/unmapped_hic.{name}.{readsblast}_unmapped_reads_vs_nt_25cul1_1e25.megablast.organisms.txt" 
      params:
        scripts_dir = scripts_dir,
        outd = "{directory}/blast/",
        readforblast = config['HiC']['reads_for_blast'],
        blastdb = config['HiC']['blastdb']
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_screening.{name}.{readsblast}_reads.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_screening.{name}.{readsblast}_reads.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_screening.{name}.{readsblast}_reads.benchmark.txt"
      threads:  config['HiC']['blast_cores']

    use rule pairtools_processing_parse from eval_workflow with:
      input: 
        mapped = lambda wildcards: hic_bams[wildcards.name],
        alength = lambda wildcards: asslength[wildcards.name]
      output:
        psamo = "{directory}/{name}_mq{mq}_tmp/full_hic.{mq}.{name}.pairsam.gz"
      wildcard_constraints:
        mq="\d+",
      params:
        scripts_dir = scripts_dir,
        mq = '{mq}',
        name = '{name}',
        outd = '{directory}/pairtools_out',
        tmpd = '{directory}/{name}_mq{mq}_tmp',
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools_parse.mq{mq}.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools_parse.mq{mq}.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_pairtools_parse.benchmark.mq{mq}.{name}.txt"
      threads: config['Parameters']['pairtools_parse_cores']

    use rule pairtools_processing_sort from eval_workflow with:
      input: 
        psamo = "{directory}/{name}_mq{mq}_tmp/full_hic.{mq}.{name}.pairsam.gz"
      output:
        spsamo = "{directory}/{name}_mq{mq}_tmp/full_hic.{mq}.{name}.sorted.pairsam.gz"
      wildcard_constraints:
        mq="\d+",
      params:
        scripts_dir = scripts_dir,
        mq = '{mq}',
        name = '{name}',
        outd = '{directory}/pairtools_out',
        tmpd = '{directory}/{name}_mq{mq}_tmp',
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools_sort.mq{mq}.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools_sort.mq{mq}.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_pairtools_sort.benchmark.mq{mq}.{name}.txt"
      threads: config['Parameters']['pairtools_sort_cores']

    use rule pairtools_processing_dedup from eval_workflow with:
      input: 
        spsamo = "{directory}/{name}_mq{mq}_tmp/full_hic.{mq}.{name}.sorted.pairsam.gz"
      output:
        spsamdedup = "{directory}/pairtools_out/full_hic.{mq}.{name}.pairsam.dedupped.gz",
        stats = "{directory}/pairtools_out/stats.mq{mq}.{name}.txt",
      wildcard_constraints:
        mq="\d+",
      params:
        scripts_dir = scripts_dir,
        mq = '{mq}',
        name = '{name}',
        outd = '{directory}/pairtools_out',
        tmpd = '{directory}/{name}_mq{mq}_tmp',
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools_dedup.mq{mq}.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools_dedup.mq{mq}.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_pairtools_dedup.benchmark.mq{mq}.{name}.txt"
      threads: config['Parameters']['pairtools_dedup_cores']

    use rule pairtools_processing_split from eval_workflow with:
      input: 
        spsamdedup = "{directory}/pairtools_out/full_hic.{mq}.{name}.pairsam.dedupped.gz",
      output:
        mappedptsort = "{directory}/pairtools_out/mapped.PT.mq{mq}.{name}.name_sorted.bam",
        mappedpt = "{directory}/pairtools_out/mapped.PT.mq{mq}.{name}.bam",
      wildcard_constraints:
        mq="\d+",
      params:
        scripts_dir = scripts_dir,
        mq = '{mq}',
        name = '{name}',
        outd = '{directory}/pairtools_out',
        tmpd = '{directory}/{name}_mq{mq}_tmp',
        rmcmd = "rm -r {directory}/{name}_mq{mq}_tmp;" if keepfiles == False
               else ""
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools_split.mq{mq}.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_pairtools_split.mq{mq}.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_pairtools_split.benchmark.mq{mq}.{name}.txt"
      threads: config['Parameters']['pairtools_split_cores']

    use rule qc_statistics from eval_workflow with: 
      input:
        statmq = "{directory}/pairtools_out/stats.mq{mq}.{name}.txt",
        mapbam = "{directory}/pairtools_out/mapped.PT.mq{mq}.{name}.bam"
      output:
        libstats = "{directory}/HiC_Final_LibraryStats_mq{mq}.{name}.txt" if config['HiC']['deepseq'] 
                 else "{directory}/HiC_QC_LibraryStats_mq{mq}.{name}.txt",
      wildcard_constraints:
        mq="\d+",
      params:
        scripts_dir = scripts_dir,
        outd = '{directory}/pairtools_out',
        assemblylength = config['HiC']['qc_assemblylen'],
        deepseq = config['HiC']['deepseq'],
        pslibstats = "{directory}/HiC_QC_LibraryStats_extrapolated_mq{mq}.{name}.txt",
        add_preseq_opts = config["HiC"]["add_preseq_opts"]
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_qc_stats.mq{mq}.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_qc_stats.mq{mq}.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_qc_stats.benchmark.mq{mq}.{name}.txt"
      threads: 2

  elif config['HiC']['get_pretext']:
    use rule filter_chromap from eval_workflow with:
      input:
        bam = "{directory}/mappings/{name}.CM.mq0.sorted.bam"
      output:
        filtered = "{directory}/mappings/{name}.CM.mq{mq}.sorted.bam"
      params:
        mq = "{mq}",
        mem = "4G",
      wildcard_constraints:
        mq="(?!0)\d+"  
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_filter_CM.mq{mq}.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_filter_CM.mq{mq}.{name}.err",
      benchmark:
        "{directory}/logs/" + str(date) + ".j%j.rule_filter_CM.benchmark.mq{mq}.{name}.txt",
      threads: config["Parameters"]["BWA_cores"]

    use rule get_extension_gaps from eval_workflow with:
      input:
        sla = lambda wildcards: hic_assemblies[wildcards.name],
      output:
        gaps_bed = "{directory}/{name}.gaps.bed",
        gaps = "{directory}/{name}.gaps.bg"
      params:
        scripts_dir = scripts_dir,
      log:
        "{directory}/logs/" + str(date) + ".j%j.rule_gaps.{name}.out",
        "{directory}/logs/" + str(date) + ".j%j.rule_gaps.{name}.err"
      benchmark:
        "{directory}/logs/" + str(date) + ".rule_gaps.{name}.benchmark.txt"
      threads:  1

    if len(minimap2) > 0:
      use rule get_extension_cov from eval_workflow with:
        input:
          sbam = lambda wildcards: pretext_lrmap[wildcards.name],
        output: 
          ontcov = "{directory}/{name}.LRcoverage.bg"
        log:
          "{directory}/logs/" + str(date) + ".j%j.rule_cov_bed.{name}.out",
          "{directory}/logs/" + str(date) + ".j%j.rule_cov_bed.{name}.err"
        benchmark:
          "{directory}/logs/" + str(date) + ".rule_cov_bed.{name}.benchmark.txt"
        threads:  1
    
  use rule generate_diploid from eval_workflow with:
    input:
      hap1=lambda wildcards: hap1_files[wildcards.directory +"/" + wildcards.name],
      hap2=lambda wildcards: hap2_files[wildcards.directory +"/" + wildcards.name]
    output:
      diploid= "{directory}/{name}.diploid.fa"
    log:
      "{directory}/logs/" + str(date) + ".j%j.rule_get_diploid.{name}.out",
      "{directory}/logs/" + str(date) + ".j%j.rule_get_diploid.{name}.err",
    benchmark:
      "{directory}/logs/" + str(date) + ".rule_get_diploid.{name}.benchmark.txt",

#2- Run evaluations

use rule tidk_search from eval_workflow with:
  input:
    sla = lambda wildcards: in_files[wildcards.directory + "/" + wildcards.name],
  output:
    tel = "{directory}/telomeres/{name}." + config["Parameters"]["telo_repeat"] + "_telo.bg"
  params:
    outd = "{directory}/telomeres/",
    teloseq = config["Parameters"]["telo_repeat"],
    outname = "{name}." + config["Parameters"]["telo_repeat"]
  log:
    "{directory}/logs/" + str(date) + ".j%j.rule_telo_search.{name}.out",
    "{directory}/logs/" + str(date) + ".j%j.rule_telo_search.{name}.err"
  benchmark:
    "{directory}/logs/" + str(date) + ".rule_telo_search.{name}.benchmark.txt"
  threads:  2

use rule get_stats_gfa from eval_workflow with:
  input:
    #assembly_fa =  lambda wildcards: in_files[eval_dir + wildcards.dir + "/" + wildcards.buscobase],
    assembly_gfa = lambda wildcards: gfa_files[eval_dir + wildcards.dir + "/" + wildcards.buscobase]
                   if re.search("hifiasm", wildcards.dir)
                   else in_files[eval_dir + wildcards.dir + "/" + wildcards.buscobase]
  output: 
    # nseries = report(eval_dir + "{dir}/stats/{buscobase}.nseries.txt",
    #           caption="../report/stats.rst",
    #           category = "Evaluate assemblies",
    #           subcategory = "{dir}"),
    stats = report (eval_dir + "{dir}/stats/{buscobase}.stats.txt",
              caption="../report/stats.rst",
              category = "Evaluate assemblies",
              subcategory = "{dir}"),
  params:
    outbase = "{buscobase}",
    scripts_dir = scripts_dir,
    params = lambda wildcards: " --discover-paths " if re.search("hifiasm", wildcards.dir)
                                  else " " 
  log:
    eval_dir + "{dir}/logs/" + str(date) + ".j%j.get_stats.{buscobase}.out",
    eval_dir + "{dir}/logs/" + str(date) + ".j%j.get_stats.{buscobase}.err"
  benchmark:  
    eval_dir + "{dir}/logs/" + str(date) + ".get_stats.{buscobase}.benchmark.out",
  conda:
    '../envs/gfastastats1.3.10.yaml'
  threads: 1

if config["Finalize"]["BUSCO lineage"] != None:
  use rule run_busco from eval_workflow with:
    input:
      assembly = lambda wildcards: in_files[eval_dir + wildcards.dir + "/" + wildcards.buscobase],
      lineage = config["Finalize"]["BUSCO lineage"]
    output:
      summary = report(eval_dir + "{dir}/busco/{buscobase}.{buscodb}.short_summary.txt",
                caption="../report/busco.rst",
                category = "Evaluate assemblies",
                subcategory = "{dir}"),
      full = eval_dir + "{dir}/busco/{buscobase}.{buscodb}.full_table.tsv",
    params:
      out_path = lambda wildcards: evals_dir[eval_dir + wildcards.dir +"/" + wildcards.buscobase],
      odb = os.path.basename(config["Finalize"]["BUSCO lineage"]),
      buscobase = lambda wildcards:  wildcards.buscobase,
      rmcmd =  lambda wildcards: "echo 'Removing BUSCO run dir'; rm -r " + \
             evals_dir[eval_dir + wildcards.dir +"/" + wildcards.buscobase] + wildcards.buscobase + ";" \
             if keepfiles == False else "" 
    log:
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.busco.{buscobase}.{buscodb}.out",
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.busco.{buscobase}.{buscodb}.err",
    benchmark:
      eval_dir + "{dir}/logs/" + str(date) + ".busco.{buscobase}.{buscodb}.benchmark.txt"
    conda:
      '../envs/busco6.0.0.yaml'
    threads: config["Parameters"]["busco_cores"]

if  config["Finalize"]["Merqury db"] != None:
  use rule run_merqury from eval_workflow with:
    input:
      meryl_db = config["Finalize"]["Merqury db"],
      assembly = lambda wildcards: in_files[eval_dir + wildcards.dir +"/" + wildcards.merqbase],
    output:
      completeness = report (eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.completeness.stats",
                    caption="../report/merqury.rst",
                     category = "Evaluate assemblies",
                     subcategory = "{dir}"),
      qv = report(eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.qv",
           caption="../report/merqury.rst",
           category = "Evaluate assemblies",
           subcategory = "{dir}"),
      hist = eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.hist",
      plots = report ([eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.ln.png", eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.fl.png", eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.spectra-cn.st.png"],
              caption="../report/merqury.rst",
              category = "Evaluate assemblies",
              subcategory = "{dir}"),
      false_dups = report(eval_dir + "{dir}/merqury/{merqbase}/{merqbase}.false_duplications.txt",
           caption="../report/merqury.rst",
           category = "Evaluate assemblies",
           subcategory = "{dir}"),
    params:
      out_pref = "{merqbase}",
      directory= eval_dir + "{dir}/merqury/{merqbase}",
      additional_plot_opts = config["Finalize"]["Merqury plot opts"]
    log:
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.merqury.{merqbase}.out",
      eval_dir + "{dir}/logs/" + str(date) + ".j%j.merqury.{merqbase}.err"
    benchmark:
      eval_dir + "{dir}/logs/" + str(date) + ".merqury.{merqbase}.benchmark.txt"
    conda:
      "../envs/merqury1.3.yaml"
    threads: 
      config["Finalize"]["Meryl threads"]

#Run final job of the pipeline
rmcmd = ""
ass_cleand = []
if keepfiles == False:
  rmcmd = "echo 'Pipeline has been completed succesffully, we are now going to delete temporary files:';"
  if config["Inputs"]["processed_illumina"] != None:
    t = config["Inputs"]["processed_illumina"]
    illumina_list = config["Wildcards"]["illumina_wildcards"].split(',')
    if os.path.exists(config["Inputs"]["processed_illumina"] + illumina_list[0] + ".1_val_1.fq.gz"):
      rmcmd += "echo 'Deleting fastqs in " + t + "'; rm " + t + "*.gz;"
  if config["Finalize"]["Merqury db"] != None:
    t = os.path.dirname(config["Finalize"]["Merqury db"]) + "/tmp_meryl/"
    if (os.path.exists(t)):
      rmcmd += "echo 'Deleting " + t + "'; rm -r " + t + ";"
  if diploid_fasta:
     assemblies.append(diploid_fasta)
  for i in assemblies:
    rundir = os.path.dirname(i) + "/"
    if rundir not in ass_cleand:
      ass_cleand.append(rundir) 
      if os.path.exists(rundir + "mappings") and not i in pretext_in:
        rmcmd += "echo 'Deleting mappings in " + rundir + "mappings'; rm -r " + rundir +  "mappings;"
      elif i in pretext_in:
        if os.path.exists(rundir + "in_pretext/pairtools_out"):
          rmcmd += "echo 'Deleting coverage bedgraphs in " + rundir + "in_pretext'; for i in  `find " + rundir + "in_pretext -name '*LRcoverage.bg'`; do rm $i; done;"
          rmcmd += "echo 'Deleting files in " + rundir + "in_pretext/pairtools_out'; rm -r " + rundir +  "in_pretext/pairtools_out;"
        if os.path.exists(rundir + "out_pretext/pairtools_out"):
          rmcmd += "echo 'Deleting coverage bedgraphs in " + rundir + "out_pretext'; for i in  `find " + rundir + "out_pretext -name '*LRcoverage.bg'`; do rm $i; done;"
          rmcmd += "echo 'Deleting files in " + rundir + "out_pretext/pairtools_out'; rm -r " + rundir +  "out_pretext/pairtools_out;"
        if os.path.exists(rundir+"mappings"):
          rmcmd += "echo 'Deleting hic alignments in " + rundir + "mappings'; for i in  `find " + rundir + "mappings -name '*CM*'`; do rm $i; done;"
      if os.path.exists(rundir + "hypo/aux"):
        rmcmd += "echo 'Deleting " + rundir + "hypo/aux;'; rm -r " + rundir + "/hypo/aux;"
  for i in reads_loc:
    if "_bam.fastq" in i:
      if os.path.exists(reads_loc[i]):
        rmcmd += "echo 'Deleting " + reads_loc[i] + ";'; rm -r " + reads_loc[i] + ";"  

use rule finalize from eval_workflow with:
  input:
    assembly = assemblies,
    buscos = BuscoSummaries,
    stats= StatsFiles,
    merqs=MerquryQV,
    telos=telomeres,
    pretext = lambda wildcards: pretext_files if config["HiC"]["deepseq"] == True else [],
    tpf = lambda wildcards: tpf_files if config["HiC"]["deepseq"] == True else []
  output:
    output = report(config["Outputs"]["stats_out"],
             caption="../report/final.rst",
             category = "Final Table"),
  params:
    scripts_dir = scripts_dir,
    rmcmd = rmcmd,
  log:
    logs_dir + str(date) + ".j%j.finalize.out",
    logs_dir + str(date) + ".j%j.finalize.err"
  benchmark:
    logs_dir + str(date) + ".finalize.benchmark.txt"
  conda:
    '../envs/ass_base.yaml'
  threads: 1
  