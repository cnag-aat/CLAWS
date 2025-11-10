#!/usr/bin/env python3
import os
import json
import argparse
import sys
import re
from configparser import ConfigParser
nextdenovo_config = ConfigParser()

##Author: Jessica Gomez-Garrido
##CNAG
##email: jessica.gomez@cnag.eu

def get_wildcards(dir, wildcards, ext):
  for r, d, f in os.walk(dir):
    for file in f:
      if file.endswith(ext):
      #if re.search(ext, file):
        a = file.replace(ext,'')
        if wildcards == None:
          wildcards = a
        else:
          wildcards += "," + a
  return wildcards

def require_ont():
  if args.ONT_filtered:
    args.ONT_filtered = os.path.abspath(args.ONT_filtered)
    args.preprocess_lr = os.path.dirname(args.ONT_filtered) + "/"
  else:
    if args.preprocess_lr == None:
      args.preprocess_lr = args.pipeline_workdir + "s" + args.preprocess_lr_step + "_p01.1_Preprocess_LR/"
    else:
      args.preprocess_lr = os.path.abspath(args.preprocess_lr) + "/" 
    args.ONT_filtered = args.preprocess_lr + "ont.filtlong.fastq.gz"

  if args.ONT_dir == None and args.ONT_reads == None and not os.path.exists(args.ONT_filtered):
    parser.print_help()
    print ("ONT reads are needed")
    sys.exit(-1)
  if args.ONT_dir:
    args.ONT_dir = os.path.abspath(args.ONT_dir) + "/" 
    if not os.path.exists(args.ONT_dir):
      parser.print_help()
      print (args.ONT_dir + " not found")
      sys.exit(-1)
    elif args.ONT_wildcards == None:
      args.ONT_wildcards = get_wildcards(args.ONT_dir, args.ONT_wildcards, '.fastq.gz')

  if args.ONT_reads:
    args.ONT_reads = os.path.abspath(args.ONT_reads) 
    if not os.path.exists(args.ONT_reads):
      if args.ONT_dir == None:
        parser.print_help()
        print (args.ONT_reads + " not found")
        sys.exit(-1)

def require_pb():
  if args.hifi_dir == None and args.hifi_reads == None:
    parser.print_help()
    print ("Pacbio Hifi reads are needed")
    sys.exit(-1)
  if args.hifi_dir:
    args.hifi_dir = os.path.abspath(args.hifi_dir) + "/" 
    if not os.path.exists(args.hifi_dir):
      parser.print_help()
      print (args.hifi_dir + " not found")
      sys.exit(-1)
    elif args.hifi_wildcards == None:
      args.hifi_wildcards = get_wildcards(args.hifi_dir, args.hifi_wildcards, '.fastq.gz')
    
  if args.preprocess_lr == None:
    args.preprocess_lr = args.pipeline_workdir + "s" + args.preprocess_lr_step + "_p01.1_Preprocess_LR/"
  else:
    args.preprocess_lr = os.path.abspath(args.preprocess_lr) + "/" 

  if args.hifi_reads:
    args.hifi_reads = os.path.abspath(args.hifi_reads) 
    if not os.path.exists(args.hifi_reads):
      if args.hifi_dir == None:
        parser.print_help()
        print (args.hifi_reads + " not found")
        sys.exit(-1)

def require_lr():
  if re.search("nano", args.lr_type):
     require_ont()
  else:
     require_pb()

def require_illumina():
  args.r10X_reads = {}
  if args.illumina_dir == None and args.pe1 == None and args.pe2==None and args.r10X==None and args.processed_illumina == None and len(args.raw_10X) == 0 and args.processed_10X == None:
    parser.print_help()
    print ("The illumina reads are needed")
    sys.exit(-1)
  elif args.pe1 or args.pe2:
    if args.pe1:
      args.pe1 = os.path.abspath(args.pe1) 
    else: 
      parser.print_help()
      print ("Both illumina pairs are needed, missing pe1")
      sys.exit(-1)
    if args.pe2:
      args.pe2 = os.path.abspath(args.pe2)      
    else: 
      parser.print_help()
      print ("Both illumina pairs are needed, missing pe2")
      sys.exit(-1)
    
    if not os.path.exists(args.pe1) or not os.path.exists(args.pe2):
      if args.illumina_dir == None and args.processed_illumina == None:
        parser.print_help()
        print ("Illumina reads are not found and are needed")
        sys.exit(-1)

  elif args.r10X:
    args.r10X = os.path.abspath(args.r10X) 
    if not os.path.exists(args.r10X):
      if len(args.raw_10X) == 0 and args.processed_10X == None:
        parser.print_help()
        print ("Illumina reads are not found and are needed")
        sys.exit(-1)
  
  elif args.illumina_dir != None:
    args.illumina_dir = os.path.abspath(args.illumina_dir) + "/" 
    if args.processed_illumina == None:
      args.processed_illumina =  args.pipeline_workdir + "s" + args.preprocess_illumina_step + "_p01.1_preprocess_illumina_reads/trim/"
    else: 
      args.processed_illumina = os.path.abspath(args.processed_illumina) + "/"
    if not os.path.exists(args.illumina_dir):
      parser.print_help()
      print (args.illumina_dir + " not found")
      sys.exit(-1)
    elif args.illumina_wildcards == None:
      args.illumina_wildcards = get_wildcards(args.illumina_dir, args.illumina_wildcards, '1.fastq.gz')
  
  elif args.processed_illumina != None:
    args.processed_illumina = os.path.abspath(args.processed_illumina) + "/"
    if not os.path.exists(args.processed_illumina):
      parser.print_help()
      print (args.processed_illumina + " not found and missing raw illumina directory")
      sys.exit(-1) 
    elif args.illumina_wildcards == None:
      args.illumina_wildcards = get_wildcards(args.processed_illumina, args.illumina_wildcards, '1.fastq.gz')
  else:
    if len(args.raw_10X):
      for my_dict in args.raw_10X:
        for key in my_dict:
          args.r10X_reads[os.path.abspath(key)] = my_dict[key]
          if args.r10X_wildcards:
            args.r10X_wildcards += ","
            args.r10X_wildcards += my_dict[key]
          if not os.path.exists(os.path.abspath(key)):
            parser.print_help()
            print (os.path.abspath(key) + " not found")
            sys.exit(-1)
          elif key == None:
            parser.print_help()
            print ('If you want to process the 10X reads with longranger, you need to provide the 10X basenames together with the directory')
            sys.exit(-1)    
      if args.processed_10X == None:
                args.processed_10X = args.pipeline_workdir + "s" + self.preprocess_10X_step + "_p01.1_preprocess_10X_linkedreads"
    
    if args.processed_10X != None:
      args.processed_10X = os.path.abspath(args.processed_10X) + "/"
      if args.r10X_wildcards == None:
        if not os.path.exists(args.processed_10X):
          parser.print_help()
          print (args.processed_10X + " not found")
          sys.exit(-1)
        else:
          args.r10X_wildcards = get_wildcards(args.processed_10X, args.r10X_wildcards, '.barcoded.fastq.gz') 

class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
      
    def __init__(self):
        """Class constructor"""
        #GENERAL PARAMETERS
        self.configFile = "assembly.config"                                                       #Name of the json configuration file to be created
        self.specFile = "assembly.spec"                                                           #Name of the spec file to be created
        self.ndconfFile = "nextdenovo.config"                                                     #Name of the nextdenovo config file to be created
        self.keepintermediate = False                                                             #Set this to True if you do not want intermediate files to be removed
        self.lr_type  = "nano-hq"                                                                 #Type of long reads (options are flye read-type options)
        self.base_name = None                                                                     #Base name for the project
        self.species = None                                                                       #Name of the species to be assembled
        self.genome_size = None         						                                              #Estimated genome size
        self.ploidy = 2       
        self.telo_string = None
        self.run_flye = True       
        self.run_hifiasm = True
        self.run_nextdenovo = False   
        self.nextpolish_ont_rounds = 0                                                            #Number of rounds for running Nexpolish with ONT
        self.nextpolish_ill_rounds = 0                                                            #Number of rounds for running Nexpolish with illumina
        self.hypo_rounds = 1                                                                      #Number of rounds for running hypo
        self.run_purgedups = True 
        self.run_yahs = True
        self.run_tigmint = False
        self.run_kraken2 = False
        self.run_smudgeplot = True
        self.genomescope_additional = " -m -1 "  
        self.preprocess_lr_step = "02.1"                                                         #Step directory for preprocessing ont
        self.preprocess_10X_step = "02.2"
        self.preprocess_illumina_step = "02.2"
        self.preprocess_hic_step = "02.3"
        self.flye_step = "03.1"                                                                   #Step direcotory for running flye
        self.nextdenovo_step = "03.2"                                                             #Step direcotory for running nextdenovo
        self.hifiasm_step = "03.3"                                                                #Step direcotory for running hifiasm
        self.concat_cores = 4                                                                     #Number of threads to concatenate the reads and to run filtlong
        self.minimap2_cores = 32                                                                  #Number of threads to run the alignment with minimap2
        self.bwa_cores = 16                                                                       #Number of threads to run the alignment for the nextpolish step
        self.nextpolish_cores = 24                                                                #Number of threads to run the nextpolish step
        self.hypo_cores = 24                                                                      #Number of threads to tun the hypo step
        self.pairtools_parse_cores = 32
        self.pairtools_sort_cores = 16
        self.pairtools_dedup_cores = 8
        self.pairtools_split_cores = 16
        self.busco_cores = 32                                                                      #Number of threads to tun the BUSCO    
        self.longranger_cores = 16                                                                 #Number of threads to run longranger   
        self.longranger_path = "/scratch/project/devel/aateam/src/10X/longranger-2.2.2" 

        #ALL SPEC PARAMETERS
        self.all_qos = "test"
        self.all_time = "00:05:00"
        self.all_queue = "general"

        #INPUT PARAMETERS
        self.scripts_dir = os.path.dirname(sys.argv[0]) + "/../scripts/"                          #Directory with the different scripts for the pipeline
        self.ONT_dir = None                                                                       #Directory with the ONT reads, give this option if you don't have a single fastq file with all the reads
        self.hifi_dir = None
        self.illumina_dir = None                                                                  #Directory with the illumina raw reads, give this option if you don't have a single fastq file with all the reads
        self.hic_dir = None
        self.raw_10X = {}                                                                         #List with 10X raw read directories, it has to be the mkfastq dir. You must specify as well the sampleIDs from this run. Example: 
        self.ONT_reads = None                                                                     #File with all the ONT reads  
        self.hifi_reads = None                                                                    #File with all the HiFi reads
        self.pe1 = None                                                                           #File with the illumina paired-end fastqs, already trimeed, pair 1
        self.pe2 = None                                                                           #File with the illumina paired-end fastqs, already trimmed, pair 2
        self.r10X = None                                                                          #File with barcoded 10X reads in fastq.gz format, concatenated
        self.ONT_filtered = None                                                                  #File with the ONT reads after running filtlong         
        self.processed_10X = None                                                                 #Directory to Processed 10X reads, already there or to be produced by the pipeline
        self.processed_illumina = None                                                            #Directory to Processed illumina reads, already there or to be produced by the pipeline
        self.assembly_in = {}                                                                     #List of input assemblies that need to be polished but are not assembled by the pipeline
        self.assemblies = {}
        self.postpolish_assemblies = {}                                                           #List of input assemblies for which postpolishing steps need to be run but are not produced by the pipeline
        self.final_assemblies = {}
        self.curated_assemblies = {}
        self.r10X_reads = {}

        #OUTPUT PARAMETERS
        self.pipeline_workdir = os.getcwd() + "/"                                                 #Base directory for the pipeline run
        self.preprocess_lr = "s" + self.preprocess_lr_step + "_p01.1_Preprocess_LR"                    #Directory to process the ONT reads
        self.concat_hic_dir = "s" + self.preprocess_hic_step + "_p01.1_Concat_HiC"
        self.flye_dir = "s" + self.flye_step + "_p" + self.preprocess_lr_step + "_flye/"         #Directory to run flye 
        self.nextdenovo_dir =  "s" + self.nextdenovo_step + "_p" + self.preprocess_lr_step + "_nextdenovo/"         #Directory to run Nextdenovo 
        self.hifiasm_dir = "s" + self.hifiasm_step + "_p" + self.preprocess_lr_step + "_hifiasm/"         #Directory to run hifiasm 
        self.flye_out = self.flye_dir + "fl.asm.fasta"
        self.nextdenovo_out = self.nextdenovo_dir + "nd.asm.fasta"
        self.hifiasm_out = [self.hifiasm_dir + "hfsm.asm.bp.fa", self.hifiasm_dir + "hfsm.asm.bp.hap1.fa", self.hifiasm_dir + "hfsm.asm.bp.hap2.p_ctg.fa"]
        self.polish_flye_dir = "s04.1_p" + self.flye_step + "_polishing/"                          #Base directory to run polishing pipeline in flye assembly
        self.polish_nextdenovo_dir = "s04.2_p" + self.nextdenovo_step + "_polishing/"              #Base directory to run polishing pipeline in nextdenovo assembly  
        self.eval_dir = "evaluations/"                                                             #Base directory for the evaluations
        self.stats_out = None   
        self.hic_qc_dir = "hic_qc/"    

        #CONCAT READS SPEC PARAMETERS
        self.concat_reads_qos = "normal"
        self.concat_reads_time = "12:00:00"
        self.concat_reads_queue = "general"
        self.concat_reads_mem = "32G"

        #NANOPLOT SPEC PARAMETERS
        self.nanoplot_qos = "normal"
        self.nanoplot_time = "12:00:00"
        self.nanoplot_queue = "general"
        self.nanoplot_mem = "32G"

        #KRAKEN PARAMETERS
        self.kraken2_db = None  
        self.kraken2_kmers = None
        self.kraken2_threads = 16
        self.additional_kraken2_opts = ""

        #KRAKEN SPEC PARAMETERS
        self.kraken2_qos = "vlong"
        self.kraken2_time = "48:00:00"
        self.kraken2_queue = "general"
        self.kraken2_mem = "10G"

        #FILTLONG PARAMETERS
        self.filtlong_minlen = "1000"
        self.filtlong_min_mean_q = "80"
        self.filtlong_opts = None

        #FILTLONG SPEC PARAMETERS
        self.filtlong_qos = "normal"
        self.filtlong_time = "12:00:00"
        self.filtlong_queue = "general"
        self.filtlong_mem = "32G"

        #BAM2FASTQ SPEC PARAMETERS
        self.bam2fq_qos = "normal"
        self.bam2fq_time = "6:00:00"
        self.bam2fq_queue = "general"
        self.bam2fq_mem = "50G"

        #TRIMGALORE PARAMETERS
        self.trim_galore_opts = "--max_n 0 --gzip -q 20 --paired --retain_unpaired"
        self.Trim_Illumina_cores = 8                                                              #Number of threads to run the trim Illumina step

        #TRIMGALORE SPEC PARAMETERS
        self.trimgalore_qos = "normal"
        self.trimgalore_time = "6:00:00"
        self.trimgalore_queue = "general"
        self.trimgalore_mem = "50G"

        #LONGRANGER SPEC PARAMETERS
        self.longranger_qos = "normal"
        self.longranger_time = "12:00:00"
        self.longranger_queue = "general"
        self.longranger_mem = "50G"

        #BUILD MERYL SPEC PARAMETERS
        self.build_meryl_qos = "normal"
        self.build_meryl_time = "12:00:00"
        self.build_meryl_queue = "general"
        self.build_meryl_mem = "50G"

        #CONCAT MERYL SPEC PARAMETERS
        self.concat_meryl_qos = "normal"
        self.concat_meryl_time = "6:00:00"
        self.concat_meryl_queue = "general"
        self.concat_meryl_mem = "10G"

        #SMUDGEPLOT SPEC PARAMETERS
        self.smudgeplot_qos = "marathon_assembly"
        self.smudgeplot_time = "7-00:00:00"
        self.smudgeplot_queue = "general"
        self.smudgeplot_mem = "750G"

        #GENOMESCOPE2 SPEC PARAMETERS
        self.genomescope_qos = "short"
        self.genomescope_time = "1:00:00"
        self.genomescope_queue = "general"
        self.genomescope_mem = "100G"

        #FLYE PARAMETERS
        self.flye_cores = 128	                                                                  #Number of threads to run Flye
        self.flye_pol_it = 2				      			                  #Number of polishing iterations to use with FLYE
        self.other_flye_opts = " --scaffold "                                                     #include here genome size in pipeline											

        #FLYE SPEC PARAMETERS
        self.flye_qos = "marathon"
        self.flye_time = "100:00:00"
        self.flye_queue = "general"
        self.flye_mem = "950G"

        #HIFIASM PARAMETERS
        self.hifiasm_cores = 50	                                                                  #Number of threads to run Flye
        self.other_hifiasm_opts = " --ont "                                                     #include here genome size in pipeline		
        self.hifiasm_ext_purge = False
        self.hifiasm_phasing = False

        #HIFIASM SPEC PARAMETERS
        self.hifiasm_qos = "vlong"
        self.hifiasm_time = "48:00:00"
        self.hifiasm_queue = "general"
        self.hifiasm_mem = "500G"

        #NEXTDENOVO PARAMETERS
        self.nextdenovo_cores = 2	                                                        #Number of threads to run nextdenovo        
        self.nextdenovo_type = "slurm"
        self.nextdenovo_task = "all"
        self.nextdenovo_rewrite = "yes"
        self.nextdenovo_parallel_jobs = 50
        self.nextdenovo_minreadlen = "1k"
        self.nextdenovo_seeddepth = 45
        self.nextdenovo_seedcutoff = 0
        self.nextdenovo_blocksize = "1g"
        self.nextdenovo_pa_correction = 100
        self.nextdenovo_minimap_raw = "-t 30"
        self.nextdenovo_sort = "-m 400g -t 20"
        self.nextdenovo_correction_opts = "-p 30 -dbuf"                      
        self.nextdenovo_minimap_cns = "-t 30 "
        self.nextdenovo_minimap_map = "-t 30 --no-kalloc"              
        self.nextdenovo_nextgraph_opt = "-a 1"

        #NEXTDENOVO SPEC PARAMETERS
        self.nextdenovo_qos = "eternal"
        self.nextdenovo_time = "480:00:00"
        self.nextdenovo_queue = "general"
        self.nextdenovo_mem = "10G"

        #MINIMAP2 SPEC PARAMETERS
        self.minimap_qos = "normal"
        self.minimap_time = "12:00:00"
        self.minimap_queue = "general"
        self.minimap_mem = "500G"
        
        #BWA SPEC PARAMETERS
        self.bwa_qos = "normal"
        self.bwa_time = "6:00:00"
        self.bwa_queue = "general"
        self.bwa_mem = "100G"

        #HYPO PARAMETERS
        self.ill_cov = 0                                                                          #Approximate short read coverage for hypo
        self.hypo_processes = 6                                                                   #Number of contigs to be processed in parallel by hypo
        self.hypo_lr = True                                                                       #Set this to true if you want to run hypo with both long and short reads
        self.hypo_opts = None                                                                     #Extra options to run Hypo 

        #HYPO SPEC PARAMETERS
        self.hypo_qos = "normal"
        self.hypo_time = "6:00:00"
        self.hypo_queue = "general"
        self.hypo_mem = "250G"

        #NEXTPOLISH LR SPEC PARAMETERS
        self.nextpolish_lr_qos = "normal"
        self.nextpolish_lr_time = "6:00:00"
        self.nextpolish_lr_queue = "general"
        self.nextpolish_lr_mem = "100G"
        
        #NEXTPOLISH SR SPEC PARAMETERS
        self.nextpolish_sr_qos = "normal"
        self.nextpolish_sr_time = "6:00:00"
        self.nextpolish_sr_queue = "general"
        self.nextpolish_sr_mem = "200G"

        #PURGEDUPS PARAMETERS
        self.purgedups_cores = 8 
        self.calcuts_opts = None                                                                  #Adjusted values to run calcuts for purgedups

        #PURGEDUPS SPEC PARAMETERS
        self.purgedups_qos = "normal"
        self.purgedups_time = "1:00:00"
        self.purgedups_queue = "general"
        self.purgedups_mem = "200G"

        #10X SCAFFOLDING PARAMETERS
        self.tigmint_opts = None                                                                  #Additional option to give to the 10X scaffolding step
        self.tigmint_cores = 12

        #10X SCAFFOLDING SPEC PARAMETERS
        self.tigmint_qos = "long"
        self.tigmint_time = "24:00:00"
        self.tigmint_queue = "general"
        self.tigmint_mem = "100G"   

        #HiC PARAMETERS
        self.hic_deepseq = True     
        self.get_pretext = True     
        self.subsample_hic = False
        self.add_preseq_opts = " "
        self.sort_pretext = " nosort "                                                           #How to sort pretext file (you can say --nosort or change the --sortby param)
        self.yahs_cores = 48
        self.yahs_mq = 10
        self.yahs_contig_ec = False
        self.yahs_opts = " "
        self.assembly_qc = None   
        self.hic_map_opts = " -5SP -T0 "                                                               #Path to the assembly to be used perfom the QC of the HiC reads
        self.mq = [0,10]                                                                          #Mapping qualities to use for producing the outputs
        self.hic_qc_assemblylen = ""                                                              #Length of the assembly to be used for hic_qc
        self.hic_readsblast = 100                                                                     #Number of unmapped hic reads to blast
        self.blast_cores = 8 
        self.blastdb = "/scratch_isilon/groups/assembly/data/blastdbs"                          #Database to use for running blast against the unmapped hic reads 

        #ASSEMBLY PREPARE SPEC PARAMETERS
        self.ass_prepare_qos = "short"
        self.ass_prepare_time = "2:00:00"
        self.ass_prepare_queue = "general"
        self.ass_prepare_mem = "150G"        

        #MAP HIC SPEC PARAMETERS
        self.map_hic_qos = "normal"
        self.map_hic_time = "12:00:00"
        self.map_hic_queue = "general"
        self.map_hic_mem = "100G"  

        #PAIRTOOLS SPEC PARAMETERS
        self.pairtools_qos = "normal"
        self.pairtools_time = "12:00:00"
        self.pairtools_queue = "general"
        self.pairtools_mem = "200G"  

        #QCSTATS SPEC PARAMETERS
        self.qcstats_qos = "short"
        self.qcstats_time = "3:00:00"
        self.qcstats_queue = "general"
        self.qcstats_mem = "50G" 

        #BLAST SPEC PARAMETERS
        self.blast_qos = "short"
        self.blast_time = "3:00:00"
        self.blast_queue = "general"
        self.blast_mem = "50G" 

        #YAHS SPEC PARAMETERS
        self.yahs_qos = "normal"
        self.yahs_time = "12:00:00"
        self.yahs_queue = "general"
        self.yahs_mem = "50G" 

        #PRETEXT SPEC PARAMETERS
        self.pretext_qos = "normal"
        self.pretext_time = "10:00:00"
        self.pretext_queue = "general"
        self.pretext_mem = "200G" 

        #TPF SPEC PARAMETERS
        self.tpf_qos = "short"
        self.tpf_time = "1:00:00"
        self.tpf_queue = "general"
        self.tpf_mem = "50G" 

        #TELOMERE_EXT SPEC PARAMETERS
        self.telext_qos = "normal"
        self.telext_time = "10:00:00"
        self.telext_queue = "general"
        self.telext_mem = "100G"

       #STATS SPEC PARAMETERS
        self.stats_qos = "vshort"
        self.stats_time = "01:00:00"
        self.stats_queue = "general"  
        self.stats_mem = "16G"
        
        #BUSCO SPEC PARAMETERS
        self.busco_qos = "short"
        self.busco_time = "6:00:00"
        self.busco_queue = "general"  
        self.busco_mem = "50G"

        #MERQURY SPEC PARAMETERS
        self.merq_qos = "normal"
        self.merq_time = "3:00:00"
        self.merq_queue = "general"
        self.merq_mem = "100G"         

        #FINALIZE PARAMETERS
        self.final_evals = True                                                                  #Set this to true if you want evaluations to be run on each of the final assemblies     
        self.busco_lineage = None                                                                 #Busco lineage to be used
        self.merqury_db = None
        self.merqury_plot_opts = None
        self.meryl_k = None
        self.meryl_threads = 4
        self.meryl_reads = "ont illumina"

        #FINALIZE SPEC PARAMETERS
        self.fin_qos = "short"
        self.fin_time = "2:00:00"
        self.fin_queue = "general"
        self.fin_mem = "5G" 

        #WILDCARDS
        self.ONT_wildcards = None
        self.hifi_wildcards = None
        self.illumina_wildcards = None
        self.r10X_wildcards = None                                                                #For raw 10X we need to give this argument, for processed 10X reads, the pipeline can obtain it
        self.hic_wildcards = None
   
###
        #DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.allSpecParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.concatreadsSpecParameters = {}
        self.subsampleHICSpecParameters = {}
        self.nanoplotSpecParameters = {}
        self.kraken2Parameters = {}
        self.kraken2SpecParameters = {}
        self.filtlongParameters = {}
        self.filtlongSpecParameters = {}
        self.bam2fqSpecParameters = {}
        self.trimgaloreParameters = {}
        self.trimgaloreSpecParameters = {}
        self.longrangerSpecParameters = {}
        self.buildmeryl1SpecParameters = {}
        self.buildmerylSpecParameters = {}
        self.concatmerylSpecParameters = {}
        self.genomescopeSpecParameters = {}
        self.smudgeplotSpecParameters = {}
        self.flyeParameters = {}
        self.flyeSpecParameters = {}
        self.hifiasmParameters = {}
        self.hifiasmSpecParameters = {}
        self.nextdenovoParameters = {}
        self.nextdenovoSpecParameters = {}
        self.nextdenovoConfig = {}
        self.minimapSpecParameters = {}
        self.bwaSpecParameters= {}
        self.hypoParameters = {}
        self.hypoSpecParameters = {}
        self.nextpolishlrSpecParameters = {}
        self.nextpolishsrSpecParameters = {}
        self.purgedupsParameters = {}
        self.purgedupsSpecParameters = {}
        self.scaffold10XParameters = {}
        self.scaffold10XSpecParameters = {}
        self.hicParameters = {}
        self.assprepSpecParameters = {}
        self.mapHicSpecParameters = {}
        self.pairtoolsParseSpecParameters = {}
        self.pairtoolsSortSpecParameters = {}
        self.pairtoolsDedupSpecParameters = {}
        self.pairtoolsSplitSpecParameters = {}
        self.blastSpecParameters = {}
        self.qcstatsSpecParameters = {}
        self.yahsSpecParameters = {}
        self.pretextSpecParameters = {}
        self.tpfSpecParameters = {}      
        self.diploidSpecParameters = {}
        self.epretextSpecParameters = {}
        self.tidkSpecParameters = {}
        self.telextSpecParameters = {}
        self.gapsSpecParameters = {}
        self.ontbgSpecParameters = {}
        self.statsSpecParameters = {}
        self.buscoSpecParameters = {}
        self.merqSpecParameters = {}
        self.finalizeParameters = {}
        self.finSpecParameters = {}
        self.wildcardParameters = {}

####

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_filtlong(parser)
        self.register_trimgalore(parser)
        self.register_kraken2(parser)
        self.register_flye(parser)
        self.register_hifiasm(parser)
        self.register_nextdenovo(parser)
        self.register_hypo(parser)
        self.register_purgedups(parser)
        self.register_scaffold10X(parser)        
        self.register_hic(parser)
        self.register_finalize(parser)
        self.register_wildcards(parser)

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('General Parameters')
        general_group.add_argument('--configFile', dest="configFile", metavar="configFile", help='Configuration JSON to be generated. Default %s' % self.configFile)
        general_group.add_argument('--specFile', dest="specFile", metavar="specFile", help='Cluster specifications JSON  fileto be generated. Default %s' % self.specFile)
        general_group.add_argument('--ndconfFile', dest="ndconfFile", metavar="ndconfFile", help='Name pf the nextdenovo config file. Default %s' % self.ndconfFile)
        general_group.add_argument('--keep-intermediate', dest="keepintermediate", action="store_true", help='Set this to True if you do not want intermediate files to be removed. Default %s' % self.keepintermediate)
        general_group.add_argument('--lr-type', dest="lr_type", metavar="lr_type", default=self.lr_type, choices=['pacbio-raw', 'pacbio-corr', 'pacbio-hifi', 'nano-raw', 'nano-corr', 'nano-hq'],  help='Type of long reads (options are: pacbio-raw, pacbio-corr, pacbio-hifi, nano-raw, nano-corr, nano-hq). Default %s' % self.lr_type)
        general_group.add_argument('--basename', dest="base_name", metavar="base_name", help='Base name for the project. Default %s' % self.base_name)
        general_group.add_argument('--species', dest="species", metavar="species", help='Name of the species to be assembled. Default %s' % self.species)
        general_group.add_argument('--genome-size', dest="genome_size", metavar="genome_size", help='Approximate genome size. Example: 615m or 2.6g. Default %s' % self.genome_size)
        general_group.add_argument('--ploidy', type = int, dest="ploidy", metavar="ploidy", default=self.ploidy, help='Expected ploidy. Default %s' % self.ploidy) 
        general_group.add_argument('--telo', dest="telo_string", metavar="telo_string", help='Expected telomere string. Default %s' % self.telo_string)
        general_group.add_argument('--no-flye', dest="run_flye", action="store_false", help='Give this option if you do not want to run Flye.')
        general_group.add_argument('--no-hifiasm', dest="run_hifiasm", action="store_false", help='Give this option if you do not want to run Hifiasm.')
        general_group.add_argument('--run-nextdenovo', dest="run_nextdenovo", action="store_true", help='Give this option if you do want to run Nextdenovo.')
        general_group.add_argument('--nextpolish-ont-rounds', type = int, dest="nextpolish_ont_rounds", metavar="nextpolish_ont_rounds", default=self.nextpolish_ont_rounds, help='Number of rounds to run the Nextpolish with ONT step. Default %s' % self.nextpolish_ont_rounds)
        general_group.add_argument('--nextpolish-ill-rounds', type = int, dest="nextpolish_ill_rounds", metavar="nextpolish_ill_rounds", default=self.nextpolish_ill_rounds, help='Number of rounds to run the Nextpolish with illumina step. Default %s' % self.nextpolish_ill_rounds)
        general_group.add_argument('--hypo-rounds', type = int, dest="hypo_rounds", metavar="hypo_rounds", default=self.hypo_rounds, help='Number of rounds to run the Hypostep. Default %s' % self.hypo_rounds)
        general_group.add_argument('--no-purgedups', dest="run_purgedups", action="store_false", help='Give this option if you do not want to run Purgedups on the Flye and Nextdenovo assemblies.')
        general_group.add_argument('--no-yahs', dest="run_yahs", action="store_false", help='Give this option if you do not want to run yahs.')
        general_group.add_argument('--no-smudgeplot', dest="run_smudgeplot", action="store_false", help='Give this option if you do not want to run smudgeplot.')
        general_group.add_argument('--run-tigmint', dest="run_tigmint", action="store_true", help='Give this option if you want to run the scaffolding with 10X reads step.')
        general_group.add_argument('--run-kraken2', dest="run_kraken2", action="store_true", help='Give this option if you want to run Kraken2 on the input reads.')
        general_group.add_argument('--genomescope-opts', dest="genomescope_additional", metavar="genomescope_additional", default=self.genomescope_additional, help='Additional options to run Genomescope2 with. Default %s' % self.genomescope_additional)
        general_group.add_argument('--preprocess-lr-step', dest="preprocess_lr_step", default=self.preprocess_lr_step, help='Step for preprocessing long-reads. Default %s' % self.preprocess_lr_step)
        general_group.add_argument('--preprocess-10X-step', dest="preprocess_10X_step", default=self.preprocess_10X_step, help='Step for preprocessing 10X reads. Default %s' % self.preprocess_10X_step)
        general_group.add_argument('--preprocess-illumina-step', dest="preprocess_illumina_step", default=self.preprocess_illumina_step, help='Step for preprocessing illumina reads. Default %s' % self.preprocess_illumina_step)
        general_group.add_argument('--preprocess-hic-step', dest="preprocess_hic_step", default=self.preprocess_hic_step, help='Step for preprocessing hic reads. Default %s' % self.preprocess_hic_step)
        general_group.add_argument('--flye-step', dest="flye_step", default=self.flye_step, help='Step for running flye. Default %s' % self.flye_step)
        general_group.add_argument('--hifiasm-step', dest="hifiasm_step", default=self.hifiasm_step, help='Step for running hifiasm. Default %s' % self.hifiasm_step)        
        general_group.add_argument('--nextdenovo-step', dest="nextdenovo_step", default=self.nextdenovo_step, help='Step for running nextdenovo. Default %s' % self.nextdenovo_step)
        general_group.add_argument('--concat-cores', type = int, dest="concat_cores", metavar="concat_cores", default=self.concat_cores, help='Number of threads to concatenate reads and to run filtlong. Default %s' % self.concat_cores)
        general_group.add_argument('--minimap2-cores', type = int, dest="minimap2_cores", metavar="minimap2_cores", default=self.minimap2_cores, help='Number of threads to run the alignment with minimap2. Default %s' % self.minimap2_cores)
        general_group.add_argument('--bwa-cores', type = int, dest="bwa_cores", metavar="bwa_cores", default=self.bwa_cores, help='Number of threads to run the alignments with BWA-Mem2. Default %s' % self.bwa_cores)
        general_group.add_argument('--hypo-cores', type = int, dest="hypo_cores", metavar="hypo_cores", default=self.hypo_cores, help='Number of threads to run the hypo step. Default %s' % self.hypo_cores)
        general_group.add_argument('--nextpolish-cores', type = int, dest="nextpolish_cores", metavar="nextpolish_cores", default=self.nextpolish_cores, help='Number of threads to run the nextpolish step. Default %s' % self.nextpolish_cores)
        general_group.add_argument('--pairtools-parse-cores', type = int, dest="pairtools_parse_cores", metavar="pairtools_parse_cores", default=self.pairtools_parse_cores, help='Number of threads to run the pairtools parse step. Default %s' % self.pairtools_parse_cores)
        general_group.add_argument('--pairtools-sort-cores', type = int, dest="pairtools_sort_cores", metavar="pairtools_sort_cores", default=self.pairtools_sort_cores, help='Number of threads to run the pairtools sort step. Default %s' % self.pairtools_sort_cores)
        general_group.add_argument('--pairtools-dedup-cores', type = int, dest="pairtools_dedup_cores", metavar="pairtools_dedup_cores", default=self.pairtools_dedup_cores, help='Number of threads to run the pairtools dedup step. Default %s' % self.pairtools_dedup_cores)
        general_group.add_argument('--pairtools-split-cores', type = int, dest="pairtools_split_cores", metavar="pairtools_split_cores", default=self.pairtools_split_cores, help='Number of threads to run the pairtools split step. Default %s' % self.pairtools_split_cores)
        general_group.add_argument('--busco-cores', type = int, dest="busco_cores", metavar="busco_cores", default=self.busco_cores, help='Number of threads to run BUSCO. Default %s' % self.busco_cores)
        general_group.add_argument('--longranger-cores', type = int, dest="longranger_cores", metavar="longranger_cores", default=self.longranger_cores, help='Number of threads to run longranger. Default %s' % self.longranger_cores)
        general_group.add_argument('--longranger-path', dest="longranger_path", metavar="longranger_path", default=self.longranger_path, help='Path to longranger executable. Default %s' % self.longranger_path)

    def register_input(self, parser):
        """Register all input parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Inputs')
        input_group.add_argument('--scripts-dir', dest="scripts_dir", help='Directory with the different scripts for the pipeline. Default %s' % self.scripts_dir)
        input_group.add_argument('--ont-dir', dest="ONT_dir", help='Directory where the ONT fastqs are stored. Default %s' % self.ONT_dir)
        input_group.add_argument('--hifi-dir', dest="hifi_dir", help='Directory where the hifi reads are stored. Default %s' % self.hifi_dir)
        input_group.add_argument('--illumina-dir', dest="illumina_dir", help='Directory where the raw illumina fastqs are stored. Default %s' % self.illumina_dir)
        input_group.add_argument('--hic-dir', dest="hic_dir", help='Directory where the HiC fastqs are stored. Default %s' % self.hic_dir)
        input_group.add_argument('--raw-10X', dest="raw_10X",  nargs="+", type=json.loads, default=self.raw_10X, help='Dictionary with 10X raw read directories, it has to be the mkfastq dir. You must specify as well the sampleIDs from this run. Example: \'{\"mkfastq-dir\":\"sample1,sample2,sample3\"}\'...')
        input_group.add_argument('--ont-reads', dest="ONT_reads", help='File with all the ONT reads. It can either be in fastq, fasta or bam format Default %s' % self.ONT_reads)
        input_group.add_argument('--hifi-reads', dest="hifi_reads", help='File with all the HiFi reads It can be either in fastq, fasta or bam format. Default %s' % self.hifi_reads)
        input_group.add_argument('--pe1', dest="pe1", help='File with the illumina paired-end fastqs, already trimmed,  pair 1.')
        input_group.add_argument('--pe2', dest="pe2", help='File with the illumina paired-end fastqs, already trimmed, pair 2.')
        input_group.add_argument('--10X', dest="r10X", help='File with barcoded 10X reads in fastq.gz format, concatenated.')
        input_group.add_argument('--processed-illumina', dest="processed_illumina", help='Directory to Processed illumina reads. Already there or to be produced by the pipeline.')
        input_group.add_argument('--processed-10X', dest="processed_10X", help='Directory to Processed 10X reads. Already there or to be produced by the pipeline.')
        input_group.add_argument('--ont-filt', dest="ONT_filtered", help='File with the ONT reads after running filtlong on them. Default %s' % self.ONT_filtered)
        input_group.add_argument('--assembly-in', dest="assembly_in", nargs="+", type=json.loads, default=self.assembly_in, help='Dictionary with assemblies that need to be polished but not assembled and directory where they should be polished. Example: \'{\"assembly1\":\"polishing_dir1\"}\' \'{\"assembly2\"=\"polishing_dir2\"}\' ...')
        input_group.add_argument('--postpolish-assemblies', dest="postpolish_assemblies", nargs="+", type=json.loads, default=self.postpolish_assemblies, help='Dictionary with assemblies for which postpolishing steps need to be run but that are not assembled and base step for the directory where the first postpolishing step should be run. Example: \'{\"assembly1\":\"s04.1_p03.1\"}\' \'{\"assembly2\":\"s04.2_p03.2\"}\' ...')
        input_group.add_argument('--curated-assemblies', dest="curated_assemblies", nargs="+", type=json.loads, default=self.curated_assemblies, help='Dictionary with assemblies that have already been curated and directory where read alignment should be run. Evaluations and read alignment will be performed. Example: \'{\"assembly1\":\"s04.1_p03.1\"}\' \'{\"assembly2\":\"s04.2_p03.2\"}\' ...')

    def register_output(self, parser):
        """Register all output parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        output_group = parser.add_argument_group('Outputs')
        output_group.add_argument('--pipeline-workdir', dest="pipeline_workdir", help='Base directory for the pipeline run. Default %s' % self.pipeline_workdir)
        output_group.add_argument('--preprocess-lr', dest="preprocess_lr",  help='Directory to process the long-reads. Default %s' % self.preprocess_lr)
        output_group.add_argument('--concat-hic-dir', dest="concat_hic_dir",  help='Directory to concatenate the HiC reads. Default %s' % self.concat_hic_dir)
        output_group.add_argument('--flye-dir', dest="flye_dir",  help='Directory to run flye. Default %s' % self.flye_dir)
        output_group.add_argument('--nextdenovo-dir', dest="nextdenovo_dir",  help='Directory to run nextdenovo. Default %s' % self.nextdenovo_dir)
        output_group.add_argument('--hifiasm-dir', dest="hifiasm_dir",  help='Directory to run hifiasm. Default %s' % self.hifiasm_dir)
        output_group.add_argument('--flye-polishing-dir', dest="polish_flye_dir",  help='Directory to polish the flye assembly. Default %s' % self.polish_flye_dir)
        output_group.add_argument('--nextdenovo-polishing-dir', dest="polish_nextdenovo_dir",  help='Directory to run nextdenovo. Default %s' % self.polish_nextdenovo_dir)
        output_group.add_argument('--eval-dir', dest="eval_dir", metavar="eval_dir",  help='Base directory for the evaluations. Default %s' %self.eval_dir)
        output_group.add_argument('--stats-out', dest="stats_out", metavar="stats_out",  help='Path to the file with the final statistics.')
        output_group.add_argument('--hic-qc-dir', dest="hic_qc_dir", metavar="hic_qc_dir",  help='Directory to run the hic_qc. Default %s' %self.hic_qc_dir)

    def register_filtlong(self, parser):
        """Register all filtlong parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        filtlong_group = parser.add_argument_group('Filtlong')
        filtlong_group.add_argument('--filtlong-minlen', dest="filtlong_minlen", metavar="filtlong_minlen", default=self.filtlong_minlen, type = int, help='Minimum read length to use with Filtlong. Default %s' % self.filtlong_minlen)
        filtlong_group.add_argument('--filtlong-min-mean-q', dest="filtlong_min_mean_q", metavar="filtlong_min_mean_q", default=self.filtlong_min_mean_q, type = int, help='Minimum mean quality to use with Filtlong. Default %s' % self.filtlong_min_mean_q)
        filtlong_group.add_argument('--filtlong-opts', dest="filtlong_opts", metavar="filtlong_opts", default=self.filtlong_opts, help='Extra options to run Filtlong (eg. -t 4000000000)')

    def register_trimgalore(self, parser):
        """Register all trimgalore parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        trimgalore_group = parser.add_argument_group('Trim_Galore')
        trimgalore_group.add_argument('--trim-galore-opts', dest="trim_galore_opts", metavar="trim_galore_opts", default=self.trim_galore_opts, help='Optional parameters for the rule trim_galore. Default %s' % self.trim_galore_opts)
        trimgalore_group.add_argument('--trim-Illumina-cores', type = int, dest="Trim_Illumina_cores", metavar="Trim_Illumina_cores", default=self.Trim_Illumina_cores, help='Number of threads to run the Illumina trimming step. Default %s' % self.Trim_Illumina_cores)

    def register_kraken2(self, parser):
        """Register all Kraken2 parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        kraken2_group = parser.add_argument_group('Kraken2')
        kraken2_group.add_argument('--kraken2-db', dest="kraken2_db", metavar="kraken2_db", default=self.kraken2_db, help='Database to be used for running Kraken2. Default %s' % self.kraken2_db)
        kraken2_group.add_argument('--kraken2-kmer', dest="kraken2_kmers", metavar="kraken2_kmers", default=self.kraken2_kmers, help='Database to be used for running Kraken2. Default %s' % self.kraken2_kmers)
        kraken2_group.add_argument('--kraken2-opts', dest="additional_kraken2_opts", metavar="additional_kraken2_opts", default=self.additional_kraken2_opts, help='Optional parameters for the rule Kraken2. Default %s' % self.additional_kraken2_opts)
        kraken2_group.add_argument('--kraken2-cores', type = int, dest="kraken2_threads", metavar="kraken2_threads", default=self.kraken2_threads, help='Number of threads to run the Kraken2 step. Default %s' % self.kraken2_threads)

    def register_flye(self, parser):
        """Register all flye parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        flye_group = parser.add_argument_group('Flye')
        flye_group.add_argument('--flye-cores', dest="flye_cores", metavar="flye_cores", default=self.flye_cores, type = int, help='Number of threads to run FLYE. Default %s' % self.flye_cores)
        flye_group.add_argument('--flye-polishing-iterations', dest="flye_pol_it", metavar="flye_pol_it", default=self.flye_pol_it, type = int, help='Number of polishing iterations to use with FLYE. Default %s' % self.flye_pol_it)
        flye_group.add_argument('--other-flye-opts', dest="other_flye_opts", metavar="other_flye_opts", default=self.other_flye_opts, help='Additional options to run Flye. Default %s' % self.other_flye_opts)

    def register_hifiasm(self, parser):
        """Register all hifiasm parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        hifiasm_group = parser.add_argument_group('Hifiasm')
        hifiasm_group.add_argument('--hifiasm-cores', dest="hifiasm_cores", metavar="hifiasm", default=self.hifiasm_cores, type = int, help='Number of threads to run Hifiasm. Default %s' % self.hifiasm_cores)
        hifiasm_group.add_argument('--other-hifiasm-opts', dest="other_hifiasm_opts", metavar="other_hifiasm_opts", default=self.other_hifiasm_opts, help='Additional options to run Hifiasm. Default %s' % self.other_hifiasm_opts)
        hifiasm_group.add_argument('--purge-hifiasm', dest="hifiasm_ext_purge", default = self.hifiasm_ext_purge, action="store_true", help='Give this option if you want to run purgedups externally on the hifiasm output.')
        hifiasm_group.add_argument('--phase-hifiasm', dest="hifiasm_phasing", default = self.hifiasm_phasing, action="store_true", help='Give this option if you want to phase with hic reads the hifiasm assembly')

    def register_nextdenovo(self, parser):
        """Register all nextdenovo parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        nextdenovo_group = parser.add_argument_group('Nextdenovo')
        nextdenovo_group.add_argument('--nextdenovo-cores', dest="nextdenovo_cores", metavar="nextdenovo_cores", default=self.nextdenovo_cores, type = int, help='Number of threads to run nextdenovo. Default %s' % self.nextdenovo_cores)
        nextdenovo_group.add_argument('--nextdenovo-jobtype', dest="nextdenovo_type", metavar="nextdenovo_type", choices=['local', 'slurm'], default=self.nextdenovo_type, help='Job_type for nextdenovo. Default %s' % self.nextdenovo_type)
        nextdenovo_group.add_argument('--nextdenovo-task', dest="nextdenovo_task", metavar="nextdenovo_task", choices=['all', 'correct', 'assemble'], default=self.nextdenovo_task, help='Task need to run. Default %s' % self.nextdenovo_task)
        nextdenovo_group.add_argument('--nextdenovo-rewrite', dest="nextdenovo_rewrite", metavar="nextdenovo_rewrite", choices=['yes', 'no'], default=self.nextdenovo_rewrite, help='Overwrite existing directory. Default %s' % self.nextdenovo_rewrite)
        nextdenovo_group.add_argument('--nextdenovo-parallel_jobs', dest="nextdenovo_parallel_jobs", metavar="nextdenovo_parallel_jobs", default=self.nextdenovo_parallel_jobs, type = int, help='Number of tasks used to run in parallel. Default %s' % self.nextdenovo_parallel_jobs)
        nextdenovo_group.add_argument('--nextdenovo-minreadlen', dest="nextdenovo_minreadlen", metavar="nextdenovo_minreadlen", default=self.nextdenovo_minreadlen, help='Filter reads with length < minreadlen. Default %s' % self.nextdenovo_minreadlen)
        nextdenovo_group.add_argument('--nextdenovo-seeddepth', dest="nextdenovo_seeddepth", metavar="nextdenovo_seeddepth", default=self.nextdenovo_seeddepth, type = int, help='Expected seed depth, used to calculate seed_cutoff, co-use with genome_size, you can try to set it 30-45 to get a better assembly result. Default %s' % self.nextdenovo_seeddepth)
        nextdenovo_group.add_argument('--nextdenovo-seedcutoff', dest="nextdenovo_seedcutoff", metavar="nextdenovo_seedcutoff", default=self.nextdenovo_seedcutoff, type = int, help='Minimum seed length, <=0 means calculate it automatically using bin/seq_stat. Default %s' % self.nextdenovo_seedcutoff)
        nextdenovo_group.add_argument('--nextdenovo-blocksize', dest="nextdenovo_blocksize", metavar="nextdenovo_blocksize", default=self.nextdenovo_blocksize, help='Block size for parallel running, split non-seed reads into small files, the maximum size of each file is blocksize. Default %s' % self.nextdenovo_blocksize)
        nextdenovo_group.add_argument('--nextdenovo-pa-correction ', dest="nextdenovo_pa_correction", metavar="nextdenovo_pa_correction", default=self.nextdenovo_pa_correction, type = int, help='number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage, overwrite parallel_jobs only for this step. Default %s' % self.nextdenovo_pa_correction )
        nextdenovo_group.add_argument('--nextdenovo-minimap_raw', dest="nextdenovo_minimap_raw", metavar="nextdenovo_minimap_raw", default=self.nextdenovo_minimap_raw , help='minimap2 options, used to find overlaps between raw reads, see minimap2-nd for details. Default %s' % self.nextdenovo_minimap_raw )
        nextdenovo_group.add_argument('--nextdenovo-minimap_cns', dest="nextdenovo_minimap_cns", metavar="nextdenovo_minimap_cns", default=self.nextdenovo_minimap_cns , help='minimap2 options, used to find overlaps between corrected reads. Default %s' % self.nextdenovo_minimap_cns )
        nextdenovo_group.add_argument('--nextdenovo-minimap_map', dest="nextdenovo_minimap_map", metavar="nextdenovo_minimap_map", default=self.nextdenovo_minimap_map , help='minimap2 options, used to map reads back to the assembly. Default %s' % self.nextdenovo_minimap_map)
        nextdenovo_group.add_argument('--nextdenovo-sort', dest="nextdenovo_sort", metavar="nextdenovo_sort", default=self.nextdenovo_sort, help='sort options, see ovl_sort for details. Default %s' % self.nextdenovo_sort)
        nextdenovo_group.add_argument('--nextdenovo-correction_opts', dest="nextdenovo_correction_opts", metavar="nextdenovo_correction_opts", default=self.nextdenovo_correction_opts, help='Correction options. Default %s' % self.nextdenovo_correction_opts)
        nextdenovo_group.add_argument('--nextdenovo-nextgraph_opt', dest="nextdenovo_nextgraph_opt", metavar="nextdenovo_nextgraph_opt", default=self.nextdenovo_nextgraph_opt, help='nextgraph options, see nextgraph for details. Default %s' % self.nextdenovo_nextgraph_opt)

    def register_hypo(self, parser):
        """Register all hypo parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        hypo_group = parser.add_argument_group('Hypo')
        hypo_group.add_argument('--sr-cov', dest="ill_cov", metavar="ill_cov", default= self.ill_cov, type=int, help='Approximate short read coverage for hypo Default %s' % self.ill_cov)
        hypo_group.add_argument('--hypo-proc', dest="hypo_processes", metavar="hypo_processes", default=self.hypo_processes, type = int, help='Number of contigs to be processed in parallel by HyPo. Default %s' % self.hypo_processes)
        hypo_group.add_argument('--hypo-no-lr', dest="hypo_lr", default=self.hypo_lr, action= "store_false", help='Set this to false if you don¡t want to run hypo with long reads. Default %s' % self.hypo_lr)
        hypo_group.add_argument('--hypo-opts', dest="hypo_opts", metavar="hypo_opts", default=self.hypo_opts, help='Additional options to run Hypo. Default %s' % self.hypo_opts)

    def register_purgedups(self, parser):
        """Register all purgedups parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        purgedups_group = parser.add_argument_group('Purge_dups')
        purgedups_group.add_argument('--purgedups-cores', type = int, dest="purgedups_cores", metavar="purgedups_cores", default=self.purgedups_cores, help='Number of threads to run purgedups. Default %s' % self.purgedups_cores)
        purgedups_group.add_argument('--purgedups-calcuts-opts', dest="calcuts_opts", metavar="calcuts_opts", default = self.calcuts_opts, help='Adjusted values to run calcuts for purgedups. Default %s' % self.calcuts_opts)

    def register_scaffold10X(self, parser):
        """Register all 10X scaffolding parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        scaffold10X_group = parser.add_argument_group('Scaffold_with_10X')
        scaffold10X_group.add_argument('--tigmint-cores', type = int, dest="tigmint_cores", metavar="tigmint_cores", default=self.tigmint_cores, help='Number of threads to run the 10X scaffolding step. Default %s' % self.tigmint_cores)
        scaffold10X_group.add_argument('--tigmint-opts', dest="tigmint_opts", metavar="tigmint_opts", default = self.tigmint_opts, help='Adjusted values to run the scaffolding with 10X reads.  Default %s' % self.tigmint_opts)

    def register_hic(self, parser):
        """Register all the HiC related parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        hic_group = parser.add_argument_group('HiC')
        hic_group.add_argument('--hic-qc', dest="hic_deepseq", action="store_false", help='Give this option if only QC of the HiC data needs to be done.')        
        hic_group.add_argument('--subsample-hic', dest="subsample_hic", action="store_true", help='Give this option if you want to subsample the hic data for qc.')        
        hic_group.add_argument('--add-preseq-opts', dest="add_preseq_opts", default=self.add_preseq_opts, help='Additional options to give for preseq etrapolation. E.g. -D. Default: %s' %self.add_preseq_opts)        
        hic_group.add_argument('--no-pretext', dest="get_pretext", action="store_false", help='Give this option if you do not want to generate the pretext file') 
        hic_group.add_argument('--sort-pretext', dest="sort_pretext", default=self.sort_pretext, help='Specify how to sort the pretext (eg. --nosort or --sortby something. Default: %s' %self.sort_pretext)        
        hic_group.add_argument('--assembly-qc', dest="assembly_qc", metavar = 'assembly_qc', help='Path to the assembly to be used perfom the QC of the HiC reads.')   
        hic_group.add_argument('--yahs-cores', dest="yahs_cores", metavar = 'yahs_cores', default = self.yahs_cores, help='Number of threads to run YAHS. Default %s' %self.yahs_cores)   
        hic_group.add_argument('--yahs-mq', dest="yahs_mq", metavar = 'yahs_mq', default = self.yahs_mq, help='Mapping quality to use when running yahs.Default %s' %self.yahs_mq)   
        hic_group.add_argument('--yahs-contig-ec', dest="yahs_contig_ec", action="store_true", help='Give this option if you want to allow yahs perform contig breaks.')   
        hic_group.add_argument('--yahs-opts', dest="yahs_opts", metavar = 'yahs_opts', default = self.yahs_opts, help='Additional options to give to YAHS.Default %s' %self.yahs_opts)   
        hic_group.add_argument('--hic-map-opts', dest="hic_map_opts", metavar = 'hic_map_opts', default = self.hic_map_opts, help='Options to use with bwa mem when aligning the HiC reads. Deafault %s' %self.hic_map_opts)   
        hic_group.add_argument('--mq', dest="mq", metavar = 'mq', default = self.mq, nargs = "+", help='Mapping qualities to use for processing the hic mappings. Default %s' %self.mq)   
        hic_group.add_argument('--hic-qc-assemblylen', dest="hic_qc_assemblylen", metavar = 'hic_qc_assemblylen', default = self.hic_qc_assemblylen, help='Lentgh of the assembly to be used for HiC QC')   
        hic_group.add_argument('--blast-cores', dest="blast_cores", metavar = 'blast_cores', default = self.blast_cores, help='Number of threads to run blast with the HiC unmapped reads.Default %s' %self.blast_cores)   
        hic_group.add_argument('--hic-blastdb', dest="blastdb", metavar = 'blastdb', default = self.blastdb, help='BLAST Database to use to classify the hic unmapped reads. Default %s' %self.blastdb)   
        hic_group.add_argument('--hic-readsblast', dest="hic_readsblast", metavar = 'hic_readsblast', default = self.hic_readsblast, help='Number of unmapped hic reads to classify with blast. Default %s' %self.hic_readsblast)   

    def register_finalize(self, parser):
        """Register all finalize parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        finalize_group = parser.add_argument_group('Finalize')
        finalize_group.add_argument('--no-final-evals', dest="final_evals", action="store_false", help='If specified, do not run evaluations on final assemblies. Default %s' % self.final_evals)
        finalize_group.add_argument('--busco-lin', dest="busco_lineage", metavar="busco_lineage", help='Path to the busco lineage to be used.')
        finalize_group.add_argument('--merqury-db', dest="merqury_db", metavar="merqury_db", help='Meryl database. Default %s' % self.merqury_db)
        finalize_group.add_argument('--merqury-plot-opts', dest="merqury_plot_opts", metavar="merqury_plot_opts", help='Meryl database. Default %s' % self.merqury_plot_opts)
        finalize_group.add_argument('--meryl-k', dest="meryl_k", metavar="meryl_k", type = int, help='Merqury plot additional options, for example \" -m 200 -n 6000|". Default %s' % self.merqury_plot_opts)
        finalize_group.add_argument('--meryl-threads', dest="meryl_threads", metavar="meryl_threads", type = int, default = self.meryl_threads, help='Number of threads to run meryl and merqury. Default %s' % self.meryl_threads)
        finalize_group.add_argument('--meryl-reads', dest="meryl_reads", metavar="meryl_reads", choices = ["illumina", "ont", "hifi"], nargs = "+",  help='Type of reads to be used to build the meryldb. Default %s' % self.meryl_reads)

    def register_wildcards(self, parser):
        """Register all wildcards parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        wildcards_group = parser.add_argument_group('Wildcards')
        wildcards_group.add_argument('--ont-list', dest="ONT_wildcards", metavar="ONT_wildcards", help='List with basename of the ONT fastqs that will be used. Default %s' % self.ONT_wildcards)
        wildcards_group.add_argument('--hifi-list', dest="hifi_wildcards", metavar="hifi_wildcards", help='List with basename of the ONT fastqs that will be used. Default %s' % self.hifi_wildcards)
        wildcards_group.add_argument('--illumina-list', dest="illumina_wildcards", metavar="illumina_wildcards", help='List with basename of the illumina fastqs. Default %s' % self.illumina_wildcards)
        wildcards_group.add_argument('--r10X-list', dest="r10X_wildcards", metavar="r10X_wildcards", help='List with basename of the raw 10X fastqs. Default %s' % self.r10X_wildcards)
        wildcards_group.add_argument('--hic-list', dest="hic_wildcards", metavar="hic_wildcards", help='List with basename of the raw hic fastqs. Default %s' % self.hic_wildcards)
 
####

    def check_parameters(self,args):
        """Check parameters consistency
            
        args -- set of parsed arguments"""

        if args.configFile != None:
          args.configFile = os.path.abspath(args.configFile) 
        else:
          args.configFile = os.path.abspath(self.configFile) 

        if args.specFile != None:
          args.specFile = os.path.abspath(args.specFile) 
        else:
          args.specFile = os.path.abspath(self.specFile) 

        if args.ndconfFile != None:
          args.ndconfFile = os.path.abspath(args.ndconfFile) 
        else:
          args.ndconfFile = os.path.abspath(self.ndconfFile) 

        if args.base_name == None:
          parser.print_help()
          print ("You need to provide a base name for the project")
          sys.exit(-1)   

        if args.species == None:
          parser.print_help()
          print ("You need to provide the species name")
          sys.exit(-1)   

        if args.genome_size == None:
          parser.print_help()
          print ("You need to provide a genome size estimate")
          sys.exit(-1)
        else:
          if re.search("m", args.genome_size):
            gsize = float(re.split("m", args.genome_size)[0])
          elif re.search("g", args.genome_size):
            gsize = float(re.split("g", args.genome_size)[0]) * 1000
          else:
            parser.print_help()
            print ("Provided genome size estimate is in unknown format, please check suggestions and give and appropriate value")
            sys.exit(-1)
          print ("Genome size is " + str(gsize) + " megabases")

        print("Warning: Long-read type is set to " + args.lr_type + ", make sure this is your data type.")

        if args.scripts_dir:
          args.scripts_dir = os.path.abspath(args.scripts_dir) + "/"
        else:
          args.scripts_dir = os.path.abspath(self.scripts_dir) + "/"
        if not os.path.exists(args.scripts_dir):
          print (args.scripts_dir + " not found")

        if args.run_yahs == True and args.hic_deepseq == False:
          parser.print_help()
          print ("Running yahs is not compatible with hic-qc, please select the appropriate option")
          sys.exit(-1)

        if args.pipeline_workdir != None:
          args.pipeline_workdir = os.path.abspath(args.pipeline_workdir) + "/"
        else:
          args.pipeline_workdir = os.getcwd() + "/"

        if args.eval_dir:
          args.eval_dir = os.path.abspath(args.eval_dir) + "/"
        else:
          args.eval_dir = args.pipeline_workdir  + self.eval_dir

        if args.stats_out == None:
          args.stats_out = args.eval_dir + args.base_name + ".stats_summary.txt"
        else:
          args.stats_out = os.path.abspath(args.stats_out)

        if args.hic_qc_dir:
          args.hic_qc_dir = os.path.abspath(args.hic_qc_dir) + "/"
        else:
          args.hic_qc_dir= args.pipeline_workdir  + self.hic_qc_dir

        args.assemblies={}
        if len(args.assembly_in):
          for my_dict in args.assembly_in:
            for key in my_dict:
              args.assemblies[os.path.abspath(key)] = os.path.abspath(my_dict[key]) + "/"

        args.assemblies_cur = {}
        if len(args.postpolish_assemblies):
          for my_dict in args.postpolish_assemblies:
            for key in my_dict:
              args.assemblies_cur[os.path.abspath(key)] = my_dict[key]

        args.final_assemblies = {}
        if len(args.curated_assemblies):
          for my_dict in args.curated_assemblies:
            for key in my_dict:
              args.final_assemblies[os.path.abspath(key)] = os.path.abspath(my_dict[key])

        if args.lr_type == "pacbio-hifi":
           args.hypo_rounds = 0
           args.nextpolish_lr = 0
           args.nextpolish_sr = 0

        if args.run_flye == False and args.run_nextdenovo == False:
           args.hypo_rounds = 0
        
        if args.run_flye == True or args.run_nextdenovo == True or args.run_hifiasm == True or args.nextpolish_ont_rounds > 0 or args.run_purgedups == True:
          require_lr()
        elif args.hypo_rounds > 0 and args.hypo_lr == True:
          require_lr()
        
        args.r10X_reads = {}

        if args.nextpolish_ill_rounds > 0 or args.hypo_rounds >0  or args.run_tigmint == True:
          if not re.search("hifi", args.lr_type):
            require_illumina()
          
            if args.r10X_reads: 
              if args.longranger_path:
                args.longranger_path = os.path.abspath(args.longranger_path)
              else:
                args.longranger_path =  os.path.abspath(self.longranger_path)
              if not os.path.exists(args.longranger_path):
                print (args.longranger_path + " not found")

        if args.flye_dir == None:
          args.flye_dir = args.pipeline_workdir + "s" + args.flye_step + "_p" + args.preprocess_lr_step + "_flye/"
        else:
          args.flye_dir = os.path.abspath(args.flye_dir) + "/"  
        args.flye_out = args.flye_dir + args.base_name + ".fl.asm.fasta"

        if args.nextdenovo_dir == None:
          args.nextdenovo_dir = args.pipeline_workdir + "s" + args.nextdenovo_step + "_p" + args.preprocess_lr_step + "_nextdenovo/"
        else:
          args.nextdenovo_dir = os.path.abspath(args.nextdenovo_dir) + "/" 
        args.nextdenovo_out = args.nextdenovo_dir + args.base_name + ".nd.asm.fasta"

        if args.hifiasm_dir == None:
          args.hifiasm_dir = args.pipeline_workdir + "s" + args.hifiasm_step + "_p" + args.preprocess_lr_step + "_hifiasm/"
        else:
          args.hifiasm_dir = os.path.abspath(args.hifiasm_dir) + "/"  
        args.hifiasm_out = []

        if args.hifiasm_phasing == True:
          args.hifiasm_out.append(args.hifiasm_dir + args.base_name + ".hfsm.asm.hic.fa")
          for n in range(1, args.ploidy+1):
            args.hifiasm_out.append(args.hifiasm_dir + args.base_name + ".hfsm.asm.hic.hap" + str(n) + ".fa")
        else:
          args.hifiasm_out.append(args.hifiasm_dir + args.base_name + ".hfsm.asm.bp.fa")
          for n in range(1,3):
            args.hifiasm_out.append(args.hifiasm_dir + args.base_name + ".hfsm.asm.bp.hap" + str(n) + ".fa")
            
        if args.nextdenovo_dir == None:
          args.nextdenovo_dir = args.pipeline_workdir + "s" + args.nextdenovo_step + "_p" + args.preprocess_ont_step + "_nextdenovo/"
        else:
          args.nextdenovo_dir = os.path.abspath(args.nextdenovo_dir) + "/" 
        args.nextdenovo_out = args.nextdenovo_dir + "nd.asm.fasta"

        if args.merqury_db:
          args.merqury_db = os.path.abspath(args.merqury_db)
          if args.meryl_reads == None:
              if args.lr_type == "pacbio-hifi":
                args.meryl_reads = "hifi"
              else:
                args.meryl_reads=self.meryl_reads 
          if not os.path.exists(args.merqury_db):
            if "ont" in args.meryl_reads:
              require_ont()
            if "illumina" in args.meryl_reads:
              require_illumina()
            if "hifi" in args.meryl_reads:
              require_pb()

            if args.meryl_k == None:
              parser.print_help()
              print (args.merqury_db + " not found, the pipeline will create it but you need to provide the kmer value")
              sys.exit(-1)             
            else:
              print (args.merqury_db + " not found, the pipeline will create it")  
        elif args.final_evals == True:
          print ("WARNING: If you want to run merqury, please provide a merqury_db path, if it does not exist, the pipeline will create it")

        if args.busco_lineage:
          args.busco_lineage = os.path.abspath(args.busco_lineage)
          if not os.path.exists(args.busco_lineage):
            print (args.busco_lineage + " not found")
        elif args.final_evals == True:
          print ("WARNING: busco lineage is needed if you want to run Busco")

        if args.run_tigmint == True and args.r10X == None:  
          if args.processed_10X == None:
            parser.print_help()
            print ("10X reads are required to run the 10X scaffolding step")
            sys.exit(-1) 

        if args.run_kraken2 == True:
          if args.kraken2_db == None or args.kraken2_kmers == None:
            parser.print_help()
            print ("Please, specify a kraken2 database and a kmer distribution file if you want to run Kraken2 and Bracken on the reads. Otherwise, do not use the --run-kraken option.")
            sys.exit(-1)
          else:
            args.kraken2_db = os.path.abspath(args.kraken2_db)
            args.kraken2_kmers = os.path.abspath(args.kraken2_kmers)
            if not os.path.exists(args.kraken2_db):
              parser.print_help()
              print ("\n" + args.kraken2_db + " should exist and it does not.")
              sys.exit(-1)
            if not os.path.exists(args.kraken2_kmers):
              parser.print_help()
              print ("\n" + args.kraken2_kmers + " should exist and it does not.")
              sys.exit(-1)

        if not re.search("-g ",args.other_flye_opts) and not re.search("--genome-size ", args.other_flye_opts):
          args.other_flye_opts += " -g " + args.genome_size + " "

        if args.ploidy != 2 and not re.search("-n-hap", args.other_hifiasm_opts):
          args.other_hifiasm_opts += " --n-hap " + str(args.ploidy) + " "
        if args.lr_type != "pacbio-hifi":
          if not re.search("--ont ",args.other_hifiasm_opts):
            args.other_hifiasm_opts += " --ont "
        else:
          if re.search("--ont", args.other_hifiasm_opts):
            args.other_hifiasm_opts =  args.other_hifiasm_opts.replace("--ont", "")

        if args.nextpolish_ill_rounds > 0 or args.hypo_rounds >0 or args.nextpolish_ont_rounds:
          if args.run_flye == True:
            if args.polish_flye_dir != None:
              args.polish_flye_dir = os.path.abspath(args.polish_flye_dir)
            else:
              step = float(args.flye_step) + 1
              args.polish_flye_dir = args.pipeline_workdir + "s0" + str(step) + "_p" + args.flye_step + "_polishing/"
            args.assemblies[args.flye_out] = args.polish_flye_dir
          if args.run_nextdenovo== True:
            if args.polish_nextdenovo_dir != None:
              args.polish_nextdenovo_dir = os.path.abspath(args.polish_nextdenovo_dir)
            else:
              step = float(args.nextdenovo_step) + 1
              args.polish_nextdenovo_dir = args.pipeline_workdir + "s0" + str(step) + "_p" + args.nextdenovo_step + "_polishing/"
            args.assemblies[args.nextdenovo_out] = args.polish_nextdenovo_dir
        elif args.run_purgedups == True:
          if args.run_flye == True:
            args.assemblies[args.flye_out] = args.flye_dir
          if args.run_nextdenovo == True:
            args.assemblies[args.nextdenovo_out] = args.nextdenovo_dir

        paths = 0
        if args.run_purgedups == True:      
          if len(args.assemblies) > 0:
            pol_bases = {}
            base_tmp = ""
            if  args.hypo_rounds >0:
              pol_bases["hypo"] = "hyp" + str(args.hypo_rounds)
            if args.nextpolish_ont_rounds > 0:
              base_tmp+= "npo" + str(args.nextpolish_ont_rounds)
            if args.nextpolish_ill_rounds > 0:
              if base_tmp != "":
                base_tmp += "."
              base_tmp += "npi" + str(args.nextpolish_ill_rounds)
            if base_tmp != "":
              pol_bases["nextpolish"] = base_tmp
            
            for m in args.assemblies:
              bpol = os.path.splitext(os.path.basename(m))[0]
              path = args.assemblies[m]
              pstep = path.split('/')[-2].split('_')[0]
                
              nstep = pstep.replace('s','')
                #cstep = float(nstep) + 1 + paths
              if paths != 0:
                paths -= 0.1
              if pol_bases:
                for p in pol_bases:
                  cstep = round(float(nstep) + 1 + paths,1)
                  args.assemblies_cur[args.assemblies[m] + p + "/" + bpol + "." +  pol_bases[p] + ".fasta"] = "s0" + str(cstep) + "_p" + nstep
                  paths += 0.1
              else:
                cstep = round(float(nstep) + 1 + paths,1)
                args.assemblies_cur[args.assemblies[m] + bpol + ".fasta"] = "s0" + str(cstep) + "_p" + nstep
                paths += 0.1

        if args.run_hifiasm and len(args.hifiasm_out):
          if args.hifiasm_ext_purge == True or args.run_yahs == True:
            if paths != 0:
              paths -= 0.1
            for m in args.hifiasm_out:
              path = os.path.basename(os.path.dirname(m))
              pstep = path.split('_')[0]
              nstep = pstep.replace('s','')
              cstep = round(float(nstep) + 1 + paths,1)
              args.assemblies_cur[m] = "s0" + str(cstep) + "_p" + nstep
              paths += 0.1
        
        if args.hifi_dir:
           require_pb()

        if args.ONT_dir:
           require_ont()

        if args.hic_dir:
          args.hic_dir = os.path.abspath(args.hic_dir) + "/"          
          args.hic_wildcards = get_wildcards(args.hic_dir, args.hic_wildcards, '.1.fastq.gz')
          if args.concat_hic_dir:
            args.concat_hic_dir = os.path.abspath(args.concat_hic_dir) + "/"      
          else:
            args.concat_hic_dir = args.pipeline_workdir + "s" + args.preprocess_hic_step + "_p01.1_Concat_HiC/"   
         
          if not args.assembly_qc and args.hic_deepseq == False:
            parser.print_help()
            print ("Please, specify an assembly to use as input for the hic qc.")
            sys.exit(-1)             
          elif args.assembly_qc:
            args.assembly_qc = os.path.abspath(args.assembly_qc)
            args.map_hic_qos = 'vshort'
            args.map_hic_time = '1:00:00'
            args.map_hic_mem = '50G'
            args.pairtools_qos = 'vshort'
            args.pairtools_time = '1:00:00'
            args.pairtools_mem = '50G'
            args.qcstats_qos = 'vshort'
            args.qcstats_time = '1:00:00'
            args.qcstats_mem = '50G'

            if args.get_pretext == True:
              print ("Warning. Turning get_pretext parameter to false, since it makes little sense to run this step with the hic QC data.")
              args.get_pretext = False

            if not args.hic_qc_assemblylen:
              parser.print_help()
              print ("Please, specify the total length of " + args.assembly_qc + " to perform the HiC QC step.")
              sys.exit(-1)    

        args.concat_reads_qos =  self.concat_reads_qos
        args.concat_reads_time = self.concat_reads_time 
        args.concat_reads_queue = self.concat_reads_queue
        args.concat_reads_mem =  self.concat_reads_mem 

        args.filtlong_qos =  self.filtlong_qos
        args.filtlong_time = self.filtlong_time 
        args.filtlong_queue = self.filtlong_queue 
        args.filtlong_mem =  self.filtlong_mem   

        args.bam2fq_qos =  self.bam2fq_qos
        args.bam2fq_time = self.bam2fq_time 
        args.bam2fq_queue = self.bam2fq_queue 
        args.bam2fq_mem =  self.bam2fq_mem   

        args.trimgalore_qos = self.trimgalore_qos
        args.trimgalore_time = self.trimgalore_time
        args.trimgalore_queue = self.trimgalore_queue
        args.trimgalore_mem = self.trimgalore_mem

        args.longranger_qos =  self.longranger_qos
        args.longranger_time = self.longranger_time 
        args.longranger_queue = self.longranger_queue   
        args.longranger_mem = self.longranger_mem  

        args.nanoplot_qos =  self.nanoplot_qos
        args.nanoplot_time = self.nanoplot_time 
        args.nanoplot_queue = self.nanoplot_queue
        args.nanoplot_mem =  self.nanoplot_mem     

        args.kraken2_qos = self.kraken2_qos
        args.kraken2_time = self.kraken2_time 
        args.kraken2_queue = self.kraken2_queue
        args.kraken2_mem = self.kraken2_mem
        
        args.build_meryl_qos =  self.build_meryl_qos
        args.build_meryl_time = self.build_meryl_time 
        args.build_meryl_queue = self.build_meryl_queue
        args.build_meryl_mem = self.build_meryl_mem

        args.concat_meryl_qos =  self.concat_meryl_qos
        args.concat_meryl_time = self.concat_meryl_time 
        args.concat_meryl_queue = self.concat_meryl_queue
        args.concat_meryl_mem = self.concat_meryl_mem

        args.genomescope_qos =  self.genomescope_qos
        args.genomescope_time = self.genomescope_time 
        args.genomescope_queue = self.genomescope_queue
        args.genomescope_mem = self.genomescope_mem

        args.smudgeplot_qos =  self.smudgeplot_qos
        args.smudgeplot_time = self.smudgeplot_time 
        args.smudgeplot_queue = self.smudgeplot_queue
        args.smudgeplot_mem = self.smudgeplot_mem

        args.flye_qos =  self.flye_qos
        args.flye_time = self.flye_time 
        args.flye_queue = self.flye_queue
        args.flye_mem =  self.flye_mem      

        args.hifiasm_qos =  self.hifiasm_qos
        args.hifiasm_time = self.hifiasm_time 
        args.hifiasm_queue = self.hifiasm_queue
        args.hifiasm_mem =  self.hifiasm_mem  

        args.nextdenovo_qos =  self.nextdenovo_qos
        args.nextdenovo_time = self.nextdenovo_time 
        args.nextdenovo_queue = self.nextdenovo_queue
        args.nextdenovo_mem =  self.nextdenovo_mem

        if args.nextdenovo_type == "local":
          args.nextdenovo_mem =  "900G"
          args.nextdenovo_cores = 128
          args.nextdenovo_qos = "marathon"
          args.nextdenovo_time = "100:00:00"

        args.minimap_qos =  self.minimap_qos
        args.minimap_time = self.minimap_time 
        args.minimap_queue = self.minimap_queue
        args.minimap_mem = self.minimap_mem

        args.bwa_qos =  self.bwa_qos
        args.bwa_time = self.bwa_time 
        args.bwa_queue = self.bwa_queue
        args.bwa_mem = self.bwa_mem

        args.hypo_qos =  self.hypo_qos
        args.hypo_time = self.hypo_time 
        args.hypo_queue = self.hypo_queue
        args.hypo_mem = self.hypo_mem

        args.nextpolish_lr_qos =  self.nextpolish_lr_qos
        args.nextpolish_lr_time = self.nextpolish_lr_time 
        args.nextpolish_lr_queue = self.nextpolish_lr_queue
        args.nextpolish_lr_mem = self.nextpolish_lr_mem
        
        args.nextpolish_sr_qos =  self.nextpolish_sr_qos
        args.nextpolish_sr_time = self.nextpolish_sr_time 
        args.nextpolish_sr_queue = self.nextpolish_sr_queue
        args.nextpolish_sr_mem = self.nextpolish_sr_mem

        args.purgedups_qos =  self.purgedups_qos
        args.purgedups_time = self.purgedups_time 
        args.purgedups_queue = self.purgedups_queue
        args.purgedups_mem = self.purgedups_mem

        args.tigmint_qos =  self.tigmint_qos
        args.tigmint_time = self.tigmint_time 
        args.tigmint_queue = self.tigmint_queue
        args.tigmint_mem = self.tigmint_mem

        args.ass_prepare_qos =  self.ass_prepare_qos
        args.ass_prepare_time = self.ass_prepare_time 
        args.ass_prepare_queue = self.ass_prepare_queue
        args.ass_prepare_mem = self.ass_prepare_mem

        args.map_hic_qos =  self.map_hic_qos
        args.map_hic_time = self.map_hic_time 
        args.map_hic_queue = self.map_hic_queue
        args.map_hic_mem = self.map_hic_mem

        args.pairtools_qos =  self.pairtools_qos
        args.pairtools_time = self.pairtools_time 
        args.pairtools_queue = self.pairtools_queue
        args.pairtools_mem = self.pairtools_mem

        args.blast_qos =  self.blast_qos
        args.blast_time = self.blast_time 
        args.blast_queue = self.blast_queue
        args.blast_mem = self.blast_mem

        args.qcstats_qos =  self.qcstats_qos
        args.qcstats_time = self.qcstats_time 
        args.qcstats_queue = self.qcstats_queue
        args.qcstats_mem = self.qcstats_mem

        args.yahs_qos =  self.yahs_qos
        args.yahs_time = self.yahs_time 
        args.yahs_queue = self.yahs_queue
        args.yahs_mem = self.yahs_mem

        args.pretext_qos =  self.pretext_qos
        args.pretext_time = self.pretext_time 
        args.pretext_queue = self.pretext_queue
        args.pretext_mem = self.pretext_mem

        args.telext_qos =  self.telext_qos
        args.telext_time = self.telext_time 
        args.telext_queue = self.telext_queue
        args.telext_mem = self.telext_mem

        args.tpf_qos =  self.tpf_qos
        args.tpf_time = self.tpf_time 
        args.tpf_queue = self.tpf_queue
        args.tpf_mem = self.tpf_mem
  
        args.stats_qos =  self.stats_qos
        args.stats_time = self.stats_time 
        args.stats_queue = self.stats_queue
        args.stats_mem = self.stats_mem

        args.busco_qos =  self.busco_qos
        args.busco_time = self.busco_time 
        args.busco_queue = self.busco_queue
        args.busco_mem = self.busco_mem

        args.merq_qos =  self.merq_qos
        args.merq_time = self.merq_time 
        args.merq_queue = self.merq_queue
        args.merq_mem = self.merq_mem

        args.fin_qos =  self.fin_qos
        args.fin_time = self.fin_time 
        args.fin_queue = self.fin_queue
        args.fin_mem = self.fin_mem

        if gsize > 500:
          args.java_opts = "Xmx150g"
          args.flye_cores = 100
          args.flye_mem = "700G"

        if gsize > 1000:
          args.concat_cores = 16
          args.flye_qos = "marathon"
          args.flye_time = "150:00:00"
          if args.nextdenovo_type == "local":
            args.nextdenovo_time = "150:00:00"
          args.nanoplot_time = "48:00:00"
          args.nanoplot_queue = "vlong"
          args.bwa_time = "24:00:00"
          args.bwa_qos = "long"
          args.bwa_mem = "300G"
          args.bwa_cores = 48
          args.minimap_time = "24:00:00"
          args.minimap_qos = "long"
          args.minimap2_cores = 80
          args.hypo_time = "12:00:00"
          args.hypo_cores = 32
          args.busco_time = "24:00:00"
          args.busco_qos = "long"
          args.busco_cores = 64
          args.merq_time = "6:00:00"
          args.ass_prepare_time = "6:00:00"
          args.ass_prepare_mem = "250G"
          args.map_hic_qos = 'long'
          args.map_hic_time = "24:00:00"
          args.map_hic_mem = "300G"
          args.yahs_qos = 'long'
          args.yahs_time = "24:00:00"
          args.yahs_mem = "100G"
          args.pairtools_mem = "500G"
          args.pairtools_parse_cores = 128
          args.pairtools_sort_cores = 128
          args.pairtools_dedup_cores = 4
          args.pairtools_split_cores = 16
          args.pairtools_time = "48:00:00"
          args.pairtools_qos = "vlong"
          args.qcstats_qos = 'normal'
          args.qcstats_time = "12:00:00"
          args.qcstats_mem = "100G"

        if gsize > 3000:
          args.ass_prepare_mem = "500G"
          args.flye_qos = "eternal"
          args.flye_time = "500:00:00"
          args.hifiasm_qos = "marathon_assembly"
          args.hifiasm_time = "150:00:00"
          args.hifiasm_mem = "950G"
          args.hifiasm_cores = 128
          args.minimap_time = "48:00:00"
          args.minimap_qos = "vlong"
          args.hypo_processes = 2
          args.hypo_mem = "500G"
          args.hypo_qos = "long"
          args.hypo_time = "24:00:00"
          args.stats_time = "1:00:00"
          args.stats_qos = "vshort"
          args.busco_cores = 100
          args.busco_time = "48:00:00"
          args.busco_qos = "vlong"
          args.busco_mem = "150G"
          args.purgedups_time = "12:00:00"
          args.pairtools_time = "124:00:00"
          args.pairtools_qos = "marathon"

          if args.nextdenovo_type == "local":
            args.nextdenovo_time = "500:00:00"
           
        if gsize <= 100:
          args.busco_time = "2:00:00"
          args.merq_time = "2:00:00"
          args.fin_time = "0:30:00"
          if args.nextdenovo_type == "local":
            args.nextdenovo_qos = "long"
            args.nextdenovo_time = "24:00:00"
            args.nextdenovo_mem = "100G"
            args.nextdenovo_cores = 50
          args.flye_cores = 20
          args.flye_mem = "100G"
          args.hifiasm_qos="normal"
          args.hifiasm_time = "12:00:00"

###

    def storeGeneralParameters(self,args):
        """Updates general parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.generalParameters["configFile"] = args.configFile
        self.generalParameters["specFile"] = args.specFile
        self.generalParameters["ndconfFile"] = args.ndconfFile
        self.generalParameters["keep_intermediate"] = args.keepintermediate
        self.generalParameters["lr_type"] = args.lr_type
        self.generalParameters["base_name"] = args.base_name
        self.generalParameters["species"] = args.species
        self.generalParameters["genome_size"] = args.genome_size
        self.generalParameters["ploidy"] = args.ploidy
        self.generalParameters["telo_repeat"] = args.telo_string
        self.generalParameters["run_flye"] = args.run_flye
        self.generalParameters["run_hifiasm"] = args.run_hifiasm
        self.generalParameters["run_nextdenovo"] = args.run_nextdenovo
        self.generalParameters["nextpolish_ont_rounds"] = args.nextpolish_ont_rounds
        self.generalParameters["nextpolish_ill_rounds"] = args.nextpolish_ill_rounds
        self.generalParameters["hypo_rounds"] = args.hypo_rounds
        self.generalParameters["run_purgedups"] = args.run_purgedups
        self.generalParameters["run_yahs"] = args.run_yahs
        self.generalParameters["run_tigmint"] = args.run_tigmint
        self.generalParameters["run_kraken2"] = args.run_kraken2
        self.generalParameters["run_smudgeplot"] = args.run_smudgeplot
        self.generalParameters["genomescope_additional_options"] = args.genomescope_additional
        self.generalParameters["preprocess_lr_step"] = args.preprocess_lr_step
        self.generalParameters["preprocess_illumina_step"] = args.preprocess_illumina_step
        self.generalParameters["preprocess_10X_step"] = args.preprocess_10X_step
        self.generalParameters["preprocess_hic_step"] = args.preprocess_hic_step
        self.generalParameters["flye_step"] = args.flye_step
        self.generalParameters["hifiasm_step"] = args.hifiasm_step
        self.generalParameters["nextdenovo_step"] = args.nextdenovo_step
        self.generalParameters["concat_cores"] = args.concat_cores
        self.generalParameters["minimap2_cores"] = args.minimap2_cores
        self.generalParameters["BWA_cores"] = args.bwa_cores
        self.generalParameters["hypo_cores"] = args.hypo_cores
        self.generalParameters["nextpolish_cores"] = args.nextpolish_cores
        self.generalParameters["pairtools_parse_cores"] = args.pairtools_parse_cores
        self.generalParameters["pairtools_sort_cores"] = args.pairtools_sort_cores
        self.generalParameters["pairtools_dedup_cores"] = args.pairtools_dedup_cores
        self.generalParameters["pairtools_split_cores"] = args.pairtools_split_cores
        self.generalParameters["busco_cores"] = args.busco_cores
        self.generalParameters["longranger_cores"] = args.longranger_cores
        self.generalParameters["longranger_path"] = args.longranger_path
        self.allParameters["Parameters"] = self.generalParameters

    def storeallSpecParameters(self,args):
        """Updates rule all cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.allSpecParameters["name"] = "{rule}_" + args.base_name + "_assembly_pipeline"
        self.allSpecParameters["qos"] = self.all_qos
        self.allSpecParameters["time"] = self.all_time
        self.allSpecParameters["queue"] = self.all_queue
        self.allParameters ["all"] = self.allSpecParameters

    def storeInputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.inputParameters["scripts_dir"] = args.scripts_dir
        self.inputParameters["ONT_dir"] = args.ONT_dir
        self.inputParameters["Hifi_dir"] = args.hifi_dir
        self.inputParameters["illumina_dir"] = args.illumina_dir
        self.inputParameters["raw_10X"] = args.r10X_reads
        self.inputParameters["HiC_dir"] = args.hic_dir
        self.inputParameters["ONT_reads"] = args.ONT_reads
        self.inputParameters["hifi_reads"] = args.hifi_reads
        self.inputParameters["ILLUMINA_pe1"] = args.pe1
        self.inputParameters["ILLUMINA_pe2"] = args.pe2
        self.inputParameters["ILLUMINA_10X"] = args.r10X
        self.inputParameters["processed_illumina"] = args.processed_illumina
        self.inputParameters["ONT_filtered"] = args.ONT_filtered
        self.inputParameters["processed_10X"] = args.processed_10X
        self.inputParameters["Assemblies for polishing"] = args.assemblies
        self.inputParameters["Assemblies for postpolishing"] = args.assemblies_cur
        self.inputParameters["Curated Assemblies"] = args.final_assemblies
        self.allParameters ["Inputs"] = self.inputParameters

    def storeOutputParameters(self,args):
        """Updates output parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.outputParameters["base_dir"] = args.pipeline_workdir 
        self.outputParameters["preprocess_lr"] = args.preprocess_lr
        self.outputParameters["concat_HiC_dir"] = args.concat_hic_dir
        self.outputParameters["flye_dir"] = args.flye_dir
        self.outputParameters["nextdenovo_dir"] = args.nextdenovo_dir
        self.outputParameters["hifiasm_dir"] = args.hifiasm_dir
        self.outputParameters["flye_out"] = args.flye_out
        self.outputParameters["hifiasm_out"] = args.hifiasm_out
        self.outputParameters["nextdenovo_out"] = args.nextdenovo_out
        self.outputParameters["eval_dir"] = args.eval_dir
        self.outputParameters["stats_out"] = args.stats_out
        self.outputParameters["hic_qc_dir"] = args.hic_qc_dir
        self.allParameters ["Outputs"] = self.outputParameters

    def storeconcatreadsSpecParameters(self,args):
        """Updates concat reads cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.concatreadsSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.ext}"
        self.concatreadsSpecParameters["qos"] = args.concat_reads_qos
        self.concatreadsSpecParameters["time"] = args.concat_reads_time
        self.concatreadsSpecParameters["queue"] = args.concat_reads_queue
        self.concatreadsSpecParameters["mem"] = args.concat_reads_mem
        self.allParameters ["concat_reads"] = self.concatreadsSpecParameters

    def storesubsampleHICSpecParameters(self,args):
        """Updates subsample hic reads cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.subsampleHICSpecParameters["name"] = "{rule}_" + args.base_name
        self.subsampleHICSpecParameters["qos"] = args.concat_reads_qos
        self.subsampleHICSpecParameters["time"] = args.concat_reads_time
        self.subsampleHICSpecParameters["queue"] = args.concat_reads_queue
        self.subsampleHICSpecParameters["mem"] = args.concat_reads_mem
        self.allParameters ["subsample_hic"] = self.subsampleHICSpecParameters

    def storeFiltlongParameters(self,args):
        """Updates filtlong parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.filtlongParameters["Filtlong minlen"] = args.filtlong_minlen
        self.filtlongParameters["Filtlong min_mean_q"] = args.filtlong_min_mean_q
        self.filtlongParameters["options"] = args.filtlong_opts
        self.allParameters ["Filtlong"] = self.filtlongParameters

    def storefiltlongSpecParameters(self,args):
        """Updates filtlong cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.filtlongSpecParameters["name"] = "{rule}_" + args.base_name + "_s" + args.preprocess_lr_step 
        self.filtlongSpecParameters["qos"] = args.filtlong_qos
        self.filtlongSpecParameters["time"] = args.filtlong_time
        self.filtlongSpecParameters["mem"] = args.filtlong_mem
        self.filtlongSpecParameters["queue"] = args.filtlong_queue
        self.allParameters ["filtlong"] = self.filtlongSpecParameters

    def storebam2fqSpecParameters(self,args):
        """Updates bam2fastq cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.bam2fqSpecParameters["name"] = "{rule}_" + args.base_name + "_s" + args.preprocess_lr_step 
        self.bam2fqSpecParameters["qos"] = args.bam2fq_qos
        self.bam2fqSpecParameters["time"] = args.bam2fq_time
        self.bam2fqSpecParameters["mem"] = args.bam2fq_mem
        self.bam2fqSpecParameters["queue"] = args.bam2fq_queue
        self.allParameters ["bam2fastq"] = self.bam2fqSpecParameters

    def storeTrimgaloreParameters(self,args):
        """Updates the Trim_Galore parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.trimgaloreParameters["options"] = args.trim_galore_opts
        self.trimgaloreParameters["Trim_Illumina_cores"] = args.Trim_Illumina_cores
        self.allParameters ["Trim_Galore"] = self.trimgaloreParameters

    def storetrimgaloreSpecParameters(self,args):
        """Updates trimgalore cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.trimgaloreSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.file}"
        self.trimgaloreSpecParameters["qos"] = args.trimgalore_qos
        self.trimgaloreSpecParameters["time"] = args.trimgalore_time
        self.trimgaloreSpecParameters["queue"] = args.trimgalore_queue
        self.trimgaloreSpecParameters["mem"] = args.trimgalore_mem
        self.allParameters ["trim_galore"] = self.trimgaloreSpecParameters

    def storelongrangerSpecParameters(self,args):
        """Updates longranger cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.longrangerSpecParameters["name"] = "{rule}_" + args.base_name + "_{base}_{wildcards.bname}"
        self.longrangerSpecParameters["qos"] = args.longranger_qos
        self.longrangerSpecParameters["time"] = args.longranger_time
        self.longrangerSpecParameters["queue"] = args.longranger_queue
        self.longrangerSpecParameters["mem"] = args.longranger_mem
        self.allParameters ["long_ranger"] = self.longrangerSpecParameters

    def storenanoplotSpecParameters(self,args):
        """Updates Nanoplot cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nanoplotSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.prefix}"
        self.nanoplotSpecParameters["qos"] = args.nanoplot_qos
        self.nanoplotSpecParameters["time"] = args.nanoplot_time
        self.nanoplotSpecParameters["queue"] = args.nanoplot_queue
        self.nanoplotSpecParameters["mem"] = args.nanoplot_mem
        self.allParameters ["nanoplot"] = self.nanoplotSpecParameters

    def storekraken2Parameters(self,args):
        """Updates Kraken2 parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.kraken2Parameters["database"] = args.kraken2_db
        self.kraken2Parameters["kmer_dist"] = args.kraken2_kmers
        self.kraken2Parameters["threads"] = args.kraken2_threads
        self.kraken2Parameters["additional_opts"] = args.additional_kraken2_opts
        self.allParameters ["Kraken2"] = self.kraken2Parameters

    def storekraken2SpecParameters(self,args):
        """Updates Kraken2 cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.kraken2SpecParameters["name"] = "{rule}_" + args.basename + "_{params.prefix}"
        self.kraken2SpecParameters["qos"] = args.kraken2_qos
        self.kraken2SpecParameters["time"] = args.kraken2_time
        self.kraken2SpecParameters["queue"] = args.kraken2_queue
        self.kraken2SpecParameters["mem"] = args.kraken2_mem
        self.allParameters ["Kraken2"] = self.kraken2SpecParameters

    def storebuildmerylSpecParameters(self,args):
        """Updates build meryl cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.buildmerylSpecParameters["name"] = "{rule}_" + args.base_name + "_{base}_{wildcards.db}"
        self.buildmerylSpecParameters["qos"] = args.build_meryl_qos
        self.buildmerylSpecParameters["time"] = args.build_meryl_time
        self.buildmerylSpecParameters["queue"] = args.build_meryl_queue
        self.buildmerylSpecParameters["mem"] = args.build_meryl_mem
        self.allParameters ["build_meryl_db"] = self.buildmerylSpecParameters

    def storeconcatmerylSpecParameters(self,args):
        """Updates concat meryl cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.concatmerylSpecParameters["name"] = "{rule}_" + args.base_name + "_{base}"
        self.concatmerylSpecParameters["qos"] = args.concat_meryl_qos
        self.concatmerylSpecParameters["time"] = args.concat_meryl_time
        self.concatmerylSpecParameters["queue"] = args.concat_meryl_queue
        self.concatmerylSpecParameters["mem"] = args.concat_meryl_mem
        self.allParameters ["concat_meryl"] = self.concatmerylSpecParameters

    def storebuildmeryl1SpecParameters(self,args):
        """Updates concat meryl cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.buildmeryl1SpecParameters["name"] = "{rule}_" + args.base_name + "_{base}"
        self.buildmeryl1SpecParameters["qos"] = args.concat_meryl_qos
        self.buildmeryl1SpecParameters["time"] = args.concat_meryl_time
        self.buildmeryl1SpecParameters["queue"] = args.concat_meryl_queue
        self.buildmeryl1SpecParameters["mem"] = args.concat_meryl_mem
        self.allParameters ["build_meryl"] = self.buildmeryl1SpecParameters

  
    def storesmudgeplotSpecParameters(self,args):
        """Updates smudgeplot cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.smudgeplotSpecParameters["name"] = "{rule}_" + args.base_name + "_{base}"
        self.smudgeplotSpecParameters["qos"] = args.smudgeplot_qos
        self.smudgeplotSpecParameters["time"] = args.smudgeplot_time
        self.smudgeplotSpecParameters["queue"] = args.smudgeplot_queue
        self.smudgeplotSpecParameters["mem"] = args.smudgeplot_mem
        self.allParameters ["smudgeplot"] = self.smudgeplotSpecParameters

    def storegenomescopeSpecParameters(self,args):
        """Updates genomescope cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.genomescopeSpecParameters["name"] = "{rule}_" + args.base_name + "_{base}"
        self.genomescopeSpecParameters["qos"] = args.genomescope_qos
        self.genomescopeSpecParameters["time"] = args.genomescope_time
        self.genomescopeSpecParameters["queue"] = args.genomescope_queue
        self.genomescopeSpecParameters["mem"] = args.genomescope_mem
        self.allParameters ["genomescope2"] = self.genomescopeSpecParameters

    def storeFlyeParameters(self,args):
        """Updates flye parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.flyeParameters["Flye cores"] = args.flye_cores
        self.flyeParameters["Flye polishing iterations"] = args.flye_pol_it
        self.flyeParameters["options"] = args.other_flye_opts
        self.allParameters ["Flye"] = self.flyeParameters

    def storeflyeSpecParameters(self,args):
        """Updates flye cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.flyeSpecParameters["name"] = "{rule}_" + args.base_name + "_s" + args.flye_step 
        self.flyeSpecParameters["qos"] = args.flye_qos
        self.flyeSpecParameters["time"] = args.flye_time
        self.flyeSpecParameters["queue"] = args.flye_queue
        self.flyeSpecParameters["mem"] = args.flye_mem
        self.allParameters ["flye"] = self.flyeSpecParameters

    def storeHifiasmParameters(self,args):
        """Updates hifiasm parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.hifiasmParameters["Hifiasm cores"] = args.hifiasm_cores
        self.hifiasmParameters["options"] = args.other_hifiasm_opts
        self.hifiasmParameters["external_purging"] = args.hifiasm_ext_purge
        self.hifiasmParameters["phase"] = args.hifiasm_phasing
        self.allParameters ["Hifiasm"] = self.hifiasmParameters

    def storeHifiasmSpecParameters(self,args):
        """Updates hifiasm cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.hifiasmSpecParameters["name"] = "{rule}_" + args.base_name + "_s" + args.hifiasm_step 
        self.hifiasmSpecParameters["qos"] = args.hifiasm_qos
        self.hifiasmSpecParameters["time"] = args.hifiasm_time
        self.hifiasmSpecParameters["queue"] = args.hifiasm_queue
        self.hifiasmSpecParameters["mem"] = args.hifiasm_mem
        self.allParameters ["hifiasm"] = self.hifiasmSpecParameters

    def storeNextdenovoParameters(self,args):
        """Updates nextdenovo parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nextdenovoParameters["Nextdenovo cores"] = args.nextdenovo_cores
        self.allParameters ["Nextdenovo"] = self.nextdenovoParameters

    def storenextdenovoSpecParameters(self,args):
        """Updates nextdenovo cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nextdenovoSpecParameters["name"] = "{rule}_" + args.base_name + "_s" + args.nextdenovo_step 
        self.nextdenovoSpecParameters["qos"] = args.nextdenovo_qos
        self.nextdenovoSpecParameters["time"] = args.nextdenovo_time
        self.nextdenovoSpecParameters["queue"] = args.nextdenovo_queue
        self.nextdenovoSpecParameters["mem"] = args.nextdenovo_mem
        self.allParameters ["nextdenovo"] = self.nextdenovoSpecParameters

    def storestatsSpecParameters(self,args):
        """Updates stats cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.statsSpecParameters["name"] = "{rule}_" + args.base_name + ".{params.outbase}"
        self.statsSpecParameters["qos"] = args.stats_qos
        self.statsSpecParameters["time"] = args.stats_time
        self.statsSpecParameters["queue"] = args.stats_queue
        self.statsSpecParameters["mem"] = args.stats_mem
        self.allParameters ["get_stats_gfa"] = self.statsSpecParameters

    def storebuscoSpecParameters(self,args):
        """Updates busco cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.buscoSpecParameters["name"] = "{rule}_" + args.base_name + ".{params.buscobase}"
        self.buscoSpecParameters["qos"] = args.busco_qos
        self.buscoSpecParameters["time"] = args.busco_time
        self.buscoSpecParameters["queue"] = args.busco_queue
        self.buscoSpecParameters["mem"] = args.busco_mem
        self.allParameters ["run_busco"] = self.buscoSpecParameters

    def storemerqurySpecParameters(self,args):
        """Updates merqury cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.merqSpecParameters["name"] = "{rule}_" + args.base_name + ".{wildcards.merqbase}"
        self.merqSpecParameters["qos"] = args.merq_qos
        self.merqSpecParameters["time"] = args.merq_time
        self.merqSpecParameters["queue"] = args.merq_queue
        self.merqSpecParameters["mem"] = args.merq_mem
        self.allParameters ["run_merqury"] = self.merqSpecParameters

    def storeminimapSpecParameters(self,args):
        """Updates minimap2 cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.minimapSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}_{wildcards.ext}" 
        self.minimapSpecParameters["qos"] = args.minimap_qos
        self.minimapSpecParameters["time"] = args.minimap_time
        self.minimapSpecParameters["queue"] = args.minimap_queue
        self.minimapSpecParameters["mem"] = args.minimap_mem
        self.allParameters ["align_lr"] = self.minimapSpecParameters

    def storebwaSpecParameters(self,args):
        """Updates BWA cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.bwaSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}" 
        self.bwaSpecParameters["qos"] = args.bwa_qos
        self.bwaSpecParameters["time"] = args.bwa_time
        self.bwaSpecParameters["queue"] = args.bwa_queue
        self.bwaSpecParameters["mem"] = args.bwa_mem
        self.allParameters ["align_illumina"] = self.bwaSpecParameters

    def storeHypoParameters(self,args):
        """Updates hypo parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.hypoParameters["illumina coverage"] = args.ill_cov
        self.hypoParameters["processes"] = args.hypo_processes
        self.hypoParameters["long_reads"] = args.hypo_lr
        self.hypoParameters["options"] = args.hypo_opts
        self.allParameters ["Hypo"] = self.hypoParameters

    def storehypoSpecParameters(self,args):
        """Updates hypo cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.hypoSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.hypoSpecParameters["qos"] = args.hypo_qos
        self.hypoSpecParameters["time"] = args.hypo_time
        self.hypoSpecParameters["queue"] = args.hypo_queue
        self.hypoSpecParameters["mem"] = args.hypo_mem
        self.allParameters ["hypo"] = self.hypoSpecParameters

    def storenextpolishlrSpecParameters(self,args):
        """Updates nextpolish lr cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nextpolishlrSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.nextpolishlrSpecParameters["qos"] = args.nextpolish_lr_qos
        self.nextpolishlrSpecParameters["time"] = args.nextpolish_lr_time
        self.nextpolishlrSpecParameters["queue"] = args.nextpolish_lr_queue
        self.nextpolishlrSpecParameters["mem"] = args.nextpolish_lr_mem
        self.allParameters ["nextpolish_lr"] = self.nextpolishlrSpecParameters

    def storenextpolishsrSpecParameters(self,args):
        """Updates nextpolish sr cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.nextpolishsrSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base}.{wildcards.param}"
        self.nextpolishsrSpecParameters["qos"] = args.nextpolish_sr_qos
        self.nextpolishsrSpecParameters["time"] = args.nextpolish_sr_time
        self.nextpolishsrSpecParameters["queue"] = args.nextpolish_sr_queue
        self.nextpolishsrSpecParameters["mem"] = args.nextpolish_sr_mem
        self.allParameters ["nextpolish_sr"] = self.nextpolishsrSpecParameters

    def storePurgedupsParameters(self,args):
        """Updates purge_dups parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.purgedupsParameters["purgedups_cores"] = args.purgedups_cores
        self.purgedupsParameters["calcuts_options"] = args.calcuts_opts
        self.allParameters ["Purge_dups"] = self.purgedupsParameters

    def storepurgedupsSpecParameters(self,args):
        """Updates purgedups cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.purgedupsSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base_in}"
        self.purgedupsSpecParameters["qos"] = args.purgedups_qos
        self.purgedupsSpecParameters["time"] = args.purgedups_time
        self.purgedupsSpecParameters["queue"] = args.purgedups_queue
        self.purgedupsSpecParameters["mem"] = args.purgedups_mem
        self.allParameters ["purge_dups"] = self.purgedupsSpecParameters

    def storescaffold10XParameters(self,args):
        """Updates 10X scaffolding parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.scaffold10XParameters["tigmint_cores"] = args.tigmint_cores
        self.scaffold10XParameters["tigmint_options"] = args.tigmint_opts
        self.allParameters ["scaffolding_10X"] = self.scaffold10XParameters

    def storescaffold10XSpecParameters(self,args):
        """Updates 10X scaffolding cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.scaffold10XSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.base_in}"
        self.scaffold10XSpecParameters["qos"] = args.tigmint_qos
        self.scaffold10XSpecParameters["time"] = args.tigmint_time
        self.scaffold10XSpecParameters["queue"] = args.tigmint_queue
        self.scaffold10XSpecParameters["mem"] = args.tigmint_mem
        self.allParameters ["scaffolding_10X"] = self.scaffold10XSpecParameters

    def storehicParameters(self,args):
        """Updates HiC related parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.hicParameters["deepseq"] = args.hic_deepseq
        self.hicParameters["subsample"] = args.subsample_hic
        self.hicParameters["add_preseq_opts"] = args.add_preseq_opts
        self.hicParameters["get_pretext"] = args.get_pretext
        self.hicParameters["sort_pretext"] = args.sort_pretext
        self.hicParameters["yahs_cores"] = args.yahs_cores
        self.hicParameters["yahs_mq"] = args.yahs_mq
        self.hicParameters["yahs_opts"] = args.yahs_opts
        self.hicParameters["yahs_contig_ec"] = args.yahs_contig_ec
        self.hicParameters["assembly_qc"] = args.assembly_qc
        self.hicParameters["qc_assemblylen"] = args.hic_qc_assemblylen
        self.hicParameters["align_opts"] = args.hic_map_opts
        self.hicParameters["MQ"] = args.mq
        self.hicParameters["reads_for_blast"] = args.hic_readsblast
        self.hicParameters["blastdb"] = args.blastdb
        self.hicParameters["blast_cores"] = args.blast_cores
        self.allParameters ["HiC"] = self.hicParameters  

    def storeassprepSpecParameters(self,args):
        """Updates assembly prepare cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.assprepSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}"
        self.assprepSpecParameters["qos"] = args.ass_prepare_qos
        self.assprepSpecParameters["time"] = args.ass_prepare_time
        self.assprepSpecParameters["queue"] = args.ass_prepare_queue
        self.assprepSpecParameters["mem"] = args.ass_prepare_mem
        self.allParameters ["assembly_prepare"] = self.assprepSpecParameters

    def storemapHicSpecParameters(self,args):
        """Updates map HiC cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.mapHicSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}"
        self.mapHicSpecParameters["qos"] = args.map_hic_qos
        self.mapHicSpecParameters["time"] = args.map_hic_time
        self.mapHicSpecParameters["queue"] = args.map_hic_queue
        self.mapHicSpecParameters["mem"] = args.map_hic_mem
        self.allParameters ["align_hic"] = self.mapHicSpecParameters

    def storepairtoolsParseSpecParameters(self,args):
        """Updates pairtools cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pairtoolsParseSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}_mq{wildcards.mq}"
        self.pairtoolsParseSpecParameters["qos"] = args.pairtools_qos
        self.pairtoolsParseSpecParameters["time"] = args.pairtools_time
        self.pairtoolsParseSpecParameters["queue"] = args.pairtools_queue
        self.pairtoolsParseSpecParameters["mem"] = args.pairtools_mem
        self.allParameters ["pairtools_processing_parse"] = self.pairtoolsParseSpecParameters

    def storepairtoolsSortSpecParameters(self,args):
        """Updates pairtools cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pairtoolsSortSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}_mq{wildcards.mq}"
        self.pairtoolsSortSpecParameters["qos"] = args.pairtools_qos
        self.pairtoolsSortSpecParameters["time"] = args.pairtools_time
        self.pairtoolsSortSpecParameters["queue"] = args.pairtools_queue
        self.pairtoolsSortSpecParameters["mem"] = args.pairtools_mem
        self.allParameters ["pairtools_processing_sort"] = self.pairtoolsSortSpecParameters

    def storepairtoolsDedupSpecParameters(self,args):
        """Updates pairtools cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pairtoolsDedupSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}_mq{wildcards.mq}"
        self.pairtoolsDedupSpecParameters["qos"] = args.pairtools_qos
        self.pairtoolsDedupSpecParameters["time"] = args.pairtools_time
        self.pairtoolsDedupSpecParameters["queue"] = args.pairtools_queue
        self.pairtoolsDedupSpecParameters["mem"] = args.pairtools_mem
        self.allParameters ["pairtools_processing_dedup"] = self.pairtoolsDedupSpecParameters

    def storepairtoolsSplitSpecParameters(self,args):
        """Updates pairtools cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pairtoolsSplitSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}_mq{wildcards.mq}"
        self.pairtoolsSplitSpecParameters["qos"] = args.pairtools_qos
        self.pairtoolsSplitSpecParameters["time"] = args.pairtools_time
        self.pairtoolsSplitSpecParameters["queue"] = args.pairtools_queue
        self.pairtoolsSplitSpecParameters["mem"] = args.pairtools_mem
        self.allParameters ["pairtools_processing_split"] = self.pairtoolsSplitSpecParameters

    def storeblastSpecParameters(self,args):
        """Updates blast cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.blastSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}"
        self.blastSpecParameters["qos"] = args.blast_qos
        self.blastSpecParameters["time"] = args.blast_time
        self.blastSpecParameters["queue"] = args.blast_queue
        self.blastSpecParameters["mem"] = args.blast_mem
        self.allParameters ["read_screening"] = self.blastSpecParameters

    def storeyahsSpecParameters(self,args):
        """Updates yahs cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.yahsSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}"
        self.yahsSpecParameters["qos"] = args.yahs_qos
        self.yahsSpecParameters["time"] = args.yahs_time
        self.yahsSpecParameters["queue"] = args.yahs_queue
        self.yahsSpecParameters["mem"] = args.yahs_mem
        self.allParameters ["run_yahs"] = self.yahsSpecParameters

    def storegapsSpecParameters(self,args):
        """Updates the get gaps extension cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.gapsSpecParameters["name"] = "{rule}_" + args.base_name
        self.gapsSpecParameters["qos"] = args.stats_qos
        self.gapsSpecParameters["time"] = args.stats_time
        self.gapsSpecParameters["queue"] = args.stats_queue
        self.gapsSpecParameters["mem"] = args.stats_mem
        self.allParameters ["get_extension_gaps"] = self.gapsSpecParameters

    def storeontbgSpecParameters(self,args):
        """Updates the get ONT extension cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.ontbgSpecParameters["name"] = "{rule}_" + args.base_name
        self.ontbgSpecParameters["qos"] = args.telext_qos
        self.ontbgSpecParameters["time"] = args.telext_time
        self.ontbgSpecParameters["queue"] = args.telext_queue
        self.ontbgSpecParameters["mem"] = args.telext_mem
        self.allParameters ["get_extension_cov"] = self.ontbgSpecParameters

    def storetidkSpecParameters(self,args):
        """Updates the run tidk search cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.tidkSpecParameters["name"] = "{rule}_" + args.base_name
        self.tidkSpecParameters["qos"] = args.telext_qos
        self.tidkSpecParameters["time"] = args.telext_time
        self.tidkSpecParameters["queue"] = args.telext_queue
        self.tidkSpecParameters["mem"] = args.telext_mem
        self.allParameters ["tidk_search"] = self.tidkSpecParameters

    def storetelextSpecParameters(self,args):
        """Updates the get telomere extension cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.telextSpecParameters["name"] = "{rule}_" + args.base_name
        self.telextSpecParameters["qos"] = args.telext_qos
        self.telextSpecParameters["time"] = args.telext_time
        self.telextSpecParameters["queue"] = args.telext_queue
        self.telextSpecParameters["mem"] = args.telext_mem
        self.allParameters ["get_extension_telomeres"] = self.telextSpecParameters

    def storepretextSpecParameters(self,args):
        """Updates pretext cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pretextSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}_{wildcards.mq}"
        self.pretextSpecParameters["qos"] = args.pretext_qos
        self.pretextSpecParameters["time"] = args.pretext_time
        self.pretextSpecParameters["queue"] = args.pretext_queue
        self.pretextSpecParameters["mem"] = args.pretext_mem
        self.allParameters ["generate_pretext"] = self.pretextSpecParameters

    def storetpfSpecParameters(self,args):
        """Updates tpf cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.tpfSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}"
        self.tpfSpecParameters["qos"] = args.tpf_qos
        self.tpfSpecParameters["time"] = args.tpf_time
        self.tpfSpecParameters["queue"] = args.tpf_queue
        self.tpfSpecParameters["mem"] = args.tpf_mem
        self.allParameters ["get_tpf"] = self.tpfSpecParameters

    def storediploidSpecParameters(self,args):
        """Updates get diploid cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.diploidSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}"
        self.diploidSpecParameters["qos"] = args.tpf_qos
        self.diploidSpecParameters["time"] = args.tpf_time
        self.diploidSpecParameters["queue"] = args.tpf_queue
        self.diploidSpecParameters["mem"] = args.tpf_mem
        self.allParameters ["generate_diploid"] = self.diploidSpecParameters

    def storeepretextSpecParameters(self,args):
        """Updates pretext cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.epretextSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}_{wildcards.mq}"
        self.epretextSpecParameters["qos"] = args.pretext_qos
        self.epretextSpecParameters["time"] = args.pretext_time
        self.epretextSpecParameters["queue"] = args.pretext_queue
        self.epretextSpecParameters["mem"] = args.pretext_mem
        self.allParameters ["add_extensions_pretext"] = self.epretextSpecParameters

    def storeqcstatsSpecParameters(self,args):
        """Updates qcstats cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.qcstatsSpecParameters["name"] = "{rule}_" + args.base_name + "_{wildcards.name}_{wildcards.mq}"
        self.qcstatsSpecParameters["qos"] = args.qcstats_qos
        self.qcstatsSpecParameters["time"] = args.qcstats_time
        self.qcstatsSpecParameters["queue"] = args.qcstats_queue
        self.qcstatsSpecParameters["mem"] = args.qcstats_mem
        self.allParameters ["qc_statistics"] = self.qcstatsSpecParameters

    def storeFinalizeParameters(self,args):
        """Updates finalize parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.finalizeParameters["final Evaluations"] = args.final_evals
        self.finalizeParameters["BUSCO lineage"] = args.busco_lineage
        self.finalizeParameters["Merqury db"] = args.merqury_db
        self.finalizeParameters["Merqury plot opts"] = args.merqury_plot_opts
        self.finalizeParameters["Meryl K"] = args.meryl_k
        self.finalizeParameters["Meryl threads"] = args.meryl_threads
        self.finalizeParameters["Meryl reads"] = args.meryl_reads
        self.allParameters ["Finalize"] = self.finalizeParameters

    def storefinalizeSpecParameters(self,args):
        """Updates finalize cluster spec parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.finSpecParameters["name"] = "{rule}_" + args.base_name
        self.finSpecParameters["qos"] = args.fin_qos
        self.finSpecParameters["time"] = args.fin_time
        self.finSpecParameters["queue"] = args.fin_queue
        self.finSpecParameters["mem"] = args.fin_mem
        self.allParameters ["finalize"] = self.finSpecParameters

    def storeWildcardParameters(self,args):
        """Updates wildcard parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.wildcardParameters["ONT_wildcards"] = args.ONT_wildcards
        self.wildcardParameters["hifi_wildcards"] = args.hifi_wildcards
        self.wildcardParameters["illumina_wildcards"] = args.illumina_wildcards
        self.wildcardParameters["10X_wildcards"] = args.r10X_wildcards
        self.wildcardParameters["HiC_wildcards"] = args.hic_wildcards
        self.allParameters ["Wildcards"] = self.wildcardParameters

    def create_nextdenovo_config(self, args):
        nextdenovo_config.add_section('General')
        nextdenovo_config.set('General', 'job_type', args.nextdenovo_type)
        nextdenovo_config.set('General', 'task', args.nextdenovo_task)
        nextdenovo_config.set('General', 'rewrite', args.nextdenovo_rewrite)
        nextdenovo_config.set('General', 'parallel_jobs', str(args.nextdenovo_parallel_jobs))

        if args.keepintermediate == True:
          nextdenovo_config.set('General', 'deltmp', 'no')
        else:
          nextdenovo_config.set('General', 'deltmp', 'yes')

        if re.search("nano", args.lr_type):
          nextdenovo_config.set('General', 'read_type', 'ont')
          if args.lr_type == "nano-raw":
            nextdenovo_config.set('General', 'input_type', 'raw')
          else:
            nextdenovo_config.set('General', 'input_type', 'corrected')
        elif re.search("hifi", args.lr_type):
          nextdenovo_config.set('General', 'read_type', 'hifi')
        else:
          print ("Need to select proper read type for running Nextdenovo")

        nextdenovo_config.set('General',  'workdir', args.nextdenovo_dir)
        nextdenovo_config.set('General',  'input_fofn', args.nextdenovo_dir + 'long_reads.fofn')

        nextdenovo_config.add_section('correct_option')
        nextdenovo_config.set('correct_option',  'read_cutoff', args.nextdenovo_minreadlen)
        nextdenovo_config.set('correct_option',  'genome_size', args.genome_size)
        nextdenovo_config.set('correct_option',  'seed_depth', str(args.nextdenovo_seeddepth))
        nextdenovo_config.set('correct_option',  'seed_cutoff', str(args.nextdenovo_seedcutoff))
        nextdenovo_config.set('correct_option',  'blocksize', args.nextdenovo_blocksize)
        nextdenovo_config.set('correct_option',  'pa_correction', str(args.nextdenovo_pa_correction))
        nextdenovo_config.set('correct_option',  'minimap2_options_raw', args.nextdenovo_minimap_raw)
        nextdenovo_config.set('correct_option',  'sort_options', args.nextdenovo_sort)
        nextdenovo_config.set('correct_option',  'correction_options', args.nextdenovo_correction_opts)

        nextdenovo_config.add_section('assemble_option')
        nextdenovo_config.set('assemble_option',  'minimap2_options_cns', args.nextdenovo_minimap_cns)
        nextdenovo_config.set('assemble_option',  'minimap2_options_map', args.nextdenovo_minimap_map)
        nextdenovo_config.set('assemble_option',  'nextgraph_options', args.nextdenovo_nextgraph_opt)


#1.Create object class Configuration File
configManager = CreateConfigurationFile()
specManager = CreateConfigurationFile()
NDConfManager = CreateConfigurationFile()

#2.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="create_configuration_file",
                description="Create a configuration json file for the assembly pipeline."
                )     
#2.1 Updates arguments and parsing
configManager.register_parameter(parser)

args = parser.parse_args()

#2.2 Check Parameters
configManager.check_parameters(args)

#3. store arguments to super map structure
configManager.storeGeneralParameters(args)
configManager.storeInputParameters(args)
configManager.storeOutputParameters(args)
configManager.storeFiltlongParameters(args)
configManager.storeTrimgaloreParameters(args)
configManager.storekraken2Parameters(args)
configManager.storeFlyeParameters(args)
configManager.storeHifiasmParameters(args)
configManager.storeNextdenovoParameters(args)
configManager.storeHypoParameters(args)
configManager.storePurgedupsParameters(args)
configManager.storescaffold10XParameters(args)
configManager.storehicParameters(args)
configManager.storeFinalizeParameters(args)
configManager.storeWildcardParameters(args)

specManager.storeallSpecParameters(args)

if args.ONT_wildcards != None or args.hifi_wildcards or args.illumina_wildcards != None or args.r10X_wildcards or args.hic_wildcards:
  specManager.storeconcatreadsSpecParameters(args)
if args.subsample_hic:
   specManager.storesubsampleHICSpecParameters(args)
if args.ONT_wildcards != None or args.ONT_reads != None:
  if not os.path.exists(args.ONT_filtered):
    specManager.storefiltlongSpecParameters(args)
specManager.storenanoplotSpecParameters(args)
if args.illumina_dir != None:
  specManager.storetrimgaloreSpecParameters(args)
if args.r10X_wildcards != None:
  specManager.storelongrangerSpecParameters(args)

if args.hifi_reads or args.ONT_reads:
  specManager.storebam2fqSpecParameters(args)
if args.merqury_db:
  if not os.path.exists(args.merqury_db):
    specManager.storebuildmeryl1SpecParameters(args)
    specManager.storebuildmerylSpecParameters(args)
    specManager.storeconcatmerylSpecParameters(args)
  if args.run_smudgeplot:
    specManager.storesmudgeplotSpecParameters(args)
  specManager.storegenomescopeSpecParameters(args)
  specManager.storemerqurySpecParameters(args)

if args.final_evals:
  specManager.storestatsSpecParameters(args)
  specManager.storebuscoSpecParameters(args)
  specManager.storetidkSpecParameters(args)

if args.run_kraken2 == True:
  specManager.storekraken2SpecParameters(args) 

if args.run_flye == True:
  specManager.storeflyeSpecParameters(args)
if args.run_hifiasm == True:
  specManager.storeHifiasmSpecParameters(args)

if args.run_nextdenovo == True:
  specManager.storenextdenovoSpecParameters(args)
  with open(args.ndconfFile, 'w') as ndconf:
    configManager.create_nextdenovo_config(args)
    nextdenovo_config.write(ndconf)
if args.nextpolish_ill_rounds > 0 or args.hypo_rounds >0:
  specManager.storebwaSpecParameters(args)
if args.nextpolish_ont_rounds > 0 or args.hypo_rounds >0 or args.run_purgedups == True or args.get_pretext:
  specManager.storeminimapSpecParameters(args)
if args.hypo_rounds > 0:
  specManager.storehypoSpecParameters(args)
if args.nextpolish_ont_rounds > 0:
  specManager.storenextpolishlrSpecParameters(args)
if args.nextpolish_ill_rounds > 0:
  specManager.storenextpolishsrSpecParameters(args)

if args.run_purgedups == True or args.hifiasm_ext_purge:
  specManager.storepurgedupsSpecParameters(args)
if args.run_tigmint == True:
  specManager.storescaffold10XSpecParameters(args)
if args.hic_dir:
  specManager.storeassprepSpecParameters(args)
  specManager.storemapHicSpecParameters(args)
  specManager.storepairtoolsParseSpecParameters(args)
  specManager.storepairtoolsSortSpecParameters(args)
  specManager.storepairtoolsDedupSpecParameters(args)
  specManager.storepairtoolsSplitSpecParameters(args)
  specManager.storeqcstatsSpecParameters(args)
  specManager.storeblastSpecParameters(args)
  if args.run_yahs == True:
    specManager.storeyahsSpecParameters(args)
  if args.get_pretext == True:
    specManager.storepretextSpecParameters(args)
    specManager.storetpfSpecParameters(args) 
    specManager.storegapsSpecParameters(args)
    specManager.storeontbgSpecParameters(args)
    specManager.storetelextSpecParameters(args)
    specManager.storeepretextSpecParameters(args)
    specManager.storediploidSpecParameters(args)
specManager.storefinalizeSpecParameters(args)
    
#4. Store JSON file
with open(args.configFile, 'w') as of:
    json.dump(configManager.allParameters, of, indent=2)
with open(args.specFile, 'w') as of:
    json.dump(specManager.allParameters, of, indent=2)
