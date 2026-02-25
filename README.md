# CLAWS (CNAG's Long-read Assembly Workflow in Snakemake)
 Snakemake Pipeline used for de novo genome assembly @CNAG. It has been developed for Snakemake v6.0.5.

It accepts Oxford Nanopore Technologies (ONT) reads, PacBio HFi reads, illumina paired-end data, illumina 10X data and Hi-C reads. It does the preprocessing of the reads, assembly, polishing, purge_dups, scaffolding, different evaluation steps and generation of pretext files for curation. Default behavior is to preprocess the reads, assemble with (1. Flye + hypo + purgedups + yahs) and (2. Hifiasm + yahs), evaluate the resulting assemblies with BUSCO, MERQURY and gfastats and produce high and low resolution pretext maps with mq10 and mq0.
It needs a config file and a spec file (json file with instructions on which resources should slurm use for each of the jobs). Both files are created by the script "create_config_assembly.py" that is located in the bin directory. To check all the options accepted by the script, do:

```
bin/create_config_assembly.py -h
```

Once the 2 config files are produced, the pipeline can be launched using snakemake like this:

``snakemake --notemp -j 999 --snakefile assembly_pipeline.smk --configfile assembly.config --is --cluster-conf assembly.spec --use-conda --use-envmodules``

If you are using an HPC cluster, please check how should you run snakemake to launch the jobs to the cluster. 

Most of the tools used will be installed via conda using the environments of the "envs" directory after providing the "--use-conda" option to snakemake. However, a few optional tools cannot be installed via conda and will have to be available in your PATH, or as a module in the cluster. Those tools are:

- NextDenovo/2.5.0
- NextPolish/1.4.1

If you want to assemble with NextDenovo or polish with NextPolish you'll need to install them first, but if you just want to run CLAWS with the default assemblers and polishers, do not worry about this step. 

# How to provide input data:

There are several ways of providing the reads. All types of data are optional and depending on your inputs and steps selection, you can personalize your run. 

### 1- ONT reads

1.1 Using the option ``--ont-dir {DIR}`` in create_config_assembly.py.

If you do so, it will look for all the files in the directory that end in '.fastq.gz' and will add the basenames to "ONT_wildcards". These wildcards will be processed by the pipeline that will:

- Concatenate all the files into a single file

- Run filtlong with the default or specified parameters. 

- Use the resulting file for assembly, polishing and/or purging, as well as for building the meryldb.

You can also specify the basenames of the files that you want to use with the ``--ont-list `` option. In this case, the pipeline will use the wildcards that you're providing instead of merging all the files in the directory.

1.2 Using the option ```--ont-reads {FILE}``` in create_config_assembly.py.

File with all the ONT reads. It can either be in fastq, fasta or bam format. It will be filtered by filtlong and the resulting file will be used for assembly, polishing, purging and building a kmer database.

### 2-Pacbio HiFi data

2.1 Using the option ``--hifi-dir {DIR}`` in create_config_assembly.py.

If you do so, it will look for all the files in the directory that end in '.fastq.gz' and will add the basenames to "ONT_wildcards". These wildcards will be processed by the pipeline that will:

- Concatenate all the files into a single file

- Use the resulting file for assembly, polishing and/or purging, as well as for building the meryldb.

You can also specify the basenames of the files that you want to use with the ``--hifi-list `` option. In this case, the pipeline will use the wildcards that you're providing instead of merging all the files in the directory.

1.2 Using the option ```--hifi-reads {FILE}``` in create_config_assembly.py.

File with all the HiFi reads. It can be either in fastq, fasta or bam format. It will be used for assembly, purging and building a meryldb.

### 3- Illumina short-read data

CLAWS can optionally take Illumina reads and use them for building meryl_dbs and/or polishing. 

3.1 Using the ``--illumina-dir {DIR}`` option, that will look for all the files in the directory that end in '.1.fastq.gz' and will add the basenames to "illumina_wildcards". These wildcards will be processed by the pipeline that will: 

- Trim adaptors with Trimgalore

- Concatenate all the trimmed *.1.fastq.gz and the *2.fastq.gz into one file per pair. 

- The resulting reads will be used for building meryldbs and polishing. 

3.2 Using the ``--processed-illumina`` option. If the directory exists and contains files, the pipeline will look for all the files in the directory that end in '.1.fastq.gz' and will add the basenames to "illumina_wildcards". These wildcards will be processed by the pipeline that will:

- Concatenate all the trimmed *.1.fastq.gz and the *2.fastq.gz in one file per pair. 

- The resulting reads will be used for building meryldbs and polishing. 

3.3 Using the ``--pe1 {FILE} and --pe2 {FILE}`` options. That will consider that these are the paired files containing all the illumina reads ready to be used and will build meryldbs and polish with them.

### 4- Hi-C data

4.1 Using the ``--hic-dir {DIR}`` option, that will look for all the files in the directory that end in '.1.fastq.gz' and will add the basenames to "illumina_wildcards". These wildcards will be concatenated into a single file that can be used by the pipeline for phasing, scaffolding and contact map generation. 

### 5-Illumina 10X-linked data

5.1 Using the  ```--raw-10X {DIR:list}``` option. 

Dictionary with 10X raw read directories, it has to be the mkfastq dir. You must specify as well the sampleIDs from this run. Example: '{"mkfastq-                        dir":"sample1,sample2,sample3"}'...

It will take each basename in the list to get the fastqs from the corresponding directory and run longranger on each sample. Afterwards, it will build meryldbs for each "barcoded" file. Finally, it will concatenate all the meryldbs and "barcoded" files. Resulting "barcoded" file will be used for polishing. 

5.2 Using the ``--processed-10X {DIR}`` parameter. 

This directory can already be there or be produced by the pipeline as described in step 2.1. Once all the "barcoded" fastq files are there, meryldbs will be built for each "barcoded" file.  Finally, it will concatenate all the meryldbs and "barcoded" files. Resulting "barcoded" file will be used for polishing. 

5.3 Using the ``--10X`` option. 

The argument to this is the path to the concatenated ".barcoded" file that needs to be used for polishing. If the pre-concatenated files are not given, meryldbs will be directly generated with this file, but it may run out of memory. 

### 6- Input assemblies

If you want to polish an already assembled assembly, you can give it to the pipeline by using the option ``--assembly-in ASSEMBLY_IN [ASSEMBLY_IN ...]
                        Dictionary with assemblies that need to be polished but not assembled and directory where they should
                        be polished. Example: '{"assembly1":"polishing_dir1"}' '{"assembly2"="polishing_dir2"}' ...``
			
If you want to start CLAWS after polishing on an already existing assembly, you can give it to the pipeline by using the option ``--postpolish-assemblies POSTPOLISH_ASSEMBLIES [POSTPOLISH_ASSEMBLIES ...]
                        Dictionary with assemblies for which postpolishing steps need to be run but that are not assembled and
                        base step for the directory where the first postpolishing step should be run. Example:
                        '{"assembly1":"s04.1_p03.1"}' '{"assembly2"="s04.2_p03.2"}' ...``

To evaluate and produce the final pretext file on a curated assembly, use ``--curated-assemblies CURATED_ASSEMBLIES [CURATED_ASSEMBLIES ...]
                        Dictionary with assemblies that have already been curated. Evaluations and read alignment will be perforder. Example:
                        '{"assembly1":"s04.1_p03.1"}' '{"assembly2":"s04.2_p03.2"}' ...``



# Description of implemented rules

1- Preprocessing:
	
- **Read concatenation:**

``zcat {input.fastqs} | pigz -p {threads} -c  > {output.final_fastq}``
	
- **Longranger for 10X reads**: it uses the Longranger version installed in the path specified in the configfile

``longranger basic --id={params.sample} --sample={params.sample} --fastqs={input.mkfastq_dir} --localcores={threads}``

- **Trimgalore:** By default it gives the ``--max_n 0 --gzip -q 20 --paired --retain_unpaired`` options, but it can be changed with the ``--trim-galore-opts `` argument. 

- **Filtlong:** it uses the Filtlong version installed in the path specified in the configfile. By default it gives the min_length and min_mean_q parameters, but extra parameters can be added with the ``--filtlong-opts`` option.
	
- **Build meryldb**: it uses the merqury conda environment specified in the configfile. It takes as argument the `--meryl-k` value that needs to be estimated first for the genome size. It can run either on the illumina reads, the ont reads or both, default behaviour is both. 
	
- **Align ONT (Minimap2):** it aligns the reads using minimap2 and outputs the alignment either in bam or in paf.gz formats. It uses the minimap2 conda environment specified in the configfile.

- **Align Illumina (BWA-MEM):** it aligns the reads with BWA-mem and outputs a bam file.

- **Align Hi-C (Chromap):** it aligns the reads with Chromap and outputs a bam file.

2- Assembly

- **Hifiasm (default)**. It is run by default, if you don't want the pipeline to run it, you can give `--no-hifiasm` option when creating the config. It uses the conda environment specified in the config. It can run in phasing mode if the option "--phase-hifiasm" is given. Extra options can be provided with the `--other-hifiasm-opts`. You need to provide a telomeric motif  If you want purgedups to be run on the output, please give the "purge-hifiasm" option. 

- **Flye (default)**. It is run by default, if you don't want the pipeline to run it, you can give `--no-flye` option when creating the config. It uses the conda environment specified in the config. By default it is set to 2 polishing iterations and gives the genome-size estimate that has been given when creating the config. Extra options can be provided with the `--flye-opts`.
	
- **Nextdenovo (if ``run-nextdenovo``):** It uses the cluster module specified in the config. If nextdenovo option is turned on, the create_config script will also create the nextdenovo config file. Check the create_config help to see which options can be modified on it. 

3- Polishing

- **Hypo (default):** It is the polisher that the pipeline uses by default, it can be turned off specifying ``--hypo-rounds 0`` when creating the config. If selected, the reads will be aligned in previous rules and then hypo will be run, it requires illumina data. It uses the conda environment specified in the config. It only runs for ont assemblies and it doesn't run after hifiasm. 
	
- **Nextpolish ont (if turned on):** to run nextpolish with ONT reads, specify ``--nextpolish-ont-rounds`` and the number of rounds you want to run of it. 

- **Nextpolish illumina (if turned on):** to run nextpolish with ONT reads, specify ``--nextpolish-ill-rounds`` and the number of rounds you want to run of it. 

4- Post-assembly

- **Purge_dups (by default):** select ``--no-purgedups`` if you don't want to run it. If no manual cutoffs are given, it'll run purgedups with automatic cutoffs and then will rerun it selecting the mean cutoff as 0.75\*cov. It uses the version installed in the cluster module specified in the config. It doesn't run after hifiasm unless "--purge-hifiasm" option is given. 

- **Yahs (by default):** select ``--no-yahs``if you do not want to run it. By default it uses mq10 and the no-contig-ec option, this can be changed by giving to the config the options "yahs-mq" and "yahs-contig-ec". "--yahs-opts" to change any other options.

5- Evaluations
	
- **Merqury:** It runs on each assembly produced by the pipeline. Hap1 and hap2 files are evaluated in the same run. 
	
- **Busco:** It uses the conda environment specified in the config as well as the parameters specified, you need to provide the lineage directory. 
	
- **gfastastas:**  Gfastastats with the nstar option runs on every assembly produced or given as input to CLAWS.

# Examples 

We will show a few examples on how to create the config files and perform a snakemake dry-run that will show which commands will be performed. When you have your case ready, you can remove the "-np" option and let snakemake execute the pipeline

1- Assemble a genome with ONT + Hi-C using default options and providing the phasing option to hifiasm 

````
CLAWS/bin/create_config_assembly.py  --configFile configs/qqChaOliv.config --specFile configs/qqChaOliv.spec --basename qqChaOliv --species "Chaetopelma olivaceum"  --genome-size 5g --ont-dir reads/ont/ --busco-lin busco_downloads/lineages/arachnida_odb12/ --illumina-dir reads/illumina/ --merqury-db qqChaOliv.meryl --meryl-k 20  --filtlong-opts " --target_bases 300000000000" --hic-dir reads/hic/   --telo ACCCCG  --phase-hifiasm
Genome size is 5000.0 megabases
Warning: Long-read type is set to nano-hq, make sure this is your data type.
qqChaOliv.meryl not found, the pipeline will create it

snakemake --notemp -j 999 --snakefile CLAWS/bin/assembly_pipeline.smk --configfile configs/qqChaOliv.config  --cluster-conf configs/qqChaOliv.spec   --cluster "python3 Snakemake-CNAG/sbatch-cnag.py {dependencies}" --is  --use-conda --use-envmodules --conda-prefix conda/assembly_pipeline_envs/ -np 
Building DAG of jobs...
Job counts:
	count	jobs
	17	add_extensions_pretext
	9	align_hic_chromap
	1	align_illumina
	11	align_lr
	1	all
	4	assembly_prepare
	3	build_meryl_db
	1	concat_meryl
	2	concat_reads
	8	filter_chromap
	1	filtlong
	1	finalize
	1	flye
	1	generate_diploid
	17	generate_pretext
	1	genomescope2
	9	get_extension_cov
	9	get_extension_gaps
	10	get_stats_gfa
	5	get_tpf
	1	hifiasm
	1	hypo
	2	nanoplot
	1	purge_dups
	10	run_busco
	8	run_merqury
	4	run_yahs
	1	smudgeplot
	11	tidk_search
	1	trim_galore
	152
````

This will make a meryl-db from both Illumina and ONT reads. The kmer histogram produced will be used to run genomescope and smudgeplot. 
It will assemble the ONT reads with 2 different paths: 1. Flye + hypo + purge_dups + yahs and 2. Hifiasm with phasing option + yahs. 
PretextMaps will be generated for all pre and post-scaffolding assemblies, with mapping qualities 10 and 0 and in low and High resolution mode. 
Every assembly produced will be evaluated with gfastastats, busco and merqury. 

2- Assemble a genome with HiFi + Hi-C and purge the Hifiasm produced assemblies. 

````
CLAWS/bin/create_config_assembly.py --configFile configs/ilEumAren.config --specFile configs/ilEumAren.spec --basename ilEumAren --species "Eumannia arenbergeri" --genome-size 1g --merqury-db ilEumAren.meryl --meryl-k 19 --hifi-reads reads/hifi/reads.bam --busco-lin busco_downloads/lineages/lepidoptera_odb12/ --lr-type pacbio-hifi --hic-dir reads/hic/ --telo TTAGG --purge-hifiasm
Genome size is 1000.0 megabases
Warning: Long-read type is set to pacbio-hifi, make sure this is your data type.
ilEumAren.meryl not found, the pipeline will create it

snakemake --notemp -j 999 --snakefile CLAWS/bin/assembly_pipeline.smk --configfile configs/ilEumAren.config   --cluster-conf configs/ilEumAren.spec  --cluster "python3 Snakemake-CNAG/sbatch-cnag.py {dependencies}" --is  --use-conda --use-envmodules --conda-prefix conda/assembly_pipeline_envs/ -np --reason
Building DAG of jobs...
Job counts:
	count	jobs
	17	add_extensions_pretext
	9	align_hic_chromap
	13	align_lr
	1	all
	4	assembly_prepare
	1	bam2fastq
	1	build_meryl
	8	filter_chromap
	1	finalize
	1	flye
	1	generate_diploid
	17	generate_pretext
	1	genomescope2
	9	get_extension_cov
	9	get_extension_gaps
	12	get_stats_gfa
	5	get_tpf
	1	hifiasm
	1	nanoplot
	4	purge_dups
	12	run_busco
	9	run_merqury
	4	run_yahs
	1	smudgeplot
	13	tidk_search
	155

````
This will make a meryl-db from the HiFi reads. The kmer histogram produced will be used to run genomescope and smudgeplot. 
It will assemble the HiFi reads with 2 different paths: 1. Flye + purge_dups + yahs and 2. Hifiasm + purge_dups + yahs. 
PretextMaps will be generated for all pre and post-scaffolding assemblies, with mapping qualities 10 and 0 and in low and High resolution mode. 
Every assembly produced will be evaluated with gfastastats, busco and merqury. 

3- Evaluate and produce pretextmaps for an already curated diploid assembly

```
CLAWS/bin/create_config_assembly.py --configFile configs/ilEumAren.CM.config --specFile configs/ilEumAren.CM.spec --basename ilEumAren --species "Eumannia arenbergeri" --genome-size 1g --merqury-db ilEumAren.meryl --hifi-reads reads/hifi/reads.bam --busco-lin busco_downloads/lineages/lepidoptera_odb12/ --lr-type pacbio-hifi --hic-dir reads/hic/ --telo TTAGG  --no-flye --no-hifiasm --no-purgedups --hypo-rounds 0 --no-yahs  --no-smudgeplot --curated-assemblies '{"assemblies/ilEumAren1.1.hap2.fa":"assemblies", "assemblies/ilEumAren1.1.hap2.fa":"assemblies"}' --mq0
Genome size is 1000.0 megabases
Warning: Long-read type is set to pacbio-hifi, make sure this is your data type.

 snakemake --notemp -j 999 --snakefile CLAWS/bin/assembly_pipeline.smk  --configfile configs/ilEumAren.CM.config  --cluster-conf configs/ilEumAren.CM.spec --cluster "python3 Snakemake-CNAG/sbatch-cnag.py {dependencies}" --is  --use-conda --use-envmodules --conda-prefix conda/assembly_pipeline_envs/  $PWD/evaluations/ilEumAren.stats_summary.txt -np --reason
Building DAG of jobs...
Job counts:
	count	jobs
	1	add_extensions_pretext
	1	align_hic_chromap
	1	align_lr
	1	bam2fastq
	1	build_meryl
	1	finalize
	1	generate_pretext
	1	get_extension_cov
	1	get_extension_gaps
	1	get_stats_gfa
	1	get_tpf
	1	run_busco
	2	run_merqury
	1	tidk_search
	15
```

This will evaluate the curated assemblies with gfastastats, busco and merqury. It will also produce low and hogh resolution pretextfiles with mapping quality 0.

# Description of all options
```
CLAWS/bin/create_config_assembly.py
usage: create_configuration_file [-h] [--configFile configFile] [--specFile specFile] [--ndconfFile ndconfFile] [--keep-intermediate] [--lr-type lr_type] [--basename base_name] [--species species]
                                 [--genome-size genome_size] [--ploidy ploidy] [--telo telo_string] [--no-flye] [--no-hifiasm] [--run-nextdenovo] [--nextpolish-ont-rounds nextpolish_ont_rounds]
                                 [--nextpolish-ill-rounds nextpolish_ill_rounds] [--hypo-rounds hypo_rounds] [--no-purgedups] [--no-yahs] [--no-smudgeplot] [--run-tigmint] [--run-kraken2]
                                 [--genomescope-opts genomescope_additional] [--preprocess-lr-step PREPROCESS_LR_STEP] [--preprocess-10X-step PREPROCESS_10X_STEP]
                                 [--preprocess-illumina-step PREPROCESS_ILLUMINA_STEP] [--preprocess-hic-step PREPROCESS_HIC_STEP] [--flye-step FLYE_STEP] [--hifiasm-step HIFIASM_STEP]
                                 [--nextdenovo-step NEXTDENOVO_STEP] [--concat-cores concat_cores] [--minimap2-cores minimap2_cores] [--bwa-cores bwa_cores] [--hypo-cores hypo_cores]
                                 [--nextpolish-cores nextpolish_cores] [--pairtools-parse-cores pairtools_parse_cores] [--pairtools-sort-cores pairtools_sort_cores]
                                 [--pairtools-dedup-cores pairtools_dedup_cores] [--pairtools-split-cores pairtools_split_cores] [--busco-cores busco_cores] [--longranger-cores longranger_cores]
                                 [--longranger-path longranger_path] [--scripts-dir SCRIPTS_DIR] [--ont-dir ONT_DIR] [--hifi-dir HIFI_DIR] [--illumina-dir ILLUMINA_DIR] [--hic-dir HIC_DIR]
                                 [--raw-10X RAW_10X [RAW_10X ...]] [--ont-reads ONT_READS] [--hifi-reads HIFI_READS] [--pe1 PE1] [--pe2 PE2] [--10X R10X] [--processed-illumina PROCESSED_ILLUMINA]
                                 [--processed-10X PROCESSED_10X] [--ont-filt ONT_FILTERED] [--assembly-in ASSEMBLY_IN [ASSEMBLY_IN ...]]
                                 [--postpolish-assemblies POSTPOLISH_ASSEMBLIES [POSTPOLISH_ASSEMBLIES ...]] [--curated-assemblies CURATED_ASSEMBLIES [CURATED_ASSEMBLIES ...]]
                                 [--pipeline-workdir PIPELINE_WORKDIR] [--preprocess-lr PREPROCESS_LR] [--concat-hic-dir CONCAT_HIC_DIR] [--flye-dir FLYE_DIR] [--nextdenovo-dir NEXTDENOVO_DIR]
                                 [--hifiasm-dir HIFIASM_DIR] [--flye-polishing-dir POLISH_FLYE_DIR] [--nextdenovo-polishing-dir POLISH_NEXTDENOVO_DIR] [--eval-dir eval_dir] [--stats-out stats_out]
                                 [--hic-qc-dir hic_qc_dir] [--filtlong-minlen filtlong_minlen] [--filtlong-min-mean-q filtlong_min_mean_q] [--filtlong-opts filtlong_opts]
                                 [--trim-galore-opts trim_galore_opts] [--trim-Illumina-cores Trim_Illumina_cores] [--kraken2-db kraken2_db] [--kraken2-kmer kraken2_kmers]
                                 [--kraken2-opts additional_kraken2_opts] [--kraken2-cores kraken2_threads] [--flye-cores flye_cores] [--flye-polishing-iterations flye_pol_it]
                                 [--other-flye-opts other_flye_opts] [--hifiasm-cores hifiasm] [--other-hifiasm-opts other_hifiasm_opts] [--purge-hifiasm] [--phase-hifiasm]
                                 [--nextdenovo-cores nextdenovo_cores] [--nextdenovo-jobtype nextdenovo_type] [--nextdenovo-task nextdenovo_task] [--nextdenovo-rewrite nextdenovo_rewrite]
                                 [--nextdenovo-parallel_jobs nextdenovo_parallel_jobs] [--nextdenovo-minreadlen nextdenovo_minreadlen] [--nextdenovo-seeddepth nextdenovo_seeddepth]
                                 [--nextdenovo-seedcutoff nextdenovo_seedcutoff] [--nextdenovo-blocksize nextdenovo_blocksize] [--nextdenovo-pa-correction  nextdenovo_pa_correction]
                                 [--nextdenovo-minimap_raw nextdenovo_minimap_raw] [--nextdenovo-minimap_cns nextdenovo_minimap_cns] [--nextdenovo-minimap_map nextdenovo_minimap_map]
                                 [--nextdenovo-sort nextdenovo_sort] [--nextdenovo-correction_opts nextdenovo_correction_opts] [--nextdenovo-nextgraph_opt nextdenovo_nextgraph_opt] [--sr-cov ill_cov]
                                 [--hypo-proc hypo_processes] [--hypo-no-lr] [--hypo-opts hypo_opts] [--purgedups-cores purgedups_cores] [--purgedups-calcuts-opts calcuts_opts]
                                 [--tigmint-cores tigmint_cores] [--tigmint-opts tigmint_opts] [--hic-qc] [--subsample-hic] [--add-preseq-opts ADD_PRESEQ_OPTS] [--no-pretext] [--sort-pretext SORT_PRETEXT]
                                 [--assembly-qc assembly_qc] [--yahs-cores yahs_cores] [--yahs-mq yahs_mq] [--yahs-contig-ec] [--yahs-opts yahs_opts] [--hic-map-opts hic_map_opts] [--mq mq [mq ...]]
                                 [--hic-qc-assemblylen hic_qc_assemblylen] [--blast-cores blast_cores] [--hic-blastdb blastdb] [--hic-readsblast hic_readsblast] [--no-final-evals]
                                 [--busco-lin busco_lineage] [--merqury-db merqury_db] [--merqury-plot-opts merqury_plot_opts] [--meryl-k meryl_k] [--meryl-threads meryl_threads]
                                 [--meryl-reads meryl_reads [meryl_reads ...]] [--ont-list ONT_wildcards] [--hifi-list hifi_wildcards] [--illumina-list illumina_wildcards] [--r10X-list r10X_wildcards]
                                 [--hic-list hic_wildcards]

Create a configuration json file for the assembly pipeline.

options:
  -h, --help            show this help message and exit

General Parameters:
  --configFile configFile
                        Configuration JSON to be generated. Default assembly.config
  --specFile specFile   Cluster specifications JSON fileto be generated. Default assembly.spec
  --ndconfFile ndconfFile
                        Name pf the nextdenovo config file. Default nextdenovo.config
  --keep-intermediate   Set this to True if you do not want intermediate files to be removed. Default False
  --lr-type lr_type     Type of long reads (options are: pacbio-raw, pacbio-corr, pacbio-hifi, nano-raw, nano-corr, nano-hq). Default nano-hq
  --basename base_name  Base name for the project. Default None
  --species species     Name of the species to be assembled. Default None
  --genome-size genome_size
                        Approximate genome size. Example: 615m or 2.6g. Default None
  --ploidy ploidy       Expected ploidy. Default 2
  --telo telo_string    Expected telomere string. Default None
  --no-flye             Give this option if you do not want to run Flye.
  --no-hifiasm          Give this option if you do not want to run Hifiasm.
  --run-nextdenovo      Give this option if you do want to run Nextdenovo.
  --nextpolish-ont-rounds nextpolish_ont_rounds
                        Number of rounds to run the Nextpolish with ONT step. Default 0
  --nextpolish-ill-rounds nextpolish_ill_rounds
                        Number of rounds to run the Nextpolish with illumina step. Default 0
  --hypo-rounds hypo_rounds
                        Number of rounds to run the Hypostep. Default 1
  --no-purgedups        Give this option if you do not want to run Purgedups on the Flye and Nextdenovo assemblies.
  --no-yahs             Give this option if you do not want to run yahs.
  --no-smudgeplot       Give this option if you do not want to run smudgeplot.
  --run-tigmint         Give this option if you want to run the scaffolding with 10X reads step.
  --run-kraken2         Give this option if you want to run Kraken2 on the input reads.
  --genomescope-opts genomescope_additional
                        Additional options to run Genomescope2 with. Default -m -1
  --preprocess-lr-step PREPROCESS_LR_STEP
                        Step for preprocessing long-reads. Default 02.1
  --preprocess-10X-step PREPROCESS_10X_STEP
                        Step for preprocessing 10X reads. Default 02.2
  --preprocess-illumina-step PREPROCESS_ILLUMINA_STEP
                        Step for preprocessing illumina reads. Default 02.2
  --preprocess-hic-step PREPROCESS_HIC_STEP
                        Step for preprocessing hic reads. Default 02.3
  --flye-step FLYE_STEP
                        Step for running flye. Default 03.1
  --hifiasm-step HIFIASM_STEP
                        Step for running hifiasm. Default 03.3
  --nextdenovo-step NEXTDENOVO_STEP
                        Step for running nextdenovo. Default 03.2
  --concat-cores concat_cores
                        Number of threads to concatenate reads and to run filtlong. Default 4
  --minimap2-cores minimap2_cores
                        Number of threads to run the alignment with minimap2. Default 32
  --bwa-cores bwa_cores
                        Number of threads to run the alignments with BWA-Mem2. Default 16
  --hypo-cores hypo_cores
                        Number of threads to run the hypo step. Default 24
  --nextpolish-cores nextpolish_cores
                        Number of threads to run the nextpolish step. Default 24
  --pairtools-parse-cores pairtools_parse_cores
                        Number of threads to run the pairtools parse step. Default 32
  --pairtools-sort-cores pairtools_sort_cores
                        Number of threads to run the pairtools sort step. Default 16
  --pairtools-dedup-cores pairtools_dedup_cores
                        Number of threads to run the pairtools dedup step. Default 8
  --pairtools-split-cores pairtools_split_cores
                        Number of threads to run the pairtools split step. Default 16
  --busco-cores busco_cores
                        Number of threads to run BUSCO. Default 32
  --longranger-cores longranger_cores
                        Number of threads to run longranger. Default 16
  --longranger-path longranger_path
                        Path to longranger executable. Default /scratch/project/devel/aateam/src/10X/longranger-2.2.2

Inputs:
  --scripts-dir SCRIPTS_DIR
                        Directory with the different scripts for the pipeline. Default /software/assembly/pipelines/Assembly_pipeline/CLAWS/bin/../scripts/
  --ont-dir ONT_DIR     Directory where the ONT fastqs are stored. Default None
  --hifi-dir HIFI_DIR   Directory where the hifi reads are stored. Default None
  --illumina-dir ILLUMINA_DIR
                        Directory where the raw illumina fastqs are stored. Default None
  --hic-dir HIC_DIR     Directory where the HiC fastqs are stored. Default None
  --raw-10X RAW_10X [RAW_10X ...]
                        Dictionary with 10X raw read directories, it has to be the mkfastq dir. You must specify as well the sampleIDs from this run. Example: '{"mkfastq-
                        dir":"sample1,sample2,sample3"}'...
  --ont-reads ONT_READS
                        File with all the ONT reads. It can either be in fastq, fasta or bam format Default None
  --hifi-reads HIFI_READS
                        File with all the HiFi reads It can be either in fastq, fasta or bam format. Default None
  --pe1 PE1             File with the illumina paired-end fastqs, already trimmed, pair 1.
  --pe2 PE2             File with the illumina paired-end fastqs, already trimmed, pair 2.
  --10X R10X            File with barcoded 10X reads in fastq.gz format, concatenated.
  --processed-illumina PROCESSED_ILLUMINA
                        Directory to Processed illumina reads. Already there or to be produced by the pipeline.
  --processed-10X PROCESSED_10X
                        Directory to Processed 10X reads. Already there or to be produced by the pipeline.
  --ont-filt ONT_FILTERED
                        File with the ONT reads after running filtlong on them. Default None
  --assembly-in ASSEMBLY_IN [ASSEMBLY_IN ...]
                        Dictionary with assemblies that need to be polished but not assembled and directory where they should be polished. Example: '{"assembly1":"polishing_dir1"}'
                        '{"assembly2"="polishing_dir2"}' ...
  --postpolish-assemblies POSTPOLISH_ASSEMBLIES [POSTPOLISH_ASSEMBLIES ...]
                        Dictionary with assemblies for which postpolishing steps need to be run but that are not assembled and base step for the directory where the first postpolishing step should be run.
                        Example: '{"assembly1":"s04.1_p03.1"}' '{"assembly2":"s04.2_p03.2"}' ...
  --curated-assemblies CURATED_ASSEMBLIES [CURATED_ASSEMBLIES ...]
                        Dictionary with assemblies that have already been curated and directory where read alignment should be run. Evaluations and read alignment will be performed. Example:
                        '{"assembly1":"s04.1_p03.1"}' '{"assembly2":"s04.2_p03.2"}' ...

Outputs:
  --pipeline-workdir PIPELINE_WORKDIR
                        Base directory for the pipeline run. Default /scratch_isilon/groups/assembly/jgomez/Annotation_AAT_pipeline/
  --preprocess-lr PREPROCESS_LR
                        Directory to process the long-reads. Default s02.1_p01.1_Preprocess_LR
  --concat-hic-dir CONCAT_HIC_DIR
                        Directory to concatenate the HiC reads. Default s02.3_p01.1_Concat_HiC
  --flye-dir FLYE_DIR   Directory to run flye. Default s03.1_p02.1_flye/
  --nextdenovo-dir NEXTDENOVO_DIR
                        Directory to run nextdenovo. Default s03.2_p02.1_nextdenovo/
  --hifiasm-dir HIFIASM_DIR
                        Directory to run hifiasm. Default s03.3_p02.1_hifiasm/
  --flye-polishing-dir POLISH_FLYE_DIR
                        Directory to polish the flye assembly. Default s04.1_p03.1_polishing/
  --nextdenovo-polishing-dir POLISH_NEXTDENOVO_DIR
                        Directory to run nextdenovo. Default s04.2_p03.2_polishing/
  --eval-dir eval_dir   Base directory for the evaluations. Default evaluations/
  --stats-out stats_out
                        Path to the file with the final statistics.
  --hic-qc-dir hic_qc_dir
                        Directory to run the hic_qc. Default hic_qc/

Filtlong:
  --filtlong-minlen filtlong_minlen
                        Minimum read length to use with Filtlong. Default 1000
  --filtlong-min-mean-q filtlong_min_mean_q
                        Minimum mean quality to use with Filtlong. Default 80
  --filtlong-opts filtlong_opts
                        Extra options to run Filtlong (eg. -t 4000000000)

Trim_Galore:
  --trim-galore-opts trim_galore_opts
                        Optional parameters for the rule trim_galore. Default --max_n 0 --gzip -q 20 --paired --retain_unpaired
  --trim-Illumina-cores Trim_Illumina_cores
                        Number of threads to run the Illumina trimming step. Default 8

Kraken2:
  --kraken2-db kraken2_db
                        Database to be used for running Kraken2. Default None
  --kraken2-kmer kraken2_kmers
                        Database to be used for running Kraken2. Default None
  --kraken2-opts additional_kraken2_opts
                        Optional parameters for the rule Kraken2. Default
  --kraken2-cores kraken2_threads
                        Number of threads to run the Kraken2 step. Default 16

Flye:
  --flye-cores flye_cores
                        Number of threads to run FLYE. Default 128
  --flye-polishing-iterations flye_pol_it
                        Number of polishing iterations to use with FLYE. Default 2
  --other-flye-opts other_flye_opts
                        Additional options to run Flye. Default --scaffold

Hifiasm:
  --hifiasm-cores hifiasm
                        Number of threads to run Hifiasm. Default 50
  --other-hifiasm-opts other_hifiasm_opts
                        Additional options to run Hifiasm. Default --ont
  --purge-hifiasm       Give this option if you want to run purgedups externally on the hifiasm output.
  --phase-hifiasm       Give this option if you want to phase with hic reads the hifiasm assembly

Nextdenovo:
  --nextdenovo-cores nextdenovo_cores
                        Number of threads to run nextdenovo. Default 2
  --nextdenovo-jobtype nextdenovo_type
                        Job_type for nextdenovo. Default slurm
  --nextdenovo-task nextdenovo_task
                        Task need to run. Default all
  --nextdenovo-rewrite nextdenovo_rewrite
                        Overwrite existing directory. Default yes
  --nextdenovo-parallel_jobs nextdenovo_parallel_jobs
                        Number of tasks used to run in parallel. Default 50
  --nextdenovo-minreadlen nextdenovo_minreadlen
                        Filter reads with length < minreadlen. Default 1k
  --nextdenovo-seeddepth nextdenovo_seeddepth
                        Expected seed depth, used to calculate seed_cutoff, co-use with genome_size, you can try to set it 30-45 to get a better assembly result. Default 45
  --nextdenovo-seedcutoff nextdenovo_seedcutoff
                        Minimum seed length, <=0 means calculate it automatically using bin/seq_stat. Default 0
  --nextdenovo-blocksize nextdenovo_blocksize
                        Block size for parallel running, split non-seed reads into small files, the maximum size of each file is blocksize. Default 1g
  --nextdenovo-pa-correction  nextdenovo_pa_correction
                        number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage, overwrite parallel_jobs only for this step. Default 100
  --nextdenovo-minimap_raw nextdenovo_minimap_raw
                        minimap2 options, used to find overlaps between raw reads, see minimap2-nd for details. Default -t 30
  --nextdenovo-minimap_cns nextdenovo_minimap_cns
                        minimap2 options, used to find overlaps between corrected reads. Default -t 30
  --nextdenovo-minimap_map nextdenovo_minimap_map
                        minimap2 options, used to map reads back to the assembly. Default -t 30 --no-kalloc
  --nextdenovo-sort nextdenovo_sort
                        sort options, see ovl_sort for details. Default -m 400g -t 20
  --nextdenovo-correction_opts nextdenovo_correction_opts
                        Correction options. Default -p 30 -dbuf
  --nextdenovo-nextgraph_opt nextdenovo_nextgraph_opt
                        nextgraph options, see nextgraph for details. Default -a 1

Hypo:
  --sr-cov ill_cov      Approximate short read coverage for hypo Default 0
  --hypo-proc hypo_processes
                        Number of contigs to be processed in parallel by HyPo. Default 6
  --hypo-no-lr          Set this to false if you don¡t want to run hypo with long reads. Default True
  --hypo-opts hypo_opts
                        Additional options to run Hypo. Default None

Purge_dups:
  --purgedups-cores purgedups_cores
                        Number of threads to run purgedups. Default 8
  --purgedups-calcuts-opts calcuts_opts
                        Adjusted values to run calcuts for purgedups. Default None

Scaffold_with_10X:
  --tigmint-cores tigmint_cores
                        Number of threads to run the 10X scaffolding step. Default 12
  --tigmint-opts tigmint_opts
                        Adjusted values to run the scaffolding with 10X reads. Default None

HiC:
  --hic-qc              Give this option if only QC of the HiC data needs to be done.
  --subsample-hic       Give this option if you want to subsample the hic data for qc.
  --add-preseq-opts ADD_PRESEQ_OPTS
                        Additional options to give for preseq etrapolation. E.g. -D. Default:
  --no-pretext          Give this option if you do not want to generate the pretext file
  --sort-pretext SORT_PRETEXT
                        Specify how to sort the pretext (eg. --nosort or --sortby something. Default: nosort
  --assembly-qc assembly_qc
                        Path to the assembly to be used perfom the QC of the HiC reads.
  --yahs-cores yahs_cores
                        Number of threads to run YAHS. Default 48
  --yahs-mq yahs_mq     Mapping quality to use when running yahs.Default 10
  --yahs-contig-ec      Give this option if you want to allow yahs perform contig breaks.
  --yahs-opts yahs_opts
                        Additional options to give to YAHS.Default
  --hic-map-opts hic_map_opts
                        Options to use with bwa mem when aligning the HiC reads. Deafault -5SP -T0
  --mq mq [mq ...]      Mapping qualities to use for processing the hic mappings. Default [0, 10]
  --hic-qc-assemblylen hic_qc_assemblylen
                        Lentgh of the assembly to be used for HiC QC
  --blast-cores blast_cores
                        Number of threads to run blast with the HiC unmapped reads.Default 8
  --hic-blastdb blastdb
                        BLAST Database to use to classify the hic unmapped reads. Default /scratch_isilon/groups/assembly/data/blastdbs
  --hic-readsblast hic_readsblast
                        Number of unmapped hic reads to classify with blast. Default 100

Finalize:
  --no-final-evals      If specified, do not run evaluations on final assemblies. Default True
  --busco-lin busco_lineage
                        Path to the busco lineage to be used.
  --merqury-db merqury_db
                        Meryl database. Default None
  --merqury-plot-opts merqury_plot_opts
                        Meryl database. Default None
  --meryl-k meryl_k     Merqury plot additional options, for example " -m 200 -n 6000|". Default None
  --meryl-threads meryl_threads
                        Number of threads to run meryl and merqury. Default 4
  --meryl-reads meryl_reads [meryl_reads ...]
                        Type of reads to be used to build the meryldb. Default ont illumina

Wildcards:
  --ont-list ONT_wildcards
                        List with basename of the ONT fastqs that will be used. Default None
  --hifi-list hifi_wildcards
                        List with basename of the ONT fastqs that will be used. Default None
  --illumina-list illumina_wildcards
                        List with basename of the illumina fastqs. Default None
  --r10X-list r10X_wildcards
                        List with basename of the raw 10X fastqs. Default None
  --hic-list hic_wildcards
                        List with basename of the raw hic fastqs. Default None
```

# Changes made to v3.2:

1. Pairtools has been replaced by Chromap, improving a lot the efficiency of the Hi-C mapping and processing. 

# Changes made to v3.0:

1. Assembly nomenclature

	Nomenclature: hypo is now hyp

	Add mq to yahs name

	Basename has been added as prefix for the assemblies

2. Behaviour changes:
   
	New options to provide Hifi reads have been implemented: 

  		``--hifi-reads:`` file with all the HiFi reads It can be either in fastq, fasta or bam format.

  		``--hifi-dir:`` directory where the hifi reads are stored. In this case, the files need to be in fastq format.  

	ONT reads can now also be in .bam format 

	Illumina reads suffix can now be "_1.fastq.gz" and ".R1.fastq.gz" or ".1.fastq.gz" (only the latter used to be possible in previous versions of CLAWS) 

	Meryl_dbs can now be built on "hifi" data if specified with "--meryl-reads" and/or with " –lr-type pacbio-hifi" option. 

	Fasta-stats has been replaced by gfastats 

	Polishing and filtlong have been turned off if long-read type is hifi reads.  

	Hifiasm expects only 2 haps output if no phasing and as many haps as given ploidy if phasing. 

	Busco version has been updated to v6.0.0 and odb12 databases 

	Change mq defaults to 0, 10 for pretext and 10 for yahs 

	New "add_preseq_opts" option 

	New subsample hic rule, possible to turn it on for hic_qc with "--subsample-hic" option 

	Generation of plot for hic qc has been added 

	``--split-prefix`` option has been added to minimap2 for genomes larger than 4G 

	.csi indexes are now made instead of .bai 

	Option "-no-contig-ec" is now default for yahs, breaking can be activated with the new "--yahs-contig-ec" option 

	Bug that was not properly running merqury in haps mode for postassembly steps has been fixed.  

	PretextGraph has now been updated to version 0.0.9 

	New files have been added to the cleaning step: 

		- Pairtools_out directory 

		- Hic alignments 

		- Hifi temporary fastq when bam is given 

		- Long-read alignments against scaffolded assemblies are now kept 

	HighrRes pretext files are now generated 

	Tidk has been updated from v0.2.0 to 0.2.65 

	Telomere string can now be given as option to the pipeline and it will be used both for hifiasm and for running tidk find on every assembly, in the evaluations directory.  

# Changes made to v2.3: 

1. Assembly nomenclature 

	Previous flye.assembly --> fl.asm 

	Previous nextdenovo.assembly --> nd.asm 

	Previous hypo --> hp 

	Previous nextpolish_ont --> npo 

	Previous nextpolish_ill --> npi 

	Previous purged --> pgd 

	Previous yahs_scaffolds_final --> yhs_scffs 

2. Versions

	Trim_galore version has been updated from 0.6.7 to 0.6.10 

	Flye version has been updated from 2.9.1 to 2.9.5 

	Busco version has been updated from 5.4.0 to 5.5.0 

3. Behaviour changes 

    Genomescope default options have changed, from " -m 10000 " to " -m –1 " 

	Added --no-smudgeplot option to skip running smudgeplot if desired 

	Hifiasm has been included 

	If lr_type = pacbio-hifi, hifiasm will run without ont option.  

	Purgedups can be optionally run on Hifiasm assemblies (option –purge-hifiasm). By default it will only scaffold the hifiasm assemblies 

	Merqury will run hap1 and hap2 together 

	Busco output filenames contain now the name of the db used.  

	Removed the generation of tmp file in rule align_lr 

	Align_hic rule is now using bwa-mem2 

	Added optimized thread options to pairtools rules (thanks to Francisco)

	Sort option has been added to generate pretext rule (--sort-pretext, default is "nosort") 


# Changes made to v2.2: 

1. General: 

	Now default read_type is nano-hq 

2. Rule trim_galore: 

	"--max_n 0" has been added to the default behaviour of "--trim-galore-opts" 

3. Meryl: 

	New option "--meryl-reads" has been added to the config. Default is "Illumina ont" to build the meryl database using both type of reads, it can be changed to one or the other 

4. Merqury: 

	Option "--merqury-plot-opts" has been added to config file. It can be used to modify the x and y axis maximum values (eg. --merqury-plot-opts " -m 200 -n 6000") 

5. Genomescope: 

	"-m 10000" is now part of the default behavior of "--genomescope-opts" 

6. Hic_statistics: 

	This is now running for each assembly and mq for which a pretext file is generated 

7. Assembly inputs for different steps: 

	a. "--assembly-in" to start after assembly step (eg. Evaluation, polishing, purging and scaffolding) 

	b. "--postpolish-assemblies" to start after polishing step (eg. Evaluation, purging and scaffolding) 

	c. "--curated-assemblies" to start after scaffolding step (eg. Evaluation and pretext generation) 
