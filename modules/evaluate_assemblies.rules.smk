from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')
scripts_dir = ""
keepfiles = False

rule finalize:
  input:
    assembly = "assembly.fasta",
    buscos = "busco_short_summary.txt",
    stats = "assembly.stats.txt",
    merqs= "assembly.qv",
    pretext = "assembly.pretext"
  output:
    output = "stats.txt",
  params:
    scripts_dir = "../scripts/",
    rmcmd = "",
  conda:
    '../envs/ass_base.yaml'
  threads: 1
  shell:
    "{params.scripts_dir}get_final_tbl_mult.py -s {input.stats} -b {input.buscos} -m {input.merqs} > {output.output};"
    "{params.rmcmd}"
    "echo 'Pipeline completed';"

# rule get_report:
#   input:
#     stats = "stats,txt",
#     config = "assembly_pipeline.config"
#   output:
#     report = "{base}.report.zip"
#   conda:
#     '../envs/ass_base.yaml'
#   threads: 1
#   shell:
#     "snakemake --snakefile /software/assembly/pipelines/Assembly_pipeline/v2.0/assembly_pipeline/bin/assembly_pipeline.smk --configfile {input.config} --report {output.report};"

rule get_stats:
  input:
    assembly = "assembly.fasta"
  output:
    nseries =  "assembly.nseries.txt",
    stats = "assembly.stats.txt",
    gaps = "assembly.gaps.txt"
  params:
    outbase = "assembly",
    scripts_dir = scripts_dir
  log:
    "logs/" + str(date) + ".j%j.assembly_stats.out",
    "logs/" + str(date) + ".j%j.assembly_stats.err"
  benchmark:
    "logs/" + str(date) + ".assembly_stats.benchmark.txt"
  conda:
    "../envs/ass_base.yaml"
  threads: 1
  shell:
    "{params.scripts_dir}fastalength {input.assembly} | {params.scripts_dir}Nseries.pl > {output.nseries};"
    "{params.scripts_dir}fasta-stats.py -f {input.assembly} -s {output.stats} -r {output.gaps};"
    "sleep 4m"

rule get_stats_gfa:
  input:
    #assembly_fa = "assembly.fa",
    assembly_gfa = "assembly.gfa"
  output:
   # nseries =  "assembly.nseries.txt",
    stats = "assembly.stats.txt",
  params:
    outbase = "assembly",
    scripts_dir = scripts_dir,
    params = " --discover-paths "
  log:
    "logs/" + str(date) + ".j%j.assembly_stats.out",
    "logs/" + str(date) + ".j%j.assembly_stats.err"
  benchmark:
    "logs/" + str(date) + ".assembly_stats.benchmark.txt"
  conda:
    '../envs/gfastastats1.3.10.yaml'
  threads: 1
  shell:
    #"{params.scripts_dir}fastalength {input.assembly_fa} | {params.scripts_dir}Nseries.pl > {output.nseries};"
    "gfastats {params.params} --nstar-report {input.assembly_gfa} > {output.stats};"
    "sleep 4m"

rule run_merqury:
  input:
    meryl_db = os.getcwd() + "/database.meryl",
    assembly = [os.getcwd() + "/assembly.fasta"]
  output:
    completeness = os.getcwd() + "/merqury_run/assembly.completeness.stats",
    qv= os.getcwd() + "/merqury_run/assembly.qv",
    hist = os.getcwd() + "/merqury_run/assembly.assembly.spectra-cn.hist",
    false_dups = os.getcwd() + "/merqury_run/assembly.false_duplications.txt",
    plots = [os.getcwd() + "/merqury_run/assembly.spectra-cn.ln.png", os.getcwd() + "/merqury_run/assembly.spectra-cn.fl.png", os.getcwd() + "/merqury_run/assembly.spectra-cn.st.png"]
  params:
    out_pref = "assembly",
    directory = os.getcwd() + "/merqury_run",
    additional_plot_opts = " "
  log:
    "logs/" + date + ".merqury.out",
    "logs/" + date + ".merqury.err",
  conda:
    "../envs/merqury1.3.yaml"
  threads: 4
  shell:
    "mkdir -p {params.directory}; cd {params.directory};"
    "merqury.sh {input.meryl_db} {input.assembly} {params.out_pref};"
    "if [ ! -f {output.hist} ]; then ln -s *{params.out_pref}.spectra-cn.hist {output.hist}; ln -s *{params.out_pref}.spectra-cn.ln.png {output.plots[0]}; ln -s *{params.out_pref}.spectra-cn.fl.png {output.plots[1]};  ln -s *{params.out_pref}.spectra-cn.st.png {output.plots[2]};  fi;"
    "$CONDA_PREFIX/share/merqury/eval/false_duplications.sh {output.hist} > {output.false_dups};"
    "for i in {input.assembly}; do b=`basename $i \".fa\"`;$CONDA_PREFIX/share/merqury/eval/false_duplications.sh {params.directory}/{params.out_pref}.$b.spectra-cn.hist > {params.directory}/$b.false_duplications.txt; done"

rule run_busco:
  input:
    assembly = os.getcwd() + "/assembly.fasta",
    lineage = "/scratch/groups/assembly/shared/databases/busco_v5/busco_v5_lineages/lineages/vertebrata_odb10"
  output:
    summary = os.getcwd() + "/busco/short_summary.txt",
    full = os.getcwd() + "/busco/full_table.tsv",
  params:
    out_path = os.getcwd() + "/busco",
    odb = "vertebrata_odb10",
    buscobase = "assembly",
    rmcmd = "echo 'Removing BUSCO run dir: busco'; \
            rm -r busco;" if keepfiles == True else "" 
  conda:
    '../envs/busco6.0.0.yaml'
  threads: 8
  shell:
    "cd {params.out_path};" 
    "busco -i {input.assembly} -c {threads} -m genome --offline -o {params.buscobase} --out_path {params.out_path} -l {input.lineage};"
    "mv {params.out_path}{params.buscobase}/run_{params.odb}/short_summary.txt {output.summary};"
    "mv {params.out_path}{params.buscobase}/run_{params.odb}/full_table.tsv {output.full};"
    "{params.rmcmd}"
  
rule align_lr:
  input:
    genome = os.getcwd() + "/assembly.fasta",
    reads = os.getcwd() + "/reads.ont.fastq.gz"
  output:
    mapping = "minimap2.bam"
  params:
    align_opts = "ax",
    type = "map-ont",
    split = "",
    compress_cmd = lambda wildcards: "samtools view -Sb - | " \
                   "samtools sort -@ " + str(config["Parameters"]["minimap2_cores"]) +" -o " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext +";" +\
                   "samtools index -c " + wildcards.directory + "/mappings/" + wildcards.name + "_" + wildcards.ext 
  conda:
    "../envs/minimap2.24.yaml"
  threads: 4
  shell:
    "minimap2 -{params.align_opts} {params.type} {params.split} -t {threads} {input.genome} {input.reads} | {params.compress_cmd};"

rule align_illumina:
  input:
    genome = os.getcwd() + "/assembly.fasta",
    reads = [os.getcwd() + "/reads.illumina.1.fastq.gz", os.getcwd() + "/reads.illumina.2.fastq.gz"]
  output:
    mapping = "illumina.bam",
    stats = "illumina.stats.txt"
  params:
    options = ""
  conda:
    "../envs/bwa-mem2.2.1.yaml"
  threads: 4
  shell:
    "bwa-mem2 index {input.genome};"
    "bwa-mem2 mem -Y {params.options} -t {threads} {input.genome} {input.reads} | samtools view -Sb - | samtools sort -@ {threads} -o {output.mapping} -;"
    "samtools index -c {output.mapping};"
    "samtools stats -@ {threads} {output.mapping} > {output.stats};"

rule align_hic:
  input:
    ass = "assembly.fa",
    bwt = "assembly.fa.bwt.2bit.64",
    read1 = "HiC.R1.fastq.gz",
    read2 = "HiC.R2.fastq.gz"
  output:
    mapped = "full_hic.bam",
    unmapped = "unmapped_hic.bam" 
  params:
    outd = "mappings",
    name = "assembly",
    options = "-5SP -T0"
  conda:
    "../envs/bwa-mem2.2.1.yaml"
  shell:
    """
    mkdir -p {params.outd};
    #samtools does not require more than 1 thread
    bwa-mem2 mem {params.options} -t {threads} {input.ass} {input.read1} {input.read2} | samtools view -bS -@ 1 - | \
    tee {output.mapped} | samtools view -@ 1 -b -f 4 \
    > {output.unmapped}
    """

rule pairtools_processing:
  input: 
    mapped = "full_hic.bam",
    alength = "assembly.genome"
  output:
    stats = "stats.mq40.txt",
    mappedptsort = "mapped.PT.mq40.name_sorted.bam",
    mappedpt = "mapped.PT.mq40.bam"
  params:
    scripts_dir = "../scripts/",
    mq = 40,
    outd = "pairtools_out",
    tmpd = "tmp",
    name = 'assembly',
    rmcmd = "rm -r {params.tmpd};"
  conda:
    "../envs/dovetail_tools_v2_R.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH";
    mkdir -p {params.tmpd};
    mkdir -p {params.outd};

    samtools view -h -@ {threads} {input.mapped} | \
    pairtools parse --min-mapq {params.mq} --walks-policy 5unique --max-inter-align-gap 30 --nproc-in {threads} --nproc-out {threads} \
    --chroms-path {input.alength} -o {params.tmpd}/full_hic.{params.mq}.{params.name}.pairsam.gz;

    pairtools sort --tmpdir={params.tmpd} --nproc {threads} \
    -o {params.tmpd}/full_hic.{params.mq}.{params.name}.sorted.pairsam.gz {params.tmpd}/full_hic.{params.mq}.{params.name}.pairsam.gz

    pairtools dedup --nproc-in {threads} --nproc-out {threads} --mark-dups \
    --output-stats {output.stats} \
    --output {params.outd}/full_hic.{params.mq}.{params.name}.pairsam.dedupped.gz {params.tmpd}/full_hic.{params.mq}.{params.name}.sorted.pairsam.gz

    pairtools split --nproc-in {threads} --nproc-out {threads} {params.outd}/full_hic.{params.mq}.{params.name}.pairsam.dedupped.gz \
    --output-pairs {params.outd}/mapped.mq{params.mq}.{params.name}.pairs \
    --output-sam {params.tmpd}/mapped.PT.mq{params.mq}.{params.name}.sam 
    
    samtools view -bS -@ {threads} {params.tmpd}/mapped.PT.mq{params.mq}.{params.name}.sam | samtools sort -@ {threads} \
    -o {output.mappedpt}; 

    samtools index -c -@ {threads} {output.mappedpt};

    samtools sort -n -@ {threads} {output.mappedpt} -o {output.mappedptsort}

    echo pairtools done, now running {params.rmcmd}
    {params.rmcmd}
    """

rule pairtools_processing_parse:
  input: 
    mapped = "full_hic.bam",
    alength = "assembly.genome"
  output:
    psamo = "mq.pairsam.gz",
  params:
    scripts_dir = "../scripts/", 
    mq = 40,
    outd = "pairtools_out",
    tmpd = "tmp",
    name = 'assembly',
  conda:
    "../envs/dovetail_tools_v2_R.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH";
    mkdir -p {params.tmpd};
    mkdir -p {params.outd};
    threads_half=$(({threads}/2));
    samtools view -h -@ 1 {input.mapped} | \
    pairtools parse --min-mapq {params.mq} --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $threads_half --nproc-out $threads_half \
    --chroms-path {input.alength} -o {output.psamo};

    """

rule pairtools_processing_sort:
  input: 
    psamo = "mq.pairsam.gz"
  output:
    spsamo = "mq.sorted.pairsam.gz"
  params:
    scripts_dir = "../scripts/", 
    mq = 40,
    outd = "pairtools_out",
    tmpd = "tmp",  
    name = 'assembly',
  conda:
    "../envs/dovetail_tools_v2_R.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH";

    pairtools sort --tmpdir={params.tmpd} --nproc {threads} \
    -o {output.spsamo} {input.psamo}
    """

rule pairtools_processing_dedup:
  input: 
    spsamo = "mq.sorted.pairsam.gz"
  output:
    spsamdedup = "mq.pairsam.dedupped.gz",
    stats = "stats.mq40.txt",
  params:
    scripts_dir = "../scripts/", 
    mq = 40,
    outd = "pairtools_out",
    tmpd = "tmp", 
    name = 'assembly',
  conda:
    "../envs/dovetail_tools_v2_R.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH";
    threads_half=$(({threads}/2));
    pairtools dedup --nproc-in $threads_half --nproc-out $threads_half --mark-dups \
    --output-stats {output.stats} \
    --output {output.spsamdedup} {input.spsamo}
 
    """

rule pairtools_processing_split:
  input: 
    spsamdedup = "mq.pairsam.dedupped.gz"
  output:
    mappedptsort = "mapped.PT.mq40.name_sorted.bam",
    mappedpt = "mapped.PT.mq40.bam"
  params:
    scripts_dir = "../scripts/", 
    mq = 40,
    outd = "pairtools_out",
    tmpd = "tmp",   
    name = 'assembly',
  conda:
    "../envs/dovetail_tools_v2_R.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH";
    threads_half=$(({threads}/2));
    cd {params.outd};
    pairtools split --nproc-in $threads_half --nproc-out $threads_half {input.spsamdedup} \
    --output-pairs {params.outd}/mapped.mq{params.mq}.{params.name}.pairs \
    --output-sam {params.tmpd}/mapped.PT.mq{params.mq}.{params.name}.sam 
    
    cis_pairs_dist.v03.R {params.outd}/mapped.mq{params.mq}.{params.name}.pairs ../{params.name}.mq{params.mq}
   
    samtools view -bS {params.tmpd}/mapped.PT.mq{params.mq}.{params.name}.sam > {params.tmpd}/tmp.mapped.PT.mq{params.mq}.{params.name}.bam 
   
    samtools sort -@ {threads} -m 20G {params.tmpd}/tmp.mapped.PT.mq{params.mq}.{params.name}.bam -o {output.mappedpt}; 

    samtools index -c -@ {threads} {output.mappedpt};

    samtools sort -n -@ {threads} -m 20G {output.mappedpt} -o {output.mappedptsort}
    
    echo pairtools done, now running {params.rmcmd}
    {params.rmcmd}

    """

rule qc_statistics:
  input:
    statmq = "stats.mq40.txt",
    mapbam = "mapped.PT.mq40.bam",
  output:
    libstats = "HiC_QC_LibraryStats_mq40.txt",
  params:
    scripts_dir = "../scripts/",
    assemblylength = "",
    deepseq = False,
    outd = "pairtools_out",
    pslibstats = "HiC_QC_LibraryStats_extrapolated_mq40.txt",
    add_preseq_opts = ""
  conda:
    "../envs/dovetail_tools_v2_R.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH"; 
    # QC Library statistics
    get_qc.py -p {input.statmq} > {output.libstats};
    sleep 5m
  # preseq extrapolation and QC Library statistics
    if [ {params.deepseq} == "False" ]; then  
    preseq lc_extrap -B -P -e 2.1e9 -s 1e8 {params.add_preseq_opts}  -seg_len {params.assemblylength} \
    -o {params.outd}/mapped.PT.mq{wildcards.mq}.preseq {input.mapbam};
    get_qc_preseq.py -p {input.statmq} -d {params.outd}/mapped.PT.mq{wildcards.mq}.preseq \
    > {params.pslibstats};
    fi
    """

rule read_screening:
  input:
    umapped = "unmapped_hic.bam",
  output:
    ufasta = "unmapped_hic.fasta",
    subufasta = "unmapped_hic.100_reads.fasta",
    blastoutbscore = "100_unmapped_reads_vs_nt_25cul1_1e25.megablast.sorted_by_bitscore.out",
    blastoutbscorethits = "100_unmapped_reads_vs_nt_25cul1_1e25.megablast.sorted_by_bitscore.tophits",
    blastoutorganisms = "100_unmapped_reads_vs_nt_25cul1_1e25.megablast.organisms.txt"
  params:
    scripts_dir = "../scripts/",
    outd = "out/blast/",
    readforblast = "",
    blastdb = "/scratch_isilon/groups/assembly/data/blastdbs",
  conda:
    "../envs/blast_samtools.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH;" 
    export BLASTDB={params.blastdb};
    mkdir -p {params.outd};

    samtools view {input.umapped} | gawk '{{print $1"\t"$10}}' | TblToFasta > {output.ufasta};
    
    FastaToTbl {output.ufasta} | gawk 'NR <= {params.readforblast} {{ print }}' - | TblToFasta > {output.subufasta};
    
    blastn -task megablast -query {output.subufasta} -db nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
    -max_target_seqs 25 -culling_limit 2 -num_threads {threads} -evalue 1e-25 | sort -k 3rn > {output.blastoutbscore};

    cut -f 3,5-6,11,14,17 {output.blastoutbscore} | head -n 5  > {output.blastoutbscorethits}

    cut -f 17-20 {output.blastoutbscore} | perl -ne '@data=split/s+/,$_; print "$data[0] $data[1]\\n";' | sort -k1,1  | uniq > {output.blastoutorganisms};

    """

rule get_extension_gaps:
  input:
    sla = "assembly.fa" ,
  output:
    gaps_bed = "gaps.bed",
    gaps = "gaps.bg"
  params:
    scripts_dir = "../scripts/"
  conda:
    '../envs/gfastastats1.3.10.yaml'
  shell:
    """
    gfastats -b gaps {input.sla} > {output.gaps_bed};
    cat {output.gaps_bed}| {params.scripts_dir}/gap_bed2bedgraph.sh > {output.gaps} ;
    sleep 4m
    """

rule get_extension_cov:
  input:
    sbam = "sorted.bam"
  output:
    ontcov = "ONTcoverage.bg"  
  conda:
    "../envs/bedtools2.30.0.yaml"
  shell:
    """
    bedtools genomecov -bga -ibam {input.sbam}  > {output.ontcov}
    """

rule tidk_search:
  input:
    sla = "assembly.fa" 
  output:
    tel = "TTAGGG_telo.bg"
  params:
    outd = "out",
    teloseq = "TTAGGG",
    outname = "TTAGGG"
  conda:
    "../envs/tidk-0.2.65.yaml"
  shell:
    """
    cd {params.outd}

    tidk search -s {params.teloseq} -o {params.outname}  -e bedgraph --dir $PWD {input.sla}
    
    ln -s {params.outname}_telomeric_repeat_windows.bedgraph {output.tel}
    
    sleep 5m
    """

rule tidk_explore:
  input:
    sla = "assembly.fa" 
  output:
    tel = "telomeres.bg"
  params:
  conda:
    "../envs/tidk-0.2.65.yaml"
  shell:
    """
    cd {params.outd}

    tidk explore -x 12 -m 5 {input.sla} > out.tidk_explore.txt
    """

rule get_extension_telomeres:
  input:
    sla = "assembly.fa" 
  output:
    tel = "telomeres.bg"
  params:
    outd = "out"
  conda:
    "../envs/tidk-0.2.0.yaml"
  shell:
    """
    #explore the telomeric repeats on the genome for motifs 5-12 long
      
    #fixed this line had tidk explore twice
    
    cd {params.outd}

    tidk explore -x 12 -m 5 --output out --dir telomeres --extension bedgraph --fasta {input.sla} \
    &> telomeric_explore.log

    #format to a proper bedgraph for pretext
    egrep -v ^id telomeres/out_telomeric_locations.bedgraph | cut -f 1-4 > {output.tel}
    sleep 3m
    """

rule generate_diploid:
  input: 
    hap1 = "assembly.hap1.fa",
    hap2 = "assembly.hap2.fa"
  output:
    diploid = "assembly.diploid.fa"
  shell:
    """
    cat {input.hap1} | sed 's/>/>H1_/' > {output.diploid}
    cat {input.hap2} | sed 's/>/>H2_/' >> {output.diploid}
    sleep 5m
    """
    