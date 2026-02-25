
from datetime import datetime
import re
import os

rule assembly_prepare:
  input: 
   "assembly_in.fa"
  output:
    glen = "assembly.genome",
    bwa = "assembly.bwt.2bit.64",
    faidx = "assembly.fai"   
  params:
    scripts_dir = "../scripts/",  
    workdir = "assembly/",
  conda:
    "../envs/bwa-mem2.2.1.yaml"
  shell:
    """
    export PATH="{params.scripts}:$PATH";
    cd {params.workdir}/; 
    fastalength {input} | gawk '{{print $2"\t"$1}}' > {output.glen};
    bwa-mem2 index {input}; 
    samtools faidx {input};
    """

rule run_yahs:
  input: 
    mappedsort = "mapped.CM.mq10.bam",
    sla = "assembly.fa",
    index = "assembly.fa.fai"
  output:
    outyahs = "yahs_scaffolds_final.fa",
    agp = "yahs_scaffolds_final.agp"
  params:
    yahsdir = "s06.1_p05.1_HiC_scaffolding/",
    contig_ec = " --no-contig-ec ",
    yahsopts = "",
    mq = 10,
    name = "assembly"
  conda:
    "../envs/yahs1.2a.2.yaml"
  shell:
    """
    mkdir -p {params.yahsdir} 
    cd {params.yahsdir}
    yahs {input.sla} {params.contig_ec} {params.yahsopts} {input.mappedsort} 
    ln -s {params.yahsdir}/yahs.out_scaffolds_final.fa {output.outyahs}
    ln -s {params.yahsdir}/yahs.out_scaffolds_final.agp {output.agp}
    """

rule generate_pretext:
  input:
    mapbam = "mapped.CM.mq10.bam"  
  output:
    pret = "assembly_mq10.pretext",
    hr_pret = "assembly_mq10.HR.pretext", 
    ptd = "snapshots/PTdonemq40.txt"
  params:
    scripts_dir = "../scripts",
    outd = "s06.1_p05.1_HiC_scaffolding",
    name = "assembly",
    mq = 40,
    sort = " nosort"
  conda:
    "../envs/pretext-suite0.0.2_plus_samtools1.6.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH;"
    mkdir -p {params.outd}/snapshots/three_wave_blue_green_yellow/
    cd {params.outd}/

    samtools view -@ {threads} -h {input.mapbam} | PretextMap -o {output.pret} \
    --sortby {params.sort} --sortorder descend --mapq {params.mq}

    samtools view -@ {threads} -h {input.mapbam} | PretextMap -o {output.hr_pret} \
    --sortby {params.sort} --sortorder descend --mapq {params.mq} --highRes   

    PretextSnapshot -m {output.pret}  --sequences "=full" \
    -o snapshots/three_wave_blue_green_yellow

    PretextSnapshot -m {output.pret} --sequences "=all" \
    -o snapshots/three_wave_blue_green_yellow

    touch {output.ptd}    
    """

rule add_extensions_pretext:
  input:
    tel = "telomeres.bg",
    gaps = "gaps.bg",
    pret = "assembly_mq.pretext",
    hr_pret = "assembly_mq.HR.pretext",
    ontcov = "ONTcoverage.bg",
    stats = "HiC_Final_LibraryStats_mq.txt"
  output:
    pretext = "assembly_mq.extensions.pretext",
    hr_pretext = "assembly_mq.extensions.HR.pretext"
  params:
    scripts_dir = "../scripts",
    outd = "s06.1_p05.1_HiC_scaffolding",
    mq = 40,
  conda:
    "../envs/PretextGraph0.0.9.yaml"
  shell:
    """
    export PATH="{params.scripts_dir}:$PATH;"
    cd {params.outd}
    cp {input.pret} {output.pretext}
    cp {input.hr_pret} {output.hr_pretext}


    if [[ -s "{input.gaps}" ]]; then
    cat {input.gaps} | PretextGraph -i {output.pretext} -n "gap"
    cat {input.gaps} | PretextGraph -i {output.hr_pretext} -n "gap"
    fi

    if [[ -s "{input.tel}" ]]; then 
    cat {input.tel} | PretextGraph -i {output.pretext} -n "telomere"
    cat {input.tel} | PretextGraph -i {output.hr_pretext} -n "telomere"
    fi

    if [[ -s "{input.ontcov}" ]]; then
    cat {input.ontcov} | PretextGraph -i {output.pretext} -n "coverage" 
    cat {input.ontcov} | PretextGraph -i {output.hr_pretext} -n "coverage" 
    fi
    """

rule get_tpf:
  input:
    fasta = "yahs.out_scaffolds_final.fa"
  output:
    tpf = "genome.yahs_scaffolded.fa.tpf"
  params:
    scripts_dir = "../scripts/",
    dir = "s06.1_p05.1_HiC_scaffolding"
  conda:
    "../envs/rapidcuration.yaml"
  shell:
    "cd {params.dir};"
    "{params.scripts_dir}rapid_split.pl -fa $(basename {input.fasta});"
    "sleep 5m;"



