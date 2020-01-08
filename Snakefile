configfile: "config.yaml"

import pandas as pd
import os
import random
df = pd.read_csv(os.path.join(config["input_dir"],config["barcodes"]), sep='\t', index_col='Sample')
SAMPLES = df.index
MONOS = df.index[df.index.str.contains("mono")]
NONMONOS=df.index[~df.index.str.contains("mono")]
flowCell = df.Flowcell[0]
lane = df.Lane[0]
projectName=str(random.randint(1,1000000)) #To ensure non overlapping path directories
tmpdir=config["tmpdir"]+"/msGBS"+projectName
param_mind=config["min-depth"]
read1=config["Read1"]
read1_sub=read1.split(".")[0]
read2=config["Read2"]
read2_sub=read2.split(".")[0]

def getParam_mind(param_mind):
    if param_mind=="default" or param_mind == "":
        mind = 2
    else:
        mind = param_mind
    return mind

param_cycle=config["cycles"]
def getParam_cycle(param_cycle):
    if param_cycle=="default" or param_cycle == "":
        cycle = 150
    else:
        cycle = param_cycle
    return cycle

include: "src/rules/ref_creation.rules"
include: "src/rules/demultiplex.rules"
include: "src/rules/mapping_stats.rules"

rule all:
    input:
        unAssembled_1=expand("{path}/output_denovo/all.joined.fastq.gz",  path=config["output_dir"]),
        merged_Assembled=expand("{path}/output_denovo/all.merged.fastq.txt.gz",  path=config["output_dir"]),
        ref=expand("{path}/output_denovo/refBlasted.fa",path=config["output_dir"]),
        log=expand("{path}/mapping/mapping_variantcalling.log",path=config["output_dir"]),
        blast_file=expand("{path}/output_blast/outputblast_kingdoms.tsv",path=config["output_dir"]),
        statscsv=expand("{path}/stats/stats.csv",path=config["output_dir"]),
        bamOut=expand("{path}/mapping/out.bam",path=config["output_dir"]),
        demulti_R1=expand("{path}/output_demultiplex/demulti_R1.fq.gz",  path=config["output_dir"]),
        demulti_R2=expand("{path}/output_demultiplex/demulti_R2.fq.gz",  path=config["output_dir"]),
        R1_sample_demulti=expand("{path}/output_demultiplex/monos/{sample}.demulti_R1.gz",  path=config["output_dir"],sample=MONOS),
        R2_sample_demulti=expand("{path}/output_demultiplex/monos/{sample}.demulti_R2.gz",  path=config["output_dir"],sample=MONOS)

