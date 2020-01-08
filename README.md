# scripts_msGBS
All msGBS scripts as described in msGBS article

# Requirements:
1.Clone pipeline from git
```
git clone https://github.com/NielsWagemaker/scripts_msGBS.git
```

2.Anaconda:

first install conda following the instructions specified at https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html 
then run the following code.
```
conda env create -f scripts/env/environment.yaml --name msGBS
conda activate msGBS
```

3.blastDB:
```
#install in directory of choice using following command:
update_blastdb.pl --decompress nt
update_blastdb.pl taxdb
tar -xzf taxdb.tar.gz
```

# Run the pipeline
1.Make barcode file

The barcode file needs at the minimum to follow the format below with each column separated by tabs.
All mono samples need to have mono in the sample name and all nonmono samples can not have mono in the sample name
```
Flowcell    Lane  Barcode_R1	Barcode_R2  Sample      ENZ_R1  ENZ_R2	Wobble_R1	Wobble_R2
H5HWHCCX2   6	  AACT	        AACT        leafmono71  PacI	  NsiI	3	        3
H5HWHCCX2   6	  GTGAGC	CCAG	    root1	PacI   	  NsiI	3	        3

```
2.Fill in config.yaml

specify locations of input reads, and the different parameters you want to use

3.Run the pipeline
```
#Run the following command in the directory that contains Snakefile
snakemake -j 12
#-j specifies the amount of cores being used when running the pipeline. 
```
