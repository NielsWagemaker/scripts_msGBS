# scripts_msGBS
All msGBS scripts as described in msGBS article

# Requirements:

1.Anaconda:

first install conda following the instructions specified at https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html 
then run the following code.
```
conda env create -f scripts/env/environment.yaml --name msGBS
conda activate msGBS
```

2.blastDB:
```
#install in directory of choice using following command:
update_blastdb.pl --decompress nt
update_blastdb.pl taxdb
tar -xzf taxdb.tar.gz
```

# Run the pipeline

1.Fill in config.yaml

specify locations of input reads, and the different parameters you want to use

2.Run the pipeline
```
snakemake -j 12
#-j specifies the amount of cores being used when running the pipeline. 
```
