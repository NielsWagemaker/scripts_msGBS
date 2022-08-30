# scripts_msGBS
All msGBS scripts as described in msGBS articles -
A:  msGBS: A new high-throughput approach to quantify the relative species abundance in root samples of multispecies plant communities, 2020, Molecular ecology Resources, https://doi.org/10.1111/1755-0998.13278
B:  msGBS: A field study, article in progress, a collaboration with Dina in 't Zandt.

The first article (A) is a method paper describing and evaluating the use of GBS with multi species plant root mixtures with the aim of quantifiying the relative contribution of the different species within the mixtures in a calibrated and non-calibrated way. The first enabling the comparison between species within samples and the second within species between samples. The second article (B, in prep) aims to apply msGBS in "natural" systems were >40 species can be present simultaniously. Here we did not apply a calibration as it is not always feasible to prepare these for highly biodiverse samples. The msGBS methods is recently also applied to Bee pollen mixtures and Diatom mixtures.

# Requirements:
1.Clone pipeline from git
```
git clone https://github.com/NielsWagemaker/scripts_msGBS.git
```

2.Anaconda:

first install conda following the instructions specified at https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html 
then run the following code.
```
conda env create -f src/env/environment.yaml --name msGBS
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
