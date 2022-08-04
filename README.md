# Allele-specific DNA-seq workflow

Nextflow workflow to process raw single-end DNA-seq data for allele-specific analysis. It takes a compressed fastq file as input and generates normalized signal tracks to identify chromosome aneuploidies or major chromosomal rearrangements.

## Author

Yuvia Alhelí PÉREZ-RICO

Affiliation: European Molecular Biology Laboratory | EMBL

## Dependencies

To use this workflow, you will need to install the following programs, the indicated version is the one that I have used to analyze my data:

- nextflow (20.04.1)
- FastQC (v0.11.8)
- Trim Galore (0.6.3)
- Bowtie2 (2.3.4.1)
- samtools (1.9)
- SNPsplit (0.3.4)
- Picard Tools (2.20.8)
- deepTools (3.5.1)

Note: I suggest to install all programs in a conda environment.

## How to use the workflow

Thank you for your interest in using the workflow!

After downloading the main script and the configuration file, you will need to prepare the following files and change the paths in the configuration file accordingly:

- An N-masked Bowtie2 index using SNPs of 2 strains.
- A compressed file with the SNPs of the 2 strains for SNPsplit.
- Bed file containing blacklist genomic regions.
- A file with the genome sequence in 2bit format.
- Bed file containing repeat annotations that will be excluded to compute the GC-bias.
- A comma-separated file indicating the sample name and path to the fastq file (an example is available in the 'docs' folder).

If you are not using SLURM, then do additional modifications to the configuration file considering the workload manager that you use. Finally, write a simple bash script that will be submitted to the cluster to activate the conda environment and start the main nextflow job:

`source /home/user/miniconda2/bin/activate /home/user/conda-envs/ASE_DNA-seq`
`nextflow run allelic_DNA-seq.nf`

Please, keep in mind that this workflow was written for the processing of mouse single-end data, therefore, you will need to do some changes to the main script before using it with data of other species, for example, changing the strain names.

## Acknowledgements

This repository is part of a project that has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 882771.

## License

Licenced under the GNU general public license.

