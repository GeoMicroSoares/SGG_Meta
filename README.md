# *anvi'o* analysis of SGG metagenomes

##Scripts used to process SGG metagenomes (Chapters 4 and 5 of my PhD thesis).
1. ```qc_trimmed_nonjob.sh``` uses ```fastp``` (https://github.com/OpenGene/fastp) to quality filter and merge Illumina paired-end reads.

2. ```merge_fastqs.sh``` simply concatenates fastq.gz files.

3. ```anvio_import_final.sh``` uses ```anvi'o``` (https://github.com/merenlab/anvio) to
	1) re-format assembled .fasta files,
	2) generate a contigs database,
	3) determine single-copy genes from trained HMMs,
	4) find COGs for and taxonomically classify contigs using ```kaiju``` (https://github.com/bioinformatics-centre/kaiju),
	5) mapping raw reads to contigs using ```bowtie2``` (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
	6) bin contigs using ```CONCOCT``` (https://github.com/BinPro/CONCOCT) and
	7) manually refine bins. Based on ```anvi'o``` tutorials available at http://merenlab.org/2016/06/22/anvio-tutorial-v2/.

4. ```anvio_pangenomics.sh``` uses ```anvi'o``` to create FeOB, SOB pangenomes and analyse them.

5. ```anvio_omnitrophica_pangenomics.sh``` uses ```ncbi-genome-download``` (https://github.com/kblin/ncbi-genome-download) and ```anvi'o``` to create *Omnitrophica*, *Cand. Methanoperedens* pangenomes and compare them to external [GTDB](https://gtdb.ecogenomic.org/) genomes using [ANIb](https://www.nature.com/articles/s41467-018-07641-9).

### Install megatools
Follow instructions on [here](https://github.com/meganz/MEGAcmd#getting-the-source).

### Set up MEGA credentials
```mega-login [EMAIL] [PASSWORD]```
### Check out MEGA cloud directories
```mega-ls SGG_metagenomes/data/anvio_sgg_mag_processing```
### Download files from MEGA cloud
```mega-get SGG_metagenomes/data/anvio_sgg_mag_processing/anvio_import_final.sh```
