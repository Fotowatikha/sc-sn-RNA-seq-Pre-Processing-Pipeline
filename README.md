# sc-sn-RNA-seq-pre-processing

Now considering that the emergence of scRNA-seq techniques provides the framework to study gene expression variability in tissue of interest, the pre-processing of high dimensional RNA-seq data also becomes an important topic. To this day, droplet-based technologies remain the technique of choice when it comes to capturing and sequencing a large number of individual cells. This method essentially barcodes single cells and tags each transcript with a unique molecular identifier (UMI) within individual droplets, which can significantly increase the throughput up to 10,000 cells per analysis. However, this technique is not free of noise, and one must carefully pre-process the data before any unbiased downstream analysis can be caried out. 

**Here, I will provide you with all the necessary information on how to use the single-cell/nuclei pre-processing pipeline, either in your own environment or on the HPC cluster.**

## About

By consider the gene count matrices as the starting point after mapping reads to the reference, the major steps in this pre-processing is to: (1) Provide the user with quality control (QC) plots to gain insight on the overall quality of the cells prior to any extensive filtering; (2) Correcting the data from cell-free ambient RNA; (3) Extensive filtering to remove droplets that are unlikely to represent intact individual cells; (4) Removal of droplets that violate the assumption for containing one cell; (5) Provide the user with QC plots to gain insight on the quality of the data post-filtering.

## Download the pipeline files from GitHub

You can either [download the files for the pipeline](https://github.com/Fotowatikha/sc-sn-RNA-seq-pre-processing/tree/main) from my Github page or download it it by using the "Code button" on the Github page and import it using Filezilla or use the line below if you have Git installed:
```sh
git clone https://github.com/Fotowatikha/sc-sn-RNA-seq-pre-processing.git
```

## Setup Conda environment

This pipeline is made to run on a Conda environment containing all the crucial packages for Python, Snakemake and R. Here, i assume that you already have Anaconda or Miniconda installed. In case you don’t have it yet, then follow the steps below:

**Anaconda/Miniconda3 installation**

1. [Download the Anaconda installer form the website](https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh) or get the latest version using the line below.
```sh
wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
```

2. Now install accordingly:
```sh
 bash "Anaconda/Miniconda3-<version>-Linux-x86_64.sh"
```

3. Check if you installed correctly by calling the list of you installed packages and tools:
```sh
conda list 
```

4. (Optional) - If you have an older version, update the conda package manager to the latest version:
```sh
conda update -n base conda
```

**Setting up the environment for the pipeline**

1. Now with Conda installed and update it is time to setup the environment that is required for the pipeline.
To install the essential packages and tools, it is important to set your Conda channels in the following order:
```sh
conda config --add channels defaults
conda config --add channels GenomeDK
conda config --add channels r
conda config --add channels bioconda
conda config --add channels conda-forge
```

2. Note that the channels you added first end up in the bottom for lower priority. This is intentional!
Make sure that the these channels are the only ones active. If you have any other channels, you should first remove them as stated below. But first check channels you currently have:
```sh
conda config --show channels
```
If there channels we do not want, then remove them:
```sh
conda config --remove channels <name of channel>
```
Then add the required channels as stated before (in the correct order!)

3. The requirements are provided in the "requirements.txt" that you downloaded form the Github page. Make a new environment and install all the packages:
```sh
conda create --name sc_sn_RNA_seq_pipeline --file requirements.txt
```

## Procedure on how to run the pipeline

**Download the latest CellRanger**

At first, we need the Cellranger tools provided by 10X genomics. You can [download](https://www.10xgenomics.com/support/software/cell-ranger/downloads) and import it with FileZilla or use the line below:
```sh
wget -O cellranger-8.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.0.tar.gz?Expires=1712696337&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=mQFiGKPIsceAJD3GxU48M3cli1ZffCc2IkmE6PV7H~eONzmfz32TIRCTiJmj9cJSj~fcbi-m04-ufoRaCVcxbp-UMMfQ6FSGNGkDrKgYmCcPIimlNtUquZIdumSsa-J0W3UhSi-TD3Jc994FH2uYJXteUeW8RnI8BvGqyk9MDqgKjYwwdbHvzWs-49n65jaaEd95Ht7PJbikQu1ta3V~YGmSxYFSU8TlRA~~KQd8YqJzU4B94HVROX79oiC83PxMcx2IbO7wwJ5rLqEpZ-N7QwehVKRvhiCDcxshaBz4TQwIT4E0MIaBnG45Hf~qM7m15cEDQbTQS~E2qgl3b7Cxsw__"
```
Unzip it using:
```sh
gzip -d cellranger-8.0.0.tar.gz
```
The files will no go to a new directory names **cellranger-8.0.0**, that contains all the cellranger dependencies.

**Running Cellranger**

The pipeline provides the option to skip any of the Cellranger tools if you wish to only carry out Quality-Control and Pre-Processing with the assumption that you already count matrices downloaded from some data base (e.g. SRA) (more on this late...). For now, we will assume that you want to run the full pipeline that includes **1.** Cellranger mkgtf, **2.** Cellranger mkref, **3.** Cellranger count and **4.** QC & Pre-processing.

1. Now go the pipeline directory (that you previously downloaded form the Github page) and open the provided **config.yaml** with your favorite text editor (e.g. Geany or Nano, or whatever you like). Now prove the full path to the installed Cellranger directory in config file as shown below:
```yaml
cellranger_directory: "/path/to/cellranger-8.0.0" 
```

The first step of the pre-processing includes the use of the CellRanger shell utility that is specifically made by 10X Genomics to handle datasets from a wide range of single-cell technologies. Here, it is expected for the user to provide a pair of FASTQ files (Read 1 - Cell barcode/UMI and Read 2 - Insert, and some index files) and a reference gene annotation from any organism of interest. If a reference annotation is not available, the pipeline enables the creation of a custom reference by using a provided reference genome sequence (FASTA file) and a gene annotations (GTF file) that is compatible with the RNA-seq aligner, STAR.

This pipeline is made to compatible with multiple samples in multiple directories. This means that may have the reads from multiple samples in the same folder. Additionally, you can also specify a the path to a folder that contains reads (folders) of multiple samples (or both). The supported examples are shown below:

**Example 1:**
```text
DATA_FOLDER
├── SAMPLE_1_READS
│   ├──  sample1_S1_L001_R1.fastq.gz
│   └──  sample1_S1_L001_RR.fastq.gz
│   └──  sample1_S1_L001_I1.fastq.gz
└── SAMPLE_2_READS
│   ├──  sample2_S1_L001_R1.fastq.gz
│   └──  sampl2e_S1_L001_RR.fastq.gz
│   └──  sample2_S1_L001_I1.fastq.gz
└── SAMPLE_3_READS
│   ├──  sample3_S1_L001_R1.fastq.gz
│   └──  sampl23_S1_L001_RR.fastq.gz
│   └──  sample3_S1_L001_I1.fastq.gz
```
**Example 2:**
```text
SAMPLE_1_2_READS
├── sample1_S1_L001_R1.fastq.gz
├── sample1_S1_L001_RR.fastq.gz
├── sample1_S1_L001_I1.fastq.gz
├── sample2_S1_L001_R1.fastq.gz
├── sample2_S1_L001_RR.fastq.gz
├── sample2_S1_L001_I1.fastq.gz
├── sample3_S1_L001_R1.fastq.gz
├── sample3_S1_L001_RR.fastq.gz
├── sample3_S1_L001_I1.fastq.gz
  ```
**Example 3:**
```text
DATA_FOLDER
├── SAMPLE_1_READS
│   ├──  sample1_S1_L001_R1.fastq.gz
│   └──  sample1_S1_L001_RR.fastq.gz
│   └──  sample1_S1_L001_I1.fastq.gz
└── SAMPLE_2_READS
│   ├──  sample2_S1_L001_R1.fastq.gz
│   └──  sampl2e_S1_L001_RR.fastq.gz
│   └──  sample2_S1_L001_I1.fastq.gz
└── SAMPLE_3_READS
│   ├──  sample3_S1_L001_R1.fastq.gz
│   └──  sampl23_S1_L001_RR.fastq.gz
│   └──  sample3_S1_L001_I1.fastq.gz
│   ├──  sample4_S1_L001_R1.fastq.gz
│   └──  sample4_S1_L001_RR.fastq.gz
│   └──  sample4_S1_L001_I1.fastq.gz
```

Even tough the examples above show that you can run many samples at the same time, keep in mind this also comes at a hefty computation cost. I general, it is more convenient to run one sample at a time (except technical replicates)

2. Now go to config.yaml and provide the the path your reads. In Example 1 and 3 (shown above) the path to your reads would look like `/path/to/DATA_FOLDER`, while in example 1 you should include the full path, including the folder where the reads are located `/path/to/DATA_FOLDER/SAMPLE_1_2_READS`. 
```yaml
reads_directory: "/path/to/DATA_FOLDER" 
```

3. With the assumption that you want to make a new reference profile, you have to go to the config.yaml and provide the path to your .gtf, .fasta (.fa), and specify that you want to run Cellranger mkgtf and Cellranger mkref. Examples are shown below:

**Select Cellranger mkref**
```yaml
make_new_ref: "yes" 
```
**Select Cellranger mkgtf**
```yaml
edit_gtf: "yes" 
```
**Select .gtf directory**
```yaml
gtf_directory: "/path/to/gtf-folder" 
```
**Select .fasta directory**
```yaml
fasta_directory: "/path/to/fasta-folder" 
```
The reference transcriptome created in this process can be used in the future for now samples. In such case leave out `make_new_ref, edit_gtf, gtf_directory, fasta_directory` in config.yaml file and only prove the path to your reads and the path the your previous reference transcriptome:
```yaml
transcriptome_directory: "/path/to/transcriptome-folder:" 
```

Since the expression assay captures transcripts by the poly-A and 3' ends, most reads will align towards that region, including the UTR. If the UTR sequence between (multiple) transcripts happens to be similar, the gene counts at an unwanted locus could become inflated, while missing the true counts for the relevant genes. Taking into account that GTF files often contain entries for non-polyA transcripts, there is a minor chance that these transcripts (and its UTRs) may overlap with actual protein-coding genes. Thus, by default, the GTF file will be manipulated to only include genes that have polyA entries. However, the option is provided to include all genes, long non-coding RNAs, antisense sequences and more... 
To ensure that we only include genes that contain polyA entries (and protein coding genes), by default in config.yaml, we set this parameter to:
```yaml
gene_biotype: "gene_biotype:protein_coding"
```
This is an addition to the `edit_gtf:` parameter and you can leave this out change it according to your need. More options are provided [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-3p-references).

For the STAR read mapping, we include alignments to intronic regions as this greatly improves the number of detected genes in single-nuclei data (due to abundant pre-mRNA) and could potentially improve the detection of cell populations that have low expression of genes (e.g. neutrophils and mast cells). To if you do not want to include intronic reads, then you can go to the config.yaml and adjust the following:
```yaml
cellranger_counts_options: "--include-introns=false"
```
You can choose to include other options as shown [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct).

Finally, CellRanger returns two count matrices (filtered and raw) and by the default option, the filtered matrix will be used for the pre-processing as this contains only cells (droplets) that have at least 500 transcripts (unique molecular identifiers - UMIs). The use of the raw matrix remains optional.

