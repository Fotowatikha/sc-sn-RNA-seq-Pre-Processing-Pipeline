# sn/sn-RNA-seq Pre-Processing Pipeline

Now considering that the emergence of scRNA-seq techniques provides the framework to study gene expression variability in tissue of interest, the pre-processing of high dimensional RNA-seq data also becomes an important topic. To this day, droplet-based technologies remain the technique of choice when it comes to capturing and sequencing a large number of individual cells. This method essentially barcodes single cells and tags each transcript with a unique molecular identifier (UMI) within individual droplets, which can significantly increase the throughput up to 10,000 cells per analysis. However, this technique is not free of noise, and one must carefully pre-process the data before any unbiased downstream analysis can be caried out. 

**Here, I will provide you with all the necessary information on how to use the single-cell/nuclei pre-processing pipeline, either in your own environment or on the HPC cluster.**

## About

The first step of the pre-processing pipeline includes the use of the CellRanger shell utility that is specifically made by 10X Genomics to handle datasets from a wide range of single-cell technologies. In order to generate your single cell feature counts for the library, we utilize [**cellranger count**](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count) that takes your reads (FASTQ files) and performs alignment, filtering, barcode counting and UMI counting. CellRanger already provides pre-built human, mouse, and barnyard (human & mouse) [reference packages](https://www.10xgenomics.com/support/software/cell-ranger/downloads) for read alignment and gene expression quantification. If you do not have a reference transcriptome, the pipeline enables the creation of a custom reference by using a provided reference genome sequence (FASTA file) and a gene annotations (GTF file) that is compatible with the RNA-seq aligner, STAR. Here, we utilize [**cellranger mkgtf**](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-mr) and [**cellranger mkref**](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-3p-references#) to generate the reference transcriptome for your specie of interest.

By consider the gene count matrices as the starting point after mapping reads to the reference, the major steps in this pre-processing is to: **(1)** Provide you with quality control (QC) plots to gain insight on the overall quality of the cells prior to any extensive filtering; **(2)** Correcting the data from cell-free ambient RNA (using [**SoupX**](https://github.com/constantAmateur/SoupX)); **(3)** Extensive filtering to remove droplets that are unlikely to represent intact individual cells (using [**DoubeletFinder**](https://github.com/chris-mcginnis-ucsf/DoubletFinder)); **(4)** Removal of droplets that violate the assumption for containing one cell; **(5)** Provide you with QC plots to gain insight on the quality of the data post-filtering.

A schematic image of the steps in the pipeline is shown below:

<img width="704" alt="Pipeline - Fig  1" src="https://github.com/Fotowatikha/sc-sn-RNA-seq-pre-processing/assets/157910396/d98fb83c-9fe7-4f96-9e55-42ac07274cbd">

## Download the pipeline files from GitHub

You can either [download the files for the pipeline](https://github.com/Fotowatikha/sc-sn-RNA-seq-pre-processing/tree/main) from my Github page or download it it by using the "Code button" on the Github page and import it using Filezilla or use the line below if you have Git installed:
```sh
git clone https://github.com/Fotowatikha/sc-sn-RNA-seq-pre-processing.git
```

## Setup Conda environment

This pipeline is made to run on a Conda environment containing all the crucial packages for Python, Snakemake and R. Here, i assume that you already have Anaconda or Miniconda installed. In case you don’t have it yet, then follow the steps below:

**Anaconda/Miniconda3 installation:**

**1.** [Download the Anaconda installer form the website](https://repo.anaconda.com/archive/Anaconda3-2024.02-1-Linux-x86_64.sh) or get the latest version using the line below.
```sh
wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
```

**2.** Now install accordingly:
```sh
 bash "Anaconda/Miniconda3-<version>-Linux-x86_64.sh"
```

**3.** Check if you installed correctly by calling the list of you installed packages and tools:
```sh
conda list 
```

**4.** (Optional) - If you have an older version, update the conda package manager to the latest version:
```sh
conda update -n base conda
```

**Setting up the environment for the pipeline:**

**1.** Now with Conda installed and update it is time to setup the environment that is required for the pipeline.
To install the essential packages and tools, it is important to set your Conda channels in the following order:
```sh
conda config --add channels defaults
conda config --add channels GenomeDK
conda config --add channels r
conda config --add channels bioconda
conda config --add channels conda-forge
```

**2.** Note that the channels you added first end up in the bottom for lower priority. This is intentional!
Make sure that the these channels are the only ones active. If you have any other channels, you should first remove them as stated below. But first check channels you currently have:
```sh
conda config --show channels
```
If there channels we do not want, then remove them:
```sh
conda config --remove channels <name of channel>
```
Then add the required channels as stated before (in the correct order!)

**3.** The requirements are provided in the **"requirements.txt"** that you downloaded form the Github page. Make a new environment and install all the packages:
```sh
conda create --name sc_sn_RNA_seq_pipeline --file requirements.txt
```

## Procedure on how to run the pipeline

**Download the latest CellRanger software:**

At first, we need the Cellranger tools provided by 10X genomics. You can [download](https://www.10xgenomics.com/support/software/cell-ranger/downloads) and import it with FileZilla or use the line below:
```sh
wget -O cellranger-8.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.0.tar.gz?Expires=1712696337&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=mQFiGKPIsceAJD3GxU48M3cli1ZffCc2IkmE6PV7H~eONzmfz32TIRCTiJmj9cJSj~fcbi-m04-ufoRaCVcxbp-UMMfQ6FSGNGkDrKgYmCcPIimlNtUquZIdumSsa-J0W3UhSi-TD3Jc994FH2uYJXteUeW8RnI8BvGqyk9MDqgKjYwwdbHvzWs-49n65jaaEd95Ht7PJbikQu1ta3V~YGmSxYFSU8TlRA~~KQd8YqJzU4B94HVROX79oiC83PxMcx2IbO7wwJ5rLqEpZ-N7QwehVKRvhiCDcxshaBz4TQwIT4E0MIaBnG45Hf~qM7m15cEDQbTQS~E2qgl3b7Cxsw__"
```
Unzip it using:
```sh
gzip -d cellranger-8.0.0.tar.gz
```
The files will no go to a new directory named **cellranger-8.0.0**, that contains all the CellRanger dependencies.

**Setup Cellranger:**

The pipeline provides the option to skip any of the Cellranger tools if you wish to only carry out QC and Pre-Processing with the assumption that you already have count matrices downloaded from some data base (e.g. SRA) (more on this late...). For now, we will assume that you want to run the full pipeline that includes **1.** Cellranger mkgtf, **2.** Cellranger mkref, **3.** Cellranger count and **4.** QC & Pre-processing.

**1.** Now go the pipeline directory (that you previously downloaded form the Github page) and open the provided **config.yaml** with your favorite text editor (e.g. Geany or Nano, or whatever you like). Provide the full path to the installed Cellranger directory in the config file as shown below:
```yaml
cellranger_directory: "/path/to/cellranger-8.0.0" 
```

For CellRanger to work, it is expected for you to provides a pair of FASTQ files (Read 1 - Cell barcode/UMI and Read 2 - Insert, Index files) and a reference gene annotation from any organism of interest. If a reference annotation is not available, the pipeline enables the creation of a custom reference by using a provided reference genome sequence (FASTA file) and a gene annotations (GTF file) that is compatible with the RNA-seq aligner, STAR. You can download these from [Ensembl](https://useast.ensembl.org/index.html).

This pipeline is made to compatible with multiple samples in multiple directories. This means that you may have the reads from multiple samples in the same folder. Additionally, you can also specify a the path to a parent-directory that contains multiple samples (folders). 

The supported examples are shown below:

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
│   └──  sample3_S1_L001_RR.fastq.gz
│   └──  sample3_S1_L001_I1.fastq.gz
│   ├──  *sample4_S1_L001_R1.fastq.gz*
│   └──  *sample4_S1_L001_RR.fastq.gz*
│   └──  *sample4_S1_L001_I1.fastq.gz*
```
**Example 3:**
```text
SAMPLE_1_2_READS
├── sample1_S1_L001_R1.fastq.gz
├── sample1_S1_L001_RR.fastq.gz
├── sample1_S1_L001_I1.fastq.gz
├── *sample2_S1_L001_R1.fastq.gz*
├── *sample2_S1_L001_RR.fastq.gz*
├── *sample2_S1_L001_I1.fastq.gz*
├── **sample3_S1_L001_R1.fastq.gz**
├── **sample3_S1_L001_RR.fastq.gz**
├── **sample3_S1_L001_I1.fastq.gz**
  ```


Even tough the examples above show that you can run many samples at the same time, keep in mind this also comes at a hefty computation cost. In general, it is more convenient to run one sample at a time (except technical replicates)

**2.** Now go to config.yaml and provide the path your reads (FASTQs). In Example 1 and 2 (shown above) the path to your reads would look like `/path/to/DATA_FOLDER`. In example 3 you should include the full path, including the folder where the reads are located: `/path/to/DATA_FOLDER/SAMPLE_1_2_READS`. 
```yaml
reads_directory: "/path/to/DATA_FOLDER" 
```

**3.** Now select a directory where you want your output files to be generated. By selecting such directory, all your results from multiple samples will be stored here, including all raw gene counts matrices, corrected gene counts matrices and HTML outputs (containing images for QC). Go to the config.yaml and provide the an outs directory. If you do not provide such directorty, your current working directory will be automatically selected (not recommended, unless you run it from the pipeline folder where the Snakefile is located):
```yaml
out_directory: "/path/to/output-folder" 
```
**DO NOT ATTEMPT TO RUN THE PIPELINE WITHIN THE out_directory OR reads_directory. Due to a bug in "CellRanger count" some temporary folders will not be deleted and result in a RuntimeError** 

**4.** With the assumption that you want to make a new reference profile, you have to go to the config.yaml and provide the path to your .gtf, .fasta (.fa), and specify that you want to run Cellranger mkgtf and Cellranger mkref. If you do not want to edit your .gtf file, then you can leave out `edit_gtf` Examples are shown below:
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
The reference transcriptome created in this process can be used in the future for new samples (form the same specie of course...). In such case leave out `make_new_ref, edit_gtf, gtf_directory, fasta_directory` in config.yaml file and only prove the path to your reads and the path the your previous reference transcriptome:
```yaml
transcriptome_directory: "/path/to/transcriptome-folder:" 
```

**5.** Since the expression assay captures transcripts by the poly-A and 3' ends, most reads will align towards that region, including the UTR. If the UTR sequence between (multiple) transcripts happens to be similar, the gene counts at an unwanted locus could become inflated, while missing the true counts for the relevant genes. Taking into account that GTF files often contain entries for non-polyA transcripts, there is a minor chance that these transcripts (and its UTRs) may overlap with actual protein-coding genes. Thus, by default, the GTF file will be manipulated to only include genes that have polyA entries. However, the option is provided to include all genes, long non-coding RNAs, antisense sequences and more... 
To ensure that we only include genes that contain polyA entries (and protein coding genes), by default in config.yaml, we set this parameter to:
```yaml
gene_biotype: "gene_biotype:protein_coding"
```
This is an addition to **Cellranger mkgtf** that select by `edit_gtf:` parameter. You can also choose to leave this out or change it according to your need. More options are provided [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-3p-references).

**6.** For the STAR read mapping, we include alignments to intronic regions as this greatly improves the number of detected genes in single-nuclei data (due to abundant pre-mRNA) and could potentially improve the detection of cell populations that have low expression of genes (e.g. neutrophils and mast cells). If you do not want to include intronic reads, then you can go to the config.yaml and adjust the following:
```yaml
cellranger_counts_options: "--include-introns=false"
```
You can choose to include other options as shown [here](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct).


**Setup the Quality-Control and Pre-Processing:**

In the next step of the pipeline, a QC will be carried out. To do so, we leverage the gene counts matrix output provided by CellRanger (filtered_feature_bc_matrix and raw_feature_bc_matrix) and implement it in the pipeline using the [**Seurat R.Package**](https://satijalab.org/seurat/) where all the processed information will be stored a list of vectors. Furthermore, we will use the same counts matrices to apply the pre-processing and filter. After this step, you will be provided with a new cleaned-up counts matrices free of bias and an HTML document (Rmd format) with all the necessary graphs to provide insight in the quality of the processed data.

Here, you can go to the config.yaml and change many parameters to your liking. However, our default parameters have been thoroughly tested on both high and low complexity data. So it is recommended to try the default parameters first. Alright, now lets go trough all the parameters you can change...

**1.** The QC will present to you the distribution of mitochondrial (MT) reads per cell. In this step, we automatically import the MT genes for your specie from Ensembl, but you must correctly specify the name of your specie in the config.yaml as shown in the line below. Now if you run the full pipeline (which includes CellRagner), we will parse the name of your species from the previously provided FASTA file. If you choose to only run the QC and Pre-Processing step of the pipeline (more on this later), then you probably have not specified a path to your FASTA genome assembly (you still can!), and thus you must fill in the the name of the specie yourself. Make sure to fill it in correctly according to the [Ensebml Biomart Species](https://useast.ensembl.org/info/about/species.html) list. An example of the pig is shown below:
```yaml
animal: "Sus scrofa"
```

**2.** Next, you are provided with the option to either use the filtered counts matrix **(filtered_feature_bc_matrix, filtered by 500 UMIs/cell)** or raw counts matrix **(raw_feature_bc_matrix, includes all droplets)**. Be default, we will leverage the filtered matrix and use the raw only to estimate the background contamination fraction (ambient RNA). If you wish to choose for the raw counts to detect cell types that contain less than 500 transcripts per cell, the pipeline will automatically remove all droplets that contain less than 100 genes as an initial filtering step before the pre-processing is carried out. Now if you want to to use the raw matrix only, then go the config.yaml and adjust the following (not recommended, unless you think you are missing some leukocytes, like Neutrophils):
```yaml
use.raw.only: "" 
```
We even allow you to deviate from the 100 genes/cell threshold if you choose to work with the raw matrix (be thoughtful about this, do you actually expect to find cells that express <100 genes? Or are you in love with semi-empty droplets?). If you want to change this parameters, do it as shown below:
```yaml
min.feature: "" 
```
You can also choose to remove genes that are not expressed in a given number of cells. If you do this, then then correction of ambient RNA will fail as the original raw matrix will have different number of genes when compared to your filter one. If you still choose to do so, then make sure to specify it in config.yaml as shown below and to skip on ambient RNA correction (more on this soon...)
```yaml
min.cell: "" 
```

**3.** On important step before clustering the cells is to select a set of highly variable genes that help us to distinguish between the cells. By selecting highly variable genes (based on standardized variance), we ensure that house-keeping genes with similar expressions across different cells do not interfere with the unsupervised clustering. There is no definitive answer on how many genes one should include for clustering. Thus it is recommended to carefully consider the right number or by trying a range of 2000-5000 variable genes. Generally, low-complexity data (i.e. data with similar cells) will have lower number of variable genes when compared to complex tissue. By default, we will select 2500 highly variable genes as this shown to work perfect with both high and semi-low complexity data. If you choose to deviate from our default, then this parameter in the config.yaml as shown in the line below:
```yaml
nfeatures: "" 
```
Another parameter than impact the total number of the clustering is the resolution of the Louvain algorithm that joins to initial clusters post KNN by joining the groups interactively until it fully converges. We have found that a resolution of 0.8 does a decent job by not under or over clustering high and semi-low complexity data. Since ambient RNA correction partially relies on finding differential expressed genes between clusters to estimate number transcripts that are contributing to the calculated background contamination (form empty droplet in the raw matrix), by changing the Louvain resolution, you can make sure that sufficient number of clusters are formed in your ultra-complex data. Only do do this if you have high percentage of ambient RNA in your cells in combination with low-complex data but sufficient number of reads per cell. If you already have very few read per cell, it might the result of RNA breakdown that contributed to higher ambient RNA concentrations in the first place. You can chang this parameter in the config.yaml as shown below:
```yaml
umap.resolution: "" 
```

**4.** In case you suspect that your ambient RNA correction results are biased due to low complexity data, then feel free to leave out this step as shown in the line below. If you have previously selected to remove genes that were not expressed in a given number of cells (`min.cell:`), the you must skip ambient RNA correction in the config.yaml as shown below:
```yaml
low.complexity: "yes" 
```
When you do correct for ambient RNA, we offer the option to show expression of marker candidates that were used for the estimation of the contamination fraction and their level of correction. By default, top 3 markers are used for plotting, but you are provided with the option to choose anything between 0-100. change this in the config.yaml as shown below:
```yaml 
n.markers.genes: "" 
```
You can also choose to show the most-affected genes and their level of correction. As the most-affected genes are not expressed by many cells, these plots might not be very representative. This option is tuned off by default, but you may turn it on and select the number plots.
```yaml
n.highly_affected.genes: ""
```
Additionally, we over the option to also include your genes of interest. To do, you must go to the pipeline folder and located the gene.txt file inside the markers folder (`pipeline/markers/gene.txt`). Open the **gene.txt** with your favorite text editor and put in the name of your genes of interest as shown in the image below (make sure to save it accordingly in the same folder). Make sure the gene names are the same as found in your counts matrix. If the names are not the same, you will get a runtime errro!

<img width="582" alt="Scherm­afbeelding 2024-04-10 om 06 32 54" src="https://github.com/Fotowatikha/sn-sn-RNA-seq-Pre-Processing-Pipeline/assets/157910396/fdaa4d6a-559b-4f4c-b111-bb5c3f2b53fd">

**5.** During the extensive filtering, by default, we automatically filter out cells with a percentage of MT reads above 10%, with gene and transcript counts below the 2% quantile as a lower limit. No filtering is applied on the higher limit as these cells could represent doublets that we want to classify by DoubletFinder instead of removing them by force. However, if you have prior knowledge on the dataset, filtering can also be applied by user-defined cutoff values for the number of genes, transcripts and percentage of MT reads. Additionally, you may also adjust the upper and lower quantiles for the parameters mentioned above. It is important to consider that by removing cells in the upper limit by user-defined cutoff values or quantiles, there is a good chance that true doublets are being removed before identified in an algorithmic fashion. You can adjust the these parameters in the config.yaml as shown below:
```yaml
# By numbers

featureLOW: ""
countLOW: ""
featureHIGH: ""
countHIGH: ""

# By quantile (between 0-100&):

feature.quantile.low: "" 
UMI.quantile.low: ""
feature.quantile.high: ""
UMI.quantile.high: ""

# MT percentage (between 0-100%):

mt.Con: ""
```
Note that we have protected the parameters above against nonsensical values... For instance, you cannot filter out more genes than the maximum number of genes in your data. If you do so, the filtering based on your given cutoff value will be skipped. Likewise, you cannot filter a given MT percentage above the maximum in your data. In such case, the parameter will change to default of 10%. The same is true if the lower quantiles are larger then upper the quantiles during filtering. In such case, the quantiles will be set to default.
Any mistakes in these steps will not raise warnings and the code will continue with default parameters, so be careful (or get ready for some unexpected results!).

**5.** The final step in the includes filtering out doublets (multiplets) from your data. If you wish to skip this procedure, then go to the config.yaml and change to parameters as shown below. Even if you choose to skip Doubletfinder, we will still run it to show you the quality of the data, but the true doublets will not be filtered out from the counts matrix.
You can also change the number of CPU cores used during the doublet detection. Since the QC and Pre-Processing step runs on the local cluster (even if you submit it as a job), the number of cores are set to 1 to ensure that you do not compete for resources. You can choose anything between 1-8. We have set to maximum to 8 so that you do not slow down the cluster for your colleagues.
```yaml
skip.DoubletFinder: "yes"

CPU.cores: ""
```
