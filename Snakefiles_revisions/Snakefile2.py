import sys
import os
import glob
from snakemake.utils import makedirs
from pathlib import Path
pipeline = "sc-sn-RNA-seq-Pre-Processing-Pipeline"

### GENERATE "SLURM_LOGS" FOLDER IN CURRENT DIRECTORY FOR ###
makedirs("SLURM_LOGS") # for slurm log output



### SET THE PATH OF THE SNAKEFILE AND CONFIG.YAML/JSON ###
"""
# We specify the name and path config.yaml/json file with --configfile(s), and here we will
  import the config file(s) in the Snakefile.
# Similarly, we specify the path the Snakefile with --snakefile.
"""
from snakemake.exceptions import WorkflowError

# ---- Config hygiene helpers (robust against empty/mistyped config values) ----
def _clean(v):
    """Normalize config values to stripped strings; treat None as empty."""
    if v is None:
        return ""
    # Keep non-strings (e.g. int/bool) but stringify consistently
    return str(v).strip()

def _as_yesno(v, key, default=""):
    """Accept yes/no/true/false/1/0; return "yes" or "no" or "" (if default="")."""
    s = _clean(v).lower()
    if s == "" and default in ("", "yes", "no"):
        return default
    if s in ("y","yes","true","t","1","on"):
        return "yes"
    if s in ("n","no","false","f","0","off"):
        return "no"
    raise WorkflowError(f'Config "{key}" must be yes/no (or true/false/1/0). Got: {v!r}')

def _require_nonempty(v, key, hint=""):
    s = _clean(v)
    if s == "":
        msg = f'Missing required config "{key}".'
        if hint:
            msg += " " + hint
        raise WorkflowError(msg)
    return s

# Ensure we have *some* config loaded (i.e. user passed --configfile or used a profile that loads it)
if not config:
    raise WorkflowError(
        "No config was loaded. Run Snakemake with e.g.:\n"
        "  snakemake --snakefile Snakefile.py --configfile config.yaml\n"
        "(or use a profile that sets --configfile).")

# Make missing keys harmless (avoid KeyError), but validate later depending on mode
EXPECTED_KEYS = [
    "out_directory","reads_directory","sample_names",
    "cellranger_directory","transcriptome_directory",
    "make_new_ref","edit_gtf","gtf_directory","fasta_directory","gene_biotype",
    "cellranger_counts_options","animal",
    "use.local.cache.MT","use.raw.only","skip.DoubletFinder","CPU.cores",
    "min.cell","min.feature","nfeatures","umap.resolution","low.complexity",
    "n.highly_affected.genes","n.markers.genes","canonical.markers",
    "featureLOW","countLOW","featureHIGH","countHIGH","mt.Con",
    "feature.quantile.low","UMI.quantile.low","feature.quantile.high","UMI.quantile.high",
]
for k in EXPECTED_KEYS:
    config.setdefault(k, "")
for k in list(config.keys()):
    config[k] = _clean(config[k])

# ---------------------------------------------------------------------------



### LOAD IN USER-DEFINED OR CURRENT WORK/OUTPUT DIRECTORY ###
"""
# This directory will be used to generate all the output from Cellranger mkgtf,
  mkref, counts, quality control and pre-processing of counts.
# If the user chooses to not run any of the Cellranger tools, but only the 
  quality control and pre-procecing, then the out directory will become the path to the
  new sampes
"""
OUT_DIRECTORY = config["out_directory"]
if OUT_DIRECTORY == "":
    OUT_DIRECTORY = os.getcwd()



### LOAD IN REQUIREMENTS FOR CELLRANGER ###
"""
# All the essential requirements for Cellranger mkgtf, mkref and counts 
  are imported from the config.yaml.
# The requirements include the parameters to run Cellranger counts and optionally 
  generate a new reference with/withour filtered gtf. file.
# The "SAMPLE_NAMES" determines the output folder for Cellranger counts. This is 
  automatically detected based on the names/sample of raw reads. If "SAMPLE_NAMES" is 
  given by the user, then it is assumded that CellRranger outpud already exits and 
  the results will be used for QC and pre-processing step.
"""
# Cellranger software path:
CELLRANGER_DIRECTORY = config["cellranger_directory"] # determines path to seoftware. if not found then no Cellranger tools will be skipped

# Collect sample names based on name of raw reads:
READS_DIRECTORY = config["reads_directory"] # set directory of the raw reads
SAMPLE_NAMES = config["sample_names"]

# ---- Basic mode validation (helps catch config mistakes early) ----
READS_DIRECTORY = _clean(READS_DIRECTORY)
SAMPLE_NAMES = _clean(SAMPLE_NAMES)
if READS_DIRECTORY != "" and SAMPLE_NAMES != "":
    raise WorkflowError(
        'You set BOTH "reads_directory" and "sample_names". Choose one:\n'
        '  • reads_directory: run Cell Ranger count from FASTQs\n'
        '  • sample_names: QC/preprocess existing Cell Ranger outputs')
if READS_DIRECTORY == "" and SAMPLE_NAMES == "":
    raise WorkflowError(
        'You must set either "reads_directory" (FASTQs) OR "sample_names" (existing outputs).\n'
        'In your config, both are empty.')
if READS_DIRECTORY != "" and not os.path.isdir(READS_DIRECTORY):
    raise WorkflowError(f'"reads_directory" does not exist or is not a directory: {READS_DIRECTORY}')
if SAMPLE_NAMES != "" and not os.path.isdir(SAMPLE_NAMES):
    raise WorkflowError(f'"sample_names" does not exist or is not a directory: {SAMPLE_NAMES}')

# If FASTQs mode: Cell Ranger path must be provided
if READS_DIRECTORY != "" and CELLRANGER_DIRECTORY == "":
    raise WorkflowError(
        'FASTQ mode detected (reads_directory is set) but "cellranger_directory" is empty.\n'
        'Set "cellranger_directory" to the folder that contains the cellranger binary.')
if SAMPLE_NAMES == "" : # if SAMPLE_NAMES not given, then it is assumed that a new one must be made based on the name of the reads (which should represent your sample name)
    SN, SNUM, SEQ_LANE, = glob_wildcards(os.path.join(READS_DIRECTORY, "{sn}_S{snum}_L{lane}_R1_001.fastq.gz"))

    # split by "/" if user specifies a nested directory
    for index_sn, sn in enumerate(SN):
        if '/' in sn:
            SN[index_sn] = sn.split("/")[1]

    SN_SNUM_LIST = [] # To obtain new directory names if reads from multiple samples are in the same folder
    for i_sn in range(len(SN)):
        SN_SNUM_LIST.append(f"{SN[i_sn]}_SAMP_{SNUM[i_sn]}")
        SN_SNUM_LIST = list(set(SN_SNUM_LIST))
    
    SAMPLE_NAMES_LIST = [] # To obtain read names that can be used as sample name for Cellranger counnts --id=
    for i_snum in SN_SNUM_LIST:
        read_name = i_snum.split("_SAMP_")
        SAMPLE_NAMES_LIST.append(i_snum.split("_SAMP_")[0]) # take first element, which is read name
    # Put the SAMPLE NAMES AND NUMBER IN A DICTIONARY. (If we keep in list, then wildcards will try all combinations)
    SAMP_DICT = dict(zip(SN_SNUM_LIST, SAMPLE_NAMES_LIST)) # Values become keys and viceversa
    #SAMP_DICT = dict(zip(SAMP_DICT.values(), SAMP_DICT.keys())) # reverse it 



# Cellranger mkref params:
MAKE_NEW_REF = _as_yesno(config["make_new_ref"], "make_new_ref", default="no")
EDIT_GTF = "" # We dont import this parameter from the config until MAKE_NEW_REF is seen... so we give it empty value
if MAKE_NEW_REF == "yes":
    EDIT_GTF = _as_yesno(config["edit_gtf"], "edit_gtf", default="no")
    GENE_BIOTYPE = config["gene_biotype"] # Cellranger mkgtf options
    GTF_DIRECTORY = config["gtf_directory"]
    if GTF_DIRECTORY != "":
        GTF_FILE, = glob_wildcards(os.path.join(GTF_DIRECTORY, "{prefix}.gtf")) # parse name of file from directory (prefix of <prefix>.gtf)
        GTF_FILTERED = os.path.join(GTF_DIRECTORY, f"{GTF_FILE[0]}.filtered.gtf")
        GTF_ORIGINAL = os.path.join(GTF_DIRECTORY, f"{GTF_FILE[0]}.gtf")

# Additional Cellranger mkgtf params:
# We keep fasta directory out of mkref if statement as this file can be used to extract ensembl species for MT genes
FASTA_DIRECTORY = config["fasta_directory"]
if FASTA_DIRECTORY != "":
    FASTA_FILE = glob.glob(f"{FASTA_DIRECTORY}/*")[0].split("/")[-1] # parse name of file from directory
    FASTA_ORIGINAL = os.path.join(FASTA_DIRECTORY, FASTA_FILE) # parse full path and file
    
    # Import specie name
    with open(FASTA_ORIGINAL, 'r') as fa_file:
        fa_list =  fa_file.readlines()
        if " ".join(fa_list[0].split(" ")[0:1]).startswith(">1"): # Ensembl .fa files start with chromosome numbers (often chr1 as starting point)
            SPECIES = FASTA_FILE.split(sep= ".", maxsplit=1)[0] # only parse the specie_genus string from (Ensembl) file name
        else:
            SPECIES = " ".join(fa_list[0].split(" ")[1:3]).replace(" ", "_") # Parse specie name from fasta header and substitute "_" between specie_genus

    TRANSCRIPTOME_MKREF = os.path.join(OUT_DIRECTORY, f"{SPECIES}_mkref_transcriptome") # transcriptome output from mkref
    #TRANSCRIPTOME_MKREF = directory(expand("{samp_outs_dir}/{transcriptome_out_folder}", samp_outs_dir = OUT_DIRECTORY, transcriptome_out_folder = f"{SPECIES}_mkref_transcriptome")), # outdir = os.path.join(OUT_DIRECTORY, f"{SPECIES}_mkref_transcriptome")
# else, we already have a reference transcriptome from previous run, and thus this path will be used instead (if given by user)
TRANSCRIPTOME = config["transcriptome_directory"] # from user

# Cellranger counts options and outputs:
CELLRANGER_COUNTS_OPTIONS = config["cellranger_counts_options"]
# Accept options with or without leading dashes
if CELLRANGER_COUNTS_OPTIONS != "" and not CELLRANGER_COUNTS_OPTIONS.lstrip().startswith("-"):
    CELLRANGER_COUNTS_OPTIONS = "--" + CELLRANGER_COUNTS_OPTIONS
cellranger_counts_output_results = ["web_summary.html"] # not used!
cellranger_counts_matrix_folders = ["filtered_feature_bc_matrix", "raw_feature_bc_matrix"]



### LOAD IN REQUIREMENTS FOR QC/PRE-PROCESSING Rscript ###
"""
# All the essential requirements for qc and pre-processing are imported 
  from the config.yaml.
# All requirements are can be adjusted by user and will imact results.
# The pre-processing has been thoroughly benchmarked and tested on variety of single-cell 
  and single-nuclei data and by default the user does not have to change 
  anything as the most optimized settings have already been selected for you. 
  And some are even hardcoded (e.g. community merging algorithm) considering 
  that their effect on the final results is still in debate.
# In futute interations, some parameters will be added as user-adjustable if we learn more about then.
# The parameters "SAMPLE_NAMES" determines the output folder for external samples.
  This parameter should be the path to external samples that have been downloaded from
  external database. In each samples (folder) there should be at least another folder that countains the
  barcodes.tsv.gz, barcodes.tsv.gz, barcodes.tsv.gz. This folder schould contain the charecters 
  "raw" or "filtered". (or both).
  These external samples the will be used for QC and pre-processing step. 
"""
if SAMPLE_NAMES != "" and READS_DIRECTORY == "": # if path to new samples given (and no reads available), then it is assumed to carry out QC
    OUT_DIRECTORY = SAMPLE_NAMES # set out directory to the directroy of the external samples
    SN, = glob_wildcards(os.path.join(SAMPLE_NAMES, "{sample}.gz")) # find all folder to .gz
    
    MATRIX_LIST = []
    SAMP_LIST = []
    # split by "/" if user specifies a nested directory
    for index_sn, sn in enumerate(SN[::3]): # loop in steps of 3 (because there are three .gz files)
        if '/' in sn:
            SAMP_LIST.append(sn.split("/")[0]) # get the sample folder name
            MATRIX_LIST.append(("/".join(sn.split("/")[1:-1]))) # get absolute path to raw/filtered the counts with the sample folder as the start
    # Make a nested list
    COMB = [list(i) for i in zip(SAMP_LIST, MATRIX_LIST)]
    
    # Make a nested list of unique samples and unique path of counts matrices and put in dict
    SSAMP_DICT = {}
    for samps in COMB:
        if samps[0] not in SSAMP_DICT: # if first element of nested list not in dict
            SSAMP_DICT[samps[0]] = [samps[0]] # then append it as key in dict (the key also becaomes the value)
        SSAMP_DICT[samps[0]].append(samps[1]) # with duplicates gone, now append the second element (path to counts matrix) to the key (sample)
    NEW_COMB = list(SSAMP_DICT.values()) # take the key value and put in list

    # Make dictionary of nested list with first element (sample) as key and second element (path to counts) as list of values
    SAMP_DICT = {i[0]: i[1:] for i in NEW_COMB}



# Check if the required matrices exist in previous cellranger runs
"""
# These two parameters below are greyed out, but we will define them as wildcard params in the QC and pre-processing rule. 
# These two parameters determine how the preprocessing is caried, including ambient RNA removal and extensive filtering.
# This is based on the availability of the raw counts matrix, filtered counts matrix or both.
"""
#COUNTS_PATH_FILTERED = lambda wildcards: next((path_counts for path_counts in SAMP_DICT[wildcards.samples] if 'filtered' in path_counts), "")
#COUNTS_PATH_RAW = lambda wildcards: next((path_counts for path_counts in SAMP_DICT[wildcards.samples] if 'raw' in path_counts), "")

# Name of the graphs per sample
"""
# These parameters below is greyed out, but we will define them as wildcard params in the QC and pre-processing rule.
# This parameter will become the title of the graphs when running the qc and pre-processing.
# In the local version of the Rscript, this parameter has also other funtions and thus greyed out here.
# If the names are too long it can make the graphs look a little bit sloppy. (We can try to parse a piece of the string instead)
"""
#GRAPH_NAMES = lambda wildcards: wildcards.samples

# Animal names Ensembl to extact MT genes:
ANIMAL = config["animal"]
if (FASTA_DIRECTORY == "") or (FASTA_DIRECTORY != "" and ANIMAL != ""):
    ANIMAL = config["animal"]
elif FASTA_DIRECTORY != "" and ANIMAL == "":
    ANIMAL = SPECIES.replace("_", " ")

# Use local MT chance of specie from Conda envrioment
USE_LOCAL_CACHE_MT = config["use.local.cache.MT"].lower()

# Use raw matrix only: 
USE_RAW_ONLY = config["use.raw.only"].lower()
# Initial filtering:
MIN_CELL = config["min.cell"]
MIN_FEATURE = config["min.feature"]

# n highly variable features
nFEATURES = config["nfeatures"]

# UMAP resolution
UMAP_RESOLUTION = config["umap.resolution"]

# For SoupX:
LOW_COMPLEXITY = config["low.complexity"]
n_HIGHLY_AFFECTED_GENES  = config["n.highly_affected.genes"]
n_MARKER_GENES = config["n.markers.genes"]
CANONICAL_MARKERS = os.path.join(workflow.basedir, "markers/gene.txt") # patch to where the pipeline is called (workflow.basedir)

# Extensive filtering:
FEATURE_LOW = config["featureLOW"]
COUNT_LOW = config["countLOW"]
FEATURE_HIGH = config["featureHIGH"]
COUNT_HIGH = config["countHIGH"]
MT_CON = config["mt.Con"]
# By quantile
FEATURE_QUANTILE_LOW = config["feature.quantile.low"]
UMI_QUANTILE_LOW = config["UMI.quantile.low"]
FEATURE_QUANTILE_HIGH = config["feature.quantile.high"]
UMI_QUANTILE_HIGH = config["UMI.quantile.high"]

# DoubletFinder:
CPU_CORES = config["CPU.cores"]
SKIP_DOUBLETFINDER = config["skip.DoubletFinder"]



### FUNTION REQUIREMENTS FOR RULES ###
"""
# Some rules can work with multiple inputs, and one of them is Cellranger mkref. 
  This rules can either use a original uneddited .gtf file from Ensembl 
  or filtered by the user based on specific biotypes (e.g. proteins coding).
  The same is true for the type of reference transcriptome used
# The funtions below determine which input to use based on the outcome of 
  previous rules or user defined paramters.
# One of the functions is specifically made to raise erros is the user is completely clueless
"""
# ERROR funtions:
def raise_error():
    """
    # This funtion returns erros if the some conbination of user parameters ere not set correctly
    # For instance, you will get an error if Cellranger counts want to run, 
      but you have not provided any transcriptome either from your own path or from Cellranger mkref
    """
    if CELLRANGER_DIRECTORY == "" and (READS_DIRECTORY != "" or MAKE_NEW_REF == "yes" or EDIT_GTF == "yes" or TRANSCRIPTOME != ""):
        raise WorkflowError('It apears that you want to run Cellranger tools, but forgot '
        'to import the path to the Cellranger software package. Go back '
        'to the config file and take a close look. If "cellranger_directory" is not given, while other '
        'related parameters are specified (e.g. "reads_directory", "transcriptome_directory", '
        '"make_new_ref", "edit_gtf"), then this error will occur.')
    if (READS_DIRECTORY != "" and (MAKE_NEW_REF == "" and TRANSCRIPTOME == "")):
        raise WorkflowError('It apears that you have do not have a transcriptome from Cellranger counts ' 
        'by either skipping Cellranger mkref ("make_new_ref") or a personal path to a 10X supported ' 
        'transcriptome ("transcriptome_directory"). This only happens if the path to your reads is given ("reads_directory"). '
        'In this case, the pipeline thinks that you want to run Cellranger counts '
        'So please go back to your config file and make sure to specify '
        'a transcriptome or remove the path to your raw reads to make the '
        'pipeline think that you do not want to run any of the Cellranger tools.')
    if MAKE_NEW_REF == "yes" and (GTF_DIRECTORY == "" or FASTA_DIRECTORY == ""):
        raise WorkflowError('It seems you are trying to make a new reference transcriptome '
        'using Cellranker mkref, bot forgot to specify the path to your .gtf or .fa files ("gtf_directory", "fasta_directory") '
        'Go back ro the config file and try again')
    if (SAMPLE_NAMES == "" and READS_DIRECTORY == "") or (SAMPLE_NAMES != "" and READS_DIRECTORY != ""):
        raise WorkflowError('It apears that you have no sample names ("sample_names") and other samples '
        'cannot be detected as the directory to raw reads is also lacking ("reads_directory"). '
        'Please go back to the config file and make sure that you specify a driectory '
        'to your reads if you want to run the full pipeline. Or specify the name of your '
        'previous Cellranger count folders if you only want to run qc/pre-processing again. '
        'This error also occurs when you specify both sample names and derectory to your reads. '
        'In such case, the pipeline doesn not know which one to use...')
    if ANIMAL == "":
        raise WorkflowError('It apears that you have not specified the Ensembl name of the specie you are working with. '
        'This means that the name of your specie cannot be parsed dua to lack of a fasta or user defined parameter. '
        'The specie is important for the QC and Pre-processinf as this information is used to extraxt the '
        'MT genes from the Ensembl database. Go back to the config file and specify either the name of your specie ("animal") '
        'or specify a path to the the fasta file ("fasta_directory")')
raise_error()

# Cellranfer mkref function for .gtf:
def cellranger_mkref_with_or_without_edited_GTF():
    """
    # This function return the original or edited .gtf file based on user parameter
    # If EDIT_GTF is not given by user, then the original GTF file will be used for mkref
    """
    if EDIT_GTF == "yes":
        return(GTF_FILTERED)
    elif EDIT_GTF == "" and MAKE_NEW_REF == "yes":
        return(GTF_ORIGINAL)

# Cellranfer mkref function for transcriptome:
def cellranger_mkref_new_or_old_transcriptome():
    """
    # This function return the user-available transcriptome or new one to be made 
      with Cellranker mkref file based on user parameter
    # This ensures that the rules follow the correct order
    """
    if MAKE_NEW_REF == "yes":
        return(TRANSCRIPTOME_MKREF)
    elif MAKE_NEW_REF == "" and TRANSCRIPTOME != "":
        return(TRANSCRIPTOME)

# QC and preprocessing when smaples are given by user:
def manual_samples_QC():
    """
    # This function returns an already user defined parent directory to the sample folders to
      ensures that the only qc and pre-proccesing runs when sample_sames are given instead of reads.
    """
    if SAMPLE_NAMES != "":
        return(OUT_DIRECTORY)



### RUN QC AND PRE-PROCESSING AND   RENAMING THE KNIT LOCALLY ###
localrules: QC_and_Pre_processing, rename_knits,  # if renaming the knit locally after the cluster operation, feel free to remove it from localrules



### OUTPUT RULES ###
"""
# The fist rules establishes the output of the pipeline. (No idea why snakemake works like that...?)
# The input of 'rule all' will be set as the output of the QC and pre-processing R script
# Global variables here can be used as wildcards in the rules
"""
rule all:
    input:
        expand(
            os.path.join(OUT_DIRECTORY, "{samples}", "HTML_outs", "QC-Preprocessing_{samples}.html"),
            samples = SAMP_DICT.keys(),
        ),




### CELLRANGER RULES ###
"""
# Here we implement the Cellranger mkgtf, mkref and counts rules.
# Cellranger mkgtf will only run if "EDIT_GTF" is set on "yes".
# Cellranger mkref will use the edited/fitered .gtf file if the EDIT_GTF" 
# is set on "yes" (determined by previous funtion)
# Cellranger counts can also use path to users own transcriptome and can run without mkref and mkgtf
"""
# Running  Cellranger mkgtf, mkref and counts on user demand:
# Run Cellranger mkgtf:
if CELLRANGER_DIRECTORY != "" and MAKE_NEW_REF == "yes":
    rule cellranger_mkgtf:
        input:
            GTF_ORIGINAL
        output:
            GTF_FILTERED
        message:
            'Rule {rule} processing'
        params:
            mkgtf_options = f"--attribute={GENE_BIOTYPE}",
            cellranger_directory = CELLRANGER_DIRECTORY,
        shell:
            """
            {params.cellranger_directory}/cellranger mkgtf \
            {input} {output} {params.mkgtf_options}
            """
# Run Cellranger mkref:
    rule cellranger_mkref:
        input:
            fa = FASTA_ORIGINAL,
            gtf = cellranger_mkref_with_or_without_edited_GTF(),
        output:
            directory(expand(TRANSCRIPTOME_MKREF))
        message:
            'Rule {rule} processing'
        params:
            cellranger_directory = CELLRANGER_DIRECTORY,
            animal_genome = f"{SPECIES}_mkref_transcriptome",
            move_to = OUT_DIRECTORY
        shell:
            """
            {params.cellranger_directory}/cellranger mkref \
            --genome={params.animal_genome} \
            --fasta={input.fa} \
            --genes={input.gtf} \
            --nthreads=8
            
            mv -n {params.animal_genome}/ {params.move_to}
            """
# Run Cellranger counts:
if CELLRANGER_DIRECTORY != "" and READS_DIRECTORY != "" :
    rule cellranger_count:
        input:
            reads = READS_DIRECTORY,
            transcriptome = cellranger_mkref_new_or_old_transcriptome(),
        output:
            html = os.path.join(OUT_DIRECTORY, "{samples}", "HTML_outs", "COUNTS-QC_{samples}.html"),
            filtered_matrix = directory(os.path.join(OUT_DIRECTORY, "{samples}", "outs", "filtered_feature_bc_matrix")),
            raw_matrix = directory(os.path.join(OUT_DIRECTORY, "{samples}", "outs", "raw_feature_bc_matrix")),
        message:
            'Rule {rule} processing'
        params:
            cellranger_directory = CELLRANGER_DIRECTORY,
            options = CELLRANGER_COUNTS_OPTIONS,
            previous_run_dir = lambda wildcards: os.path.join(OUT_DIRECTORY, wildcards.samples),
            sample_numbers = lambda wildcards: wildcards.samples,
            sample_names = lambda wildcards: SAMP_DICT[wildcards.samples],
            move_to = OUT_DIRECTORY,
        shell:
            """
            target="{params.move_to}/{params.sample_numbers}"
            if [[ -z "{params.sample_numbers}" || -z "{params.move_to}" || "{params.move_to}" == "/" ]]; then
                echo "Refusing to delete: $target"; exit 1;
            fi
            rm -rf "$target"
            
            {params.cellranger_directory}/cellranger count \
            --id={params.sample_numbers} \
            --transcriptome={input.transcriptome} \
            --fastqs={input.reads} \
            --sample={params.sample_names} \
            --create-bam=false \
            {params.options}
            
            mkdir -p {params.sample_numbers}/HTML_outs
            
            mv {params.sample_numbers}/outs/web_summary.html {params.sample_numbers}/HTML_outs/COUNTS-QC_{params.sample_numbers}.html
            
            mv -n {params.sample_numbers}/ {params.move_to}
            """

# mv -n {params.sample_numbers}/ {params.move_to}
# mkdir -p {params.move_to}/{params.sample_numbers}
# mv -n {params.sample_numbers}/* {params.move_to}/{params.sample_numbers}/
# rmdir --ignore-fail-on-non-empty {params.sample_numbers}


### QC and pre-processing RULES ###
"""
# For every 10X sample, quality control and cleanup will be carried out. 
# The results of the quality control are provided in a HTML document (in the output directory).
# The cleaned up counts matrices will be saved in the cellranger output folder as "new-outs"
  in the same 10X format (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz).
"""
# Run QC when sample names were determined by reads and when results are from recent cellranger run
if SAMPLE_NAMES == "" :
    rule QC_and_Pre_processing:
        input:
            rules.cellranger_count.output,
        output:
            os.path.join(OUT_DIRECTORY, "{samples}", "scRNA.nb.html"),
        message:
            'Rule {rule} processing'
        params:
            parent_dir = OUT_DIRECTORY,
            sample_names = lambda wildcards: wildcards.samples,
            counts_path_filtered = "outs/filtered_feature_bc_matrix",
            counts_path_raw = "outs/raw_feature_bc_matrix",
            graph_names = lambda wildcards: wildcards.samples,
            animal = ANIMAL,
            use_local_cache_MT = USE_LOCAL_CACHE_MT,
            use_raw_only = USE_RAW_ONLY,
            MIN_cell = MIN_CELL,
            MIN_feature = MIN_FEATURE,
            nfeatures = nFEATURES,
            umap_resolution = UMAP_RESOLUTION,
            low_complexity = LOW_COMPLEXITY,
            n_highly_affected_genes = n_HIGHLY_AFFECTED_GENES,
            n_markers_genes = n_MARKER_GENES,
            canonical_markers = CANONICAL_MARKERS,
            feature_LOW = FEATURE_LOW,
            count_LOW = COUNT_LOW,
            feature_HIGH = FEATURE_HIGH,
            count_HIGH = COUNT_HIGH,
            mt_Con = MT_CON,
            feature_quantile_low = FEATURE_QUANTILE_LOW,
            UMI_quantile_low = UMI_QUANTILE_LOW,
            feature_quantile_high = FEATURE_QUANTILE_HIGH,
            UMI_quantile_high = UMI_QUANTILE_HIGH,
            CPU_cores = CPU_CORES,
            skip_DoubletFinder = SKIP_DOUBLETFINDER,
            knits = lambda wildcards: os.path.join(OUT_DIRECTORY, wildcards.samples),
        script:
            os.path.join(workflow.basedir, "script/scRNA.Rmd")

# Run QC when sample names were determined by user (previous samples)
"""
# This rules is a standalone version that can run on any sample from 10X that as at
  least a "filtered_feature_bc_matrix" or "raw_feature_bc_matrix". 
# At least one of the two matrices are required for a succesfull outcome or the Rscript will fail.
# However the lack of the raw matrix will ensure that ambient RNA correction and extensive filtering is not carried out.
"""
if SAMPLE_NAMES != "" :
    rule QC_and_Pre_processing:
        input:
            manual_samples_QC(),
        output:
            os.path.join(OUT_DIRECTORY, "{samples}", "scRNA.nb.html"),
        message:
            'Rule {rule} processing'
        params:
            parent_dir = OUT_DIRECTORY,
            sample_names = lambda wildcards: wildcards.samples,
            counts_path_filtered = lambda wildcards: next((path_counts for path_counts in SAMP_DICT[wildcards.samples] if 'filtered' in path_counts), ""),
            counts_path_raw = lambda wildcards: next((path_counts for path_counts in SAMP_DICT[wildcards.samples] if 'raw' in path_counts), ""),
            graph_names = lambda wildcards: wildcards.samples,
            animal = ANIMAL,
            use_local_cache_MT = USE_LOCAL_CACHE_MT,
            use_raw_only = USE_RAW_ONLY,
            MIN_cell = MIN_CELL,
            MIN_feature = MIN_FEATURE,
            nfeatures = nFEATURES,
            umap_resolution = UMAP_RESOLUTION,
            low_complexity = LOW_COMPLEXITY,
            n_highly_affected_genes = n_HIGHLY_AFFECTED_GENES,
            n_markers_genes = n_MARKER_GENES,
            canonical_markers = CANONICAL_MARKERS,
            feature_LOW = FEATURE_LOW,
            count_LOW = COUNT_LOW,
            feature_HIGH = FEATURE_HIGH,
            count_HIGH = COUNT_HIGH,
            mt_Con = MT_CON,
            feature_quantile_low = FEATURE_QUANTILE_LOW,
            UMI_quantile_low = UMI_QUANTILE_LOW,
            feature_quantile_high = FEATURE_QUANTILE_HIGH,
            UMI_quantile_high = UMI_QUANTILE_HIGH,
            CPU_cores = CPU_CORES,
            skip_DoubletFinder = SKIP_DOUBLETFINDER,
            knits = lambda wildcards: os.path.join(OUT_DIRECTORY, f"{wildcards.samples}"),
        script:
            os.path.join(workflow.basedir, "script/scRNA.Rmd")

# Generate final output for QC and pre-processing
"""
# A final rule that take the output html files of the QC_and_Pre_processing rule and change the file
  to the sample name.
# The results of the quality control are provided in a HTML document (in the output directory).
"""
# Rename the html output to samplen_name.html
rule rename_knits:
    input:
        rules.QC_and_Pre_processing.output
    output:
        os.path.join(OUT_DIRECTORY, "{samples}", "HTML_outs", "QC-Preprocessing_{samples}.html")
    message:
        'Rule {rule} processing'
    params:
        path_to_knit_file = lambda wildcards: os.path.join(OUT_DIRECTORY, wildcards.samples),
        sample_numbers = lambda wildcards: wildcards.samples,
    shell:
        """
        mkdir -p {params.path_to_knit_file}/HTML_outs
        
        mv {params.path_to_knit_file}/scRNA.nb.html {params.path_to_knit_file}/HTML_outs/QC-Preprocessing_{params.sample_numbers}.html
        """