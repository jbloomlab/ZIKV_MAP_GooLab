# Mutational antigenic profiling of ZIKV E protein
Experiments by Caroline Kikawa, Jackson Barr Stuart and Leslie Goo.
Analysis by Jesse Bloom and Caroline Kikawa.

## Results
For a summary of the results, see [results/notebooks/selections_analysis.md](results/notebooks/selections_analysis.md), which is the Markdown summary of running the Jupyter notebook [selections_analysis.ipynb](selections_analysis.ipynb).

Other results are placed in [./results/](results), although not all files are tracked in the GitHub repo.

For the data and analysis used in the neutralization assays and visualizations for the paper, see the notebooks and results in [./paper_figures/](paper_figures) 

## Running analysis
First activate the [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment for the analysis. If you have prebuilt the relevant environments, you can do this just with:

    conda activate zikv_dmstools2

Otherwise, first build the *zikv_dmstools2* [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment from the [environments/environment_dmstools2.yml](environments/environment_dmstools2.yml) file, then activate it as above.

After you have activated the either [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment, simply run the Python Jupyter notebooks: [selections_analysis.ipynb](map_analysis.ipynb) or [polyclonal_analysis.ipynb](polyclonal_analysis.ipynb). On the Hutch cluster, you will first want to grab a node with 16 cores before doing this.

## Input data
The input data are in [./data/](data):

 - [./data/E.fasta](data/E.fasta): coding sequence of E protein from ZIKV MR766 strain used as parent for mutagenesis.

 - [./data/subamplicon_alignspecs.txt](data/subamplicon_alignspecs.txt): the alignment specs for the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html).

 - [./data/samplelist.csv](data/samplelist.csv): all the samples that we sequenced and the locations of the associated deep-sequencing data.
 
  - [./data/concat_fastq](data/concat_fastq): a folder containing a script, input, output and documentation for a program that will accept multiple FASTQ files for a single sample and concatenate them together.
  
   - [./data/functonal_data](data/functional_data): a folder containing [amino acid preferences](data/functional_data/unscaled_prefs.csv) and [site-wise mutational tolerance data](data/functional_data/struct_props_mut_tol.csv) published in [Sourisseau et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31511387/) 
