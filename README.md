# Mutational antigenic profiling of ZIKV E protein
Experiments by Jackson Barr Stuart and Leslie Goo.
Analysis by Jesse Bloom.

## Results
For a summary of the results, see [results/summary/map_analysis.md](results/summary/map_analysis.md), which is the Markdown summary of running the Jupyter notebook [map_analysis.ipynb](map_analysis.ipynb).

Other results are placed in [./results/](results), although not all files are tracked in the GitHub repo.

## Running analysis
First activate the [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment for the analysis.
If you are using the *BloomLab* software on the Fred Hutch computing cluster, you can do this just with:

    conda activate ZIKV_DMS

Otherwise, first build the *ZIKV_DMS* [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment from the [environment_pinned.yml](environment_pinned.yml) or [environment_unpinned.yml](environment_unpinned.yml) file (depending on whether you want fully pinned or unpinned versions), then activate it as above.

After you have activated the *ZIKV_DMS* [conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment, simply run the Python Jupyter notebook [map_analysis.ipynb](map_analysis.ipynb).

To run the notebook automatically and build the HTML summary linked to above, simply run the bash script [run_nbs.bash](run_nbs.bash).
On the Hutch cluster, you will first want to grab a node with 16 cores before doing this.

## Input data
The input data are in [./data/](data):

 - [./data/E.fasta](data/E.fasta): coding sequence of E protein from ZIKV MR766 strain used as parent for mutagenesis.

 - [./data/subamplicon_alignspecs.txt](data/subamplicon_alignspecs.txt): the alignment specs for the [barcoded subamplicon sequencing](https://jbloomlab.github.io/dms_tools2/bcsubamp.html).

 - [./data/samplelist.csv](data/samplelist.csv): all the samples that we sequenced and the locations of the associated deep-sequencing data.
