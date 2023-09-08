# ZIKV mutational antigenic profiling with broadly neutralizing antibodies

CK prepped samples for deep mutational scanning using known broadly neutralizing antibodies: EDE1-C8, EDE1-C10, and MZ4 alongside ZIKV-targeted mAbs ZKA64 and ZV-67 as positive controls. The virus library is a ZIKV MR766 E-protein virus library first described in Sourisseau et al J Virology.


**EDE1-C8** was prepared at 1800 ng/mL ("1800").

**EDE1-C10** was prepared at 300 ng/mL ("C10_300").

**MZ4** was prepared at 4800 ng/mL ("MZ4_4800").

**SiGN-3C** was prepared at two concentrations: 20 ug/mL ("SiGN_20") and 10 ug/mL ("SiGN_10").

**ZV-67** was prepared at 40 ug/mL ("ZV67_40000").

All antibody concentrations were incubated with 800,000 IU (MOI 1) of virus library and allowed to infect Vero cells for 24 hours prior to isolating total RNA - see [CK018](https://benchling.com/s/etr-YtBt2dzXr60FGcLKYF1C?m=slm-zQusZO8o9Qc9Egy4lYhe) for details. 

Viral cDNA libraries were generated for next-generation sequencing as outlined in [CK019](https://benchling.com/s/etr-Vp41lo9pFVD1ZcXYw0bM?m=slm-m2YLmBU986a39SRUKDYz).


<h1>Table of Contents<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span>
    <a href="#Deep-mutational-scanning-of-ZIKV-E-protein" data-toc-modified-id="Deep-mutational-scanning-of-ZIKV-E-protein-1">Deep mutational scanning of ZIKV E protein</a></span><ul class="toc-item"><li><span>
    <a href="#Set-up-for-analysis" data-toc-modified-id="Set-up-for-analysis-1.1">Set up for analysis</a></span></li><li><span>
    <a href="#Library-sequencing-quality-control" data-toc-modified-id="Library-sequencing-quality-control-1.2">Library sequencing quality control</a></span></li><li><span>
    <a href="#Differential-selection" data-toc-modified-id="Differential selection-1.3">Differential selection</a></span></li><li><span>
    <a href="#Figures-for-paper" data-toc-modified-id="Figures-for-paper-1.4">Figures for paper</a>

# Set up analysis


```python
# import tools
import os
import glob
import numpy as np
import pandas as pd
import scipy
from functools import reduce
from Bio import SeqIO
import dms_tools2
import dmslogo
from dms_tools2.ipython_utils import showPDF
import seaborn as sns
import matplotlib.pyplot as plt
```


```python
# set use-existing clause for dms_tools2
use_existing = 'yes'
```


```python
# ID input/output directories
samplelist = './data/samplelist.csv'
datadir = './data/'
resultsdir = './results/'
os.makedirs(resultsdir, exist_ok=True)
```


```python
# reference sequences
Erefseq = './data/E.fasta'
subamplicon_alignspecs = './data/subamplicon_alignspecs.txt'
```


```python
# read in E protein sequence data 
refseq = SeqIO.read(Erefseq, 'fasta')
E_seq = refseq.seq
E_prot = E_seq.translate()

# Bio.SeqIO documention https://biopython.org/wiki/SeqIO
```


```python
# read in sample list and add the sampleID 'name'
samples = (pd.read_csv(samplelist, index_col=False))
# pd.set_option('display.max_colwidth', None)
# code necessary for more recent versions of pandas (1.above) to set max column width

samples.insert(0, 'name', (samples['library'] + '-' + samples['selection']))
samples.style.hide_index()
```




<style  type="text/css" >
</style><table id="T_a32e663e_490a_11ee_be5e_ac1f6bca632c" ><thead>    <tr>        <th class="col_heading level0 col0" >name</th>        <th class="col_heading level0 col1" >library</th>        <th class="col_heading level0 col2" >antibody</th>        <th class="col_heading level0 col3" >selection</th>        <th class="col_heading level0 col4" >percent_infectivity</th>        <th class="col_heading level0 col5" >date</th>        <th class="col_heading level0 col6" >R1</th>    </tr></thead><tbody>
                <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow0_col0" class="data row0 col0" >lib1-C10-300</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow0_col1" class="data row0 col1" >lib1</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow0_col2" class="data row0 col2" >EDE1-C10</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow0_col3" class="data row0 col3" >C10-300</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow0_col4" class="data row0 col4" >0.11%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow0_col5" class="data row0 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow0_col6" class="data row0 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib1_C10_300_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow1_col0" class="data row1 col0" >lib1-C8-1800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow1_col1" class="data row1 col1" >lib1</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow1_col2" class="data row1 col2" >EDE1-C8</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow1_col3" class="data row1 col3" >C8-1800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow1_col4" class="data row1 col4" >0.01%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow1_col5" class="data row1 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow1_col6" class="data row1 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib1_C8_1800_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow2_col0" class="data row2 col0" >lib1-MZ4-4800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow2_col1" class="data row2 col1" >lib1</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow2_col2" class="data row2 col2" >MZ4</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow2_col3" class="data row2 col3" >MZ4-4800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow2_col4" class="data row2 col4" >0.46%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow2_col5" class="data row2 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow2_col6" class="data row2 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib1_MZ4_4800_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow3_col0" class="data row3 col0" >lib1-SIgN-20000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow3_col1" class="data row3 col1" >lib1</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow3_col2" class="data row3 col2" >SIgN-3C</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow3_col3" class="data row3 col3" >SIgN-20000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow3_col4" class="data row3 col4" >0.21%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow3_col5" class="data row3 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow3_col6" class="data row3 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib1_SIgN_20000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow4_col0" class="data row4 col0" >lib1-ZV67-40000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow4_col1" class="data row4 col1" >lib1</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow4_col2" class="data row4 col2" >ZV-67</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow4_col3" class="data row4 col3" >ZV67-40000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow4_col4" class="data row4 col4" >0.84%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow4_col5" class="data row4 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow4_col6" class="data row4 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib1_ZV67_40000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow5_col0" class="data row5 col0" >lib1-no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow5_col1" class="data row5 col1" >lib1</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow5_col2" class="data row5 col2" >no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow5_col3" class="data row5 col3" >no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow5_col4" class="data row5 col4" >100</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow5_col5" class="data row5 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow5_col6" class="data row5 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib1_no_antibody_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow6_col0" class="data row6 col0" >lib2-C10-300</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow6_col1" class="data row6 col1" >lib2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow6_col2" class="data row6 col2" >EDE1-C10</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow6_col3" class="data row6 col3" >C10-300</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow6_col4" class="data row6 col4" >0.23%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow6_col5" class="data row6 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow6_col6" class="data row6 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib2_C10_300_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow7_col0" class="data row7 col0" >lib2-C8-1800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow7_col1" class="data row7 col1" >lib2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow7_col2" class="data row7 col2" >EDE1-C8</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow7_col3" class="data row7 col3" >C8-1800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow7_col4" class="data row7 col4" >0.13%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow7_col5" class="data row7 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow7_col6" class="data row7 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib2_C8_1800_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow8_col0" class="data row8 col0" >lib2-MZ4-4800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow8_col1" class="data row8 col1" >lib2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow8_col2" class="data row8 col2" >MZ4</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow8_col3" class="data row8 col3" >MZ4-4800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow8_col4" class="data row8 col4" >0.02%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow8_col5" class="data row8 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow8_col6" class="data row8 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib2_MZ4_4800_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow9_col0" class="data row9 col0" >lib2-SIgN-10000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow9_col1" class="data row9 col1" >lib2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow9_col2" class="data row9 col2" >SIgN-3C</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow9_col3" class="data row9 col3" >SIgN-10000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow9_col4" class="data row9 col4" >0.35%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow9_col5" class="data row9 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow9_col6" class="data row9 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib2_SIgN_10000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow10_col0" class="data row10 col0" >lib2-SIgN-20000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow10_col1" class="data row10 col1" >lib2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow10_col2" class="data row10 col2" >SIgN-3C</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow10_col3" class="data row10 col3" >SIgN-20000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow10_col4" class="data row10 col4" >1.36%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow10_col5" class="data row10 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow10_col6" class="data row10 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib2_SIgN_20000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow11_col0" class="data row11 col0" >lib2-ZV67-40000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow11_col1" class="data row11 col1" >lib2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow11_col2" class="data row11 col2" >ZV-67</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow11_col3" class="data row11 col3" >ZV67-40000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow11_col4" class="data row11 col4" >0.89%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow11_col5" class="data row11 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow11_col6" class="data row11 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib2_ZV67_40000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow12_col0" class="data row12 col0" >lib2-no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow12_col1" class="data row12 col1" >lib2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow12_col2" class="data row12 col2" >no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow12_col3" class="data row12 col3" >no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow12_col4" class="data row12 col4" >80.59%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow12_col5" class="data row12 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow12_col6" class="data row12 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib2_no_antibody_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow13_col0" class="data row13 col0" >lib3-C10-300</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow13_col1" class="data row13 col1" >lib3</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow13_col2" class="data row13 col2" >EDE1-C10</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow13_col3" class="data row13 col3" >C10-300</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow13_col4" class="data row13 col4" >0.39%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow13_col5" class="data row13 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow13_col6" class="data row13 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3_C10_300_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow14_col0" class="data row14 col0" >lib3-C8-1800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow14_col1" class="data row14 col1" >lib3</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow14_col2" class="data row14 col2" >EDE1-C8</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow14_col3" class="data row14 col3" >C8-1800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow14_col4" class="data row14 col4" >0.23%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow14_col5" class="data row14 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow14_col6" class="data row14 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3_C8_1800_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow15_col0" class="data row15 col0" >lib3-MZ4-4800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow15_col1" class="data row15 col1" >lib3</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow15_col2" class="data row15 col2" >MZ4</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow15_col3" class="data row15 col3" >MZ4-4800</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow15_col4" class="data row15 col4" >0.05%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow15_col5" class="data row15 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow15_col6" class="data row15 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3_MZ4_4800_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow16_col0" class="data row16 col0" >lib3-SIgN-10000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow16_col1" class="data row16 col1" >lib3</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow16_col2" class="data row16 col2" >SIgN-3C</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow16_col3" class="data row16 col3" >SIgN-10000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow16_col4" class="data row16 col4" >1.19%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow16_col5" class="data row16 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow16_col6" class="data row16 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3_SIgN_10000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow17_col0" class="data row17 col0" >lib3-SIgN-20000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow17_col1" class="data row17 col1" >lib3</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow17_col2" class="data row17 col2" >SIgN-3C</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow17_col3" class="data row17 col3" >SIgN-20000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow17_col4" class="data row17 col4" >0.88%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow17_col5" class="data row17 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow17_col6" class="data row17 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3_SIgN_20000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow18_col0" class="data row18 col0" >lib3-ZV67-40000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow18_col1" class="data row18 col1" >lib3</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow18_col2" class="data row18 col2" >ZV-67</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow18_col3" class="data row18 col3" >ZV67-40000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow18_col4" class="data row18 col4" >0.01%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow18_col5" class="data row18 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow18_col6" class="data row18 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3_ZV67_40000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow19_col0" class="data row19 col0" >lib3-no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow19_col1" class="data row19 col1" >lib3</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow19_col2" class="data row19 col2" >no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow19_col3" class="data row19 col3" >no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow19_col4" class="data row19 col4" >209.23%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow19_col5" class="data row19 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow19_col6" class="data row19 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3_no_antibody_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow20_col0" class="data row20 col0" >lib3-v2-SIgN-10000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow20_col1" class="data row20 col1" >lib3-v2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow20_col2" class="data row20 col2" >SIgN-3C</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow20_col3" class="data row20 col3" >SIgN-10000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow20_col4" class="data row20 col4" >0.45%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow20_col5" class="data row20 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow20_col6" class="data row20 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3-v2_SIgN_10000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow21_col0" class="data row21 col0" >lib3-v2-SIgN-20000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow21_col1" class="data row21 col1" >lib3-v2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow21_col2" class="data row21 col2" >SIgN-3C</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow21_col3" class="data row21 col3" >SIgN-20000</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow21_col4" class="data row21 col4" >0.55%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow21_col5" class="data row21 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow21_col6" class="data row21 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3-v2_SIgN_20000_R1.fastq.gz</td>
            </tr>
            <tr>
                                <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow22_col0" class="data row22 col0" >lib3-v2-no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow22_col1" class="data row22 col1" >lib3-v2</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow22_col2" class="data row22 col2" >no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow22_col3" class="data row22 col3" >no-antibody</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow22_col4" class="data row22 col4" >209.23%</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow22_col5" class="data row22 col5" >220622-230109</td>
                        <td id="T_a32e663e_490a_11ee_be5e_ac1f6bca632crow22_col6" class="data row22 col6" >/fh/fast/bloom_j/computational_notebooks/ckikawa/2022/ZIKV_MAP_GooLab/data/concat_fastq/concat_fastq_out/lib3-v2_no_antibody_R1.fastq.gz</td>
            </tr>
    </tbody></table>




```python
# read in align specs
with open (subamplicon_alignspecs, 'r') as file:
    alignspecs = file.read().replace('\n','')
alignspecs
```




    '1,303,33,38 304,609,38,40 610,903,41,36 904,1200,41,37 1201,1512,36,35'




```python
# process deep sequencing data: run dms2_batch_subamp

# make folder where codon counts will go
codons_batch_subamp = os.path.join(resultsdir + 'codoncounts')
os.makedirs(codons_batch_subamp, exist_ok=True)

# tell me if sample path is missing
for c in samples['R1']:
    if c == 'NaN':
       raise ValueError('samples is missing a file path')

# create df and export csv with only name/R1 for dms2_batch_subamp   
batchfile = os.path.join(datadir + 'batchSubampIDs.csv')
samples[['name','R1']].to_csv(batchfile, index=False)

# standard trim in Bloom lab analyses is 200 bp
# from documentation: "removes low-quality nucleotides that tend to be at the end of long reads"
R1TRIM = 200
R2TRIM = 200

#needs a name to run
baka = ! dms2_batch_bcsubamp \
    --batchfile {batchfile} \
    --refseq {Erefseq} \
    --alignspecs {alignspecs} \
    --outdir {codons_batch_subamp} \
    --summaryprefix summary \
    --R1trim 200 \
    --R2trim 200 \
    --ncpus 16 \
    --use_existing {use_existing}

# ! executes following code in the command line

# add column to df.samples that includes filepath to each codoncounts.csv
samples['codoncounts'] = codons_batch_subamp + '/' + samples['name'] + '_codoncounts.csv'

# tell me how you did
print(f'dms2_batch_bcsubamp {dms_tools2.__version__} aligned sequencing data and created codon count files in {codons_batch_subamp}')
```

    dms2_batch_bcsubamp 2.6.10 aligned sequencing data and created codon count files in ./results/codoncounts


# Library sequencing quality control 
Now we'll examine some key QC aspects such as deep sequencing coverage, barcode sampling, subamplicon balancing, etc.


```python
summary = codons_batch_subamp + "/summary_"

showPDF([os.path.join(summary + "readstats.pdf"),
         os.path.join(summary + "bcstats.pdf")])
```


    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_13_0.png)
    


LEFT, ABOVE: reads per sample. Per CK019 bottleneck calculations, we required ~4.5e6 per sample. This is below target as most samples are around 3e6 reads. We knew that 2 MiSeq flowcells might not be sufficient... this may be an indication of that. We can simply submit more prep for sequencing if we decide to. 

RIGHT, ABOVE: barcodes per sample. Target was 1.5e6 barcodes per sample. Aligned barcodes misses this target, with all samples falling in the 5e5-1e6 range. 


```python
showPDF(os.path.join(summary + "readsperbc.pdf"))
```


    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_15_0.png)
    


ABOVE: reads per barcode. We need to read barcodes at least 2 times to error correct, but in many samples it looks like the majority of barcodes are being read once. This is likely due to insufficient sequencing depth.


```python
showPDF([os.path.join(summary + 'depth.pdf'),
         os.path.join(summary + 'mutfreq.pdf')])
```


    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_17_0.png)
    


LEFT, ABOVE: sequencing depth across 5 subamplicons. This looks quite uneven, especially compared to JBS083. Hopefully this resolves when we sequence more. 
RIGHT, ABOVE: mutation frequency across 5 subamplicons. Exciting to see probable selection taking place in some samples. 


```python
showPDF([os.path.join(summary + 'codonmuttypes.pdf'),
         os.path.join(summary + 'codonntchanges.pdf')])
```


    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_19_0.png)
    


LEFT, ABOVE: per-codon frequency of nonsynonymous, stop and synonymous mutations across samples. As expected, we see purging of stop codons in all my samples. 
RIGHT, ABOVE: per-codon frequency of 1-, 2-, or 3-nuc mutations across samples.


```python
showPDF([os.path.join(summary + 'singlentchanges.pdf'),
         os.path.join(summary + 'cumulmutcounts.pdf')])
```


    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_21_0.png)
    


LEFT, ABOVE: per-codon frequency of various mutation types to evaluate for oxidative damage.
RIGHT, ABOVE: fraction of mutations that occurs a given #(count) of times. The no-ab (no antibody) control contains more diversity than the antibody- or serum-selected conditions.

# Differential selection
We now will use differential selection to quantify antibody escape for each viral mutation against our panel of antibodies. 


```python
# Make base batch datafrane 
basebatch = (
    samples
    .query('selection != "no-antibody"')
    .assign(sel = lambda x: x['name'],
            group = lambda x: x['antibody'],
            mock = lambda x: x['library'] + '-no-antibody')
    .drop(columns = ['R1', 'name'])
    # .drop(none_index)
    # .drop(noAb_index)
    .assign(name = lambda x: x['library']+'-'+x['selection'])
    # .rename(columns = {'library' : 'name'})
    [['library','group', 'name', 'sel', 'mock', 'percent_infectivity']]
    .reset_index(drop = True)
    )

# make directory for individual concentrations diffsel results
diffseldir = os.path.join(resultsdir + 'diffsel')
os.makedirs(diffseldir, exist_ok = True)

# make batchfile
diffsel_batch = basebatch.rename(columns = {'group': 'antibody',
                                            # 'library': 'name',
                                            # 'name': 'group'
                                           })

# Make column labels and contents that are more interpretable and readable
diffsel_batch['grouplabel'] = (diffsel_batch['sel']
                                     .str.strip('lib1')
                                     .str.strip('lib2')
                                     .str.strip('lib3')
                                     .str.replace('-v2', '')
                                     .str.replace('-C10', 'broad antibody EDE1-C10: ')
                                     .str.replace('-C8', 'broad antibody EDE1-C8: ')
                                     .str.replace('-MZ4', 'pseudo-broad antibody MZ4: ')
                                     .str.replace('-SIgN', 'broad antibody SIgN-3C: ')
                                     .str.replace('-ZV67', 'targeted antibody ZV-67: ')
                                     .str.replace(' -', ' (')
                                     + ' ng/mL)'
                                    )
diffsel_batch['conc'] = (diffsel_batch['sel']
                         .str.split('-').str[-1]
                        ).astype(float)

diffsel_batch['group'] = (diffsel_batch['antibody'] + '-' + diffsel_batch['conc'].astype(int).astype(str))


# Sort based on custom list
sorter = ['SIgN-3C', 'EDE1-C8', 'EDE1-C10', 'MZ4', 'ZV-67']
# Custom lexicographical sort
sorterIndex = dict(zip(sorter, range(len(sorter))))
diffsel_batch['antibody_rank'] = diffsel_batch['antibody'].map(sorterIndex)
diffsel_batch.sort_values(['antibody_rank', 'name', 'conc'],
        ascending = [False, True, True], inplace = True)

# Other fiddly changes
diffsel_batch = (diffsel_batch     
                # Reorder columns
                 [[
                     'antibody', 
                   'group', 'name', 'sel', 'mock', 'percent_infectivity', 
                   'grouplabel'
                  ]]
                 .reset_index(drop=True)
                 # .rename(columns = {'group': 'name'})
                )
diffsel_batch['name'] = (diffsel_batch['name']

                                     .str.replace('lib3-v2', 'replicate3v2')
                                     .str.split('-').str[0]
                                     .str.replace('lib1', 'replicate1')
                                     .str.replace('lib2', 'replicate2')
                                     .str.replace('lib3', 'replicate3')

                                     # .str.replace('C10', 'EDE1-C10')
                                     # .str.replace('C8', 'EDE1-C8')
                                     # .str.replace('SIgN', 'SIgN-3C')
                                     # .str.replace('ZV67', 'ZV-67')
                               )




# # write individual concentration diffsel batchfile to csv
diffsel_batchfile = os.path.join(diffseldir + '/batch.csv')
diffsel_batch.to_csv(diffsel_batchfile, index = False)

diffsel_batch
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>antibody</th>
      <th>group</th>
      <th>name</th>
      <th>sel</th>
      <th>mock</th>
      <th>percent_infectivity</th>
      <th>grouplabel</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ZV-67</td>
      <td>ZV-67-40000</td>
      <td>replicate1</td>
      <td>lib1-ZV67-40000</td>
      <td>lib1-no-antibody</td>
      <td>0.84%</td>
      <td>targeted antibody ZV-67: (40000 ng/mL)</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ZV-67</td>
      <td>ZV-67-40000</td>
      <td>replicate2</td>
      <td>lib2-ZV67-40000</td>
      <td>lib2-no-antibody</td>
      <td>0.89%</td>
      <td>targeted antibody ZV-67: (40000 ng/mL)</td>
    </tr>
    <tr>
      <th>2</th>
      <td>ZV-67</td>
      <td>ZV-67-40000</td>
      <td>replicate3</td>
      <td>lib3-ZV67-40000</td>
      <td>lib3-no-antibody</td>
      <td>0.01%</td>
      <td>targeted antibody ZV-67: (40000 ng/mL)</td>
    </tr>
    <tr>
      <th>3</th>
      <td>MZ4</td>
      <td>MZ4-4800</td>
      <td>replicate1</td>
      <td>lib1-MZ4-4800</td>
      <td>lib1-no-antibody</td>
      <td>0.46%</td>
      <td>pseudo-broad antibody MZ4: (4800 ng/mL)</td>
    </tr>
    <tr>
      <th>4</th>
      <td>MZ4</td>
      <td>MZ4-4800</td>
      <td>replicate2</td>
      <td>lib2-MZ4-4800</td>
      <td>lib2-no-antibody</td>
      <td>0.02%</td>
      <td>pseudo-broad antibody MZ4: (4800 ng/mL)</td>
    </tr>
    <tr>
      <th>5</th>
      <td>MZ4</td>
      <td>MZ4-4800</td>
      <td>replicate3</td>
      <td>lib3-MZ4-4800</td>
      <td>lib3-no-antibody</td>
      <td>0.05%</td>
      <td>pseudo-broad antibody MZ4: (4800 ng/mL)</td>
    </tr>
    <tr>
      <th>6</th>
      <td>EDE1-C10</td>
      <td>EDE1-C10-300</td>
      <td>replicate1</td>
      <td>lib1-C10-300</td>
      <td>lib1-no-antibody</td>
      <td>0.11%</td>
      <td>broad antibody EDE1-C10: (300 ng/mL)</td>
    </tr>
    <tr>
      <th>7</th>
      <td>EDE1-C10</td>
      <td>EDE1-C10-300</td>
      <td>replicate2</td>
      <td>lib2-C10-300</td>
      <td>lib2-no-antibody</td>
      <td>0.23%</td>
      <td>broad antibody EDE1-C10: (300 ng/mL)</td>
    </tr>
    <tr>
      <th>8</th>
      <td>EDE1-C10</td>
      <td>EDE1-C10-300</td>
      <td>replicate3</td>
      <td>lib3-C10-300</td>
      <td>lib3-no-antibody</td>
      <td>0.39%</td>
      <td>broad antibody EDE1-C10: (300 ng/mL)</td>
    </tr>
    <tr>
      <th>9</th>
      <td>EDE1-C8</td>
      <td>EDE1-C8-1800</td>
      <td>replicate1</td>
      <td>lib1-C8-1800</td>
      <td>lib1-no-antibody</td>
      <td>0.01%</td>
      <td>broad antibody EDE1-C8: (1800 ng/mL)</td>
    </tr>
    <tr>
      <th>10</th>
      <td>EDE1-C8</td>
      <td>EDE1-C8-1800</td>
      <td>replicate2</td>
      <td>lib2-C8-1800</td>
      <td>lib2-no-antibody</td>
      <td>0.13%</td>
      <td>broad antibody EDE1-C8: (1800 ng/mL)</td>
    </tr>
    <tr>
      <th>11</th>
      <td>EDE1-C8</td>
      <td>EDE1-C8-1800</td>
      <td>replicate3</td>
      <td>lib3-C8-1800</td>
      <td>lib3-no-antibody</td>
      <td>0.23%</td>
      <td>broad antibody EDE1-C8: (1800 ng/mL)</td>
    </tr>
    <tr>
      <th>12</th>
      <td>SIgN-3C</td>
      <td>SIgN-3C-20000</td>
      <td>replicate1</td>
      <td>lib1-SIgN-20000</td>
      <td>lib1-no-antibody</td>
      <td>0.21%</td>
      <td>broad antibody SIgN-3C: (20000 ng/mL)</td>
    </tr>
    <tr>
      <th>13</th>
      <td>SIgN-3C</td>
      <td>SIgN-3C-10000</td>
      <td>replicate2</td>
      <td>lib2-SIgN-10000</td>
      <td>lib2-no-antibody</td>
      <td>0.35%</td>
      <td>broad antibody SIgN-3C: (10000 ng/mL)</td>
    </tr>
    <tr>
      <th>14</th>
      <td>SIgN-3C</td>
      <td>SIgN-3C-20000</td>
      <td>replicate2</td>
      <td>lib2-SIgN-20000</td>
      <td>lib2-no-antibody</td>
      <td>1.36%</td>
      <td>broad antibody SIgN-3C: (20000 ng/mL)</td>
    </tr>
    <tr>
      <th>15</th>
      <td>SIgN-3C</td>
      <td>SIgN-3C-10000</td>
      <td>replicate3</td>
      <td>lib3-SIgN-10000</td>
      <td>lib3-no-antibody</td>
      <td>1.19%</td>
      <td>broad antibody SIgN-3C: (10000 ng/mL)</td>
    </tr>
    <tr>
      <th>16</th>
      <td>SIgN-3C</td>
      <td>SIgN-3C-20000</td>
      <td>replicate3</td>
      <td>lib3-SIgN-20000</td>
      <td>lib3-no-antibody</td>
      <td>0.88%</td>
      <td>broad antibody SIgN-3C: (20000 ng/mL)</td>
    </tr>
    <tr>
      <th>17</th>
      <td>SIgN-3C</td>
      <td>SIgN-3C-10000</td>
      <td>replicate3v2</td>
      <td>lib3-v2-SIgN-10000</td>
      <td>lib3-v2-no-antibody</td>
      <td>0.45%</td>
      <td>broad antibody SIgN-3C: (10000 ng/mL)</td>
    </tr>
    <tr>
      <th>18</th>
      <td>SIgN-3C</td>
      <td>SIgN-3C-20000</td>
      <td>replicate3v2</td>
      <td>lib3-v2-SIgN-20000</td>
      <td>lib3-v2-no-antibody</td>
      <td>0.55%</td>
      <td>broad antibody SIgN-3C: (20000 ng/mL)</td>
    </tr>
  </tbody>
</table>
</div>




```python
# run dms2_batch_diffsel for individual concentration batchfile
log = ! dms2_batch_diffsel \
        --batchfile {diffsel_batchfile} \
        --summaryprefix summary \
        --indir {codons_batch_subamp} \
        --outdir {diffseldir} \
        --ncpus 16 \
        --use_existing {use_existing}
```


```python
# Define prefix for excess fraction surviving files to use
diffselprefix = os.path.join(diffseldir, 'summary_')

# Show OVERALL correlation plots for each of the above antibodies
for antibody in diffsel_batch['group'].unique():
    plot = [diffselprefix + antibody + '-positivesitediffselcorr.pdf']
    print(f"\nReplicate correlations for {antibody}")
    showPDF(plot)
```

    
    Replicate correlations for ZV-67-40000



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_26_1.png)
    


    
    Replicate correlations for MZ4-4800



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_26_3.png)
    


    
    Replicate correlations for EDE1-C10-300



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_26_5.png)
    


    
    Replicate correlations for EDE1-C8-1800



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_26_7.png)
    


    
    Replicate correlations for SIgN-3C-20000



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_26_9.png)
    


    
    Replicate correlations for SIgN-3C-10000



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_26_11.png)
    



```python
# Show lineplot from each antibody concentration
plot = diffselprefix + 'medianpositivediffsel.pdf'
showPDF(plot)
```


    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_27_0.png)
    



```python
# for each antibody group AND concentration in diffsel batchfile, run dms2_logoplot to show aa-level differential selection
for ab in diffsel_batch.query('group != "control-antibody"').group.unique():
    diffselfile = os.path.join(diffseldir, f'summary_{ab}-medianmutdiffsel.csv')
    baka = ! dms2_logoplot \
        --outdir {diffseldir} \
        --ncpus 16 \
        --name {ab} \
        --diffsel {diffselfile} \
        --restrictdiffsel positive \
        --nperline 101 \
        --numberevery 5 \
        --scalebar 2 "diffsel = 2" \
        --underlay yes \
        --overlay1 {diffselfile} wildtype wildtype \
        --use_existing 'no'

```


```python
# Also show median across concentrations (mostly for SIgN-3C)
diffselprefix = diffseldir + '/summary_'

# Initialize empty lists
medianfiles = []
medavgsitefiles = []
logoplots = []

# Define dictionary of sites of interest for zooms
zoom_dictionary = {}

# Iteration through antibodies; generate protein-wide logoplots
for antibody in diffsel_batch['antibody'].unique():
    print('\nGetting and plotting overall across-concentration median for {0}'.format(antibody))
    
    # list of files
    mediandiffsel_files = glob.glob('{0}*{1}-*medianmutdiffsel.csv'
            .format(diffselprefix, antibody))
    
    # Average across mutation fraction surviving
    medianmutdf = dms_tools2.diffsel.avgMutDiffSel(
            mediandiffsel_files, 'median')
    print('Here are the top 5 mutation-level effects...')
    print(medianmutdf.head(5))
    

    # Convert median mutation to median avg site fract            ion surviving
    avgsitedf = dms_tools2.diffsel.mutToSiteDiffSel(
            medianmutdf)
    print('Here are the top 10 site-level effects...')
    print(avgsitedf
          [['site', 'positive_diffsel']]
          .sort_values(by=['positive_diffsel'], ascending=False)
          .reset_index(drop=True)
          .head(10))
    
    # Write median mutation fracsurvive dataframe to csv files
    medianfile = os.path.join(diffseldir, 
            'antibody_{0}_median.csv'.format(antibody))
    medianfiles.append(medianfile)
    print("Writing across-concentration medians to {0}".format(
            medianfile))
    medianmutdf.to_csv(medianfile, index=False)
    
    # Write median average site fracsurvive dataframe to csv files
    avgsitefile = os.path.join(diffseldir,
            'antibody_{0}_median_avgsite.csv'.format(antibody))
    medavgsitefiles.append(avgsitefile)
    print("Writing across-concentration site medians to {0}".format(
            avgsitefile))
    avgsitedf.to_csv(avgsitefile, index=False)
    
    # Define antibody sites of interest and save to zoom dictionary
    top_sites = (avgsitedf.sort_values(by=['positive_diffsel'], ascending=False)['site'].head(2).to_list())
    top_muts = (medianmutdf.sort_values(by=['mutdiffsel'], ascending=False)['site'].head(2).to_list())
    
    unique_sites = list(set(top_sites + top_muts))
        
    zoom_dictionary[antibody] = unique_sites
    
    print(zoom_dictionary)
    
    # now make logo plot
    # scale bar unit is maximum effect
    scaleunit = '{0:.1g}'.format(medianmutdf['mutdiffsel'].max())
    scalelabel = '"diffsel = {0}"'.format(scaleunit)
    logoplot = os.path.join(diffseldir,
            '{0}_diffsel.pdf'.format(antibody))
    logoplots.append(logoplot)
    print("Creating logo plot {0} for {1} from {2}".format(
            logoplot, antibody, medianfile))
    log = !dms2_logoplot \
            --diffsel {medianfile} \
            --name {antibody} \
            --outdir {diffseldir} \
            --restrictdiffsel positive \
            --numberevery 5 \
            --nperline 101 \
            --colormap colorwheel \
            --underlay yes \
            --overlay1 {medianfile} wildtype wildtype \
            --scalebar {scaleunit} {scalelabel} \
            --use_existing {use_existing}

    showPDF(logoplot)
    
zoom_dictionary['EDE1-C10'] = sorted(list(set(zoom_dictionary['EDE1-C10'] + zoom_dictionary['EDE1-C8'])))
zoom_dictionary['EDE1-C8'] = sorted(list(set(zoom_dictionary['EDE1-C10'] + zoom_dictionary['EDE1-C8'])))

```

    
    Getting and plotting overall across-concentration median for ZV-67
    Here are the top 5 mutation-level effects...
       site wildtype mutation  mutdiffsel
    0   333        A        W    5.220291
    1   333        A        E    5.128572
    2   333        A        L    5.077331
    3   333        A        Q    5.061331
    4   310        A        D    4.846885
    Here are the top 10 site-level effects...
    mutation  site  positive_diffsel
    0          333         42.666457
    1          310         22.563684
    2          311         12.796853
    3          371          6.305352
    4          350          4.983127
    5          394          4.153731
    6          355          3.679925
    7          325          3.522249
    8          331          3.238815
    9          150          3.133429
    Writing across-concentration medians to ./results/diffsel/antibody_ZV-67_median.csv
    Writing across-concentration site medians to ./results/diffsel/antibody_ZV-67_median_avgsite.csv
    {'ZV-67': [333, 310]}
    Creating logo plot ./results/diffsel/ZV-67_diffsel.pdf for ZV-67 from ./results/diffsel/antibody_ZV-67_median.csv



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_29_1.png)
    


    
    Getting and plotting overall across-concentration median for MZ4
    Here are the top 5 mutation-level effects...
       site wildtype mutation  mutdiffsel
    0   301        K        H    4.217704
    1   479        T        G    4.011507
    2   301        K        Y    3.940772
    3   368        S        L    3.929273
    4   368        S        V    3.823385
    Here are the top 10 site-level effects...
    mutation  site  positive_diffsel
    0          368         30.324081
    1          301         17.964025
    2          352          4.915727
    3          353          4.324830
    4          479          4.242705
    5          348          3.985321
    6          355          3.936266
    7          320          3.908599
    8          337          3.652867
    9          306          3.374493
    Writing across-concentration medians to ./results/diffsel/antibody_MZ4_median.csv
    Writing across-concentration site medians to ./results/diffsel/antibody_MZ4_median_avgsite.csv
    {'ZV-67': [333, 310], 'MZ4': [368, 301, 479]}
    Creating logo plot ./results/diffsel/MZ4_diffsel.pdf for MZ4 from ./results/diffsel/antibody_MZ4_median.csv



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_29_3.png)
    


    
    Getting and plotting overall across-concentration median for EDE1-C10
    Here are the top 5 mutation-level effects...
       site wildtype mutation  mutdiffsel
    0    33        V        I    4.582371
    1    49        T        D    2.497055
    2   231        T        R    2.428671
    3   464        S        L    2.335385
    4   439        N        K    1.837236
    Here are the top 10 site-level effects...
    mutation  site  positive_diffsel
    0          283          5.435655
    1          315          5.035614
    2           33          4.649977
    3          439          3.173396
    4           49          2.844304
    5          231          2.629227
    6          464          2.457810
    7           12          1.954235
    8          254          1.838681
    9          345          1.761304
    Writing across-concentration medians to ./results/diffsel/antibody_EDE1-C10_median.csv
    Writing across-concentration site medians to ./results/diffsel/antibody_EDE1-C10_median_avgsite.csv
    {'ZV-67': [333, 310], 'MZ4': [368, 301, 479], 'EDE1-C10': [315, 283, 49, 33]}
    Creating logo plot ./results/diffsel/EDE1-C10_diffsel.pdf for EDE1-C10 from ./results/diffsel/antibody_EDE1-C10_median.csv



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_29_5.png)
    


    
    Getting and plotting overall across-concentration median for EDE1-C8
    Here are the top 5 mutation-level effects...
       site wildtype mutation  mutdiffsel
    0   397        T        G    3.517528
    1   439        N        K    2.968811
    2   349        M        E    2.773152
    3    55        E        V    2.677955
    4     7        S        P    2.620159
    Here are the top 10 site-level effects...
    mutation  site  positive_diffsel
    0          315          7.329833
    1            7          5.568574
    2           55          4.304360
    3          349          3.546344
    4          397          3.523397
    5          373          3.074773
    6          283          3.049768
    7          439          2.985281
    8           37          2.879558
    9           21          2.640162
    Writing across-concentration medians to ./results/diffsel/antibody_EDE1-C8_median.csv
    Writing across-concentration site medians to ./results/diffsel/antibody_EDE1-C8_median_avgsite.csv
    {'ZV-67': [333, 310], 'MZ4': [368, 301, 479], 'EDE1-C10': [315, 283, 49, 33], 'EDE1-C8': [439, 315, 397, 7]}
    Creating logo plot ./results/diffsel/EDE1-C8_diffsel.pdf for EDE1-C8 from ./results/diffsel/antibody_EDE1-C8_median.csv



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_29_7.png)
    


    
    Getting and plotting overall across-concentration median for SIgN-3C
    Here are the top 5 mutation-level effects...
       site wildtype mutation  mutdiffsel
    0     8        N        T    2.433386
    1    72        S        G    1.871529
    2    29        G        Q    1.772030
    3    49        T        D    1.742743
    4   320        E        S    1.566200
    Here are the top 10 site-level effects...
    mutation  site  positive_diffsel
    0          140          6.606207
    1          315          5.459780
    2          162          4.396402
    3          320          4.185990
    4           29          3.860917
    5          157          3.511810
    6            8          3.168485
    7           46          3.024122
    8            7          2.938824
    9           72          2.880645
    Writing across-concentration medians to ./results/diffsel/antibody_SIgN-3C_median.csv
    Writing across-concentration site medians to ./results/diffsel/antibody_SIgN-3C_median_avgsite.csv
    {'ZV-67': [333, 310], 'MZ4': [368, 301, 479], 'EDE1-C10': [315, 283, 49, 33], 'EDE1-C8': [439, 315, 397, 7], 'SIgN-3C': [8, 72, 315, 140]}
    Creating logo plot ./results/diffsel/SIgN-3C_diffsel.pdf for SIgN-3C from ./results/diffsel/antibody_SIgN-3C_median.csv



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_29_9.png)
    


# Figures for the paper
Next we will produce various plots for the paper. Namely: logoplots of functional epitopes; mutational tolerance analysis; tables of relevant data to include in the supplement. 

## Zoomed-in logoplots of regions of interest
These zoomed-in logoplots (composed of sites with largest site- and mutation-level effects) will be used to compare to structural epitopes. We will scale per-antibody.


```python
# make directory for figures
figsdir = os.path.join(resultsdir, 'figures')
os.makedirs(figsdir, exist_ok = True)
```


```python
zoom_dictionary
```




    {'ZV-67': [333, 310],
     'MZ4': [368, 301, 479],
     'EDE1-C10': [7, 33, 49, 283, 315, 397, 439],
     'EDE1-C8': [7, 33, 49, 283, 315, 397, 439],
     'SIgN-3C': [8, 72, 315, 140]}




```python
plt.clf()
# Set plot style
sns.set(style = 'white') 

# Make logoplots for each antibody
for ab in zoom_dictionary.keys():
    
    # Access sites from zoomsites dictionary
    sites = zoom_dictionary[ab]
    
    # Specific to MZ4 - rather than show one of the sites, show 182
    sites = [182 if x == 479 else x for x in sites]
    
    # Specific to EDE1 - skip 439 and 397 to save space
    # sites = [ele for ele in sites if ele != 439]
    sites = [ele for ele in sites if ele != 397]



    
    above_sites = [x + 1 for x in sites]
    below_sites = [x - 1 for x in sites]
    
    sites = sites + above_sites + below_sites

    # Get dataframe of antibody selections diffsel data
    df = pd.read_csv(diffseldir + f"/antibody_{ab}_median.csv")
    
    # Zoom in on sites identified in dictionary 
    zoomed_df = (df
                 .loc[df['site'].isin(sites)]
                )

    # Add a site label that includes the wildtype amino acid
    zoomed_df = zoomed_df.assign(site_label = lambda x: x['wildtype'] + x['site'].astype(str))
    
    # Floor selection at 0
    zoomed_df['mutdiffsel'] = zoomed_df['mutdiffsel'].mask(zoomed_df['mutdiffsel'].lt(0),0)
    
    # Draw and save the plots
    fig, axes = dmslogo.draw_logo(
        data=zoomed_df,
        x_col='site',
        letter_col='mutation',
        letter_height_col='mutdiffsel',
        addbreaks=True,
        xtick_col='site_label',
        ylabel = 'antibody escape',
        xlabel = '',
        # colorscheme = colorwheel
    )
    
    # Save to figures directory
    antibodyzoom = os.path.join(figsdir, f'{ab}_zoom.pdf')
    fig.savefig(antibodyzoom, bbox_inches='tight')
    
    print(f'showing zoomed logoplot for {ab}...')
    showPDF(antibodyzoom)
```

    showing zoomed logoplot for ZV-67...



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_33_1.png)
    


    showing zoomed logoplot for MZ4...



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_33_3.png)
    


    showing zoomed logoplot for EDE1-C10...



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_33_5.png)
    


    showing zoomed logoplot for EDE1-C8...



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_33_7.png)
    


    showing zoomed logoplot for SIgN-3C...



    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_33_9.png)
    


## Create polyclonal analysis input files 
To visualize results on 3D protein structures, we will use the b-factor reassignment function defined in polyclonal. We will use site-level summed differential selection (our metric for 'antibody escape') as b-factors. 


```python
# make reassigned bfactor results subfolder
reassignedpdbdir = os.path.join(resultsdir, 'reassignedpdb/')
os.makedirs(reassignedpdbdir, exist_ok=True)
```


```python
# Get dataframe of input files for b-factor assignment 
pointers = (diffsel_batch
                            .assign(sitediffsel = lambda x: (diffseldir + '/antibody_' + 
                                               x['antibody'] + '_median_avgsite.csv')
                                   )
                            [['antibody','sitediffsel']]
                            .drop_duplicates()
                            .reset_index(drop=True)
                           )
                
                        
pointers
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>antibody</th>
      <th>sitediffsel</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ZV-67</td>
      <td>./results/diffsel/antibody_ZV-67_median_avgsit...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>MZ4</td>
      <td>./results/diffsel/antibody_MZ4_median_avgsite.csv</td>
    </tr>
    <tr>
      <th>2</th>
      <td>EDE1-C10</td>
      <td>./results/diffsel/antibody_EDE1-C10_median_avg...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>EDE1-C8</td>
      <td>./results/diffsel/antibody_EDE1-C8_median_avgs...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>SIgN-3C</td>
      <td>./results/diffsel/antibody_SIgN-3C_median_avgs...</td>
    </tr>
  </tbody>
</table>
</div>




```python
# make base dataframe to add protein labels to
# these should contain sites, wt, mutation, condition 
# as well as the PDB file protein_chain and protein_site
# these will include mut_* values since we are interested in mapping escape mutants

# make dataframe for all conditions/concentrations with site, wt, mut and mutdiffsel values
pdbdf = pd.DataFrame()

for index, info in pointers.iterrows():
    ab = info.str.split('/t')[0][0]
    sitediffselcsv = info.str.split('/t')[1][0]
    
    pdbdf = (pdbdf.append(pd.read_csv(sitediffselcsv)
                          [['site','positive_diffsel']]
                          .sort_values(by = ['site'])
                          .reset_index(drop=True)
                          .assign(condition = ab,
                                 )
                         )
            )
          

print('Here is the diffsel dataframe we will add protein labels to...')
pdbdf
```

    Here is the diffsel dataframe we will add protein labels to...





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>site</th>
      <th>positive_diffsel</th>
      <th>condition</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>0.000000</td>
      <td>ZV-67</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>1.518964</td>
      <td>ZV-67</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>1.454714</td>
      <td>ZV-67</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>1.067576</td>
      <td>ZV-67</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>0.223613</td>
      <td>ZV-67</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>499</th>
      <td>500</td>
      <td>0.000000</td>
      <td>SIgN-3C</td>
    </tr>
    <tr>
      <th>500</th>
      <td>501</td>
      <td>0.000000</td>
      <td>SIgN-3C</td>
    </tr>
    <tr>
      <th>501</th>
      <td>502</td>
      <td>0.299438</td>
      <td>SIgN-3C</td>
    </tr>
    <tr>
      <th>502</th>
      <td>503</td>
      <td>0.592304</td>
      <td>SIgN-3C</td>
    </tr>
    <tr>
      <th>503</th>
      <td>504</td>
      <td>0.185849</td>
      <td>SIgN-3C</td>
    </tr>
  </tbody>
</table>
<p>2520 rows × 3 columns</p>
</div>




```python
# homodimer polyclonal analysis input file
homodimer_polyclonal_df = pd.DataFrame() 

# in 5ire protein structure, the three chains of E homodimer (+1 monomer) are labeled 'A', 'C', and 'E' 
# here, iteratively add the label 'A', 'C', or 'E' to each site
# the product dataframe will basically be three times the size
chainlist = ['A', 'C', 'E']
templist = []

ab_dict = {'EDE1-C10' : 'EDE1_C10',
           'EDE1-C8' : 'EDE1_C8',
           'ZKA-64' : 'ZKA_64',
           'ZV-67' : 'ZV_67',
           'SiGN-3C' : 'SiGN_3C'}
    
for chain in chainlist:
    templist.append(pdbdf
                    .assign(protein_chain = chain,
                            protein_site = lambda x: x['site'],
                            label_site = lambda x: x['site'])
                   )

homodimer_polyclonal_df = (homodimer_polyclonal_df
                           .append(templist)
                           .replace(ab_dict)
                          )
                           

# save these files and preview dataframe
medianfracsurvivecsv = os.path.join(reassignedpdbdir, 'alldiffsel.csv')
homodimer_polyclonal_df.to_csv(medianfracsurvivecsv, index=False, float_format='%.3g')

print(f"Writing CSV to {medianfracsurvivecsv}; here is a preview...")
homodimer_polyclonal_df
```

    Writing CSV to ./results/reassignedpdb/alldiffsel.csv; here is a preview...





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>site</th>
      <th>positive_diffsel</th>
      <th>condition</th>
      <th>protein_chain</th>
      <th>protein_site</th>
      <th>label_site</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>0.000000</td>
      <td>ZV_67</td>
      <td>A</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>1.518964</td>
      <td>ZV_67</td>
      <td>A</td>
      <td>2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>1.454714</td>
      <td>ZV_67</td>
      <td>A</td>
      <td>3</td>
      <td>3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>1.067576</td>
      <td>ZV_67</td>
      <td>A</td>
      <td>4</td>
      <td>4</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>0.223613</td>
      <td>ZV_67</td>
      <td>A</td>
      <td>5</td>
      <td>5</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>499</th>
      <td>500</td>
      <td>0.000000</td>
      <td>SIgN-3C</td>
      <td>E</td>
      <td>500</td>
      <td>500</td>
    </tr>
    <tr>
      <th>500</th>
      <td>501</td>
      <td>0.000000</td>
      <td>SIgN-3C</td>
      <td>E</td>
      <td>501</td>
      <td>501</td>
    </tr>
    <tr>
      <th>501</th>
      <td>502</td>
      <td>0.299438</td>
      <td>SIgN-3C</td>
      <td>E</td>
      <td>502</td>
      <td>502</td>
    </tr>
    <tr>
      <th>502</th>
      <td>503</td>
      <td>0.592304</td>
      <td>SIgN-3C</td>
      <td>E</td>
      <td>503</td>
      <td>503</td>
    </tr>
    <tr>
      <th>503</th>
      <td>504</td>
      <td>0.185849</td>
      <td>SIgN-3C</td>
      <td>E</td>
      <td>504</td>
      <td>504</td>
    </tr>
  </tbody>
</table>
<p>7560 rows × 6 columns</p>
</div>




```python
ablist = ['EDE1_C10',
          'EDE1_C8',
          'MZ4',
          'ZV_67', 
          'SiGN_3C']

for ab in ablist:
    df = homodimer_polyclonal_df.query('condition == "' + str(ab) + '"')
    csv = os.path.join(reassignedpdbdir, ab + '_diffsel.csv')
    df.to_csv(csv, index = False)
    print(f'Writing CSV for {ab} to {csv}...')
```

    Writing CSV for EDE1_C10 to ./results/reassignedpdb/EDE1_C10_diffsel.csv...
    Writing CSV for EDE1_C8 to ./results/reassignedpdb/EDE1_C8_diffsel.csv...
    Writing CSV for MZ4 to ./results/reassignedpdb/MZ4_diffsel.csv...
    Writing CSV for ZV_67 to ./results/reassignedpdb/ZV_67_diffsel.csv...
    Writing CSV for SiGN_3C to ./results/reassignedpdb/SiGN_3C_diffsel.csv...


# Functional analysis
To ask and answer questions about how different mutational tolerance across E protein (at both the individual amino acid mutation-level and at the site-level) might affect antibody escape, we will generate per-antibody scatter plots comparing different mutational tolerance/flexibility metrics to antibody escape. 

## Site-level flexibility versus site-level antibody escape
Here, we can use functional score estimates from the work performed in Sourisseau et al to compute site-level flexibility. We can compare these metrics to site-level antibody escape to see whether narrow antibodies target more flexible sites (relative to broad antibodies)


```python
# Define colors
colors = {'teal': '#3bceac',
          'red': '#DA2C38',
          'blue': '#446df6',
          'yellow': '#F6AE2D',
          'green': '#00635D',
         }
```


```python
# ID functional data directory
functionaldatadir = os.path.join(datadir, 'functional_data')

# Read in mutation-level preferences data
unscaledprefsfile = os.path.join(functionaldatadir, 'unscaled_prefs.csv')
unscaledprefs = pd.read_csv(unscaledprefsfile)
unscaledprefs.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>site</th>
      <th>A</th>
      <th>C</th>
      <th>D</th>
      <th>E</th>
      <th>F</th>
      <th>G</th>
      <th>H</th>
      <th>I</th>
      <th>K</th>
      <th>...</th>
      <th>M</th>
      <th>N</th>
      <th>P</th>
      <th>Q</th>
      <th>R</th>
      <th>S</th>
      <th>T</th>
      <th>V</th>
      <th>W</th>
      <th>Y</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>0.005438</td>
      <td>0.009065</td>
      <td>0.012163</td>
      <td>0.009180</td>
      <td>0.008687</td>
      <td>0.004282</td>
      <td>0.003339</td>
      <td>0.179643</td>
      <td>0.003888</td>
      <td>...</td>
      <td>0.279248</td>
      <td>0.008542</td>
      <td>0.001240</td>
      <td>0.011589</td>
      <td>0.001585</td>
      <td>0.003096</td>
      <td>0.316543</td>
      <td>0.112591</td>
      <td>0.016822</td>
      <td>0.008731</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>0.007809</td>
      <td>0.018669</td>
      <td>0.021458</td>
      <td>0.018187</td>
      <td>0.018150</td>
      <td>0.062636</td>
      <td>0.026099</td>
      <td>0.015878</td>
      <td>0.026583</td>
      <td>...</td>
      <td>0.021237</td>
      <td>0.022353</td>
      <td>0.006953</td>
      <td>0.019685</td>
      <td>0.579880</td>
      <td>0.023507</td>
      <td>0.012992</td>
      <td>0.012776</td>
      <td>0.052111</td>
      <td>0.024828</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>0.006306</td>
      <td>0.715622</td>
      <td>0.031189</td>
      <td>0.018934</td>
      <td>0.019171</td>
      <td>0.010131</td>
      <td>0.011007</td>
      <td>0.010956</td>
      <td>0.012948</td>
      <td>...</td>
      <td>0.019899</td>
      <td>0.012324</td>
      <td>0.004083</td>
      <td>0.013883</td>
      <td>0.010874</td>
      <td>0.016100</td>
      <td>0.007397</td>
      <td>0.006839</td>
      <td>0.015366</td>
      <td>0.052559</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>0.012582</td>
      <td>0.021175</td>
      <td>0.026134</td>
      <td>0.025855</td>
      <td>0.021090</td>
      <td>0.014896</td>
      <td>0.013630</td>
      <td>0.436679</td>
      <td>0.013627</td>
      <td>...</td>
      <td>0.074051</td>
      <td>0.043320</td>
      <td>0.006742</td>
      <td>0.020199</td>
      <td>0.005229</td>
      <td>0.014981</td>
      <td>0.027572</td>
      <td>0.153553</td>
      <td>0.033425</td>
      <td>0.018557</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>0.009296</td>
      <td>0.011905</td>
      <td>0.017206</td>
      <td>0.036911</td>
      <td>0.012321</td>
      <td>0.766760</td>
      <td>0.009187</td>
      <td>0.006683</td>
      <td>0.010623</td>
      <td>...</td>
      <td>0.026395</td>
      <td>0.010510</td>
      <td>0.004116</td>
      <td>0.012365</td>
      <td>0.012495</td>
      <td>0.006081</td>
      <td>0.004521</td>
      <td>0.013316</td>
      <td>0.014955</td>
      <td>0.011066</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 21 columns</p>
</div>




```python
# Make logoplots for same zoomsites as above

# Make base dataframe
tidy_unscaledprefs = (unscaledprefs
                        .melt(id_vars=['site'],
                              var_name='aa',
                              value_name='unscaledpref'
                             )        
                )

# Iterate through zoomsites and make plots

for ab in zoom_dictionary.keys():
    # Access sites from zoomsites dictionary
    sites = zoom_dictionary[ab]

    # Zoom in on sites identified in dictionary 
    zoomed_unscaledprefs = (tidy_unscaledprefs
                 .loc[tidy_unscaledprefs['site'].isin(sites)]
                )

    # Draw and save the plots
    fig, axes = dmslogo.draw_logo(
        data=zoomed_unscaledprefs,
        x_col='site',
        ylabel = 'mut tolerance',
        letter_col='aa',
        letter_height_col='unscaledpref',
        # colorscheme = colorwheel
    )
    
    # Save to figures directory
    prefsantibodyzoom = os.path.join(figsdir, f'prefs_{ab}_zoom.pdf')
    fig.savefig(prefsantibodyzoom, bbox_inches='tight', dpi = 'figure')
    print(f'saving zoomed preferences for {ab} at selected zoomsites...')
    # showPDF(prefsantibodyzoom)
```

    saving zoomed preferences for ZV-67 at selected zoomsites...
    saving zoomed preferences for MZ4 at selected zoomsites...
    saving zoomed preferences for EDE1-C10 at selected zoomsites...
    saving zoomed preferences for EDE1-C8 at selected zoomsites...
    saving zoomed preferences for SIgN-3C at selected zoomsites...



```python
# We can also look at the site-level
# For this, the 2019 Sourisseau paper has 2 metrics, neffective and site entropy
muteffectsfile = os.path.join(functionaldatadir, 'struct_props_mut_tol.csv')
muteffects = pd.read_csv(muteffectsfile)[['site','mutational_tolerance_measure','mutational_tolerance']].drop_duplicates()
muteffects.query('site == 1')
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>site</th>
      <th>mutational_tolerance_measure</th>
      <th>mutational_tolerance</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>entropy</td>
      <td>1.809948</td>
    </tr>
    <tr>
      <th>2004</th>
      <td>1</td>
      <td>neffective</td>
      <td>6.110127</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Separate mutational tolerance methods
muteffects_reform = (muteffects
                      .query('mutational_tolerance_measure == "entropy"')
                      .rename(columns = {'mutational_tolerance': 'entropy'})
                      .merge(muteffects.query('mutational_tolerance_measure == "neffective"')
                            .rename(columns = {'mutational_tolerance': 'neffective'}), 
                            on = ['site'])
                     [['site','entropy','neffective']]
                     )
muteffects_reform.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>site</th>
      <th>entropy</th>
      <th>neffective</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>1.809948</td>
      <td>6.110127</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>1.851241</td>
      <td>6.367718</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>1.375274</td>
      <td>3.956159</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>2.123703</td>
      <td>8.362047</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>1.190491</td>
      <td>3.288694</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Reformat dataframe
diffseldir = os.path.join(resultsdir, 'diffsel')

ab_list = diffsel_batch['antibody'].unique()


# Initialize empty dataframe for differential selection data of interest
diffsel_df = pd.DataFrame()

for ab in ab_list:
    file = glob.glob(os.path.join(diffseldir, 'antibody_' + str(ab) + "_median_avgsite.csv"))[0]
    print(f'For antibody {ab}, found median site diffsel at... {file}')
    
    # Get scratch dataframe from each antibody
    scratch_df = (pd.read_csv(file)
                  .assign(antibody = ab)
                 )  

    # Concatenate scratch and existing dataframe
    diffsel_df = pd.concat([diffsel_df, scratch_df])
    
# Rename IC50 dataframe to match variables in diffsel dataframe
df = (muteffects_reform
           .merge(diffsel_df, how = 'left', on = ['site']) # Then merge with diffsel data
          )
df.head()
```

    For antibody ZV-67, found median site diffsel at... ./results/diffsel/antibody_ZV-67_median_avgsite.csv
    For antibody MZ4, found median site diffsel at... ./results/diffsel/antibody_MZ4_median_avgsite.csv
    For antibody EDE1-C10, found median site diffsel at... ./results/diffsel/antibody_EDE1-C10_median_avgsite.csv
    For antibody EDE1-C8, found median site diffsel at... ./results/diffsel/antibody_EDE1-C8_median_avgsite.csv
    For antibody SIgN-3C, found median site diffsel at... ./results/diffsel/antibody_SIgN-3C_median_avgsite.csv





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>site</th>
      <th>entropy</th>
      <th>neffective</th>
      <th>abs_diffsel</th>
      <th>positive_diffsel</th>
      <th>negative_diffsel</th>
      <th>max_diffsel</th>
      <th>min_diffsel</th>
      <th>antibody</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>1.809948</td>
      <td>6.110127</td>
      <td>7.581099</td>
      <td>0.000000</td>
      <td>-7.581099</td>
      <td>0.000000</td>
      <td>-2.621183</td>
      <td>ZV-67</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>1.809948</td>
      <td>6.110127</td>
      <td>7.239773</td>
      <td>0.000000</td>
      <td>-7.239773</td>
      <td>0.000000</td>
      <td>-2.508192</td>
      <td>MZ4</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>1.809948</td>
      <td>6.110127</td>
      <td>6.473181</td>
      <td>0.000000</td>
      <td>-6.473181</td>
      <td>0.000000</td>
      <td>-2.140399</td>
      <td>EDE1-C10</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>1.809948</td>
      <td>6.110127</td>
      <td>5.938759</td>
      <td>0.001424</td>
      <td>-5.937335</td>
      <td>0.001424</td>
      <td>-2.731474</td>
      <td>EDE1-C8</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>1.809948</td>
      <td>6.110127</td>
      <td>2.841469</td>
      <td>0.131642</td>
      <td>-2.709827</td>
      <td>0.131573</td>
      <td>-0.978957</td>
      <td>SIgN-3C</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Using entropy, facet by antibody so we can see individal plots more easily
plt.clf()
# Set plot style
sns.set(style = 'white') 
# Set plot size
plt.figure(figsize=(6,6)) #width=8, height=4
# Set plot style
sns.set_style('whitegrid')

# Make dataframe
abdata = df.loc[df['antibody'].isin(['EDE1-C10', 'EDE1-C8', 'SIgN-3C', 'MZ4', 'ZV-67'])]

g = sns.FacetGrid(data = abdata, col = 'antibody', hue = 'antibody', palette = list(colors.values()),
                 col_wrap = 3, margin_titles=True)


g.map(sns.scatterplot, 
      'entropy', 
      'positive_diffsel',
     )

# Change titles and axes
g.set_ylabels('site antibody escape')
g.set_xlabels('site entropy')
g.set_titles(col_template="{col_name}")


# Save to figures directory
fig = os.path.join(figsdir,'site_entropy_vs_selection_faceted.pdf')
g.savefig(fig, dpi = 'figure')
showPDF(fig)
```


    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_47_0.png)
    



```python
# Repeat using neffective, facet by antibody so we can see individal plots more easily

plt.clf()
# Set plot style
sns.set(style = 'white') 
# Set plot size
plt.figure(figsize=(6,6)) #width=8, height=4
# Set plot style
sns.set_style('whitegrid')

g = sns.FacetGrid(data = abdata, col = 'antibody', hue = 'antibody', palette = list(colors.values()),
                 col_wrap = 3, margin_titles=True)
g.map(sns.scatterplot, 
      'neffective', 
      'positive_diffsel',
     )

# Change titles and axes
g.set_titles(col_template="{col_name}")
g.set_ylabels('site antibody escape')
g.set_xlabels('site neffective')


# Save to figures directory
fig = os.path.join(figsdir,'site_neffective_vs_selection_faceted.pdf')
g.savefig(fig, dpi = 'figure')
showPDF(fig)
```


    
![png](./results/notebooks/selections_analysis_files/./results/notebooks/selections_analysis_48_0.png)
    


## Mutation-level amino acid tolerance compared to mutation-level antibody escape 
Similar to the above analysis, we can see whether particular antibody escape mutations targeted by narrow antibodies are more preferred relative to those targeted by broad antibodies. 


```python
# First retrieve amino acid level diffsel
# Initialize empty dataframe for differential selection data of interest

mutdiffsel_df = pd.DataFrame()

for ab in ab_list:
    file = glob.glob(os.path.join(diffseldir, 'antibody_' + str(ab) + "_median.csv"))[0]
    print(f'For antibody {ab}, found median mutational diffsel at... {file}')
    
    # Get scratch dataframe from each antibody
    scratch_df = (pd.read_csv(file)
                  .assign(antibody = ab)
                 )  

    # Concatenate scratch and existing dataframe
    mutdiffsel_df = pd.concat([mutdiffsel_df, scratch_df])
    
mutdiffsel_df = mutdiffsel_df.assign(isite = lambda x: x.site.astype(str) + x.mutation)
```

    For antibody ZV-67, found median mutational diffsel at... ./results/diffsel/antibody_ZV-67_median.csv
    For antibody MZ4, found median mutational diffsel at... ./results/diffsel/antibody_MZ4_median.csv
    For antibody EDE1-C10, found median mutational diffsel at... ./results/diffsel/antibody_EDE1-C10_median.csv
    For antibody EDE1-C8, found median mutational diffsel at... ./results/diffsel/antibody_EDE1-C8_median.csv
    For antibody SIgN-3C, found median mutational diffsel at... ./results/diffsel/antibody_SIgN-3C_median.csv



```python
# We already have a dataframe of unmerged prefs
tidy_unscaledprefs = tidy_unscaledprefs.assign(isite = lambda x: x.site.astype(str) + x.aa)

# Rename IC50 dataframe to match variables in diffsel dataframe and merge
tidy_mutdiffsel_df = (tidy_unscaledprefs
                      .merge(mutdiffsel_df, how = 'left', on = ['isite', 'site']) # Then merge with diffsel data
                     )

tidy_mutdiffsel_df = (tidy_mutdiffsel_df
                      .assign(mutant = lambda x: x.wildtype + x.isite)
                      [['site', 'wildtype', 'mutation', 'mutant', 'unscaledpref', 'mutdiffsel', 'antibody']]
                     )

# We should also scale prefs by the WT preference
wtprefs_df = (tidy_mutdiffsel_df
              .query('wildtype == mutation')
              .sort_values(by = ['site', 'wildtype', 'mutation'])
              [['site', 'unscaledpref']]
              .drop_duplicates()
              .rename(columns = {'unscaledpref': 'wtpref'})
             )

prefs_vs_mutdiffsel_df = (wtprefs_df
           .merge(tidy_mutdiffsel_df)
          )

prefs_vs_mutdiffsel_df = prefs_vs_mutdiffsel_df.assign(scaledpref = lambda x: x.unscaledpref / x.wtpref)

# Drop WT sites
prefs_vs_mutdiffsel_df = prefs_vs_mutdiffsel_df.dropna()

# Add log-transformed scaled preferences
prefs_vs_mutdiffsel_df['logscaledpref'] = np.log2(prefs_vs_mutdiffsel_df['scaledpref'])
prefs_vs_mutdiffsel_df.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>site</th>
      <th>wtpref</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>mutant</th>
      <th>unscaledpref</th>
      <th>mutdiffsel</th>
      <th>antibody</th>
      <th>scaledpref</th>
      <th>logscaledpref</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>0.179643</td>
      <td>I</td>
      <td>A</td>
      <td>I1A</td>
      <td>0.005438</td>
      <td>-0.000814</td>
      <td>ZV-67</td>
      <td>0.030269</td>
      <td>-5.04603</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>0.179643</td>
      <td>I</td>
      <td>A</td>
      <td>I1A</td>
      <td>0.005438</td>
      <td>-0.000861</td>
      <td>MZ4</td>
      <td>0.030269</td>
      <td>-5.04603</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1</td>
      <td>0.179643</td>
      <td>I</td>
      <td>A</td>
      <td>I1A</td>
      <td>0.005438</td>
      <td>-0.000740</td>
      <td>EDE1-C10</td>
      <td>0.030269</td>
      <td>-5.04603</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1</td>
      <td>0.179643</td>
      <td>I</td>
      <td>A</td>
      <td>I1A</td>
      <td>0.005438</td>
      <td>-0.000542</td>
      <td>EDE1-C8</td>
      <td>0.030269</td>
      <td>-5.04603</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1</td>
      <td>0.179643</td>
      <td>I</td>
      <td>A</td>
      <td>I1A</td>
      <td>0.005438</td>
      <td>-0.000001</td>
      <td>SIgN-3C</td>
      <td>0.030269</td>
      <td>-5.04603</td>
    </tr>
  </tbody>
</table>
</div>




```python
# Plot scatters faceted by antibody

plt.clf()
# Set plot style
sns.set(style = 'white') 
# Set plot size
plt.figure(figsize=(6,6)) #width=8, height=4
# Set plot style
sns.set_style('whitegrid')

# Make dataframe
# abdata = df.loc[df['antibody'].isin(['EDE1-C10', 'EDE1-C8', 'SIgN-3C', 'MZ4', 'ZV-67'])]

g = sns.FacetGrid(data = prefs_vs_mutdiffsel_df, col = 'antibody', hue = 'antibody', palette = list(colors.values()),
                 col_wrap = 3)
g.map(sns.scatterplot, 
      'unscaledpref', 
      'mutdiffsel',
     )

# Change titles and axes
g.set_titles(col_template="{col_name}")
g.set_ylabels('AA antibody escape')
g.set_xlabels('unscaled AA preference')


# Save to figures directory
fig = os.path.join(figsdir,'mutationlevel_pref_vs_diffsel.pdf')
g.savefig(fig, dpi = 'figure')
showPDF(fig)
```


```python
# We should also scale prefs by the WT preference
# Plot scatters faceted by antibody
plt.clf()
# Set plot style
sns.set(style = 'white') 
# Set plot size
plt.figure(figsize=(6,6)) #width=8, height=4
# Set plot style
sns.set_style('whitegrid') 

# Make dataframe
# abdata = df.loc[df['antibody'].isin(['EDE1-C10', 'EDE1-C8', 'SIgN-3C', 'MZ4', 'ZV-67'])]

g = sns.FacetGrid(data = prefs_vs_mutdiffsel_df, col = 'antibody', hue = 'antibody', palette = list(colors.values()),
                 col_wrap = 3)
g.map(sns.scatterplot, 
      'logscaledpref', 
      'mutdiffsel',
     )

# Change titles and axes
g.set_titles(col_template="{col_name}")
g.set_ylabels('AA antibody escape')
g.set_xlabels('AA mutational tolerance')

# Save to figures directory
fig = os.path.join(figsdir,'mutationlevel_SCALEDpref_vs_diffsel.pdf')
g.savefig(fig, dpi = 'figure')
showPDF(fig)
```

# Table figures for the paper
For the supplement, it will be nice to have a table of percent infectivity and antibody concentration used for each replicate. I will generate a PDF of a subset of 'diffsel_batch', defined above. 


```python
percentinfect_table = (diffsel_batch
                       .assign(conc = lambda x: x['grouplabel'].str.split('(').str[1].str.split(')').str[0],
                               temp = lambda x: x['name'],
                               index = lambda x: x['name'].str.strip('lib')
                              )
                      [['antibody', 'conc', 'temp', 'percent_infectivity']]
                       .reset_index(drop=True)
                       
                      )

percentinfect_table['concentration'] = (percentinfect_table['conc']
                                        .str.split(' ').str[0] + ' ng/$\mu$L')

percentinfect_table['temp'] = (percentinfect_table['temp']
                                        .str.replace('-v2', 'v2')
                                             )

percentinfect_table = (percentinfect_table
                       .assign(dmslibrary = lambda x: x['temp'].str.split('-').str[0])
                       .rename(columns = {'percent_infectivity': 'percent infectivity',
                                          'dmslibrary': 'DMS library'})
                       [['antibody', 'concentration', 'DMS library', 'percent infectivity']]
                      )

percentinfect_table
```


```python
# We will render figures in other notebooks
outdir = './paper_figures/neutralizations_escape_mutants/data/'
os.makedirs(outdir, exist_ok=True)

outfile = os.path.join(outdir, 'percent_infectivity_table.csv')
percentinfect_table.to_csv(outfile,
                             index=False,) 
```

# Save a markdown version of this notebook
For ease of readability, we will this notebook in markdown format to a 'notebook' folder within the results. 


```python
notebooksdir = os.path.join(resultsdir, 'notebooks')
os.makedirs(notebooksdir, exist_ok=True)

! jupyter nbconvert selections_analysis.ipynb --to markdown --output ./results/notebooks/selections_analysis.md

```


```python

```


```python

```
