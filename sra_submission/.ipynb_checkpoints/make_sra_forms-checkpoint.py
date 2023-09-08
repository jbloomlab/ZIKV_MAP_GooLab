''' modified by Caroline Kikawa, Aug 2023, to upload Zika MR766 deep mutational scanning datasets'''

import glob
import os
import sys
import pandas as pd
import tarfile

# The experiments are listed in `samplelist.csv` 
exptsdf = (pd.read_csv('../data/samplelist.csv')
			.assign(sample_name=lambda x: x.selection + '_' + x.library))


biosampleid = 'SAMN37316695' # this biosampleid was created for the ZIKV DMS bnAb selection.
# the bioproject number is PRJNA530795 for ZIKV DMS sequencing.
bioprojectid = 'PRJNA530795' 

tar_filename = 'ZIKV_E_Ab_selection.tar'
make_tar = True # set to true to actually make the tarfile; probably want to test with *make_tar = False* at first to make sure all the files are found properly.

fnew = open('submissionform_Aug2023_ZikaDMS_bnAbSelections.tsv','w')

# header:
fnew.write('bioproject_accession\tbiosample_accession\tlibrary_ID\ttitle\tlibrary_strategy\tlibrary_source\tlibrary_selection\tlibrary_layout\tplatform\tinstrument_model\tdesign_description\tfiletype\tfilename\tfilename2\tfilename3\tfilename4\tfilename5\tfilename6\tfilename7\tfilename8\tfilename9\tfilename10\tfilename11\tfilename12\tfilename13\tfilename14\tfilename15\tfilename16\tfilename17\tfilename18\tfilename19\tfilename20\tfilename21\tfilename22\tfilename23\tfilename24\n')

# open a tar file to add all the fastq files to:
if make_tar:
    tar_out = tarfile.open(tar_filename, mode='w')

# go through samples, adding information to the submission form and adding FASTQ.gz files to the tar archive:
for sample_tup in exptsdf.itertuples():
    sample_name = sample_tup.sample_name
    print ('\nProcessing sample {0}...'.format(sample_name))

    FASTQ_directory = sample_tup.R1
    print(FASTQ_directory)

    # sample_name is "antibody" + "library" on spreadsheet
    antibody = sample_tup.antibody
    library = sample_tup.library


    if ('ZV67' in sample_name):
        title = 'MR766 E library {0} neutralized with ZV-67'.format(library)
    elif ('MZ4' in sample_name):
        title = 'MR766 E library {0} neutralized with MZ4'.format(library)
    elif ('C10' in sample_name):
        title = 'MR766 E library {0} neutralized with EDE1-C10'.format(library)
    elif ('C8' in sample_name):
        title = 'MR766 E library {0} neutralized with EDE1-C8'.format(library)
    elif ('SIgN' in sample_name):
        title = 'MR766 E library {0} neutralized with SIgN-3C'.format(library)
    elif ('no-antibody' in sample_name):
        title = 'MR766 E library {0} mock-neutralized'.format(library)


    # Each sample is a row in the spreadsheet with several filenames at the end of the row, here is everything until the filenames:
    fnew.write('{0}\t{1}\t{2}\t{3}\tAMPLICON\tVIRAL RNA\tPCR\tpaired\tILLUMINA\tIllumina MiSeq V3\t250-nt paired-end reads of Zika virus E PCR amplicons\tfastq\t'.format(bioprojectid,biosampleid,sample_name,title))

    # Now add a tab-spaced entry in this library_ID's row for each file that will be uploaded:
    files = sorted(glob.glob(FASTQ_directory)+glob.glob(FASTQ_directory.replace('R1', 'R2')))
    print ('Found these files:'), files
    
    for fname in files:
        print ('Processing {0}...'.format(fname))
        sys.stdout.flush() # forces output to be written to the terminal
        short_fname = fname[fname.rfind('/')+1:]
    
        fnew.write('{0}\t'.format(short_fname))
        fnew.flush()
    
        if make_tar:
            print ('Adding {0} to the tar archive'.format(fname))
            # use the short filename when adding to the tar instead of the full path
            tar_out.add(fname, arcname=short_fname)
    
    fnew.write('\n') # end the current sample's line in the spreadsheet.

fnew.close()

if make_tar:
    tar_out.close()
    print ('Completed.')
    print ('Contents of {0}:'.format(tar_filename))
    t = tarfile.open(tar_filename, 'r')
    for member_info in t.getmembers():
        print (member_info.name)
