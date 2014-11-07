'''
Created on Oct 31, 2014

@author: davidkavanagh
'''
import os
import re
import sys


def findMisoFiles(inputDir):
    '''
    search directory structure of miso files, and stores their paths to a dict, 
    mapping sample name to the path of the miso gene file
    '''
    myRE = re.compile(r"([A-Z]{4}_RNA_\w+).*\.miso")
    all_files = [os.path.join(root,f) for root, dir, files in os.walk(inputDir) for f in files ]
    print 'Found {0} .miso files'.format(len(all_files))
    sample_gene_miso_files = {}
    for f in all_files:
        match = myRE.search(f)
        if match:
            sampleName = match.group(1)
            sample_gene_miso_files[sampleName] = f
        sys.stdout.write('\rSearching {0} for miso files'.format(root))
        sys.stdout.flush()
    
    print('\nFound {1} miso files for {0} samples'.format(len(sample_gene_miso_files), sum([len(v) for v in sample_gene_miso_files.values()])))

    return sample_gene_miso_files

if __name__ == '__main__':
    try:
        inputDir = sys.argv[1]
        
    except IndexError:
        sys.stderr.write('No input directory specified')
        sys.exit()
    
    sample_gene_miso_files = findMisoFiles(inputDir)