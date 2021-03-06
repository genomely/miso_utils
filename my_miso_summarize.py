'''
Created on Nov 6, 2014

@author: davidkavanagh
'''
'''
Created on Oct 29, 2014

@author: davidkavanagh
'''
import ast
import os
import re
import sys
import time
import glob
import multiprocessing as mp
import numpy as np


sampleGeneSummary = {}
sampleIsoformSummary = {}

class misoGene():
    '''
    hold some gene level info. Will be used to output matrix of genes X samples, 
    with each entry being total assigned reads over all isoforms
    '''
    
    def __init__(self, geneName, chrom, strand, start, end):
        self.geneName = geneName
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.stop = end
        self.sampleData= {}
    
    def addSampleData(self, sampleName, totalAssignedReads):
        self.sampleData[sampleName] = totalAssignedReads
    

class misoIsoform():
        '''
        The idea here is to have an object represent each unique event,
        with the sampleData attribute holding the PSIs, counts, and other sample level info
        '''
        def __init__(self, event_name, transcriptID, isoformID, chrom, strand, mRNA_start, mRNA_end):
            self.event_name = event_name
            self.isoformID = isoformID
            self.transcriptID = transcriptID
            self.chrom = chrom
            self.strand = strand
            self.mRNA_start = mRNA_start
            self.mRNA_end = mRNA_end
            self.sampleData = {}
        
        def addSampleData(self,sampleName, miso_posterior_mean, miso_stdev, ci_low, ci_high,assigned_reads):
            samplesAlreadySeen = set(self.sampleData.keys())
            if sampleName not in samplesAlreadySeen:
                self.sampleData[sampleName] = {'mean':miso_posterior_mean, 'ci_low':ci_low, 'ci_high':ci_high, 'sd': miso_stdev, 'assigned_reads':assigned_reads}
            else:
                sys.stderr.write("Sample {0} data was seen twice for isoform {1} in genes {2}. Something's wrong here!".format(sampleName, self.isoformID, self.event_name))


def findMisoFiles(inputDir):
    '''
    must use a quoted string for the inputDir on the command line. 
    glob expansion happens in python and is much faster than looping through the directory structure
    '''
    files = glob.glob(inputDir)
    
    sample_gene_miso_files = {}
     
    '''
    this is a real hack, will not work if sample name scheme changes
    A better strategy would be to have an input file with two fields:
    SampleName \t Path/to/dir/with/miso/files/
    
    Then use glob.glob() on each sample path to get all miso files for each sample. 
    '''
    
    pattern = re.compile("([A-Z]{4}_RNA_[A-Z_]+_\d+)")  
    
    for f in files:
        geneName = os.path.splitext(os.path.basename(f))[0]
        sampleName = pattern.search(f).group(0)
        sample_gene_miso_files.setdefault(sampleName, []).append(f)

    sys.stdout.write('Found {0} miso files for {1} samples\n'.format(len(files),len(sample_gene_miso_files)))

    return sample_gene_miso_files

def parse_miso_line(line):
    
    fields = line[1:].rstrip().rsplit() #has a leading "#"
    
    '''
    worst header format ever!
    '''
    strFields = {}
    for i in fields: 
        entries = i.rsplit('=')
        strFields[entries[0]] = entries[1]
 
    isoforms = ast.literal_eval(strFields['isoforms'])
    
    '''
    Isoforms are of the form :
    'exon:ENST00000410086:1_exon:ENST00000410086:2_exon:ENST00000410086:3_exon:ENST00000410086:4_exon:ENST00000410086:5'
    '''
    transcriptIDs = [i.rsplit(':')[1] for i in isoforms] 
    '''
    #assuming that the transcript IDs are always the same across of the length of the isoforms (maybe not for novel stuff), 
    and that each TxID is unique for each isoform
    '''
    
    chrom = strFields['chrom']
    strand = strFields['strand']
    mRNA_starts = strFields['mRNA_starts'].rsplit(',')
    mRNA_ends = strFields['mRNA_ends'].rsplit(',')
    assigned_counts = strFields['assigned_counts'].rsplit(',')
    
    '''
    assigned_counts is sometimes shorter than isoforms - 
    the difference in length is because trailing counts of 0s are clipped.
    '''
    while len(assigned_counts) < len(isoforms):
        lastIndex = int(assigned_counts[-1].rsplit(':')[0])
        nextIndex = lastIndex + 1
        assigned_counts.append('{0}:0'.format(nextIndex))
           
    '''
    The assigned counts list is sometimes shorter than the number of isoforms
    '''
    #print '#isoforms: {0}\t #assigned counts:{1}\t {2}'.format(len(isoforms), len(assigned_counts), assigned_counts)
    
    '''
    Not using the counts anymore
    '''
    #captures the (0,0) of (0,0):100
    #counts_id = [i.rsplit(':')[0] for i in counts]
    #captures the 100 of (0,0):100
    #counts_reads = [i.rsplit(':')[1] for i in counts]
    
    #captures the 1 of 1:0
    #assigned_counts_index = [int(i.rsplit(':')[0]) for i in assigned_counts] #wasn't using this anywhere
    #captures the 0 of 1:0
    assigned_counts_reads = [int(i.rsplit(':')[1]) for i in assigned_counts]
    
    parsed_line = (transcriptIDs, isoforms, assigned_counts_reads,
                   chrom, strand, mRNA_starts, mRNA_ends )
    
    return parsed_line

def summarizeGeneMisoFiles_worker(args):
    '''
    Takes one sample and all the gene miso file for that sample and an mp queue
    Process all the genes for this one person, then return the results to the listener
    for inserting into the global objects
    '''
    sampleName = args[0]
    filePaths = args[1]
    idx = args[2]
    
    sys.stdout.write('{0} got {1} genes files from sample {2} (chunk no. {3}) \n'.format(mp.current_process().name,len(filePaths), sampleName, idx))
    
    geneSummary = {}
    isoformSummary = {}
    c = 0
    for filename in filePaths:
        geneName = os.path.splitext(os.path.basename(filename))[0]
        
        misoFile = open(filename, 'r')
        
        transcriptIDs, isoforms, assigned_counts_reads, \
        chrom, strand, mRNA_starts, mRNA_ends = parse_miso_line(misoFile.next())
        
        misoFile.next() #skip the next line it's just "sampled_psi\tlog.score"
        
        nMCMCsamplings = 2700 # this seems to be constant
        
        psis = np.zeros([nMCMCsamplings, len(isoforms)])
    
        counter = 0
        for line in misoFile:
            fields = line.rstrip().rsplit()
            psi =  [float(i) for i in fields[0].rsplit(',')]
            psis[counter,] = psi
            counter += 1
        
        miso_means = np.mean(psis, axis=0)
        miso_sds = np.std(psis, axis=0)
        miso_ci_low = np.percentile(psis, axis=0, q=2.5)
        miso_ci_high = np.percentile(psis, axis=0, q=97.5)
        
        totalAssignedReads = sum([int(i) for i in assigned_counts_reads])
              
        if strand == '+':
            most_start = min([int(i) for i in mRNA_starts])
            most_end = max([int(i) for i in mRNA_ends])
        else:
            most_start = max([int(i) for i in mRNA_starts])
            most_end = min([int(i) for i in mRNA_ends])
        '''
        bundling up the gene and isoform results for passing back to the result list. Probably not the optimal way of doing this.
        '''
        geneSummary[geneName] = {'geneName':geneName, 'chrom':chrom, 'strand':strand, 'start':most_start, 'end':most_end, 'sampleName':sampleName, 'reads':totalAssignedReads}
        isoformSummary[geneName] = (transcriptIDs, isoforms, geneName, chrom, strand, mRNA_starts, mRNA_ends, sampleName, miso_means, miso_sds, miso_ci_high, miso_ci_low, assigned_counts_reads)
        
    
    res = (sampleName, geneSummary, isoformSummary)  
    return res

def combineMappedResults(result_List):
    sampleGeneSummary = {}
    sampleIsoformSummary = {}
    
    for sampleName, geneSummary, isoformSummary in result_List:

        for gene, data in geneSummary.iteritems():
            sampleGeneSummary.setdefault(data['geneName'], misoGene(data['geneName'], data['chrom'], data['strand'], data['start'], data['end'])\
                               ).addSampleData(data['sampleName'], data['reads'])

            transcriptIDs, isoforms, geneName, chrom, strand, mRNA_starts, mRNA_ends, \
                sampleName, miso_means, miso_sds, miso_ci_high, miso_ci_low, assigned_counts_reads = isoformSummary[gene]
                 
            for i, iso in enumerate(isoforms):
                thisEventName = '{0}|{1}'.format(geneName, iso)
                sampleIsoformSummary.setdefault(thisEventName, misoIsoform(geneName, transcriptIDs[i], iso, chrom, strand, mRNA_starts[i], mRNA_ends[i])\
                                          ).addSampleData(sampleName, miso_means[i], miso_sds[i], miso_ci_low[i], miso_ci_high[i], assigned_counts_reads[i])
                                          
    return (sampleGeneSummary, sampleIsoformSummary)
                                          
def outputMatrix(filePrefix, isoObj, sampleList, dataType):
    '''
    Output different matrix
    '''
    if dataType == 'mean':
        filePostFix = '_meanPsiMatrix.txt'
    if dataType == 'sd':
        filePostFix = '_standardDevMatrix.txt'
    if dataType == 'assigned_reads':
        filePostFix = '_assignedReadsMatrix.txt'
    
    outputFileName = filePrefix + filePostFix
    sys.stdout.write('writing {1} matrix to file {0}\n'.format(outputFileName, dataType))
    outputFile = open(outputFileName, 'w')
    outputFile.write('transcriptID\t{0}\n'.format('\t'.join(sampleList)))
    for obj in isoObj.values():
        #build and output string starting with the event name and the isoform ID
        output = '{0}'.format(obj.transcriptID)
        
        samples = []
        for s in sampleList:
            try:
                samples.append(str(obj.sampleData[s][dataType]))
            except KeyError:
                samples.append('NA')
        #then add in the sample psis which were collected into a list above. Not that this was done to preseve ordering of the samples between rows
        output += '\t' + '\t'.join(samples)
        output += '\n'
        outputFile.write(output)
    outputFile.close()  
    
def outputGeneCountsMatrix(filePrefix, geneObj, sampleList):
    
    outputFileName = filePrefix + '_geneCountsMatrix.txt'
    sys.stdout.write('writing gene counts matrix to file {0}\n'.format(outputFileName))
    outputFile = open(outputFileName, 'w')
    outputFile.write('event_name\t{0}\n'.format('\t'.join(sampleList)))
    for geneID, obj in geneObj.iteritems():
        #build and output string starting with the event name and the isoform ID
        output = '{0}'.format(geneID)
        
        samples = []
        for s in sampleList:
            try:
                samples.append(str(obj.sampleData[s]))
            except KeyError:
                samples.append('NA')
        
        output += '\t' + '\t'.join(samples)
        output += '\n'
        outputFile.write(output)
    outputFile.close()
    
def outputIsoformMetaInfo(filePrefix, isoObj):
    outputFileName = filePrefix + '_isoformMetaInfo.txt'
    sys.stdout.write('writing isoform meta info to file {0}\n'.format(outputFileName))
    outputFile = open(outputFileName, 'w')
    outputFile.write('transcriptID\tisoformID\tGENE\tCHROM\tSTRAND\tSTART\tEND\n')
    for iso, obj in isoObj.iteritems():
        output = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(obj.transcriptID, obj.isoformID, obj.event_name, obj.chrom, obj.strand, obj.mRNA_start, obj.mRNA_end)
        outputFile.write(output)
    outputFile.close()
        


def slice_it(li, n):
    start = 0
    for i in xrange(n):
        stop = start + len(li[i::n])
        yield li[start:stop]
        start = stop                                         
       
def main():
    try:
        inputDir = sys.argv[1]
        
    except IndexError:
        sys.stderr.write('No input directory specified')
        sys.exit()
        
    try:
        outputFilePrefix = sys.argv[2]
        
    except IndexError:
        sys.stderr.write('No output file name specified - using "misosummary_allSamples_<type>matrix.txt')
        outputFilePrefix = 'misosummary_allSamples'
    
    try:
        nCpus = int(sys.argv[3])
    except IndexError:
        sys.stderr.write('No nCPUs given - assuming 6 (4 workers, 1 listener, 1 main)')
        nCpus = 6
    except ValueError:
        sys.stderr.write('nCPUs must be an integer - assuming 6 (4 workers, 1 listener, 1 main)')
        nCpus = 6
        
    sample_gene_miso_files = findMisoFiles(inputDir)
    sampleList = sample_gene_miso_files.keys()
    #must use Manager queue here, or will not work
    pool = mp.Pool(nCpus)

    #put listener to work first
    
    jobs = []
    
    idx = 0
    for sampleName, misoPaths in sample_gene_miso_files.iteritems():
            chunked_misoPaths = slice_it(misoPaths, (nCpus-2))
            for chunk in chunked_misoPaths:
                job = (sampleName, chunk, idx)
                jobs.append(job)
                idx += 1
     
    start = time.time()
    
    results = pool.map(summarizeGeneMisoFiles_worker, jobs)
    
    pool.close()
    pool.join()
    
    finish = time.time() - start
    
    sys.stdout.write('Finished parsing miso sample gene files\n')
    sys.stdout.write('Total parsing time was {0} seconds \n'.format(finish))
    
    sys.stdout.write('Combining results...\n')
    start = time.time()
    sampleGeneSummary, sampleIsoformSummary = combineMappedResults(results)
    finish = time.time() - start
    sys.stdout.write('Results combined. Time taken to combine: {0}\n'.format(finish))
     
    outputMatrix(outputFilePrefix, sampleIsoformSummary, sampleList, 'mean')
    outputMatrix(outputFilePrefix, sampleIsoformSummary, sampleList, 'sd')
    outputMatrix(outputFilePrefix, sampleIsoformSummary, sampleList, 'assigned_reads')
    
    outputGeneCountsMatrix(outputFilePrefix, sampleGeneSummary, sampleList)
    
    outputIsoformMetaInfo(outputFilePrefix, sampleIsoformSummary)
    
    


if __name__ == '__main__':
    main()
    sys.stdout.write('finished!\n')