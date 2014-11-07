'''
Created on Oct 29, 2014

@author: davidkavanagh
'''
import ast
import os
import re
import sys
import time

import multiprocessing as mp
import numpy as np


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
        with the sampleData attribute holding the PSIs
        '''
        def __init__(self, event_name, isoformID, chrom, strand, mRNA_start, mRNA_end):
            self.event_name = event_name
            self.isoformID = isoformID
            self.chrom = chrom
            self.strand = strand
            self.mRNA_start = mRNA_start
            self.mRNA_end = mRNA_end
            self.sampleData = {}
        
        def addSampleData(self,sampleName, miso_posterior_mean, ci_low, ci_high, miso_stdev, assigned_reads):
            samplesAlreadySeen = set(self.sampleData.keys())
            if sampleName not in samplesAlreadySeen:
                self.sampleData[sampleName] = {'mean':miso_posterior_mean, 'ci_low':ci_low, 'ci_high':ci_high, 'sd': miso_stdev, 'assigned_reads':assigned_reads}
            else:
                sys.stderr.write("Sample {0} data was seen twice for isoform {1} in genes {2}. Something's wrong here!".format(sampleName, self.isoformID, self.event_name))


def findMisoFiles(inputDir):
    '''
    search directory structure of miso files, and stores their paths to a dict, 
    mapping sample name to the path of the miso gene file
    '''
    sample_gene_miso_files = {}   
    for root, dir, files in os.walk(inputDir):
        if os.path.split(root)[1][:3] == 'chr':
            for f in files:
                if f[-7:] == '.pickle':
                    break
                elif f[-5:] == '.miso':
                    geneName = os.path.splitext(os.path.basename(f))[0]
                    sampleName = re.search("([A-Z]{4}_RNA_[A-Z_]+_\d+)", root).group(0)
                    sample_gene_miso_files.setdefault(sampleName, []).append(os.path.join(root,f))
            sys.stdout.write('Searching {0} for miso files\n'.format(root))

    
    sys.stdout.write('Found {1} miso files for {0} samples\n'.format(len(sample_gene_miso_files), sum([len(v) for v in sample_gene_miso_files.values()])))

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
    assigned_counts_index = [int(i.rsplit(':')[0]) for i in assigned_counts]
    #captures the 0 of 1:0
    assigned_counts_reads = [int(i.rsplit(':')[1]) for i in assigned_counts]
    
    parsed_line = (isoforms, assigned_counts_index, assigned_counts_reads,
                   chrom, strand, mRNA_starts, mRNA_ends )
    
    return parsed_line

def summarizeGeneMisoFiles_worker(sampleName, filePaths, q):
    '''
    Takes one sample and all the gene miso file for that sample and an mp queue
    Process all the genes for this one person, then return the results to the listener
    for inserting into the global objects
    '''
    
    sys.stdout.write('{0} got {1} genes files from sample {2}\n'.format(mp.current_process().name,len(filePaths), sampleName))
    
    geneSummary = {}
    isoformSummary = {}
    c = 0
    for filename in filePaths:
        geneName = os.path.splitext(os.path.basename(filename))[0]
        
        misoFile = open(filename, 'r')
        
        isoforms, assigned_counts_index, assigned_counts_reads, \
        chrom, strand, mRNA_starts, mRNA_ends = parse_miso_line(misoFile.next())
        
        misoFile.next() #skip the next line it's just "sampled_psi\tlog.score"
        
        nMCMCsamplings = 2700 # this seems to constant
        
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
        bundling up the gene and isoform results for passing back the listener. Probably not the optimal way of doing this.
        '''
        geneSummary[geneName] = {'geneName':geneName, 'chrom':chrom, 'strand':strand, 'start':most_start, 'end':most_end, 'sampleName':sampleName, 'reads':totalAssignedReads}
        isoformSummary[geneName] = (isoforms, geneName, chrom, strand, mRNA_starts, mRNA_ends, sampleName, miso_means, miso_sds, miso_ci_high, miso_ci_low, assigned_counts_reads)
        
        '''
        Not really necessary, but more to convince myself that it's not crashing or stalling!
        
        c += 1
        sys.stdout.write('\r{0} genes done!'.format(c))
        sys.stdout.flush()
        '''
    
    res = (sampleName, geneSummary, isoformSummary)
    q.put(res)    
    return res

def summarizeGeneMisoFile_listener(q):
    
    sampleGeneSummary = {}
    sampleIsoformSummary = {}
    
    start = time.time()
    while 1:
        m = q.get()
        if m == 'kill':
            break
        '''
        if stuff is returned from the queue - unpack the various data structures, assign the genes to global genes dict
        assign the isoforms for that gene to the global isoforms dict.
        '''
        sampleName = m[0]
        geneSummary = m[1]
        isoformSummary = m[2]

        
        for gene, data in geneSummary.iteritems():

            sampleGeneSummary.setdefault(data['geneName'], misoGene(data['geneName'], data['chrom'], data['strand'], data['start'], data['end'])\
                               ).addSampleData(data['sampleName'], data['reads'])
            
            isoforms, geneName, chrom, strand, mRNA_starts, mRNA_ends, \
            sampleName, miso_means, miso_sds, miso_ci_high, miso_ci_low, assigned_counts_reads = isoformSummary[gene]
            

            for i, iso in enumerate(isoforms):
                thisEventName = '{0}|{1}'.format(geneName, iso)
                sampleIsoformSummary.setdefault(thisEventName, misoIsoform(geneName, iso, chrom, strand, mRNA_starts[i], mRNA_ends[i])\
                                          ).addSampleData(sampleName, miso_means[i], miso_sds[i], miso_ci_low[i], miso_ci_high[i], assigned_counts_reads[i])
        sys.stdout.write('Processed {0} genes from sample {1}\n'.format(len(sampleGeneSummary), sampleName))
    
    sys.stdout.write('Finished parsing miso sample gene files\n')
    finish = time.time() - start
    sys.stdout.write('Total parsing time was {0} seconds \n'.format(finish))
    return sampleGeneSummary, sampleIsoformSummary

                                          
def outputMatrix(filePrefix, isoObj, sampleList, dataType):
    '''
    Output psis matrix
    '''
                                        
    outputFileName = filePrefix + '_{0}Matrix.txt'.format(dataType) 
    sys.stdout.write('writing {1} matrix to file {0}\n'.format(outputFileName, dataType))
    outputFile = open(outputFileName, 'w')
    outputFile.write('event_name\tisoform\t{0}\n'.format('\t'.join(sampleList)))
    for isoID, obj in isoObj.iteritems():
        #build and output string starting with the event name and the isoform ID
        output = '{0}\t{1}'.format(obj.event_name, obj.isoformID)
        
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
    manager = mp.Manager()
    q = manager.Queue()    
    pool = mp.Pool(nCpus)

    #put listener to work first
    watcher = pool.apply_async(summarizeGeneMisoFile_listener, (q,))
    
    jobs = []

    for sampleName, misoPaths in sample_gene_miso_files.iteritems():
        chunked_misoPaths = slice_it(misoPaths, (nCpus-2))
        for chunk in chunked_misoPaths:
            job = pool.apply_async(summarizeGeneMisoFiles_worker, (sampleName, chunk, q))
            jobs.append(job)
    
    for job in jobs: 
        job.get()

    #now we are done, kill the listener
    q.put('kill')
    
    pool.close()
    pool.join()
    
    sampleGeneSummary = watcher.get()[0]
    sampleIsoformSummary = watcher.get()[1]

    outputMatrix(outputFilePrefix, sampleIsoformSummary, sampleList, 'mean')
    outputMatrix(outputFilePrefix, sampleIsoformSummary, sampleList, 'sd')
    outputMatrix(outputFilePrefix, sampleIsoformSummary, sampleList, 'assigned_reads')
    
    outputGeneCountsMatrix(outputFilePrefix, sampleGeneSummary, sampleList)

if __name__ == '__main__':
    main()
    sys.stdout.write('finished!\n')