'''
Created on Oct 27, 2014

@author: davidkavanagh
'''
import os
import re
import sys


class misoEvent():
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
        
        def addSampleData(self,sampleName, miso_posterior_mean, ci_low, ci_high, disgarded_reads, informative_reads):
            samplesAlreadySeen = set(self.sampleData.keys())
            if sampleName not in samplesAlreadySeen:
                self.sampleData[sampleName] = {'mean':miso_posterior_mean, 'ci_low':ci_low, 'ci_high':ci_high, 'disgarded_reads': disgarded_reads, 'informative_reads': informative_reads}
            else:
                sys.stderr.write("Sample {0} data was seen twice for isoform {1} in genes {2}. Something's wrong here!".format(sampleName, self.isoformID, self.event_name))

def parse_miso_line(line):
    '''
    -----Miso summary files have the following fields-----
    event_name            (single value string)      
    miso_posterior_mean   (comma separated list of float)     
    ci_low                (comma separated list of float)  
    ci_high               (comma separated list of float)
    isoforms              (comma separated list of quoted strings)
    counts                (Annoyingly formatted list that does not allow easy parsing - comma separated list of the form (0,0):1,(0,1):2 - may have to use regexp to parse these!)
    assigned_counts       (comma separated list of the form 0:1,1:2) 
    chrom                 (single value string)
    strand                (single value +/-)
    mRNA_starts           (comma separated list of integers)
    mRNA_ends             (comma separated list of integers)
    '''
    fields = line.rstrip().rsplit()
    event_name = fields[0]
    miso_posterior_mean, ci_low, ci_high = fields[1].rsplit(','), fields[2].rsplit(','), fields[3].rsplit(',')  
    isoforms = fields[4].rsplit(',')
    counts = re.split(",\(", fields[5]) # counts is not used. This field would be difficult to parse
    assigned_counts = fields[6].rsplit(',') 
    chrom,  strand, mRNA_starts,  mRNA_ends = fields[7], fields[8], fields[9].rsplit(','), fields[10].rsplit(',')
    
    '''
    Annoyingly, if there are only two isoforms, Miso outputs only one value. 
    Need to catch this add the 1-x to the list
    '''
    if len(miso_posterior_mean) == 1:
        miso_posterior_mean = [float(miso_posterior_mean[0]), 1-float(miso_posterior_mean[0])]
        ci_low = [float(ci_low[0]), 1-float(ci_low[0])]
        ci_high = [float(ci_high[0]), 1-float(ci_high[0])]
        
    '''
    The assigned counts list is sometimes shorter than the number of isoforms
    '''
    
    #print '#isoforms: {0}\t #assigned counts:{1}\t {2}'.format(len(isoforms), len(assigned_counts), assigned_counts)
    
    #captures the (0,0) of (0,0):100
    counts_id = [i.rsplit(':')[0] for i in counts]
    #captures the 100 of (0,0):100
    counts_reads = [i.rsplit(':')[1] for i in counts]
    
    #captures the 1 of 1:0
    assigned_counts_index = [int(i.rsplit(':')[0]) for i in assigned_counts]
    #captures the 0 of 1:0
    assigned_counts_reads = [int(i.rsplit(':')[1]) for i in assigned_counts]
    
    parsed_line = (event_name,event_name, isoforms, 
                   chrom, strand, mRNA_starts, mRNA_ends,counts_id, counts_reads, 
                   assigned_counts_index, assigned_counts_reads, 
                   miso_posterior_mean, ci_low, ci_high)
    
    return parsed_line
            

if __name__ == '__main__':
    #take the directory of input files as a command line positional argument
    try:
        inputdir = sys.argv[1]
        
    except IndexError:
        sys.stderr.write('No input directory specified')
        sys.exit()
        
    try:
        outputFilePrefix = sys.argv[2]
        
    except IndexError:
        sys.stderr.write('No output file name specified - using "misosummary_allSamples_<type>matrix.txt')
        outputFilePrefix = 'misosummary_allSamples'
        
    
    for root, dir, files in os.path.walk(inputdir):
        for f in files:
            print f
    
      
    #add the directory path to the filename
    files = [f for f in os.listdir(inputdir)]
    filepaths = [os.path.join(inputdir, f) for f in files]
        
    all_isoforms = {}
    all_samples = []
    
    #loop through the files
    for f in filepaths:
        sampleName = os.path.splitext(os.path.basename(f))[0]
        print ('processing sample {0}'.format(sampleName))
        all_samples.append(sampleName)
        miso_summary = open(f, 'r')
        miso_summary.next() #disregard the Header
        events = {}
        for line in miso_summary:
            event_name,event_name, isoforms, \
                   chrom, strand, mRNA_starts, mRNA_ends, \
                   counts_id, counts_reads, \
                   assigned_counts_index, assigned_counts_reads, \
                   miso_posterior_mean, ci_low, ci_high = parse_miso_line(line) #use the parse function above to parse the line
            
            '''
            Assigned counts are the only isoforms with reads supporting them. 
            If we want this info then we don't get the psis of the isoforms to 
            which no reads were inferred (should be zero but they are often not).
            At the moment I am ignoring the assigned counts entirely, and outputting all psis!
            '''
            
            ''' 
            #debugging statement to show the above is true   
            if len(miso_posterior_mean) > len(assigned_counts_index):
                print '#isoforms: {0}\t#psis: {1}\t#counts: {2}\t#assigned_reads: {3}'.format(len(isoforms), len(miso_posterior_mean),len(counts_id), len(assigned_counts_index))
            '''
               
            '''
            counts is a vector howmany reads are informative for which isoforms, however, the first entry is sometimes an extra entry of reads that were disgarded.
            Want to report # disgarded reads and sum of informative reads
            if the first entry is disgarded reads then the vector is all 0s. Split them and check they sum to zero.
            '''
            disgarded_reads = 0
            informative_reads = 0

            #            first element    get rid of brackets    split zeros
            firstCombo = counts_id[0].rstrip(')').lstrip('(').rsplit(',')
            #sum and check equal to 0
            firstComboSum = sum([int(x) for x in firstCombo])
            if firstComboSum == 0:
                #save those reads to a different variable and remove from the remaining reads
                disgarded_reads = counts_reads[0]
                counts_reads = counts_reads[1:]
            #sum remaining reads
            informative_reads = sum([int(x) for x in counts_reads])
            
            #loop through the isoforms of in this row, and add psis and cis, to the sample data dictionaries for the event
            for i, iso in enumerate(isoforms):
                isoformID  = iso[1:-1]
                thisEventName = '{0}|{1}'.format(event_name, isoformID)
                    
                all_isoforms.setdefault(thisEventName, misoEvent(event_name, isoformID, chrom, strand, mRNA_starts[i], mRNA_ends[i] \
                                        )).addSampleData(sampleName, miso_posterior_mean[i], ci_low[i], ci_high[i], disgarded_reads, informative_reads)
                                          
    '''
    Output psis matrix
    '''
                                        
    outputFileName = outputFilePrefix + '_psiMatrix.txt' 
    print('writing psi matrix to file {0}'.format(outputFileName))
    outputFile = open(outputFileName, 'w')
    outputFile.write('event_name\tisoform\t{0}\n'.format('\t'.join(all_samples)))
    for isoID, isoObj in all_isoforms.iteritems():
        #build and output string starting with the event name and the isoform ID
        output = '{0}\t{1}'.format(isoObj.event_name, isoObj.isoformID)
        
        samplePsis = []
        for s in all_samples:
            try:
                samplePsis.append(str(isoObj.sampleData[s]['mean']))
            except KeyError:
                samplePsis.append('NA')
        #then add in the sample psis which were collected into a list above. Not that this was done to preseve ordering of the samples between rows
        output += '\t' + '\t'.join(samplePsis)
        output += '\n'
        outputFile.write(output)
    outputFile.close()
        
    '''
    Output  low cis matrix
    '''     
    outputFileName = outputFilePrefix + '_cilowMatrix.txt'
    print('writing cis matrix to file {0}'.format(outputFileName)) 
    outputFile = open(outputFileName, 'w')

    outputFile.write('event_name\tisoform\t{0}\n'.format('\t'.join(all_samples)))
    for isoID, isoObj in all_isoforms.iteritems():
        #build and output string starting with the event name and the isoform ID
        output = '{0}\t{1}'.format(isoObj.event_name, isoObj.isoformID)
        
        sample_cis = []
        
        for s in all_samples:
            try:
                sample_cis.append(str(isoObj.sampleData[s]['ci_low']))
            except KeyError:
                sample_cis.append('NA')
        
        output += '\t' + '\t'.join(sample_cis)
        output += '\n'
        outputFile.write(output)
        
    outputFile.close()
    
    '''
    Output  high cis matrix
    '''     
    outputFileName = outputFilePrefix + '_cihighMatrix.txt'
    print('writing cis matrix to file {0}'.format(outputFileName)) 
    outputFile = open(outputFileName, 'w')

    outputFile.write('event_name\tisoform\t{0}\n'.format('\t'.join(all_samples)))
    for isoID, isoObj in all_isoforms.iteritems():
        #build and output string starting with the event name and the isoform ID
        output = '{0}\t{1}'.format(isoObj.event_name, isoObj.isoformID)
        
        sample_cis = []
        
        for s in all_samples:
            try:
                sample_cis.append(str(isoObj.sampleData[s]['ci_high']))
            except KeyError:
                sample_cis.append('NA')
        
        output += '\t' + '\t'.join(sample_cis)
        output += '\n'
        outputFile.write(output)
        
    outputFile.close()
        
    '''
    Output disgarded and informative reads matrix
    '''     
    outputFileName = outputFilePrefix + '_readsumMatrix.txt' 
    print('writing readsum matrix to file {0}'.format(outputFileName))
    outputFile = open(outputFileName, 'w')
    
    disgard_samples = [i+'_disgard' for i in all_samples]
    inform_samples = [i+'_inform' for i in all_samples]
    
    readsum_samples = [val for pair in zip(disgard_samples,inform_samples) for val in pair]
    
    outputFile.write('event_name\tisoform\t{0}\n'.format('\t'.join(readsum_samples)))
    for isoID, isoObj in all_isoforms.iteritems():
        #build and output string starting with the event name and the isoform ID
        output = '{0}\t{1}'.format(isoObj.event_name, isoObj.isoformID)
        
        sample_disgard = []
        sample_inform = []
        
        for s in all_samples:
            try:
                sample_disgard.append(str(isoObj.sampleData[s]['disgarded_reads']))
                sample_inform.append(str(isoObj.sampleData[s]['informative_reads']))
            except KeyError:
                sample_disgard.append('NA')
                sample_inform.append('NA')
                
        sampleReadSums = [val for pair in zip(sample_disgard,sample_inform) for val in pair]
        
        output += '\t' + '\t'.join(sampleReadSums)
        output += '\n'
        
        outputFile.write(output)
    outputFile.close()
    
    print('Finished!')        
        