#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pysam
import os 
import statistics as st
from collections import defaultdict
from cigar import Cigar
import pandas as pd


def left_clipped(read):
    '''
    Function to determine is read is left clipped
    by L.santuari (sv-channels)
    '''
    
    if read.cigartuples is not None:
            if read.cigartuples[0][0] in [4,5]:
                return True
            else:
                return False
    else:
        return False
    
    
def right_clipped(read):
    '''
    Function to determine is read is left clipped
    by L.santuari (sv-channels)
    '''
    
    if read.cigartuples is not None:
            if read.cigartuples[-1][0] in [4,5]:
                return True
            else:
                return False
    else:
        return

    
def get_sa(read, both=False):
    '''
    Get supplementary Alignment
    If only one SA: returns list with four elements
    If two SAs: returns nested list (two lists within one list)
    '''
    
    if not read.is_unmapped:
        if read.has_tag('SA'):
            if len(read.get_tag('SA')) > 0:
                sa = read.get_tag('SA').split(",") # split the string
                if len(sa) < 7:
                    sa_list = [sa[0],sa[1],sa[2],sa[3]]
                    return sa_list # return the list
               
                elif len(sa) >6 and len(sa) < 12 and both:
                    sa_list = [None,None]
                    sa_list[0] =[sa[0],sa[1],sa[2],sa[3]] # append info first split alignment
                    sa_list[1] = [sa[5][-2:],sa[6],sa[7],sa[8]] # append info second split alignment
                    return sa_list # return the list

def av_distance_reads(bam):
    
    distances = []
    for read in bam.fetch():
        
        if read.is_proper_pair and (not read.is_unmapped)         and read.is_reverse != read.mate_is_reverse         and read.reference_name == read.next_reference_name:
            
            dis = abs(read.next_reference_start - read.reference_start)
            if dis < 10**3:
                
                distances.append(dis)
            
            if len(distances) > 2 * 10 ** 6:
                break
    mean = st.mean(distances)
    stdev = st.stdev(distances)
    return mean, stdev

def get_dataframe(dictionary,chr_list):
    '''
    Makes dataframe of given dictionary 
    binned on chromosome 
    binned on orientation of reads
    '''

    def get_subdf(dictionary,chr):
        if len(dictionary[chr]) == 4:
            cols = ['chr1','start','chr2','end']
        elif len(dictionary[chr]) == 8:
            cols = ['chr1 left','start left','chr2 left','end left','chr1 right','start right','chr2 right', 'end right']
        sides = ['FR','RF','FF','RR']
        for side in range(0,len(sides)):
            sides[side] = pd.DataFrame(dictionary[chr][sides[side]], columns = cols)
        df = pd.concat([sides[0],sides[1],sides[2],sides[3]])
        
        return df


    name = 'chrom{}'
    chroms = []
    for ch in chr_list:
        chroms.append(name.format(ch))   
    
    for i in range(0,len(chroms)):
        chroms[i] = get_subdf(dictionary,chr_list[i])
    total_df = pd.concat([chroms[0], chroms[1]],axis=1,keys = chr_list)
    
    return total_df

class bin_reads:
    '''
    Class to extract: discordant read pairs, split and clipped reads.
    positions of reads are returned in dictionary format
    reads are binned on the contig they map to and orientation (forward, reverse)
    discordant reads are additionally binned based on how far the distance between readpairs deviates from mean +- STDEV 
    '''
    
    
    def __init__(self,bam,chr=False,start=False,end=False):
        '''
        param bam: bamfile (pysam alignment)
        opt.param chr: contig of interest
        opt.param start = start position
        opt.param end = end position
        '''
        
        
        self.header = bam.header
        self.chr_list = list(i['SN'] for i in self.header['SQ'])
        self.orientation = ['FF','RR','FR','RF']
        if not(chr and start and end):
            self.iter = bam.fetch()
        else:
            self.iter = bam.fetch(chr,start,end)
            

    def sort_read(self,read,split=False):
        '''
        
        function to extract information about reads,
        Determines the start and end,
        Determines orientation 
        '''
        
        
        if split:
            orien = {'+':True,'-':False}
            sa_info = get_sa(read)
            if sa_info is not None: 
                if len(sa_info) == 4: # if there is only 1 SA, get_sa returns a list with the length of 4 (contig, start, orientation, cigar)
                    if left_clipped(read) and not right_clipped(read):
                        chr1 = sa_info[0]
                        start = sa_info[1]
                        chr2 = read.reference_name
                        end = read.reference_start
                    elif right_clipped(read) and not left_clipped(read):
                        chr1 = read.reference_name
                        start = read.reference_end
                        chr2 = sa_info[0]
                        end = rsa_info[1]
                    if not read.is_reverse and not orien[sa_info[2]]:
                        orien = 'FR'
                    elif read.is_reverse and orien[sa_info[2]]:
                        orien = 'RF'
                    elif read.is_reverse and not orien[sa_info[2]]:
                        orien = 'RR'
                    else:
                        orien = 'FF'
                        
                        
                    return chr1, chr2, start, end, orien
                
                elif len(sa_info) == 2: # when there are two alignments, get_sa returns a list containing two nested lists in which the contig, start, orientation and cigar are recorded
                    chr1L = sa_info[0][0]
                    startL = sa_info[0][1]
                    chr2L = read.reference_name
                    endL = read.reference_start
                    
                    chr1R = read.reference_name
                    startR = read.reference_end
                    chr2R = sa_info[1][0]
                    endR = sa_info[1][1]
                    orientations = list('orienL','orienR')
                    n = 0
                    for i in range(0,len(orientations)):
                        
                        if orien[sa_info[n][2]] and read.is_reverse:
                            i = 'FR'
                        
                        if not orien[sa_info[n][2]] and not read.is_reverse:
                            i = 'RF'
                        if orien[sa_info[n][2]] and not read.is_reverse:
                            i = 'FF'
                        if not orien[sa_info[n][2]] and read.is_reverse:
                            i = 'RR'
                    
                    return chr1L, startL, chr2L, endL, orientation[0], chr1R, startR, chr2R, endR, orientation[1]
                 
            
        elif not split:
            chr1 = read.reference_name
            chr2 = read.next_reference_name
            end = read.next_reference_start
            
            if left_clipped(read) and not right_clipped(read):
                start = read.reference_start
            elif right_clipped(read) and not left_clipped(read):
                start = read.reference_end
            elif not (left_clipped(read) and right_clipped(read)):
                start = read.reference_start
            elif left_clipped(read) and right_clipped(read):
                start = read.reference_start
                end = read.reference_end
            
            
            if not read.is_reverse and read.mate_is_reverse:
                orien = 'FR'
            elif not read.is_reverse and not read.mate_is_reverse:
                orien = 'FF'
                                                                          
            elif read.is_reverse and not read.mate_is_reverse:
                orien = 'RF'
                                                                          
            else:
                orien = 'RR'
                
            return chr1, chr2, start, end, orien
                                                                          


    def get_clipped(self):
        '''
        
        Creates dictionary with all clipped reads (which do not have a supplementary alignment)
        dictionary contains the location where the reads map back to the reference
        '''
        
        
        clipped = {ch:{i:[] for i in self.orientation} for ch in self.chr_list}
        
        for read in self.iter:
            if not read.is_unmapped:
                if left_clipped(read) or right_clipped(read):
                    if not read.has_tag('SA'):
                        chr1,chr2,start,end,orien = self.sort_read(read=read)
                        if abs(start-end) > 10:
                            clipped[chr1][orien].append([chr1,start,chr2,end])
                
                       
        return clipped
    
    def get_split(self):
        '''
        Get positions of reads with Supplementary Alignment (SA)
        returns two directories:
            one for split reads positions with position of the supplementary alignment
            one for split read pisitions and their respective mate
           '''
        
        
        split = {ch:{i:[] for i in self.orientation} for ch in self.chr_list}
        split_mate = {ch:{i:[] for i in self.orientation} for ch in self.chr_list}
        for read in self.iter:
            if not read.is_unmapped:
               
                if read.has_tag('SA'):
                     
                    sa_info = get_sa(read)
                    if sa_info is not None:
                        chr1,chr2,start,end,orien = self.sort_read(read=read,split=True)
                        if abs(start-end) > 10:
                            split[chr1][orien].append([chr1,start,chr2,end])
                    
                        chr1,chr2,start,end,orien = self.sort_read(read=read)
                        if abs(start-end) > 10:
                            split_mate[chr1][orien].append([chr1,start,chr2,end])
                
        return split, split_mate
    
    def get_both_split(self):
        '''
        Returns positions of reads that are split on both left and right and who have two supplementary alignments
        dictionary contains the coordinates: 
            SA on left side; start position read
            end position read; SA on right side
        '''
        
    
        both_split = {ch:{i:[] for i in self.orientation} for ch in self.chr_list}
        for read in self.iter:
            if not read.is_unmapped:
                if left_clipped(read) and right_clipped(read):
                    if read.has_tag('SA'):
                        sa_info = get_sa(read)
                        if len(sa_info) == 2: # if there are two SAs, sa_info returns a list with two nested lists. Otherwise it returns a list with 4 items
                            chr_1L, chr_2L, startL, endL, orienL,                             chr_1R, chr_2R, startR, endR, orienR = self.sort_read(read=read,split =True)
                            if (abs(startL-endL) and abs(startR-endR)) > 10:
                                both_split[chr_1L][orienL].append([chr_1L,startL,chr_2L,endL,                                                                chr_1R, chr_2R, startR, endR, orienR])
        return both_split
    
    def get_discordant(self,mean,stdev):
        '''
        Get positions discordant read pairs who are not clipped
        binned based on the deviation from the mean +/- STDEV 
        param mean = mean
        param stdev = standard deviation 
        '''
        mean, stdev = av_distance_reads()
        
        times = ['1XSD','2XSD','3XSD']
        discordant = {t:{ch:{i:[] for i in self.orientation} for ch in self.chr_list} for t in times}                                                         
        
        for read in self.iter:
            if not read.is_unmapped and not read.mate_is_unmapped and not right_clipped(read) and not left_clipped(read):
                dis = abs(read.next_reference_start - read.reference_start)
    
                if dis != 0:
                    chr1,chr2,start,end,orien = self.sort_read(read=read)
                    if abs(start-end) > 10:
                        if dis >= (mean+(3*stdev)) or dis <= (mean-(3*stdev)):
                            discordant['3XSD'][chr1][orien].append([chr1,start,chr2,end])
                        elif dis >= (mean+(2*stdev)) or dis <= (mean-(2*stdev)):
                            discordant['2XSD'][chr1][orien].append([chr1,start,chr2,end])
                        elif dis >= (mean+(stdev)) or dis <= (mean-(stdev)):
                            discordant['1XSD'][chr1][orien].append([chr1,start,chr2,end])       
        return discordant
    


def get_read_pos(bam,DataName):
    
    bam = pysam.AlignmentFile(bam,'rb')
    binned = bin_reads(bam)
    # get the chromosome list
    chr_list = binned.chr_list

    # Produce Dataframes for different reads
    ## clipped reads:
    clipped = binned.get_clipped()

    ## discordant reads:
    mean, stdev = av_distance_reads(bam)
    discor = binned.get_discordant(mean = mean, stdev= stdev)

    ## split reads & split reds + mate:
    split, split_mate = binned.get_split()

    ## Double split reads:
    both = binned.get_both_split()


    bam.close()


    # Generate Dataframes and save as csv:
    parent_dir  = os.getcwd()
    outdir = os.join(parent_dir,DataName)
    os.makedirs(outdir, exist_ok=True)



    ## Dataframe clipped 
    clipped_df = get_dataframe(clipped,chr_list)
    outfile = os.join(outdir,'clipped_pos.csv')
    clipped.to_csv(outfile,index=False)

    ## Dataframes discordant readpairs
    stdevs = discor.keys()
    for sd in stdevs:
        df = get_dataframe(discor[sd],chr_list)
        outfile = os.join(outdir,'discordant_{}_pos.csv'.format(sd))
        df.to_csv(outfile,index=False)
    
    ## Dataframes split reads and split + mate reads 
    split_df = get_dataframe(split,chr_list)
    split_mate_df = get_dataframe(split,chr_list)

    df_names = ['split_pos','split_mate_pos']
    ('.'join([df_names[n],'csv'])
    n=0
    for df in [split_df,split_mate_df]:
        outfile = os.join(outdir,('.'join([df_names[n],'csv']))
        df.to_csv(outfile,index=False)
        n +=1

    ## make a dataframe for read double split
    both_split = get_dataframe(both,chr_list)

    outfile = os.join(outdir,'both_split_pos.csv')
    both_split.to_csv(outfile,index=False)
                          
                          
def main():
    parser = argparse.ArgumentParser(description='Extract read positions (discordant/clipped/split)')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        help="Specify input file (BAM)")

    parser.add_argument('-DN',
                        '--DataName'
                        type=str
                        help = "Specify name data/sample")
 
 
    args = parser.parse_args()  
    
    
    get_read_pos(bam=args.bam,DataName = args.DataName)
                          
    
if __name__ == '__main__':
    main()

