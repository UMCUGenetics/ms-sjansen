#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import json
import argparse
import pysam
import os 
import statistics as st
from collections import defaultdict
from cigar import Cigar
import pandas as pd
from time import time
import logging

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

    
def get_sa(read):
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
               
                elif len(sa) >6 and len(sa) < 12:
                    sa_list = [None,None]
                    sa_list[0] =[sa[0],sa[1],sa[2],sa[3]] # append info first split alignment
                    sa_list[1] = [sa[5][-2:],sa[6],sa[7],sa[8]] # append info second split alignment
                    return sa_list # return the list

def av_distance_reads(bam):
    
    distances = []
    for read in bam.fetch():
        
        if read.is_proper_pair and (not read.is_unmapped) \
        and read.is_reverse != read.mate_is_reverse \
        and read.reference_name == read.next_reference_name:
            
            dis = abs(read.next_reference_start - read.reference_start)
            if dis < 10**3:
                
                distances.append(dis)
            
            if len(distances) > 2 * 10 ** 6:
                break
    mean = st.mean(distances)
    stdev = st.stdev(distances)
    return mean, stdev

def get_dataframe(dictionary,chr_list,bs=False):
    '''
    Makes dataframe of given dictionary 
    binned on chromosome 
    binned on orientation of reads
    '''

    def get_subdf(dictionary,chr,bs = False):
        if bs:
            cols = ['chr1_L','start_L','chr2_L','end_L','chr1_R','start_R','chr2_R', 'end_R', 'orientation', 'read type']
        else:
            cols = ['chr1','start','chr2','end','orientation', 'read type']
        sides = ['FR','RF','FF','RR']
        for side in range(0,len(sides)):
            sides[side] = pd.DataFrame(dictionary[chr][sides[side]], columns = cols)
        df = pd.concat([sides[0],sides[1],sides[2],sides[3]],keys = ['FR','RF','FF','RR'])
        
        return df


    name = 'chrom{}'
    chroms = []
    for ch in chr_list:
        chroms.append(name.format(ch))   
    
    for i in range(0,len(chroms)):
        if bs:
            chroms[i] = get_subdf(dictionary,chr_list[i],bs=True)
        else:
            chroms[i] = get_subdf(dictionary,chr_list[i])
    total_df = pd.concat([chroms[0], chroms[1]],keys = chr_list)
    total_df = total_df.dropna(axis = 'rows')
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
                    logging.info('start process split read %s' % read.qname)
                    if left_clipped(read) and not right_clipped(read):
                        chr1 = sa_info[0]
                        start = sa_info[1]
                        chr2 = read.reference_name
                        end = read.reference_start
                        type = 'SPLIT_L'
                    elif right_clipped(read) and not left_clipped(read):
                        chr1 = read.reference_name
                        start = read.reference_end
                        chr2 = sa_info[0]
                        end = sa_info[1]
                        type = 'SPLIT_R'
                    elif right_clipped(read) and left_clipped(read):
                        cigstring = read.cigarstring
                        for l in range(0,len(cigstring)):
                            if cigstring[l] in ['S','H']:
                                clip = l
                            elif cigstring[l] == 'M':
                                Match = l
                        if clip > Match:
                            chr1 = sa_info[0]
                            start = sa_info[1]
                            chr2 = read.reference_name
                            end = read.reference_start
                            type = 'SPLIT_L'
                        elif clip < Match:
                            chr1 = read.reference_name
                            start = read.reference_end
                            chr2 = sa_info[0]
                            end = sa_info[1]
                            type = 'SPLIT_R'
                        else:
                            chr1 = read.reference_name
                            start = 'NaN'
                            chr2 = read.reference_name
                            end = 'NaN'
                            type = 'NaN'
                    if not read.is_reverse and not orien[sa_info[2]]:
                        orien = 'FR'
                    elif read.is_reverse and orien[sa_info[2]]:
                        orien = 'RF'
                    elif read.is_reverse and not orien[sa_info[2]]:
                        orien = 'RR'
                    else:
                        orien = 'FF'
                    
                    logging.info('Split Read %s processed at position: %s %s' %(read.qname,read.reference_name, read.reference_start))
                        
                    return chr1, chr2, start, end, orien, type
                
                elif len(sa_info) == 2: # when there are two alignments, get_sa returns a list containing two nested lists in which the contig, start, orientation and cigar are recorded
                    #logging.info('start process two-sided split read %s' % read.qname)
                    chr1L = sa_info[0][0]
                    startL = sa_info[0][1]
                    chr2L = read.reference_name
                    endL = read.reference_start
                    
                    chr1R = read.reference_name
                    startR = read.reference_end
                    chr2R = sa_info[1][0]
                    endR = sa_info[1][1]
                    type = 'SPLIT_LR'
                    
                    orientation = ['orienL','orienR']
                    n = 0
                    for i in range(0,len(orientation)):
                        
                        if orien[sa_info[n][2]] and read.is_reverse:
                            orientation[i] = 'FR'
                        
                        if not orien[sa_info[n][2]] and not read.is_reverse:
                            orientation[i] = 'RF'
                        if orien[sa_info[n][2]] and not read.is_reverse:
                            orientation[i] = 'FF'
                        if not orien[sa_info[n][2]] and read.is_reverse:
                            orientation[i] = 'RR'
                    
                    logging.info('Two-sided split Read %s processed at position: %s %s' %(read.qname,read.reference_name, read.reference_start))

                    
                    return chr1L, startL, chr2L, endL, orientation[0], chr1R, startR, chr2R, endR, orientation[1], type
                 
            
        elif not split:
            chr1 = read.reference_name
            chr2 = read.next_reference_name
            end = read.next_reference_start
            
            if left_clipped(read) and not right_clipped(read):
                start = read.reference_start
                type = 'CLIP_R'
            elif right_clipped(read) and not left_clipped(read):
                start = read.reference_end
                type = 'CLIP_L'
            elif not (left_clipped(read) and right_clipped(read)):
                start = read.reference_start
                type = 'DISCOR'
            elif left_clipped(read) and right_clipped(read):
                start = read.reference_start
                end = read.reference_end
                type = 'CLIP_LR'
            
            
            if not read.is_reverse and read.mate_is_reverse:
                orien = 'FR'
            elif not read.is_reverse and not read.mate_is_reverse:
                orien = 'FF'
                                                                          
            elif read.is_reverse and not read.mate_is_reverse:
                orien = 'RF'
                                                                          
            else:
                orien = 'RR'
                
            return chr1, chr2, start, end, orien, type
                                                                          


    def get_clipped(self):
        '''
        
        Creates dictionary with all clipped reads (which do not have a supplementary alignment)
        dictionary contains the location where the reads map back to the reference
        '''
        
        n_r = 10 ** 6
        last_t = time()
        clipped = {ch:{i:[] for i in self.orientation} for ch in self.chr_list}

        for i, read in enumerate(self.iter):
            if not i % n_r:
                logging.info("%d clipped reads processed (%f alignments / s)" % \
                             (i, n_r / (time() - last_t)))
                last_t = time()
            if not read.is_unmapped:
                if left_clipped(read) or right_clipped(read):
                    if not read.has_tag('SA'):
                        chr1,chr2,start,end,orien, type = self.sort_read(read=read)
                        if abs(start-end) > 10:
                            clipped[chr1][orien].append([chr1,start,chr2,end,orien, type])
                
                       
        return clipped
    
    def get_split(self):
        '''
        Get positions of reads with Supplementary Alignment (SA)
        returns two directories:
            one for split reads positions with position of the supplementary alignment
            one for split read pisitions and their respective mate
           '''
        n_r = 10 ** 6
        last_t = time()
        
        split = {ch:{i:[] for i in self.orientation} for ch in self.chr_list}
        split_mate = {ch:{i:[] for i in self.orientation} for ch in self.chr_list}
        for i, read in enumerate(self.iter):
            if not i % n_r:
                logging.info("%d split reads processed (%f alignments / s)" % \
                             (i, n_r / (time() - last_t)))
                last_t = time()
            if not read.is_unmapped and read.has_tag('SA'):
                     
                sa_info = get_sa(read)
                if sa_info is not None and len(sa_info) == 4:
                    chr1, chr2, start, end, orien, type = self.sort_read(read=read,split=True)
                    if abs(int(start)-int(end)) > 10:
                        split[chr1][orien].append([chr1, start, chr2, end, orien,type])

                    chr1, chr2, start, end, orien, type = self.sort_read(read=read)
                    if abs(int(start)-int(end)) > 10:
                        split_mate[chr1][orien].append([chr1, start, chr2, end, orien, type])

        return split, split_mate
    
    def get_both_split(self):
        '''
        Returns positions of reads that are split on both left and right and who have two supplementary alignments
        dictionary contains the coordinates: 
            SA on left side; start position read
            end position read; SA on right side
        '''
        n_r = 10 ** 6
        last_t = time()
        
        both_split = {ch:{i:[] for i in self.orientation} for ch in self.chr_list}
        for i, read in enumerate(self.iter):
            if not i % n_r:
                logging.info("%d two-side split reads processed (%f alignments / s)" % \
                             (i, n_r / (time() - last_t)))
                last_t = time()
                
            if not read.is_unmapped and read.has_tag('SA'):
                sa_info = get_sa(read)
                if sa_info is None:
                    continue
                if len(sa_info) == 2:
                    chr_1L, startL, chr_2L, endL, orienL, \
                    chr_1R, chr_2R, startR, endR, orienR, type = self.sort_read(read,split =True)
                    
                    
                    if (abs(int(startL)-int(endL)) and abs(int(startR)-int(endR))) > 10:
                        both_split[chr_2L][orienL].append([chr_1L,startL,chr_2L,endL, \
                                                           chr_1R, chr_2R, startR, endR, orienL, type])
        return both_split
    
    def get_discordant(self,mean,stdev):
        '''
        Get positions discordant read pairs who are not clipped
        binned based on the deviation from the mean +/- STDEV 
        param mean = mean
        param stdev = standard deviation 
        '''
        n_r = 10 ** 6
        last_t = time()
        
        times = ['1XSD','2XSD','3XSD']
        discordant = {t:{ch:{i:[] for i in self.orientation} for ch in self.chr_list} for t in times}                                                         
        
        for i, read in enumerate(self.iter):
            if not i % n_r:
                logging.info("%d discordant reads processed (%f alignments / s)" % \
                             (i, n_r / (time() - last_t)))
                last_t = time()
            if not read.is_unmapped and not read.mate_is_unmapped and not right_clipped(read) and not left_clipped(read):
                dis = abs(read.next_reference_start - read.reference_start)
    
                if dis != 0:
                    chr1,chr2,start,end,orien,type = self.sort_read(read=read)
                    if abs(start-end) > 10:
                        if dis >= (mean+(3*stdev)) or dis <= (mean-(3*stdev)):
                            type = 'DISCOR_3XSD'
                            discordant['3XSD'][chr1][orien].append([chr1,start,chr2,end,orien,type])
                        elif dis >= (mean+(2*stdev)) or dis <= (mean-(2*stdev)):
                            type = 'DISCOR_2XSD'
                            discordant['2XSD'][chr1][orien].append([chr1,start,chr2,end,orien,type])
                        elif dis >= (mean+(stdev)) or dis <= (mean-(stdev)):
                            type = 'DISCOR_1XSD'
                            discordant['1XSD'][chr1][orien].append([chr1,start,chr2,end,orien,type])       
        return discordant
    


def get_read_pos(bam,DataName,rcsv = False,rjson = False):
    
    def get_json(df,bins):
        cols = df.columns.difference(bins)
        df2 = (df.groupby(bins)[cols] \
                      .apply(lambda x: x.to_dict('r')) \
                      .reset_index())
        return df2
    # Generate Dataframes and save as csv:
    parent_dir  = os.getcwd()
    outdir = os.path.join(parent_dir,DataName)
    os.makedirs(outdir, exist_ok=True)

    ## make logfile
    FORMAT = '%(asctime)s %(message)s'
    logfilename = os.path.join(outdir, '_'.join(['Logfile',str(time())]))   
    logging.basicConfig(format=FORMAT,
                            filename=logfilename,
                            filemode='w',
                            level=logging.INFO)
    t0 = time()
    
    bam = pysam.AlignmentFile(bam,'rb')
    binned = bin_reads(bam)
    # get the chromosome list
    chr_list = binned.chr_list

    # Produce Dataframes for different reads
    ## clipped reads:
    clipped = binned.get_clipped()

    ## discordant reads:
    mean, stdev = av_distance_reads(bam)
    binned = bin_reads(bam)

    discor = binned.get_discordant(mean = mean, stdev= stdev)

    ## split reads & split reds + mate:
    binned = bin_reads(bam)

    split, split_mate = binned.get_split()
    
    ## Double split reads:
    binned = bin_reads(bam)

    both = binned.get_both_split()
    
    bam.close()




    ## Dataframe clipped 
    clipped_df = get_dataframe(clipped,chr_list)
    if not clipped_df.empty:
        logging.info('Processed Dataframe positions: %s reads' % 'Clipped')
        print('Processed clipped reads data in dataframe')
        outfile = os.path.join(outdir,'clipped_pos.csv')

        if rcsv:
            clipped_df.to_csv(outfile,index=False)
        if rjson:
            clipped_df = to_json(clipped_df, ['chr1','read type','orientation'])
            clipped_df.to_json(outfile.replace('csv','json'),orient = 'records')
    else:
        print('Empty dataframe: Clipped reads')
    ## Dataframes discordant readpairs
    stdevs = discor.keys()
    for sd in stdevs:
        df = get_dataframe(discor[sd],chr_list)
        if not df.empty:
            logging.info('Processed Dataframe positions: %s reads, %s' % ('Discordant', sd))
            print('Processed discordant read (%s) data in dataframe' % sd)
            outfile = os.path.join(outdir,'discordant_{}_pos.csv'.format(sd))

            if rcsv:
                df.to_csv(outfile,index=False)
            if rjson:
                df = to_json(clipped_df, ['chr1','read type','orientation'])
                df.to_json(outfile.replace('csv','json'),orient = 'records')
                
        
        else:
            print('Empty dataframe: Discordant reads')
    ## Dataframes split reads and split + mate reads 
    split_df = get_dataframe(split,chr_list)
    split_mate_df = get_dataframe(split,chr_list)
    
    if not split_df.empty and not split_mate_df.empty:
        logging.info('Processed Dataframe positions: %s reads' % 'Split ')
        print('Processed split reads data in dataframe')
        outfile_s = os.path.join(outdir,'.'.join(['split_SA','csv']))
        outfile_sm = os.path.join(outdir,'.'.join(['split_mate_pos','csv']))
        
        if rcsv:
            split_df.to_csv(outfile_s, index = False)
            split_mate_df.to_csv(outfile_sm, index = False)
        if rjson:
            split_jdf = to_json(split_df, ['chr1','read type','orientation'])
            split_jdf.to_json(outfile_s.replace('csv','json'),orient = 'records')
                                  
            split_mate_df = to_json(split_mate_df, ['chr1','read type','orientation'])                      
            split_mate_df.to_json(outfile_sm.replace('csv,json'),orient = 'records')
    else:
        print('Empty dataframe: Split and Splitmate')

    ## make a dataframe for read double split
    both_split_df = get_dataframe(both,chr_list,bs=True)
    if not both_split_df.empty:
        logging.info('Processed Dataframe positions: %s reads' % 'Two-sides split')
        print('Processed two-sided split reads data in dataframe')
        
        outfile = os.path.join(outdir,'both_split_pos.csv')
                                  
        if rcsv:
            both_split_df.to_csv(outfile,index=False)
        if rjson:
            both_split_jdf = to_json(both_split_df, ['chr2_L','orientation'])
    else:
        print('Empty dataframe: Split both sides')
    
    logging.info('Time: read positions on BAM %s: %f' % \
                 (DataName, (time() - t0)))
def main():
    parser = argparse.ArgumentParser(description='Extract read positions (discordant/clipped/split)')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        help="Specify input file (BAM)")

    parser.add_argument('-DN',
                        '--DataName',
                        type=str,
                        help = "Specify name data/sample")
 
    parser.add_argument('-rcsv',
                        '--returnCSV',
                        type=bool,
                        default=False,
                        help = 'Specify whether you want the output to be in .csv format')
    parser.add_argument('-rjson',
                        '--returnJSON',
                        type=bool,
                        default=False,
                        help = 'Specify wheter you want .json formatted output')
    args = parser.parse_args()  
    
    
    get_read_pos(bam=args.bam,DataName = args.DataName, rcsv = args.returnCSV, rjson = args.returnJSON)
                          
    
if __name__ == '__main__':
       main()

