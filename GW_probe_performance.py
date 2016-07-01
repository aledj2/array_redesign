'''
Created on 29 Apr 2016
@author: ajones7
'''
import os
import MySQLdb
import numpy as np
from string import rstrip


class ISCA():

    def __init__(self):
        self.arrayfolder = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes"
        self.output_folder = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes"
        self.baylorprobes = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes\\probes_in_low_covered_genes_baylor.txt"
        self.aledv3probes = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes\\probes_in_low_covered_genes_aledv3.txt"
        
        self.aledarraydata="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes\\AledArray_data.txt"
        self.baylorarraydata="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes\\Baylor_data.txt"
        

        # dict to store the probes and the log ratios
        self.aledv3dictionary = {}
        self.baylordictionary = {}

        # Create an array to store all the files in.
        self.chosenfiles = []
        self.number_of_arrays_aled = 0
        self.number_of_arrays_baylor = 0

        self.baylorprobenames = []
        self.aledprobenames = []
        
        self.aled_log_ratios=[]
        self.baylor_log_ratios=[]

    def get_probe_list(self):
        '''This module reads the list of probes already identified as being within a gene that requires improved coverage and creates a dict with each probe as a key.'''
        # open file
        baylor_probelist = open(self.baylorprobes, 'r')
        for j in baylor_probelist:
            #splot line
            j=j.split('\t')
            #remove newline characters
            probe = str(j[1].rstrip())
            #append to a list
            self.baylorprobenames.append(probe)
        baylor_probelist.close()
        
        #repeat for a second probe list
        aledv3_probelist = open(self.aledv3probes, 'r')
        for j in aledv3_probelist:
            j=j.split('\t')
            probe = str(j[1].rstrip())
            self.aledprobenames.append(probe)
        aledv3_probelist.close()
        
        #for each list. turn into a set (remove duplicates) and add to dictionary with empty list as key
        for i in set(self.baylorprobenames):
            probename = i
            self.baylordictionary[probename] = []
        
        for i in set(self.aledprobenames):
            probename = i
            self.aledv3dictionary[probename] = []

    def read_probe_scores(self):
        '''The raw data has been taken from workbench manually and copied into a large tab delimited txt file. This needs careful handling! 
        Loop through this and if probename is in target probe put all log ratios into the dict '''
        #open log ratio file
        arraydata=open(self.aledarraydata,'r')
        for i in arraydata:
            #split and capture probename
            i=i.split('\t')
            probename=i[0]
            #coords = i[1] + ":" + i[2] + "-" + i[3]
            #probename_coord = probename + "^" + coords
            
            # if the probe is in the list count how many arrays have a log ratio present.
            # if probe failed no log ratio - assumes that there is at least one probe which hasn't failed in every array 
            if probename in self.aledv3dictionary:
                # count the highest number of columns and remove 4 columns for name, chr, start, stop
                if len(i)>self.number_of_arrays_aled:
                    self.number_of_arrays_aled=len(i)-4
                # loop through each log ratio
                for j in range(4,len(i)):
                    #if its not blank/empty append to the list in the dict value
                    if i[j]:
                        self.aledv3dictionary[probename].append(float(i[j]))
        arraydata.close()
        
        #repeat 
        arraydata2=open(self.baylorarraydata,'r')
        for i in arraydata2:
            i=i.split('\t')
            probename=i[0]
            #coords = i[1] + ":" + i[2] + "-" + i[3]
            #probename_coord = probename + "^" + coords
            if probename in self.baylordictionary:
                if len(i)>self.number_of_arrays_baylor:
                    self.number_of_arrays_baylor=len(i)-4
                
                for j in range(4,len(i)):
                    self.baylordictionary[probename].append(float(i[j]))
        arraydata2.close()

    def stats(self):
        '''calculate the spread of each probe'''
        #open output file
        baylor_output_file = open(self.output_folder + "\\baylor_probe_performance.txt", 'w+')

        # for each probe in the dictionary
        for probe in self.baylordictionary:
            #if there are more than 10 log ratios (any with less than this prob have lots of fails)
            if len(self.baylordictionary[probe]) > 10:
                #calculate averages, iqr etc
                q75, q25 = np.percentile(self.baylordictionary[probe], [75, 25])
                iqr = q75 - q25
                mean = np.mean(self.baylordictionary[probe])
                median = np.median(self.baylordictionary[probe])
                #set up a count for probes within 1st target range
                good_count = 0
                bad_count = 0
                #count number of probes within +-0.2
                for i in self.baylordictionary[probe]:
                    if -0.2 <= i <= 0.2:
                        good_count = good_count + 1
                    else:
                        bad_count = bad_count + 1
                
                #number of probes which haven't failed
                npassprobes = good_count + bad_count
                # number of failed probes
                nfailedprobes = self.number_of_arrays_baylor - npassprobes
                # how many probes (incl failed probes) that are outside the range 
                percent_outside = np.divide(float(bad_count), self.number_of_arrays_baylor) * 100
                # how many probes which did not fail that are outside the range 
                percent_of_probes_outside_excl_failedprobes = np.divide(float(bad_count), npassprobes) * 100
                
                #repeat for second range
                good_count2 = 0
                bad_count2 = 0
                for i in self.baylordictionary[probe]:
                    if -0.35 <= i <= 0.35:
                        good_count2 = good_count2 + 1
                    else:
                        bad_count2 = bad_count2 + 1
    
                npassprobes2 = good_count2 + bad_count2
    
                percent_outside2 = np.divide(
                    float(bad_count2), self.number_of_arrays_baylor) * 100
                percent_of_probes_outside2_excl_failedprobes = np.divide(
                    float(bad_count2), npassprobes) * 100
                
                #write to output file    
                baylor_output_file.write(probe + "\tnumber of probes with no score\t" + str(nfailedprobes) + "\t% outside -0.2 and 0.2 (out of "+str(self.number_of_arrays_aled) +")\t" + str(percent_outside) + "\t% outside -0.35 and 0.35(out of "+str(self.number_of_arrays_aled) +")\t"+ str(percent_outside2) + "%\t IQR \t" + str(iqr) + "\t q25 \t" + str(q25) + "\t q75 \t" + str(q75) + "\t mean \t" + str(mean) + "\t median \t" + str(median) + "\tscores\t" + str(self.baylordictionary[probe]) + "\n")
        baylor_output_file.close()
        
        #repeat for second array design
        aled_output_file = open(self.output_folder + "\\aled_probe_performance.txt", 'w+')               
        for probe in self.aledv3dictionary:
            if len(self.aledv3dictionary[probe]) > 10:
                q75, q25 = np.percentile(self.aledv3dictionary[probe], [75, 25])
                iqr = q75 - q25
                mean = np.mean(self.aledv3dictionary[probe])
                median = np.median(self.aledv3dictionary[probe])
                good_count = 0
                bad_count = 0
                for i in self.aledv3dictionary[probe]:
                    # if -0.2 <= i < 0 or  0 < i <= 0.2:
                    if -0.2 <= i <= 0.2:
                        good_count = good_count + 1
                    else:
                        bad_count = bad_count + 1
     
                npassprobes = good_count + bad_count
                nfailedprobes = self.number_of_arrays_aled - npassprobes
                percent_outside = np.divide(
                    float(bad_count), self.number_of_arrays_aled) * 100
                percent_of_probes_outside_excl_failedprobes = np.divide(
                    float(bad_count), npassprobes) * 100
     
                good_count2 = 0
                bad_count2 = 0
                for i in self.aledv3dictionary[probe]:
                    # if -0.35 <= i <= 0 or 0 < i <= 0.35:
                    if -0.35 <= i <= 0.35:
                        good_count2 = good_count2 + 1
                    else:
                        bad_count2 = bad_count2 + 1
     
                npassprobes2 = good_count2 + bad_count2
     
                percent_outside2 = np.divide(
                    float(bad_count2), self.number_of_arrays_aled) * 100
                percent_of_probes_outside2_excl_failedprobes = np.divide(
                    float(bad_count2), npassprobes) * 100
     
                # print str(probe) + "\t bad count \t"+str(bad_count)+ "\t good
                # count \t" + str(good_count)+"\t % \t"+str(percent_outside)
                aled_output_file.write(probe + "\tnumber of probes with no score\t" + str(nfailedprobes) + "\t% outside -0.2 and 0.2 (out of "+str(self.number_of_arrays_aled) +")\t" + str(percent_outside) + "\t% outside -0.35 and 0.35(out of "+str(self.number_of_arrays_aled) +")\t"+ str(percent_outside2) + "%\t IQR \t" + str(iqr) + "\t q25 \t" + str(q25) + "\t q75 \t" + str(q75) + "\t mean \t" + str(mean) + "\t median \t" + str(median) + "\tscores\t" + str(self.aledv3dictionary[probe]) + "\n")
        aled_output_file.close()
        
        
        
if __name__ == "__main__":
    a = ISCA()
    a.get_probe_list()
    a.read_probe_scores()
    a.stats()
    
