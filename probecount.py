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
        # define parameters used when connecting to database
        self.host = "localhost"
        self.port = int(3306)
        self.username = "root"
        self.passwd = "mysql"
        self.database = "array"
        
        self.gene_coords="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\filtered_bedfile.txt"
        self.gene_dict={}
        self.probe_list=[]
        self.count_dict={}
        
        self.output_file="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\probecount.txt"
        
    def count_probes_for_gene(self):
        genes=open(self.gene_coords,'r')
        for line in genes:
            splitline=line.split('\t')
            chr=splitline[0]
            start=splitline[1]
            stop=splitline[2]
            gene=splitline[3]
            self.gene_dict[gene]=[chr,start,stop]
        
    
        # read database and create a dictionary of isca probes
        db = MySQLdb.Connect(host=self.host, port=self.port,user=self.username, passwd=self.passwd, db=self.database)
        cursor = db.cursor()

        # sql statement
        probelist = "select distinct chromosome,`Start`,`Stop` from adele_probe_design"

        try:
            cursor.execute(probelist)
            probelist = cursor.fetchall()
        except MySQLdb.Error, e:
            db.rollback()
            print "fail - unable to get list of imported filenames"
            if e[0] != '###':
                raise
        finally:
            db.close()
        
        for i in probelist:
            chr=i[0]
            start=i[1]
            stop=i[2]
            if chr=='Chr23':
                chr='ChrX'
            elif chr=='Chr24':
                chr='ChrY'
            self.probe_list.append([chr,start,stop])
        
        for gene in self.gene_dict:
            gene_chr=self.gene_dict[gene][0].lower()
            gene_start=int(self.gene_dict[gene][1])
            gene_end=int(self.gene_dict[gene][2])
            #print gene_end
            count=0
            for probe in self.probe_list:
                probe_chr=probe[0].lower()
                probe_start=int(probe[1])
                probe_end=int(probe[2])
                #print gene_chr+":"+str(gene_start)+"-"+str(gene_end)
                #print probe_chr+":"+str(probe_start)+"-"+str(probe_end)
                
                if probe_chr==gene_chr and probe_end>gene_start and probe_start<gene_end:
                    #print "match"
                    count+=1
            self.count_dict[gene]=count
        
        outputfile=open(self.output_file,'w')
        for i in self.count_dict:
            #print i,self.count_dict[i]
            
            outputfile.write(i+"\t"+str(self.count_dict[i])+"\n")
        outputfile.close()
            
    

if __name__=="__main__":
    i=ISCA()
    i.count_probes_for_gene()
    #print i.gene_dict
    #print i.count_dict