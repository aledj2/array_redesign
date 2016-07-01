'''
Created on 29 Apr 2016

@author: ajones7
'''
import os
import MySQLdb
import numpy as np
from string import rstrip

class bedfile():

    def __init__(self):
        # define parameters used when connecting to database
        self.host = "localhost"
        self.port = int(3306)
        self.username = "root"
        self.passwd = "mysql"
        self.database = "array"
        
        self.gene_coords="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Guys19_Gene_20100203.txt"
        self.gene_list=[]
        self.gene_dict={}
        self.probe_list=[]
        self.count_dict={}
        
        self.output_file="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\filtered_bedfile.txt"
        
    def filtergenecoords(self):
        genes=open(self.gene_coords,'r')
        for line in genes:
            splitline=line.split('\t')
            gene=splitline[3]
            self.gene_dict[gene]=line
        
    
        # read database and create a dictionary of isca probes
        db = MySQLdb.Connect(host=self.host, port=self.port,user=self.username, passwd=self.passwd, db=self.database)
        cursor = db.cursor()

        # sql statement
        genelist = "select distinct approvedsymbol from adele_genes"

        try:
            cursor.execute(genelist)
            genelist = cursor.fetchall()
        except MySQLdb.Error, e:
            db.rollback()
            print "fail - unable to get list of imported filenames"
            if e[0] != '###':
                raise
        finally:
            db.close()
        
        for i in genelist:
            genename=i[0]
            self.gene_list.append(genename)
        
        outputfile=open(self.output_file,'w')
        for gene in self.gene_dict:
            if gene in self.gene_list:
                outputfile.write(self.gene_dict[gene])
        outputfile.close()
            

class filter_above_bedfile():
    def __init__(self):
        self.targetgenesbedfile="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\filtered_bedfile.txt"
        self.list_of_poorly_covered_genes="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\genenames_max4_probes.txt"
        self.genes_with_no_coverage="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\genes_with_no_probes_bedfile.txt"
    
    def filter(self):
        filtered=open(self.targetgenesbedfile,'r')
        genes=open(self.list_of_poorly_covered_genes,'r')
        output=open(self.genes_with_no_coverage,'w')
        
        genes2cover=[]
        for line in genes:
            splitline=line.split('\t')
            gene=splitline[0].strip()
            genes2cover.append(gene)
        
        for line2 in filtered:
            splitline2=line2.split('\t')
            gene=splitline2[3].strip()
            if gene in genes2cover:
                output.write(line2)
            

        
        filtered.close()
        genes.close()
        output.close() 

if __name__=="__main__":
    #i=bedfile()
    #i.filtergenecoords()
    #print i.gene_dict
    #print i.count_dict
    j=filter_above_bedfile()
    j.filter()