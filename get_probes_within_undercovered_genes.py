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
        
        self.gene_coords="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes\\genes_with_no_probes_bedfile.txt"
        self.gene_dict={}
        self.probe_list=[]
        self.count_dict={}
        
        #self.array_design="array_v3"
        self.array_design="baylor_array_design"
        #self.array_design="ISCA_array_design"
        
        self.output_file="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes\\probes_in_low_covered_genes_baylor.txt"
        #self.output_file="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes\\probes_in_low_covered_genes_aledv3.txt"
        #self.output_file="S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\Under_covered_genes\\probes_in_low_covered_genes_ISCA.txt"
        
        # sql statement
        self.probelist = "select distinct Probe_ID from "+self.array_design+" where Chromosome = %s and `Start` < %s and `stop` >%s"
        
        
    def count_probes_for_gene(self):
        # read database and create list of probes in adeles design 
        db = MySQLdb.Connect(host=self.host, port=self.port,user=self.username, passwd=self.passwd, db=self.database)
        cursor = db.cursor()
        list_of_adele_probes=[]
        # sql statement
        adeleprobelist = "select distinct Probe_ID from adele_probe_design"

        try:
            cursor.execute(adeleprobelist)
            adeleprobelist = cursor.fetchall()
        except MySQLdb.Error, e:
            db.rollback()
            print "fail - unable to get list of imported filenames"
            if e[0] != '###':
                raise
        finally:
            db.close()
        
        for i in adeleprobelist:
            probe_ID=i[0] 
            list_of_adele_probes.append(probe_ID)
        
        genes=open(self.gene_coords,'r')
        outputfile=open(self.output_file,'w')
        
        for line in genes:
            splitline=line.split('\t')
            chr=splitline[0]
            start=splitline[1]
            stop=splitline[2]
            gene=splitline[3]
            self.gene_dict[gene]=[chr,start,stop]
        
            # read database and get all probes within target genes
            db = MySQLdb.Connect(host=self.host, port=self.port,user=self.username, passwd=self.passwd, db=self.database)
            cursor = db.cursor()
            query="select distinct Probe_ID from "+self.array_design+" where Chromosome = '"+chr+"' and `Start` < "+stop+" and `stop` > " + start
            try:
                cursor.execute(query)
                probelist = cursor.fetchall()
            except MySQLdb.Error, e:
                db.rollback()
                print "fail - unable to get list of imported filenames"
                if e[0] != '###':
                    raise
            finally:
                db.close()

            
            #check if the probe is not already in the design.
            for i in probelist:
                #print i[0]
                probe_ID=i[0]
                if probe_ID not in list_of_adele_probes:
                    outputfile.write(gene+"\t"+probe_ID+"\n")
        outputfile.close()
        genes.close()
    
if __name__=="__main__":
    i=ISCA()
    i.count_probes_for_gene()
