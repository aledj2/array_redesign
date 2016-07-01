'''
Created on 16 Jun 2016

@author: ajones7
'''

import MySQLdb

class read_tables():
    def __init__(self):
        # define parameters used when connecting to database
        self.host = "localhost"
        self.port = int(3306)
        self.username = "root"
        self.passwd = "mysql"
        self.database = "array"
        
        self.adele_probes="adele_probe_design"
        self.probe_seq="isca_v2_probe_seq"
        
        self.adele_probe_list=[]
        self.probe_w_seq=[]
        self.allprobesiwthseq=[]
        
    def probe_list(self):
        db = MySQLdb.Connect(host=self.host, port=self.port,
                             user=self.username, passwd=self.passwd, db=self.database)
        cursor = db.cursor()

        # sql statement
        probelist = "select Probe_ID from " + self.adele_probes

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
            self.adele_probe_list.append(i[0])
    
    def get_probe_seq(self):
        db = MySQLdb.Connect(host=self.host, port=self.port,
                             user=self.username, passwd=self.passwd, db=self.database)
        cursor = db.cursor()

        # sql statement
        probelist = "select * from " + self.probe_seq + " union select probe_ID, sequence from array_v3"

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
            self.probe_w_seq.append((i[0],i[1]))
    
    def create_txt_file(self):
        output_file = open("S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\probeID_with_seq.txt", 'w+')
        count=0
        for i in self.probe_w_seq:
            if i[0] in self.adele_probe_list:
                output_file.write(i[0]+"\t"+i[1]+"\n")
                count += 1
        output_file.close()
        
        count2=0
        for k in self.probe_w_seq:
            self.allprobesiwthseq.append(k[0])
        
        for j in self.adele_probe_list:
            if j in self.allprobesiwthseq:
                count2+=1 
            else:
                print j
        
        print "len adele probes="+str(len(self.adele_probe_list))
        print "number of probes with seq:"+str(count)
        print "probes without a seq = "+str(len(self.adele_probe_list)-count2)
        
if __name__ == "__main__":
    a=read_tables()
    a.probe_list()
    #print a.adele_probe_list
    a.get_probe_seq()
    #print a.probe_w_seq
    a.create_txt_file()