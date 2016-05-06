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
        self.arrayfolder = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\GW_input_files"
        self.output_folder = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design"
        self.exclude = "S:\\Genetics_Data2\\Array\\Audits and Projects\\160421 ArrayCGH v4.0 Design\\ExcludeProbes.txt"
        # define parameters used when connecting to database
        self.host = "localhost"
        self.port = int(3306)
        self.username = "root"
        self.passwd = "mysql"
        self.database = "array"

        # table containing isca probes
        # "nrxn1" #, shox, eichler_selected_probes", "nrxn_17q11_2","baylor_genes_probes","cgh_par3","v2_41_replacements"
        self.isca_table = "isca_60k_targetted"

        # dict to store the probes and the log ratios
        self.dictionary = {}
        self.dictionary2 = {}

        # Create an array to store all the files in.
        self.chosenfiles = []
        self.number_of_arrays = 141

        self.exclude_probes = []

    def get_isca_probes(self):
        # read database and create a dictionary of isca probes
        db = MySQLdb.Connect(host=self.host, port=self.port,
                             user=self.username, passwd=self.passwd, db=self.database)
        cursor = db.cursor()

        # sql statement
        probelist = "select Probe_ID from " + self.isca_table

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

        exclude_file = open(self.exclude, 'r')

        for j in exclude_file:
            j = str(j.rstrip())
            self.exclude_probes.append(j)

        # print self.exclude_probes

        for i in probelist:
            probename = i[0]
            if probename not in self.exclude_probes:
                # print "excluded "+str(probename)
                self.dictionary[probename] = []

        # print self.dictionary

# open fe file and loop through storing the log ratio for each probe
    def read_fefiles(self):
        for afile in os.listdir(self.arrayfolder):
            print afile
            # if afile ==
            # "253174622077_S01_Guys121919_CGH_1100_Jul11_2_1_1.txt" or afile
            # =="253174622077_S01_Guys121919_CGH_1100_Jul11_2_1_2.txt":
            file2open = self.arrayfolder + "\\" + afile
            wholefile = open(file2open, 'r')
            for i, line in enumerate(wholefile):

                # splits the line on tab and appends this to a list
                splitfeatures = line.split('\t')
                #print splitfeatures
                probename = splitfeatures[0]
                coords = splitfeatures[1] + ":" + splitfeatures[2] + "-" + splitfeatures[3]
                probename_coord = probename + "^" + coords

                featurenumber = splitfeatures[4]
                logratios = []

                number_of_scores = len(splitfeatures)
                # print probename, number_of_scores
                for j in range(5, number_of_scores):
                    logratios.append(splitfeatures[j])

                if probename in self.dictionary:
                    if probename_coord in self.dictionary2:
                        print "error\t"+str(probename_coord)
                    else:
                        self.dictionary2[probename_coord] = []
                        for k in logratios:
                            k = k.rstrip()
                            #print k
                            self.dictionary2[probename_coord].append(float(k))

            # print self.dictionary

# calculate the spread of each probe
    def stats(self):
        output_file = open(
            self.output_folder + "\\GWdata_" + self.isca_table + "_stats_excluded.txt", 'w+')
        for probe in self.dictionary2:
            splitprobename = probe.split("^")
            probename = splitprobename[0]
            coord = splitprobename[1]
            #print self.dictionary2[probe]
            if len(self.dictionary2[probe]) > 10:
            
                q75, q25 = np.percentile(self.dictionary2[probe], [75, 25])
                iqr = q75 - q25
                mean = np.mean(self.dictionary2[probe])
                median = np.median(self.dictionary2[probe])
                good_count = 0
                bad_count = 0
                for i in self.dictionary2[probe]:
                    # if -0.2 <= i < 0 or  0 < i <= 0.2:
                    if -0.2 <= i <= 0.2:
                        good_count = good_count + 1
                    else:
                        bad_count = bad_count + 1
    
                npassprobes = good_count + bad_count
                nfailedprobes = self.number_of_arrays - npassprobes
                percent_outside = np.divide(
                    float(bad_count), self.number_of_arrays) * 100
                percent_of_probes_outside_excl_failedprobes = np.divide(
                    float(bad_count), npassprobes) * 100
    
                good_count2 = 0
                bad_count2 = 0
                for i in self.dictionary2[probe]:
                    # if -0.35 <= i <= 0 or 0 < i <= 0.35:
                    if -0.35 <= i <= 0.35:
                        good_count2 = good_count2 + 1
                    else:
                        bad_count2 = bad_count2 + 1
    
                npassprobes2 = good_count2 + bad_count2
    
                percent_outside2 = np.divide(
                    float(bad_count2), self.number_of_arrays) * 100
                percent_of_probes_outside2_excl_failedprobes = np.divide(
                    float(bad_count2), npassprobes) * 100
    
                # print str(probe) + "\t bad count \t"+str(bad_count)+ "\t good
                # count \t" + str(good_count)+"\t % \t"+str(percent_outside)
                output_file.write(probename + "\t" + coord + "\tnumber of probes with no score\t" + str(nfailedprobes) + "\t% outside -0.2 and 0.2 (out of 141)\t" + str(percent_outside) + "\t% outside -0.35 and 0.35(out of 141)\t" + str(percent_outside2) + "%\t IQR \t" + str(iqr) + "\t q25 \t" + str(q25) + "\t q75 \t" + str(q75) + "\t mean \t" + str(mean) + "\t median \t" + str(median) + "\tscores\t" + str(self.dictionary2[probe]) + "\n")
        output_file.close()
if __name__ == "__main__":
    a = ISCA()
    a.get_isca_probes()
    a.read_fefiles()
    a.stats()
