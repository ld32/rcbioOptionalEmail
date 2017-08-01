#!/usr/bin/env python 
import xlrd, sys

import os.path

if len(sys.argv) != 2:
        print "Usage: " + os.path.basename(__file__) + "  <sampleList.xlsx, required.>"
        print "For example, please check online for Excel table format:"
        print "https://wiki.rc.hms.harvard.edu/display/O2/Build+Folder+Structures+From+Sample+Sheet+for+Rcbio+NGS+Workflows"
        exit()

sampleSheet=sys.argv[1]

if not os.path.isfile(sampleSheet):
        print "Sample sheet not exist: " + sampleSheet
        exit()

if os.path.exists("group1"):
    print "Folder exists. please remove or rename it: group1"
    exit()

if os.path.exists("group2"):
    print "Folder exists. please remove or rename it: group2"
    exit()
    
if os.path.exists("group3"):
    print "Folder exists. please remove or rename it: group3"
    exit()    

wb = xlrd.open_workbook(os.path.join(sampleSheet))
wb.sheet_names()
sh = wb.sheet_by_index(0)

for i in xrange(sh.nrows):
        if i==0: # ignore header 
                continue
        
        print "working on row: " + str(i)
        group = sh.cell(i,0).value
                   
        if type(group) is float:        
            group=str(group).split(".")[0]
        
        if ' ' in group: 
            print "Group name can not have space: '%s'" % group
            exit()
            
        sample = sh.cell(i,1).value
        if type(sample) is float:        
            sample=str(sample).split(".")[0]
     
        if ' ' in sample: 
            print "Sample name can not have space: '%s'" % sample
            exit()
   
        library = sh.cell(i,2).value
        if type(library) is float:        
            library=str(library).split(".")[0]
        
        if ' ' in library: 
            print "Library name can not have space: '%s'" % library
            exit()
            
        path = sh.cell(i,3).value
        if type(path) is float:        
            path=str(path).split(".")[0]
        
        if ' ' in path: 
            print "Path name can not have space: '%s'" % path
            exit()
            
        read1 = os.path.expanduser(path + "/" + sh.cell(i,4).value)
        
        if ' ' in read1: 
            print "Read1 name can not have space: '%s'" % read1
            exit()
        
        if os.path.isfile(read1): 
                if not os.path.exists("group" + group + "/" + sample + "/" ): 
                        os.makedirs("group" + group + "/" + sample + "/" )
                if read1.endswith(".gz"):  
                        os.symlink(read1, "group" + group + "/" + sample + "/" + library +"_1.fq.gz" )
                elif read1.endswith(".bz2"):
                        os.symlink(read1, "group" + group + "/" + sample + "/" + library +"_1.fq.bz2" )
                else:
                        os.symlink(read1, "group" + group + "/" + sample + "/" + library +"_1.fq" )
        else: 
                print "row " + str(i) + ": read1 file not exist: " + read1 
                exit()
        if sh.cell(i,5).value != "":  # if we have read2
                read2 = os.path.expanduser(path + "/" + sh.cell(i,5).value)
                if ' ' in read2: 
                    print "Read2 name can not have space: '%s'" % read2
                    exit()
                
                if os.path.isfile(read2): 
                        if read2.endswith(".gz"):  
                             os.symlink(read2, "group" + group + "/" + sample + "/" + library +"_2.fq.gz" )
                        elif read2.endswith(".bz2"):
                             os.symlink(read2, "group" + group + "/" + sample + "/" + library +"_2.fq.bz2" )
                        else:
                             os.symlink(read2, "group" + group + "/" + sample + "/" + library +"_2.fq" )
                else: 
                        print "row " + str(i) + ": read2 file not exist: " + read2 
                        exit()
                        
if not os.path.exists("group2"):
    print "Something is wrong. Please make sure you have at least two groups."
    exit()

print "Done"
