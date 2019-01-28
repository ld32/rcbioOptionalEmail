#!/bin/sh

for i in A B; do            
    
    u=university$i.txt    
    x=abc
    #the following is added, meaning: start step1, depends on nothing, step name is "find1", want to copy u to /tmp        
    #@1,0,find1,us,,u.x     
    grep -H John $u >>  John.txt; grep -H Mike $u >>  Mike.txt        
  
    #the following is added, meaning: start step2, depends on nothing, step name is "find2", want to copy u to /tmp        
    #@2,0,find2,u,sbatch -p short -n 1 -t 50:0
    grep -H Nick $u >>  Nick.txt; grep -H Julia $u >>  Julia.txt

done 
    
#the following is added, meaning: start step3, depends on step1 and 2, step name is "merge"
#@3,1.2,merge          
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt

