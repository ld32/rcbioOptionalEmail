#!/bin/sh

for i in A B; do            
    
    u=university$i.txt   
    
    #step1, find John and Mike    
    grep -H John $u >>  John.txt; grep -H Mike $u >>  Mike.txt       
    
    #step2, find Nick and Julia
    grep -H Nick $u >>  Nick.txt; grep -H Julia $u >>  Julia.txt
    
done
    
# step3 merge          
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt

