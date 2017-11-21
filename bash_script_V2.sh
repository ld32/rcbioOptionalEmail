#!/bin/sh

#the following is added, meaning: starting a loop and use university as looping variable    
#loopStart,i	
for i in A B; do            
    
    u=university$i.txt    
    
    #the following is added, meaning: start step1, depends on nothing, step name is "find1", want to copy u to /tmp        
    #@1,0,find1,u:     
    grep -H John $u >>  John.txt; grep -H Mike $u >>  Mike.txt        
  
    #the following is added, meaning: start step2, depends on nothing, step name is "find2", want to copy u to /tmp        
    #@2,0,find2,u:
    grep -H Nick $u >>  Nick.txt; grep -H Julia $u >>  Julia.txt

#the following is added, meaning: loop end here    
#loopEnd                     
done
    
#the following is added, meaning: start step3, depends on step1 and 2, step name is "merge"
#@3,1.2,merge:           
cat John.txt Mike.txt Nick.txt Julia.txt > all.txt

# command to create testing data
# echo -e "John Paul\nMike Smith\nNick Will\nJulia Johnson\nTom Jones"  >> universityA.txt
# echo -e "John Paul\nMike Smith\nNick Will\nJulia Johnson\nTom Jones"  >> universityB.txt
