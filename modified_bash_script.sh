#!/bin/sh            

#loopStart,i
for i in `ls -d condition*`; do               #first loop           
    cd $i            
    
    #loopStart,j
    for j in `ls -d sample*`; do              #second loop
        cd $j    
      
        #loopStart,l
        for l in `ls -f *.txt`; do             #third loop
            
            #@1,0,copy1: copy original file to .copy1
            cp $l $l.copy1                     #step1
            
            #@2,1,copy2: copy .copy1 to .copy2
            cp $l.copy1 $l.copy2               #step2  
        
        #loopEnd
        done    
        
        #@3,1,merge1: merge .copy1 to .copy1.merged  
        cat *.copy1 > $j.copy1.merged          #step3
    
        #@4,2,merge2: merge .copy2 to .copy2.merged  
        cat *.copy2 > $j.copy2.merged          #step4 
        cd ..    
 
    #loopEnd
    done            
    
    #@5,3.4,mergeall: merge everything together 
    cat */*merged > $i.all                     #step5
    
    cd ..

#loopEnd    
done           
