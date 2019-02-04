#!/bin/sh

Usage="Usage: $0 \n  Note: this script will go through job name list in flag/alljobs.jid.first to see if the jobs finish successfully or not." 

[ -f flag/alljobs.jid.first ] || { echo Job list file not exist: flag/alljobs.jid.first; exit 1; }

names=`tail -n +2 flag/alljobs.jid.first | awk '{print $3}' | tr "\n" " "`

echo According to flag/alljobs.jid.first:
for name in $names; do       
   [ -f flag/$name.sh ] && [ -f flag/$name.success ] && echo Success! $name  || { echo "Failed!! $name #############" && [ -f flag/$name.out ] && echo -e "You can review the job log by:\nless flag/$name.out\n"; }
done 

echo 
echo According to flag/alljobs.jid: 
names=`tail -n +2 flag/alljobs.jid | awk '{print $3}' | tr "\n" " "`
for name in $names; do       
   [ -f flag/$name.sh ] && [ -f flag/$name.success ] && echo Success! $name  || { echo "Failed!! $name #############" && [ -f flag/$name.out ] && echo -e "You can review the job log by:\nless flag/$name.out\n"; }
done 

