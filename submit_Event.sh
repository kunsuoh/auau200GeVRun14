#! /bin/sh
for filename in $(cat fileLists/Pico.list); do
   cp ./run.csh ./run_$filename.csh
 
   echo "root4star -l -b <<EOF">>run_$filename.csh
   echo -n ".x runPicoNpeEventMaker.C(\"fileLists/PicoList/picoList_">>run_$filename.csh
   echo -n $filename>>run_$filename.csh
   echo  ".list\",\"Out/$filename\")">>run_$filename.csh
   echo "EOF">>run_$filename.csh
   qsub -hard -l h_vmem=4G -l scratchfree=500,gscratchio=1 -o Log/job_$filename.log -e Log/job_$filename.err run_$filename.csh

   mv run_$filename.csh script/. 

done   

