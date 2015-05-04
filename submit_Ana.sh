#! /bin/sh
for filename in $(cat fileLists/Analist); do
   cp ./run.csh ./run_Ana_$filename.csh
 
   echo "root4star -l -b <<EOF">>run_Ana_$filename.csh
   echo -n ".x runPicoNpeAnaMaker.C(\"">>run_Ana_$filename.csh
   echo -n $filename>>run_Ana_$filename.csh
   echo  "\")">>run_Ana_$filename.csh
   echo "EOF">>run_Ana_$filename.csh
#   qsub -hard -l h_vmem=4G -l scratchfree=500,gscratchio=1,projectio=1 -o Log/job_Ana_$filename.log -e Log/job_Ana_$filename.err run_Ana_$filename.csh
   echo ./run_Ana_$filename.csh

   mv run_Ana_$filename.csh script/. 

done   

