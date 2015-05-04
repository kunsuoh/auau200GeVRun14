#! /bin/sh
count=0
for filename in $(cat fileLists/AnaList); do
   count2=$(($count/100))
   if [ -f Out/$filename.picoNpe.root ]; then
      if [ $(($count%100)) -eq "0" ]; then 
         cp ./run.csh ./run_Ana_$count2.csh
      fi
      echo "root4star -l -b <<EOF">>run_Ana_$count2.csh
      echo -n ".x runPicoNpeAnaMaker.C(\"">>run_Ana_$count2.csh
      echo -n $filename>>run_Ana_$count2.csh
      echo  "\")">>run_Ana_$count2.csh
      echo "EOF">>run_Ana_$count2.csh
      #qsub -hard -l h_vmem=4G -l scratchfree=500,gscratchio=1,projectio=1 -o Log/job_Ana_$filename.log -e Log/job_Ana_$filename.err run_Ana_$count2.csh
      
      if [ $(($count%100)) -eq "99" ]; then 
         ./run_Ana_$count2.csh
      fi
      mv run_Ana_$count2.csh script/. 
      count=$(($count+1))
      echo $count
   else 
      echo "no input file: Out/$filename.picoNpe.root"
   fi
done 
./run_Ana_$count2.csh

