#! /bin/sh
count=0
for filename in $(cat fileLists/AnaList); do
   count2=$(($count/100))

   if [ -s $filename ]; then
      filesize=$(wc -c "Out/$filename.picoNpe.hists.root" | cut -f 1 -d ' ')
      if [ $filesize -gt 300 ]; then
         if [ $(($count%100)) -eq "0" ]; then 
            cp ./run.csh ./run_Ana_$count2.csh
            echo "Make run_Ana_$count2.csh"
         fi
         echo "root4star -l -b <<EOF">>run_Ana_$count2.csh
         echo -n ".x runPicoNpeAnaMaker.C(\"">>run_Ana_$count2.csh
         echo -n $filename>>run_Ana_$count2.csh
         echo  "\")">>run_Ana_$count2.csh
         echo "EOF">>run_Ana_$count2.csh
         if [ $(($count%100)) -eq "99" ]; then 
            echo "Submit run_Ana_$count2.csh" 
            #./run_Ana_$count2.csh
            qsub -hard -l h_vmem=4G -l scratchfree=500,gscratchio=1,projectio=1 -o Log/job_Ana_$count2.log -e Log/job_Ana_$count2.err run_Ana_$count2.csh
            mv run_Ana_$count2.csh script/. 
         fi
   
         count=$(($count+1))
         echo $count
      #else 
         #echo "no input file: Out/$filename.picoNpe.root"
      fi
   fi
done 
#./run_Ana_$count2.csh
qsub -hard -l h_vmem=4G -l scratchfree=500,gscratchio=1,projectio=1 -o Log/job_Ana_$count2.log -e Log/job_Ana_$count2.err run_Ana_$count2.csh
mv run_Ana_$count2.csh script/. 
