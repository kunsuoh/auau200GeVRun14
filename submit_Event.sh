#! /bin/sh
# 1st input parameter : DDD_NNNNNNNN

count=0
cp ./run.csh ./run_$1.csh
for filename in $(cat fileLists/Run14/AuAu/200GeV/picoLists/runs/picoList_$1.list); do
   echo $filename >> $1.list
   
   echo "root4star -l -b <<EOF">>run_$1.csh
   echo -n ".x runPicoNpeEventMaker.C(\"">>run_$1.csh
   echo -n $1.list>>run_$1.csh
   echo  "\",\"Out/$1_$count\")">>run_$1.csh
   echo "EOF">>run_$1.csh
   echo $1 $count
   count=$(($count+1))
done   
#qsub -hard -l h_vmem=4G -l scratchfree=500,gscratchio=1,projectio=1 -o Log/job_$1.log -e Log/job_$1.err run_$filename.csh
mv run_$1.csh script/
mv $1.list script/
