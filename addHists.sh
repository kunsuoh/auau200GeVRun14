#! /bin/sh
count=0
rm hists_tmp.root hists.root

for filename in $(cat fileLists/AnaList); do
	if [ -f Out/$filename.picoNpe.root ]; then
    	filesize=$(wc -c "Out/$filename.picoNpe.hists.root" | cut -f 1 -d ' ')
		if [ $filesize -gt 300 ]; then
		echo "Adding Out/$filename.picoNpe.hists.root"
	  		if [ $count -eq 0 ]; then 
				hadd hists_tmp.root Out/$filename.picoNpe.hists.root
				mv hists_tmp.root hists.root
			else 
				hadd hists_tmp.root Out/$filename.picoNpe.hists.root hists.root
				mv hists_tmp.root hists.root
			fi
			count=$(($count+1))
   		fi
   fi
done 

now="$(date '+%Y%m%d_%H%M%S')"
mv hists.root hists_$now.root
echo "Output histogram file: hists_$now.root"
