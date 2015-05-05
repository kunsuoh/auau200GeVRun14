#! /bin/sh
count=0
for filename in $(cat fileLists/AnaList); do
	if [ -f Out/$filename.picoNpe.root ]; then
    	filesize=$(wc -c "Out/$filename.picoNpe.hists.root" | cut -f 1 -d ' ')
		if [ $filesize -gt 300 ]; then
	  		if [ $count -eq 0 ]; then 
				hadd hists_tmp.root Out/$filename.picoNpe.hists.root
				mv hists_tmp.root hists.root
			else 
				hadd hists_tmp.root Out/$filename.picoNpe.hists.root hists.root
				mv hists_tmp.root hists.root
			fi
   		fi
   fi
done 
