#! /bin/sh
# 1st input parameter : run id (DDD_NNNNNNNN)


for filename in $(cat fileLists/List); do
	./submit_Event.sh $filename
done   


