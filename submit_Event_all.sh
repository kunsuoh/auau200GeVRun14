#! /bin/sh

for filename in $(cat fileLists/List); do
	./submit_Event.sh $filename
done   


