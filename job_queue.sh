#!/usr/bin/bash

while ps -p 156460 > /dev/null 
do
	sleep 10
	echo "process running"
done
echo "last command running"
