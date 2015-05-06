#!/bin/bash

#for i in {1..10000}

LIP=true
if [ `which qstat ` ]
then
    LIP=false
fi

while true ;
	do
         if [ ${LIP} ]
	 then
	     echo "running: " `qstat -u vischia | grep vischia | grep " r " | wc -l` "    total: " `qstat -u vischia | grep vischia | wc -l` 
	 else
	     echo "RUNNING: `bjobs | grep RUN | wc -l`            PENDING: `bjobs | grep PEND | wc -l`"
	 fi
	 sleep 5	
	 done
