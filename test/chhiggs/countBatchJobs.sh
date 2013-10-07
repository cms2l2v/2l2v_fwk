#!/bin/bash

#for i in {1..10000}
while true ;
	do
	 echo "RUNNING: `bjobs | grep RUN | wc -l`            PENDING: `bjobs | grep PEND | wc -l`"
		sleep 5	
	 done
