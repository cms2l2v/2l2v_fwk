#!/bin/bash

#for i in {1..10000}

USER=""
if [ $# -eq 0 ]
  then
    USER=`whoami`
  else
    USER=${1}
fi


LIP=true
if [ `which qstat ` ]
then
    LIP=false
fi

while true ;
	do
         if [ ${LIP} ]
	 then
	     if [ "${USER}" = "stalker" ]
	     then
		 echo "You are running in stalker mode, jajajajajajajajajaja"
		 echo "Vischia( run: " `qstat -u vischia | grep " r " | wc -l` ", wait: " `qstat -u vischia | grep " qw " | wc -l`  "). Cris( run: " `qstat -u cbeiraod | grep " r " | wc -l` ", wait: " `qstat -u cbeiraod | grep " qw " | wc -l` "). Bargassa( run: " `qstat -u bargassa | grep " r " | wc -l` ", wait: " `qstat -u bargassa | grep " qw " | wc -l` ")"
	     else
		 echo "running: " `qstat -u ${USER} | grep " r " | wc -l` "    waiting: " `qstat -u ${USER} | grep " qw " | wc -l` 
	     fi
	 else
	     echo "RUNNING: `bjobs | grep RUN | wc -l`            PENDING: `bjobs | grep PEND | wc -l`"
	 fi
	 sleep 5	
	 done
