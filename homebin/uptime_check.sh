#!/bin/bash

#Indefinitely sleep until machine uptime below certain threshold

#read cut-off load
co=$1
if [ "$co" = "" ];then
  co=3 #default cut-off load
fi

load=10 #initialize with big load

while [ "$load" -ge "$co" ];do
  N=$[ ( $RANDOM % 100 ) + 30 ]  #a random number between 30 and 130
  sleep $N
  load=`uptime |tr -s ' '|cut -d ' ' -f 11|cut -d '.' -f 1`
  load2=`uptime |tr -s ' '|cut -d ' ' -f 12|cut -d '.' -f 1`
  if [ "$load2" -gt "$load" ];then
    load=$load2
  fi
done
echo "load is $load. Cut-off is $co"
sleep 1s
