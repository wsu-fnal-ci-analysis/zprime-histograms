#!/bin/bash

rm -f att.txt

for file in `ls -trd1 jobs/sleep* | head -500`; do
 echo $file
 echo "rm -f $file" >> att.txt
done

for file in `ls -trd1 jobs/*.log | head -500`; do
 echo $file
 echo "rm -f $file" >> att.txt
done

bash att.txt
