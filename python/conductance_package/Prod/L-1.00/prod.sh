#!/bin/bash
for script in  0001_0100.sub 0101_0200.sub 0201_0300.sub 0301_0400.sub 0401_0500.sub 0501_0600.sub 0601_0700.sub 0701_0800.sub 0801_0900.sub 0901_1000.sub 1001_1100.sub 1101_1200.sub 1201_1300.sub 1301_1400.sub 1401_1500.sub 1501_1600.sub 1601_1700.sub 1701_1800.sub 1801_1900.sub 1901_2000.sub;do
  if [ $script == 0001_0100.sub ];then
    PBS_JOBID=$(qsub $script)
    echo "qsub $script"
  else
    PBS_JOBID=$(qsub -W depend=afterok:${PBS_JOBID} $script)
    echo "qsub -W depend=afterok:${PBS_JOBID} $run.pbs"
  fi
  sleep 3s
done
