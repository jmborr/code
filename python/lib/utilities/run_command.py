#!/usr/bin/python
import sys
from inputArgs.inputArgs import inpHand
from utilities.cluster import run_command

ih=inpHand('Usage: run_command.py',
           ' -c _R_cmd command (enclose within double quotes)',
           ' -h __help this script runs the same command on all nodes'
           )
ih.parse(locals(),sys.argv)

run_command(cmd)
