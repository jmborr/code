#!/bin/bash

mencoder $1  -oac lavc -fps 24 -ovc lavc quality=10 -lavcopts abitrate=160:vcodec=mpeg2video -of mpeg -o $1.mpg
