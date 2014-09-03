#!/bin/bash

PROJD=/projects/development/VisionElasticFit

#Find peak heights
python VisionElasticFit.py misc --kwargs "job=find peak heights,parFile=$PROJD/test/b2blb/parameters.2.dat,outFile=$PROJD/test/b2blb/peakHeights.2.dat"

#create HTML report
python VisionElasticFit.py misc --kwargs "job=create HTML report,outFile=$PROJD/test/report.html"

#Read file pixel_???.dat and find the position of the maximum for the fit
python VisionElasticFit.py misc --kwargs "job=find position of the maximum for the fit,outFile=$PROJD/test/b2blb/maximum.dat"