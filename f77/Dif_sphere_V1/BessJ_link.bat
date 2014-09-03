@echo off
g77 BessJ.for -o BessJ.exe -s -W -Wall -fno-backslash -O -Wl,--subsystem,console -lfit -lcpgplot -lpgplot -lGrWin -mwindows

