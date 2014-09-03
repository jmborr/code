@echo off
g77 check_bessel.for -o check_bessel.exe -s -W -Wall -fno-backslash -O -Wl,--subsystem,console -lfit -lcpgplot -lpgplot -lGrWin -mwindows

