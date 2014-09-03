@echo off
g77 dif_sphere_V1.for splt_pg_gl.for fun.for lspec.for -o dif_sphere_V1.exe -s -W -Wall -fno-backslash -O -Wl,--subsystem,console -lfit -lcpgplot -lpgplot -lGrWin -mwindows

