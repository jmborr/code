@echo off
g77 diffusion_sphere.for splt_pg_gl.for fun.for lspec.for -o diffusion_sphere.exe -s -W -Wall -fno-backslash -O -Wl,--subsystem,console -lfit -lcpgplot -lpgplot -lGrWin -mwindows

