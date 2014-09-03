#!/bin/bash
/bin/rm __init__.pyc _tm_score_cpp.so tm_score_cpp.pyc tm_score_cpp.py tm_score_cpp_wrap.cpp tm_score_cpp_wrap.o
swig -c++ -python -o tm_score_cpp_wrap.cpp tm_score.i
g++ -fPIC -Wno-deprecated -c pdbClasses2.cpp
gcc -fPIC -c tm_score_cpp_wrap.cpp -o tm_score_cpp_wrap.o -I/usr/include/python2.3/
g++ -shared tm_score_cpp_wrap.o pdbClasses2.o -o _tm_score_cpp.so
