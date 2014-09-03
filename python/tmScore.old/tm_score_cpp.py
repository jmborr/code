# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _tm_score_cpp

def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types



getrmsd = _tm_score_cpp.getrmsd

getCAcoords2 = _tm_score_cpp.getCAcoords2

getd0 = _tm_score_cpp.getd0

coverage = _tm_score_cpp.coverage

coverageII = _tm_score_cpp.coverageII

rmsd_cov = _tm_score_cpp.rmsd_cov

gettm = _tm_score_cpp.gettm

supCAs = _tm_score_cpp.supCAs

print_array = _tm_score_cpp.print_array

allocD2 = _tm_score_cpp.allocD2

copyD2 = _tm_score_cpp.copyD2

allocD1 = _tm_score_cpp.allocD1

allocI1 = _tm_score_cpp.allocI1

lowerNotAligned = _tm_score_cpp.lowerNotAligned

