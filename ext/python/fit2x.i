%module(directors="1") fit2x
%feature("kwargs", 1);
%feature("director");
%feature("autodoc", "docstring");

%{
// This fixes numpy int casting to std::vector,int>
// (see: https://github.com/swig/swig/issues/888)
#define SWIG_PYTHON_CAST_MODE
// This is needed for numpy as you need SWIG_FILE_WITH_INIT
#define SWIG_FILE_WITH_INIT
#include <assert.h>
#include "fits2x.h"
%}

%include "documentation.i"

#if VERBOSE_FIT2X
// Warning 509: Overloaded method ignored,
// Warning 511: Ignore overloaded functions
// Warning 401: Nothing known about base class
// Warning 453: No typemaps are defined. (Some typemaps are defined and not used)
// Warning 319: No access specifier given for base class (std::enable_shared_from_this)
// Warning 362: operator= ignored
#pragma SWIG nowarn=319,362,401,453,509,511
#endif

%include "misc_types.i"

// Python code that should be included at the begining (import, Base class, etc)
%pythoncode "../ext/python/fit2x.py"

// shared_prt
%shared_ptr(TTTR);

// for fit23, fit24, etc.
%apply (int DIM1, double* INPLACE_ARRAY1) {(int len1, double* x)}
%apply (int DIM1, short* IN_ARRAY1) {(int len2, short* fixed)}

/* Convolution and LabView interface*/
%include "lvarray.i"
%include "fsconv.i"
%include "phasor.i"

/* Fits */
%include "fit23.i"
%include "fit24.i"
%include "fit25.i"
%include "fit26.i"

/* Decay */
%include "DecayCurve.i"
%include "DecayRange.i"
%include "DecayLifetimeHandler.i"

%include "DecayModifier.i"
%include "DecayConvolution.i"
%include "DecayPattern.i"
%include "DecayPileup.i"
%include "DecayLinearization.i"
%include "DecayScale.i"

%include "DecayScore.i"

%include "fits2x.h"

%pythoncode "../ext/python/ParameterDirector.py"
%pythoncode "../ext/python/Decay.py"
