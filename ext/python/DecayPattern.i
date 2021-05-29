%{
#include "DecayPattern.h"
%}

%attribute(DecayPattern, DecayCurve*,pattern, get_pattern, set_pattern);
%attribute(DecayPattern, double, pattern_fraction, get_pattern_fraction, set_pattern_fraction);
%attribute(DecayPattern, double, constant_offset, get_constant_offset, set_constant_offset);

%include "DecayPattern.h"

