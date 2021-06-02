%{
#include "DecayModifier.h"
%}

%attribute(DecayModifier, bool, active, is_active, set_active);
%attribute(DecayModifier, DecayCurve, data, get_data, set_data);

%include "DecayModifier.h"
