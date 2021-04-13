%{
#include "../include/DecayBackground.h"
%}

%include "../include/DecayBackground.h"

%attribute(
    DecayBackground,
    double, constant_offset,
    get_constant_offset, set_constant_offset
);

