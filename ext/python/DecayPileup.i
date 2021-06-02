%{
#include "DecayPileup.h"
%}

%attribute(DecayPileup, double, instrument_dead_time, get_instrument_dead_time, set_instrument_dead_time);
%attribute(DecayPileup, double, repetition_rate, get_repetition_rate, set_repetition_rate);
%attribute(DecayPileup, DecayCurve, data, get_data, set_data);
%attributestring(DecayPileup, std::string, pile_up_model, get_pile_up_model, set_pile_up_model);

%include "DecayPileup.h"
