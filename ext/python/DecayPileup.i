%{
#include "../include/DecayPileup.h"
%}

%include "../include/DecayPileup.h"

%attribute(DecayPileup, bool, use_pile_up_correction, get_use_pile_up_correction, set_use_pile_up_correction);
%attribute(DecayPileup, double, instrument_dead_time, get_instrument_dead_time, set_instrument_dead_time);
%attribute(DecayPileup, double, repetition_rate, get_repetition_rate, set_repetition_rate);
%attributestring(DecayPileup, std::string, pile_up_model, get_pile_up_model, set_pile_up_model);
