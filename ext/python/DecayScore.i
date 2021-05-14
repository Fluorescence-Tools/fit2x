%feature("kwargs", 1);
%{
#include "../include/DecayScore.h"
%}

%include "../include/DecayScore.h"
%shared_ptr(DecayScale);

%attribute(DecayScore, int, score_range_stop, get_score_range_stop, set_score_range_stop);
%attribute(DecayScore, int, score_range_start, get_score_range_start, set_score_range_start);
%attributestring(DecayScore, std::string, score_type, get_score_type, set_score_type);
%attribute(DecayScore, double, score, evaluate);
