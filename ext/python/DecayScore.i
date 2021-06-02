%{
#include "DecayScore.h"
%}

%attributestring(DecayScore, std::string, score_type, get_score_type, set_score_type);
%attribute(DecayScore, double, score, score);
%attribute(DecayScore, DecayCurve, model, get_model, set_model);
%attribute(DecayScore, DecayCurve, data, get_data, set_data);
%attribute2(DecayScore, std::vector<double>, weighted_residuals, get_weighted_residuals);

%include "DecayScore.h"
