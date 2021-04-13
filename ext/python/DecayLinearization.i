%{
#include "../include/DecayLinearization.h"
%}

// Input arrays
%apply(double* IN_ARRAY1, int DIM1) {
    (double* linearization_table, int n_linearization_table)
}
%include "../include/DecayLinearization.h"

%attribute(DecayLinearization, bool, use_linearization, get_use_linearization, set_use_linearization);

