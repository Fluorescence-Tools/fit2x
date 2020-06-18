#include "lvarrays.h"


LVI32Array *CreateLVI32Array(size_t len) {
    auto *ts = (LVI32Array *) calloc(1, sizeof(LVI32Array));
    ts->data = (int *) calloc(len, sizeof(int));
    ts->length = len;
    return ts;
}


LVDoubleArray *CreateLVDoubleArray(size_t len) {
    auto *ts = (LVDoubleArray *) calloc(1, sizeof(LVDoubleArray));
    ts->data = (double *) calloc(len, sizeof(double));
    ts->length = len;
    return ts;
}


MParam *CreateMParam(
        double dt,
        std::vector<double> corrections,
        std::vector<double> irf,
        std::vector<double> background,
        std::vector<int> data
) {
    auto ts = (MParam *) calloc(1, sizeof(MParam));
    ts->irf = (LVDoubleArray **) calloc(1, sizeof(LVDoubleArray *));
    ts->bg = (LVDoubleArray **) calloc(1, sizeof(LVDoubleArray *));
    ts->corrections = (LVDoubleArray **) calloc(1, sizeof(LVDoubleArray *));
    ts->M = (LVDoubleArray **) calloc(1, sizeof(LVDoubleArray *));
    ts->expdata = (LVI32Array **) calloc(1, sizeof(LVI32Array *));

    int n_channels = std::max(
        {irf.size() / 2,
         background.size() / 2,
         data.size() / 2}
    );
    // all channel numbers are multiplied by two for jordi format
    *(ts->irf) = CreateLVDoubleArray((size_t) 2 * n_channels);
    *(ts->expdata) = CreateLVI32Array((size_t) 2 * n_channels);
    *(ts->bg) = CreateLVDoubleArray((size_t) 2 * n_channels);
    *(ts->corrections) = CreateLVDoubleArray(corrections.size());
    *(ts->M) = CreateLVDoubleArray(2 * n_channels);
    ts->dt = dt;

    for (int i = 0; i < irf.size(); i++) (*ts->irf)->data[i] = irf[i];
    for (int i = 0; i < background.size(); i++) (*ts->bg)->data[i] = background[i];
    for (int i = 0; i < corrections.size(); i++) (*ts->corrections)->data[i] = corrections[i];
    for (int i = 0; i < data.size(); i++) (*ts->expdata)->data[i] = data[i];

    if (irf.size() != background.size()) {
        std::cerr << "WARNING: length of background pattern and IRF differ." << std::endl;
    }
    return ts;
}
