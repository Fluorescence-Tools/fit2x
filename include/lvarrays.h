//!  Structures and functions used for LabView interface
/*!
 * fit2x was originally developed as a C backend for LabView software. Therefore,
 * the interface with fit2x uses structures that can be accessed by Labview. In
 * order to make an interfacing with Python and other languages possible there is
 * a this files defines a set of functions that facilitate the creation of the
 * LabView structures.
*/
#ifndef FIT2X_LVARRAYS_H
#define FIT2X_LVARRAYS_H

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>


typedef struct {
    int length;
    int* data;
} LVI32Array;


typedef struct {
    int length;
    double* data;
} LVDoubleArray;


typedef struct {
    LVI32Array** expdata;
    LVDoubleArray** irf;
    LVDoubleArray** bg;	// must be normalized outside!!!
    double dt;
    LVDoubleArray** corrections;
    LVDoubleArray** M;
} MParam;


/*!
 *
 * @param len
 * @return
 */
LVI32Array *CreateLVI32Array(size_t len);

/*!
 *
 * @param len
 * @return
 */
LVDoubleArray *CreateLVDoubleArray(size_t len);


MParam *CreateMParam(
    double dt=1.0,
    std::vector<double> corrections = std::vector<double>(),
    std::vector<double> irf = std::vector<double>(),
    std::vector<double> background = std::vector<double>(),
    std::vector<int> data = std::vector<int>()
);

#endif //FIT2X_LVARRAYS_H
