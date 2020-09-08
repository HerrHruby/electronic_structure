//
// Created by Ian Wu on 22/08/2020.
//

#ifndef HF_CONTRACTED_INTEGRALS_H
#define HF_CONTRACTED_INTEGRALS_H

#include "basis_function.h"
#include "integrals.h"

double S(const Basis_Function &A, const Basis_Function &B);

double T(const Basis_Function &A, const Basis_Function &B);

double V(const Basis_Function &A, const Basis_Function &B, const std::valarray<double> &centre);

double TEI(const Basis_Function &A, const Basis_Function &B, const Basis_Function &C, const Basis_Function &D);


#endif //HF_CONTRACTED_INTEGRALS_H
