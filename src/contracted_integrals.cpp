//
// Created by Ian Wu on 22/08/2020.
//

#include "../include/contracted_integrals.h"
#include <iostream>


double S(const Basis_Function &A, const Basis_Function &B) {
    double s = 0.0;
    for (int i = 0; i < A.getCoefs().size(); i++) {
        for (int j = 0; j < B.getCoefs().size(); j++) {
            s +=  A.getNorm()[i] * B.getNorm()[j] * A.getCoefs()[i] * B.getCoefs()[j] *
                    overlap(A.getExps()[i], A.getShell(), A.getOrigin(),
                    B.getExps()[j], B.getShell(), B.getOrigin());
        }
    }
    return s;
}

double T(const Basis_Function &A, const Basis_Function &B) {
    double t = 0.0;
    for (std::size_t i = 0; i < A.getCoefs().size(); i++) {
        for (std::size_t j = 0; j < B.getCoefs().size(); j++) {
            t += A.getNorm()[i] * B.getNorm()[j] * A.getCoefs()[i] * B.getCoefs()[j] *
                    kinetic(A.getExps()[i], A.getShell(), A.getOrigin(),
                         B.getExps()[j], B.getShell(), B.getOrigin());

        }
    }
    return t;
}

double V(const Basis_Function &A, const Basis_Function &B, const std::valarray<double> &centre) {
    double v = 0.0;
    for (std::size_t i = 0; i < A.getCoefs().size(); i++) {
        for (std::size_t j = 0; j < B.getCoefs().size(); j++) {
            v += A.getNorm()[i] * B.getNorm()[j] * A.getCoefs()[i] * B.getCoefs()[j] *
                 coulomb(A.getExps()[i], A.getShell(), A.getOrigin(),
                         B.getExps()[j], B.getShell(), B.getOrigin(), centre);
        }
    }
    return v;
}

double TEI(const Basis_Function &A, const Basis_Function &B, const Basis_Function &C, const Basis_Function &D) {
    long double tei = 0.0;
    for (std::size_t i = 0; i < A.getCoefs().size(); ++i) {
        for (std::size_t j = 0; j < B.getCoefs().size(); ++j) {
            for (std::size_t k = 0; k < C.getCoefs().size(); ++k) {
                for (std::size_t l = 0; l < D.getCoefs().size(); ++l) {
                    tei += A.getNorm()[i] * B.getNorm()[j] * C.getNorm()[k] * D.getNorm()[l] *
                            A.getCoefs()[i] * B.getCoefs()[j] * C.getCoefs()[k] * D.getCoefs()[l] *
                           two_electron(A.getExps()[i], A.getShell(), A.getOrigin(),
                                                                            B.getExps()[j], B.getShell(), B.getOrigin(),
                                                                            C.getExps()[k], C.getShell(), C.getOrigin(),
                                                                            D.getExps()[l], D.getShell(),
                                                                            D.getOrigin());
                }
            }
        }
    }
    return tei;
}









