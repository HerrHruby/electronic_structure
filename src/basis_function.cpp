//
// Created by Ian Wu on 21/08/2020.
//

#include "../include/basis_function.h"

#include <utility>

int fact2(int n) {
    int k = 1;
    for (int i = 3; i < n + 1; i += 2) {
        k *= i;
    }
    return k;
}

void Basis_Function::normalise() {
    int l, m, n;
    std::tie(l, m, n) = shell;
    int L = l + m + n;
    for (double i : exps) {
        double norm_val = sqrt(
                pow(2, 2.0 * L + 1.5) * pow(i, L + 1.5) / fact2(2 * l - 1) / fact2(2 * m - 1) / fact2(2 * n - 1) /
                pow(M_PI, 1.5));
        norm.push_back(norm_val);
    }
    double prefactor = pow(M_PI, 1.5) * (double)fact2(2 * l - 1) * (double)fact2(2 * m - 1) * (double)fact2(2 * n - 1) / pow(2.0, L);

    double N = 0.0;
    int num_exps = exps.size();
    for (int i = 0; i < num_exps; i++) {
        for (int j = 0; j < num_exps; j++) {
            N += norm[i] * norm[j] * coefs[i] * coefs[j] / pow(exps[i] + exps[j], L + 1.5);
        }
    }
    N *= prefactor;
    N = pow(N, -0.5);
    for (int i = 0; i < num_exps; i++) {
        coefs[i] *= N;
    }
}

Basis_Function::Basis_Function(std::valarray<double> origin, std::tuple<int, int, int> shell,
                               std::vector<double> exps, std::vector<double> coefs) :
        origin(std::move(origin)), shell(std::move(shell)), exps(std::move(exps)), coefs(std::move(coefs)) {
    normalise();
}

const std::valarray<double> &Basis_Function::getOrigin() const {
    return origin;
}

void Basis_Function::setOrigin(const std::valarray<double> &origin) {
    Basis_Function::origin = origin;
}

const std::tuple<int, int, int> &Basis_Function::getShell() const {
    return shell;
}

void Basis_Function::setShell(const std::tuple<int, int, int> &shell) {
    Basis_Function::shell = shell;
}

const std::vector<double> &Basis_Function::getExps() const {
    return exps;
}

void Basis_Function::setExps(const std::vector<double> &exps) {
    Basis_Function::exps = exps;
}

const std::vector<double> &Basis_Function::getCoefs() const {
    return coefs;
}

void Basis_Function::setCoefs(const std::vector<double> &coefs) {
    Basis_Function::coefs = coefs;
}

const std::vector<double> &Basis_Function::getNorm() const {
    return norm;
}

void Basis_Function::setNorm(const std::vector<double> &norm) {
    Basis_Function::norm = norm;
}
