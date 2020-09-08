//
// Created by Ian Wu on 21/08/2020.
//

#ifndef HF_INTEGRALS_H
#define HF_INTEGRALS_H
#include <cmath>
#include <tuple>
#include <valarray>
#include <gmp.h>

double boys(double n, double T);

std::valarray<double> gaussian_product(double a, const std::valarray<double> &A, double b, const std::valarray<double> &B);

double vec_norm(std::valarray<double> A);

double E(int i, int j, int t, double Q, double a, double b);

double overlap(double a, std::tuple<int, int, int> ang1, std::valarray<double> A, double b,
              std::tuple<int, int, int> ang2, std::valarray<double> B);

double kinetic(double a, std::tuple<int, int, int> ang1, const std::valarray<double> &A, double b,
              std::tuple<int, int, int> ang2, const std::valarray<double> &B);

double R(int t, int u, int v, int n, double p, double Dx, double Dy, double Dz, double D);

double coulomb(double a, std::tuple<int, int, int> ang1, std::valarray<double> A, double b,
              std::tuple<int, int, int> ang2, std::valarray<double> B, std::valarray<double> centre);

long double two_electron(double a, std::tuple<int, int, int> ang1, std::valarray<double> A, double b,
                   std::tuple<int, int, int> ang2, std::valarray<double> B, double c, std::tuple<int, int, int> ang3,
                   std::valarray<double> C, double d, std::tuple<int, int, int> ang4, std::valarray<double> D);

#endif //HF_INTEGRALS_H
