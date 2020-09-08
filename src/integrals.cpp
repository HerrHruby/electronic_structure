//
// Created by Ian Wu on 21/08/2020.
//

#include "../include/integrals.h"
#include "arb.h"
#include "arb_hypgeom.h"
#include <iostream>


double boys(double n, double T) {

    slong prec;
    prec = 256;
    double res_d;
    arb_t a, b, t, x, y, res;
    arb_init(a); arb_init(b); arb_init(t); arb_init(x); arb_init(y); arb_init(res);
    arb_set_d(x, 0.5);
    arb_set_d(y, 1.5);
    arb_set_d(a, n);
    arb_set_d(b, n);
    arb_set_d(t, -T);
    arb_add(a, a, x, prec);
    arb_add(b, b, y, prec);
    arb_hypgeom_1f1(res, a, b, t, 0, prec);
//    arb_printn(res, 15, 0);
    res_d = arf_get_d(arb_midref(res), ARF_RND_DOWN);
    arb_clear(x); arb_clear(y); arb_clear(a); arb_clear(b); arb_clear(t); arb_clear(res);

    return res_d/(2.0 * n + 1);
}

std::valarray<double> gaussian_product(double a,  const std::valarray<double> &A, double b, const std::valarray<double>& B) {
    return (a * A + b * B) / (a + b);
}

double vec_norm(std::valarray<double> A) {
    std::valarray<double> dotp = A * A;
    double norm = dotp.sum();
    return sqrt(norm);
}

double E(int i, int j, int t, double Q, double a, double b) {

    double p = a + b;
    double q = (a * b) / p;

    if (t < 0 || t > (i + j)) {
        return 0.0;

    } else if (i == 0 && j == 0 && t == 0) {
        return exp(-q * Q * Q);

    } else if (j == 0) {
        return (1.0 / (2.0 * p)) * E(i - 1, j, t - 1, Q, a, b) - ((q * Q) / a) * E(i - 1, j, t, Q, a, b) +
               (t + 1) * E(i - 1, j, t + 1, Q, a, b);
    } else {
        return (1.0 / (2.0 * p)) * E(i, j - 1, t - 1, Q, a, b) + ((q * Q) / b) * E(i, j - 1, t, Q, a, b) +
               (t + 1) * E(i, j - 1, t + 1, Q, a, b);
    }
}

double overlap(double a, std::tuple<int, int, int> ang1, std::valarray<double> A, double b, std::tuple<int, int, int> ang2,
              std::valarray<double> B) {

    int l1, m1, n1, l2, m2, n2;
    std::tie(l1, m1, n1) = ang1;
    std::tie(l2, m2, n2) = ang2;
    double S1 = E(l1, l2, 0, A[0] - B[0], a, b);
    double S2 = E(m1, m2, 0, A[1] - B[1], a, b);
    double S3 = E(n1, n2, 0, A[2] - B[2], a, b);

    return S1 * S2 * S3 * pow(M_PI / (a + b), 1.5);
}

double kinetic(double a, std::tuple<int, int, int> ang1, const std::valarray<double> &A, double b, std::tuple<int, int, int> ang2,
        const std::valarray<double> &B) {

    int l1, m1, n1, l2, m2, n2;
    std::tie(l1, m1, n1) = ang1;
    std::tie(l2, m2, n2) = ang2;
    double term1 = b * (double) (2 * (l2 + m2 + n2) + 3) *
                  overlap(a, std::make_tuple(l1, m1, n1), A, b, std::make_tuple(l2, m2, n2), B);
    double term2 = -2 * pow(b, 2) * (overlap(a, std::make_tuple(l1, m1, n1), A, b, std::make_tuple(l2 + 2, m2, n2), B) +
                                     overlap(a, std::make_tuple(l1, m1, n1), A, b, std::make_tuple(l2, m2 + 2, n2), B) +
                                     overlap(a, std::make_tuple(l1, m1, n1), A, b, std::make_tuple(l2, m2, n2 + 2), B));
    double term3 = -0.5 *
                  ((double) (l2 * (l2 - 1)) *
                   overlap(a, std::make_tuple(l1, m1, n1), A, b, std::make_tuple(l2 - 2, m2, n2), B) +
                   (double) (m2 * (m2 - 1)) *
                   overlap(a, std::make_tuple(l1, m1, n1), A, b, std::make_tuple(l2, m2 - 2, n2), B) +
                   (double) (n2 * (n2 - 1)) *
                   overlap(a, std::make_tuple(l1, m1, n1), A, b, std::make_tuple(l2, m2, n2 - 2), B));

    return term1 + term2 + term3;
}

double R(int t, int u, int v, int n, double p, double Dx, double Dy, double Dz, double D) {

    double T = p * D * D;

    double r = 0.0;
    if (t == 0 && u == 0 && v == 0) {
        r += pow(-2 * p, n) * boys(n, T);
    } else if (t == 0 && u == 0) {
        if (v > 1) {
            r += (v - 1) * R(t, u, v - 2, n + 1, p, Dx, Dy, Dz, D);
        }
        r += Dz * R(t, u, v - 1, n + 1, p, Dx, Dy, Dz, D);
    } else if (t == 0) {
        if (u > 1) {
            r += (u - 1) * R(t, u - 2, v, n + 1, p, Dx, Dy, Dz, D);
        }
        r += Dy * R(t, u - 1, v, n + 1, p, Dx, Dy, Dz, D);
    } else {
        if (t > 1) {
            r += (t - 1) * R(t - 2, u, v, n + 1, p, Dx, Dy, Dz, D);
        }
        r += Dx * R(t - 1, u, v, n + 1, p, Dx, Dy, Dz, D);
    }

    return r;
}

double coulomb(double a, std::tuple<int, int, int> ang1, std::valarray<double> A, double b, std::tuple<int, int, int> ang2,
              std::valarray<double> B, std::valarray<double> centre) {

    int l1, m1, n1, l2, m2, n2;
    std::tie(l1, m1, n1) = ang1;
    std::tie(l2, m2, n2) = ang2;
    double p = a + b;
    std::valarray<double> P = gaussian_product(a, A, b, B);
    double PC = vec_norm(P - centre);

    double j = 0.0;
    for (int t = 0; t < (l1 + l2 + 1); t++) {
        for (int u = 0; u < (m1 + m2 + 1); u++) {
            for (int v = 0; v < (n1 + n2 + 1); v++) {
                j += E(l1, l2, t, A[0] - B[0], a, b) * E(m1, m2, u, A[1] - B[1], a, b) *
                     E(n1, n2, v, A[2] - B[2], a, b) *
                     R(t, u, v, 0, p, P[0] - centre[0], P[1] - centre[1], P[2] - centre[2], PC);

            }
        }
    }
    j *= (2.0 * M_PI) / p;
    return j;
}


long double
two_electron(double a, std::tuple<int, int, int> ang1, std::valarray<double> A, double b, std::tuple<int, int, int> ang2,
             std::valarray<double> B, double c, std::tuple<int, int, int> ang3, std::valarray<double> C, double d,
             std::tuple<int, int, int> ang4, std::valarray<double> D) {

    int l1, m1, n1, l2, m2, n2, l3, m3, n3, l4, m4, n4;
    std::tie(l1, m1, n1) = ang1;
    std::tie(l2, m2, n2) = ang2;
    std::tie(l3, m3, n3) = ang3;
    std::tie(l4, m4, n4) = ang4;
    double p = a + b;
    double q = c + d;
    double alpha = p * q / (p + q);
    std::valarray<double> P = gaussian_product(a, A, b, B);
    std::valarray<double> Q = gaussian_product(c, C, d, D);
    double PQ = vec_norm(P - Q);

    long double z = 0.0;
    for (int t = 0; t < (l1 + l2 + 1); t++) {
        for (int u = 0; u < (m1 + m2 + 1); u++) {
            for (int v = 0; v < (n1 + n2 + 1); v++) {
                for (int w = 0; w < (l3 + l4 + 1); w++) {
                    for (int x = 0; x < (m3 + m4 + 1); x++) {
                        for (int y = 0; y < (n3 + n4 + 1); y++) {
                            z += E(l1, l2, t, A[0] - B[0], a, b) *
                                 E(m1, m2, u, A[1] - B[1], a, b) *
                                 E(n1, n2, v, A[2] - B[2], a, b) *
                                 E(l3, l4, w, C[0] - D[0], c, d) *
                                 E(m3, m4, x, C[1] - D[1], c, d) *
                                 E(n3, n4, y, C[2] - D[2], c, d) *
                                 pow(-1, w + x + y) *
                                 R(t + w, u + x, v + y, 0, alpha,
                                   P[0] - Q[0], P[1] - Q[1], P[2] - Q[2], PQ);
                        }
                    }
                }
            }
        }
    }
    z *= 2 * pow(M_PI, 2.5) / (p * q * sqrt(p + q));
    return z;
}

