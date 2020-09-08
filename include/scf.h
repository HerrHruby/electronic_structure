//
// Created by Ian Wu on 24/08/2020.
//

#ifndef HF_SCF_H
#define HF_SCF_H
#include <gmp.h>
#include <vector>
#include "basis_function.h"
#include "xyzfile_reader.h"
#include "STO3G_data.h"
#include "contracted_integrals.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <arb.h>
#include <arb_hypgeom.h>



std::vector<Basis_Function> compute_basis_list(std::tuple<std::vector<std::string>,
        std::vector<std::valarray<double>>> input_dat);

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> compute_S_mat(std::vector<Basis_Function> basis_list);

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
compute_HCore(std::vector<Basis_Function> basis_list, std::tuple<std::vector<std::string>,
        std::vector<std::valarray<double>>> input_dat);

std::vector<std::vector<std::vector<std::vector<double>>>> compute_TE_tensor(std::vector<Basis_Function> basis_list);

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
orthogonalise(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const &input_mat);

double mat_convergence(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const &M_new,
                      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const &M_old);

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> compute_P(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P,
                                                               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C,
                                                               int basis_len, int N_e);

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> update_F(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> F,
                                                              Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> HCore,
                                                              Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P,
                                                              std::vector<std::vector<std::vector<std::vector<double>>>> TE_tensor,
                                                              int basis_len);

double compute_var_E(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P,
                   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> HCore,
                   Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> F, int basis_len);

double compute_nuc_E(std::vector<std::string> atom_list, std::vector<std::valarray<double>> pos_list);

void SCF(std::vector<Basis_Function> const &basis_list, std::tuple<std::vector<std::string>,
        std::vector<std::valarray<double>>> input_dat, double convergence, int iters, int N_e);


#endif //HF_SCF_H
