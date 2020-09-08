//
// Created by Ian Wu on 22/08/2020.
//

#include "../include/scf.h"

std::vector<Basis_Function> compute_basis_list(std::tuple<std::vector<std::string>,
        std::vector<std::valarray<double>>> input_dat) {

    std::vector<Basis_Function> basis_list;
    std::vector<std::string> atom_list;
    std::vector<std::valarray<double>> pos_list;
    std::tie(atom_list, pos_list) = input_dat;
    int input_len = atom_list.size();

    for (int i = 0; i < input_len; i++) {
        std::string atom = atom_list[i];
        std::vector<std::vector<double>> coef_list = coef_table.at(atom);
        std::vector<double> coefs_1s = coef_list[0];

        basis_list.emplace_back(pos_list[i], std::tuple(0, 0, 0), coefs_1s, d1s);
        for (int j = 0; j < coefs_1s.size(); j++) {
        }

        if (coef_list.size() >= 2) {
            std::vector<double> coefs_2sp = coef_list[1];

            basis_list.emplace_back(pos_list[i], std::tuple(0, 0, 0), coefs_2sp, d2s);
            basis_list.emplace_back(pos_list[i], std::tuple(1, 0, 0), coefs_2sp, d2p);
            basis_list.emplace_back(pos_list[i], std::tuple(0, 1, 0), coefs_2sp, d2p);
            basis_list.emplace_back(pos_list[i], std::tuple(0, 0, 1), coefs_2sp, d2p);
            for (int j = 0; j < coefs_1s.size(); j++) {
            }
        }

        if (coef_list.size() >= 3) {
            std::vector<double> coefs_3sp = coef_list[2];

            basis_list.emplace_back(pos_list[i], std::tuple(0, 0, 0), coefs_3sp, d3s);
            basis_list.emplace_back(pos_list[i], std::tuple(1, 0, 0), coefs_3sp, d3p);
            basis_list.emplace_back(pos_list[i], std::tuple(0, 1, 0), coefs_3sp, d3p);
            basis_list.emplace_back(pos_list[i], std::tuple(0, 0, 1), coefs_3sp, d3p);
            for (int j = 0; j < coefs_1s.size(); j++) {
            }
        }
    }
    return basis_list;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> compute_S_mat(std::vector<Basis_Function> basis_list) {

    int basis_count = basis_list.size();
    Eigen::MatrixXd S_mat(basis_count, basis_count);
    for (int i = 0; i < basis_count; i++) {
        for (int j = 0; j < basis_count; j++) {
            double Sij = S(basis_list[i], basis_list[j]);
            S_mat(i, j) = Sij;
        }
    }
    return S_mat;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
compute_HCore(std::vector<Basis_Function> basis_list, std::tuple<std::vector<std::string>,
        std::vector<std::valarray<double>>> input_dat) {

    int basis_count = basis_list.size();

    Eigen::MatrixXd T_mat(basis_count, basis_count);
    Eigen::MatrixXd V_mat(basis_count, basis_count);

    std::vector<std::string> atom_list;
    std::vector<std::valarray<double>> pos_list;
    std::tie(atom_list, pos_list) = input_dat;

    for (int i = 0; i < basis_count; i++) {
        for (int j = 0; j < basis_count; j++) {
            double Tij = T(basis_list[i], basis_list[j]);
            T_mat(i, j) = Tij;
            double v = 0.0;
            for (int c = 0; c < atom_list.size(); c++) {
                std::string atom = atom_list[c];
                int Z = nuc_charge.at(atom);
                v += -Z * V(basis_list[i], basis_list[j], pos_list[c]);
            }
            V_mat(i, j) = v;
        }
    }
    std::cout << "T: " << T_mat << '\n' << std::endl;
    std::cout << "V: " << V_mat << '\n' << std::endl;
    return T_mat + V_mat;
}

std::vector<std::vector<std::vector<std::vector<double>>>> compute_TE_tensor(std::vector<Basis_Function> basis_list) {

    int basis_count = basis_list.size();
    std::vector<std::vector<std::vector<std::vector<double>>>>
            TE_tensor(basis_count, std::vector<std::vector<std::vector<double>>>
            (basis_count, std::vector<std::vector<double> >(basis_count, std::vector<double>(basis_count, 0))));

    int ij, kl;
    double te_tensor;
    for (int i = 0; i < basis_count; i++) {
        for (int j = 0; j < i + 1; j++) {
            ij = (i * (i + 1) / 2) + j;
            for (int k = 0; k < basis_count; k++) {
                for (int l = 0; l < k + 1; l++) {
                    kl = (k * (k + 1) / 2) + l;
                    if (ij >= kl) {
                        te_tensor = TEI(basis_list[i], basis_list[j], basis_list[k], basis_list[l]);
                        TE_tensor[i][j][k][l] = te_tensor;
                        TE_tensor[k][l][i][j] = te_tensor;
                        TE_tensor[j][i][l][k] = te_tensor;
                        TE_tensor[l][k][j][i] = te_tensor;
                        TE_tensor[j][i][k][l] = te_tensor;
                        TE_tensor[l][k][i][j] = te_tensor;
                        TE_tensor[i][j][l][k] = te_tensor;
                        TE_tensor[k][l][j][i] = te_tensor;
                    }

                }
            }
        }
    }
    return TE_tensor;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
orthogonalise(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const &input_mat) {

    int dims = input_mat.cols();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> es(input_mat);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> D = es.eigenvalues().asDiagonal();

    for (int i = 0; i < dims; i++) {
        double s = D(i, i);
        D(i, i) = pow(s, -0.5);
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P = es.eigenvectors();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X = P * (D * P.transpose());

    return X;
}

double mat_convergence(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const &M_new,
                      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const &M_old) {

    double delta = 0.0;
    int mat_size = M_new.cols();
    for (int i = 0; i < mat_size; i++) {
        for (int j = 0; j < mat_size; j++) {
            delta += pow((M_new(i, j) - M_old(i, j)), 2);
        }
    }
    delta = sqrt((delta / 4.0));
    return delta;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> compute_P(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P,
                                                               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C,
                                                               int basis_len, int N_e) {

    for (int i = 0; i < basis_len; i++) {
        for (int j = 0; j < basis_len; j++) {
            double p = 0.0;
            for (int elec = 0; elec < N_e / 2; elec++) {
                p += 2 * C(i, elec) * C(j, elec);
            }
            P(i, j) = p;
        }
    }
    return P;
}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> update_F(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> F,
                                                              Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> HCore,
                                                              Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P,
                                                              std::vector<std::vector<std::vector<std::vector<double>>>> TE_tensor,
                                                              int basis_len) {
    for (int i = 0; i < basis_len; i++) {
        for (int j = 0; j < basis_len; j++) {
            F(i, j) = HCore(i, j);
            for (int k = 0; k < basis_len; k++) {
                for (int l = 0; l < basis_len; l++) {
                    F(i, j) = F(i, j) + P(k, l) * (TE_tensor[i][j][k][l] - 0.5 * TE_tensor[i][k][j][l]);
                }
            }
        }
    }
    return F;
}

double compute_var_E(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P,
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> HCore,
                    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> F, int basis_len) {

    double energy = 0.0;
    for (int i = 0; i < basis_len; i++) {
        for (int j = 0; j < basis_len; j++) {
            energy += 0.5 * P(i, j) * (HCore(i, j) + F(i, j));
        }
    }
    return energy;
}

double compute_nuc_E(std::vector<std::string> atom_list, std::vector<std::valarray<double>> pos_list) {

    double tot_Z = 0.0;
    for (int i = 0; i < atom_list.size(); i++) {
        for (int j = 0; j < atom_list.size(); j++) {
            if (i == j) {
                continue;
            } else {
                std::string atom_a = atom_list[i];
                std::string atom_b = atom_list[j];
                int Z_a = nuc_charge.at(atom_a);
                int Z_b = nuc_charge.at(atom_b);
                double R = vec_norm(pos_list[i] - pos_list[j]);
                tot_Z += (double) (Z_a * Z_b) / R;
            }
        }
    }
    return 0.5 * tot_Z;
}

void SCF(std::vector<Basis_Function> const &basis_list, std::tuple<std::vector<std::string>,
        std::vector<std::valarray<double>>> input_dat, double convergence, int iters, int N_e) {

    std::vector<std::string> atom_list;
    std::vector<std::valarray<double>> pos_list;
    std::tie(atom_list, pos_list) = input_dat;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> S_mat = compute_S_mat(basis_list);
    std::cout << "Computing S..." << std::endl;
    std::cout << "S: " << S_mat << '\n' << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> HCore = compute_HCore(basis_list, input_dat);
    std::cout << "Computing HCore..." << std::endl;
    std::cout << "HCore: " << HCore << '\n' << std::endl;

    std::cout << "Computing Two-Electron Tensor..." << std::endl;
    std::vector<std::vector<std::vector<std::vector<double>>>> TE_tensor = compute_TE_tensor(basis_list);

    std::cout << "Orthogonalising S..." << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X = orthogonalise(S_mat);
    std::cout << "X: " << X << '\n' << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> F = HCore;

    double E;
    double delta = 0;
    int cycles = 0;
    int basis_len = basis_list.size();

    Eigen::MatrixXd P(basis_len, basis_len);
    Eigen::MatrixXd P_old(basis_len, basis_len);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C_prime;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> epsilon;

    std::cout << "Performing SCF..." << std::endl;
    while (((cycles < iters) && (abs(delta) >= convergence)) || cycles < 2) {
        std::cout << "-------------------------------------------" << '\n' << std::endl;
        std::cout << "Cycle Number: " << cycles << '\n' << std::endl;

        cycles += 1;

        F = X.transpose() * F * X;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(F);

        epsilon = es2.eigenvalues();
        C_prime = es2.eigenvectors();

        C = X * C_prime;

        if (cycles > 1) {
            P_old = P;
        }

        P = compute_P(P, C, basis_len, N_e);
        std::cout << "Density Matrix: " << '\n' << P << '\n' << std::endl;

        F = update_F(F, HCore, P, TE_tensor, basis_len);

        delta = mat_convergence(P, P_old);
        std::cout << "Density Matrix Change: " << delta << '\n' << std::endl;

        E = compute_var_E(P, HCore, F, basis_len);
        std::cout << "Variational Energy: " << E << '\n' << std::endl;
    }

    std::cout << "-------------------------------------------" << '\n' << "FINAL REPORT" << std::endl;
    if (cycles < iters) {
        std::cout << "SCF converged within " << cycles << " cycles" << '\n' << std::endl;
    } else {
        std::cout << "SCF failed to converge" << '\n' << std::endl;
    }

    double tot_Z = compute_nuc_E(atom_list, pos_list);

    std::cout << "Final Orbital Matrix: " << '\n' << C << '\n' << std::endl;
    std::cout << "Orbital Energies: " << '\n' << epsilon << '\n' << std::endl;
    std::cout << "Nuclear Repulsion Energy: " << tot_Z << '\n' << std::endl;
    std::cout << "Final Variational Energy: " << E << '\n' << std::endl;
    std::cout << "Total Energy: " << tot_Z + E << '\n' << std::endl;

    if (cycles <= iters) {
        std::cout << "Converged! Mission Success " << '\n' << std::endl;
    } else {
        std::cout << "Convergence Failure" << '\n' << std::endl;
    }

}













