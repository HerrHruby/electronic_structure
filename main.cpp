#include <iostream>
#include "include/scf.h"


int main() {

    int N_e = 26;
    int iters = 100000;
    float convergence = 0.0001;

    std::tuple <std::vector<std::string>, std::vector<std::valarray<double>>> input_tup = read_file();
    std::vector<std::string> atom_list;
    std::vector<std::valarray<double>> pos_list;
    std::tie(atom_list, pos_list) = input_tup;

    std::vector<Basis_Function> basis_list = compute_basis_list(input_tup);

    SCF(basis_list, input_tup, convergence, iters, N_e);

    return 0;


}