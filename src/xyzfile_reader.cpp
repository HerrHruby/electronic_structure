//
// Created by Ian Wu on 22/08/2020.
//

#include "../include/xyzfile_reader.h"
#include <iostream>

std::tuple <std::vector<std::string>, std::vector<std::valarray<double>>> read_file() {
    std::ifstream infile("../Molecule_Files/Ethanol.txt");

    std::string atom;
    double x, y, z;
    std::vector<std::valarray<double>> pos_list;
    std::vector<std::string> atom_list;

    while (infile >> atom >> x >> y >> z) {
        std::valarray<double> coord_arr = {x/0.52917721092, y/0.52917721092, z/0.52917721092};
        pos_list.push_back(coord_arr);
        atom_list.push_back(atom);
    }

    std::tuple <std::vector<std::string>, std::vector<std::valarray<double>>> read_tup = std::make_tuple(atom_list, pos_list);

    return read_tup;
}

