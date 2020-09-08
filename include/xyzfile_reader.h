//
// Created by Ian Wu on 22/08/2020.
//

#ifndef HF_XYZFILE_READER_H
#define HF_XYZFILE_READER_H

#include <vector>
#include <fstream>
#include <valarray>

std::tuple <std::vector<std::string>, std::vector<std::valarray<double>>> read_file();


#endif //HF_XYZFILE_READER_H
