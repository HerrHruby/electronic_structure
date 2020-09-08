//
// Created by Ian Wu on 21/08/2020.
//

#ifndef HF_BASIS_FUNCTION_H
#define HF_BASIS_FUNCTION_H

#include <utility>
#include <vector>
#include <tuple>
#include <valarray>
#include <cmath>

int fact2(int n);

class Basis_Function {

public:
    Basis_Function(std::valarray<double> origin, std::tuple<int, int, int> shell, std::vector<double> exps,
                   std::vector<double> coefs);

    void normalise();

    const std::valarray<double> &getOrigin() const;

    void setOrigin(const std::valarray<double> &origin);

    const std::tuple<int, int, int> &getShell() const;

    void setShell(const std::tuple<int, int, int> &shell);

    const std::vector<double> &getExps() const;

    void setExps(const std::vector<double> &exps);

    const std::vector<double> &getCoefs() const;

    void setCoefs(const std::vector<double> &coefs);

    const std::vector<double> &getNorm() const;

    void setNorm(const std::vector<double> &norm);


private:
    std::valarray<double> origin;
    std::tuple<int, int, int> shell;
    std::vector<double> exps;
    std::vector<double> coefs;
    std::vector<double> norm;


};


#endif //HF_BASIS_FUNCTION_H
