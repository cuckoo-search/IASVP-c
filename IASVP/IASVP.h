/**
 * Authors:
 * Rafael Arturo Trujillo Rasúa <trujillo@uci.cu>
 * Rigoberto Leander Salgado Reyes <rlsalgado2006@gmail.com>
 *
 * Copyright 2016 by Rigoberto Leander Salgado Reyes.
 *
 * This program is licensed to you under the terms of version 3 of the
 * GNU Affero General Public License. This program is distributed WITHOUT
 * ANY EXPRESS OR IMPLIED WARRANTY, INCLUDING THOSE OF NON-INFRINGEMENT,
 * MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Please refer to the
 * AGPL (http:www.gnu.org/licenses/agpl-3.0.txt) for more details.
 */

#pragma once

#include <vector>
#include <functional>
#include <algorithm>

#include <cmath>

#include <lapack.h>

#include <Utils.h>

using vdouble = std::vector<double>;
using vint = std::vector<int>;
using fn_vdouble_2_vdouble = std::function<vdouble(const vdouble &)>;

class IASVP {
private:
    const fn_vdouble_2_vdouble &matrixMaker;
    const vdouble sigma;

public:
    IASVP() = delete;

    IASVP(const IASVP &rhs) = delete;

    IASVP &operator=(const IASVP &rhs) = delete;

    IASVP(const vdouble &seed, const fn_vdouble_2_vdouble &fn);

    ~IASVP();

    const vdouble &getSigma() const;

    double RelativeError(const vdouble &seed);

    vdouble IASVPToeplitzTriInfNLES(const vdouble &seed) const;

    double FIASVPToeplitzTriInf(const vdouble &seed) const;
};

//----------------------------------------------------------------------------------------------
IASVP::IASVP(const vdouble &seed, const fn_vdouble_2_vdouble &fn) :
        matrixMaker(fn), sigma(CalcSV(seed, fn)) {
}

//----------------------------------------------------------------------------------------------
IASVP::~IASVP() {
}

//----------------------------------------------------------------------------------------------
const vdouble &IASVP::getSigma() const {
    return sigma;
}

//----------------------------------------------------------------------------------------------
double IASVP::RelativeError(const vdouble &seed) {
    return FIASVPToeplitzTriInf(seed) / NORM2(sigma);
}

//----------------------------------------------------------------------------------------------
double IASVP::FIASVPToeplitzTriInf(const vdouble &seed) const {
    auto new_sigma = CalcSV(seed, matrixMaker);
    auto it = std::cbegin(sigma);
    auto tmp = REDUCE(new_sigma, 0.0, [&it](auto acc, auto x) { return acc + pow(x - (*(it++)), 2); })

    return sqrt(tmp);
}

//----------------------------------------------------------------------------------------------
vdouble IASVP::IASVPToeplitzTriInfNLES(const vdouble &seed) const {
    auto new_sigma = CalcSV(seed, matrixMaker);
    INNER_MAP_2(new_sigma, sigma, [](auto a, auto b) { return a - b; })

    return new_sigma;
}

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
template<typename T>
class Operator;

template<typename T>
class HybridEmptyNest : public Operator<T> {
public:
    const IASVP &iasvp;
    const fn_vdouble_2_vdouble F = [this](const auto &seed) { return this->iasvp.IASVPToeplitzTriInfNLES(seed); };

    const fn_vdouble_2_vdouble Jac = std::bind(JacIASVPToeplitzTriInf, std::placeholders::_1, makeToeplitz);

    HybridEmptyNest(const IASVP &_iasvp) : iasvp(_iasvp) { }

    virtual ~HybridEmptyNest() { }

    virtual void apply(CuckooSearch<T> &cs) const override;
};

//------------------------------------------------------------
template<typename T>
void HybridEmptyNest<T>::apply(CuckooSearch<T> &cs) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    double rand = dis(gen);

    cs.shuffle();

    for (auto i = 0u; i < cs.eggs; i++) {
        if (std::isgreater(dis(gen), cs.pa)) {
            for (auto j = 0u; j < cs.nd; j++) {
                cs.newNest[i]->solution[j] = cs.nest[i]->solution[j] +
                                             rand * (cs.nest[cs.perm1[i]]->solution[j] -
                                                     cs.nest[cs.perm2[i]]->solution[j]);
            }
            int iter = 0;
            newtonBiseccionNLES(this->F, cs.newNest[i]->solution, this->Jac, 0.0000001, 0.0000001, 10, iter);
            cs.newNest[i]->evaluate();
        } else {
            *cs.newNest[i] = *cs.nest[i];
        }
    }
}
















