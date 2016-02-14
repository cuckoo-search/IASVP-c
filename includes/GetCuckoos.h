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

#include <algorithm>

#include <Utils.h>

template<typename T>
class Operator;

template<typename T>
using Nest = std::vector<std::unique_ptr<T>>;

template<typename T>
class GetCuckoos : public Operator<T> {
public:
    GetCuckoos() = default;

    virtual ~GetCuckoos() {
    }

    virtual void apply(CuckooSearch<T> &cs) const override;
};

//-------------------------------------------------------------------
template<typename T>
void GetCuckoos<T>::apply(CuckooSearch<T> &cs) const {
    double beta = 3.0 / 2.0;
    double _sigma = pow((Gamma(1.0 + beta) * sin(M_PI * beta / 2.0) /
                         (Gamma((1.0 + beta) / 2.0) * beta * pow(2.0, ((beta - 1.0) / 2.0)))),
                        (1.0 / beta));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> normal(0.0, 1.0);

    MAP(cs.nest, cs.newNest, ([&normal, _sigma, &gen, beta, &cs](const auto &x) {
        auto result = std::make_unique<T>(*x);

        for (auto j = 0u; j < x->solution.size(); j++) {
            auto u_j = normal(gen) * _sigma;
            auto v_j = normal(gen);
            auto step_j = pow(u_j / fabs(v_j), (1.0 / beta));
            auto stepsize_j = (0.01 * step_j) * (x->solution[j] - cs.nest[cs.bestNest]->solution[j]);
            result->solution[j] = x->solution[j] + stepsize_j * normal(gen);
            result->checkBounds(j);
            result->evaluate();
        }

        return result;
    }))
}

















