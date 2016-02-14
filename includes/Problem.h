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
#include <limits>
#include <algorithm>
#include <utility>
#include <cmath>

#include <Utils.h>

class Problem;

using fn_Problem_2_double = std::function<double(const Problem &)>;

using vdouble = std::vector<double>;

using fn__2_double = std::function<double()>;

using std::rel_ops::operator!=;

using std::rel_ops::operator>;

using std::rel_ops::operator<=;

using std::rel_ops::operator>=;


class Problem {
public:
    vdouble solution;
    const fn_Problem_2_double &fn;
    const fn__2_double &gen;
    double fitness;
    uint nd;
    double lb;
    double ub;

    Problem() = delete;

    Problem(const fn_Problem_2_double &fn, const fn__2_double &gen, uint nd, double lb, double ub);

    Problem(const Problem &&rhs) noexcept;

    Problem &operator=(const Problem &&rhs) noexcept;

    Problem(const Problem &rhs);

    Problem &operator=(const Problem &rhs);

    friend bool operator<(const Problem &lhs, const Problem &rhs);

    friend bool operator==(const Problem &lhs, const Problem &rhs);

    virtual ~Problem();

    virtual void evaluate();

    virtual void checkBounds(uint pos);
};

//----------------------------------------------------------------------------------------------
Problem::Problem(const fn_Problem_2_double &_fn, const fn__2_double &_gen, uint nd, double lb, double ub) :
        fn(_fn), gen(_gen) {
    this->nd = nd;
    this->lb = lb;
    this->ub = ub;
    this->solution.resize(nd);
    this->fitness = std::numeric_limits<double>::max();
    GENERATE(solution, gen)
    this->evaluate();
}

//----------------------------------------------------------------------------------------------
Problem::~Problem() {
}

//----------------------------------------------------------------------------------------------
Problem::Problem(const Problem &&rhs) noexcept : fn(rhs.fn), gen(rhs.gen) {
    nd = rhs.nd;
    lb = rhs.lb;
    ub = rhs.ub;
    solution = std::move(rhs.solution);
    fitness = rhs.fitness;
}

//----------------------------------------------------------------------------------------------
Problem::Problem(const Problem &rhs) : fn(rhs.fn), gen(rhs.gen) {
    nd = rhs.nd;
    lb = rhs.lb;
    ub = rhs.ub;
    solution = rhs.solution;
    fitness = rhs.fitness;
}

//----------------------------------------------------------------------------------------------
Problem &Problem::operator=(const Problem &&rhs) noexcept {
    if (this != &rhs) {
        nd = rhs.nd;
        lb = rhs.lb;
        ub = rhs.ub;
        solution = std::move(rhs.solution);
        fitness = rhs.fitness;
    }

    return *this;
}

//----------------------------------------------------------------------------------------------
Problem &Problem::operator=(const Problem &rhs) {
    if (this != &rhs) {
        nd = rhs.nd;
        lb = rhs.lb;
        ub = rhs.ub;
        solution = rhs.solution;
        fitness = rhs.fitness;
    }

    return *this;
}

//----------------------------------------------------------------------------------------------
bool operator<(const Problem &lhs, const Problem &rhs) {
    return std::isless(lhs.fitness, rhs.fitness);
}

//----------------------------------------------------------------------------------------------
bool operator==(const Problem &lhs, const Problem &rhs) {
    return !std::islessgreater(lhs.fitness, rhs.fitness);
}

//----------------------------------------------------------------------------------------------
void Problem::evaluate() {
    fitness = fn(*this);
}

//----------------------------------------------------------------------------------------------
void Problem::checkBounds(uint pos) {
    if (pos >= 0u && pos < solution.size()) {
        if (std::isnan(solution[pos])) {
            solution[pos] = (std::signbit(solution[pos]) ? lb : ub);
        } else {
            solution[pos] = (std::isless(solution[pos], lb) ? lb : solution[pos]);
            solution[pos] = (std::isgreater(solution[pos], ub) ? ub : solution[pos]);
        }
    }
}

//----------------------------------------------------------------------------------------------
