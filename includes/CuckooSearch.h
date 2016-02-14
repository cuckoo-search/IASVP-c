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

#include <functional>
#include <memory>
#include <vector>
#include <algorithm>

#include <Operator.h>
#include <BestNest.h>
#include <EmptyNest.h>
#include <GetCuckoos.h>


template<typename T>
using Nest = std::vector<std::unique_ptr<T>>;

template<typename T>
using Operators = std::initializer_list<std::unique_ptr<Operator<T>>>;

template<typename T>
using fn_T_2_bool = std::function<bool(const T &)>;

template<typename T>
using fn_T_2_double = std::function<double(const T &)>;

using fn__2_double = std::function<double()>;

using vint = std::vector<int>;

template<typename T>
class CuckooSearch {
public:
    const fn_T_2_double<T> &fn;
    const fn__2_double &gen;
    const fn_T_2_bool<T> &stop;

    Nest<T> nest;
    Nest<T> newNest;

    vint perm1;
    vint perm2;

    uint bestNest = 0u;
    uint niter = 0u;
    uint eggs;
    uint nd;
    double lb;
    double ub;
    float pa;

    CuckooSearch() = delete;

    CuckooSearch(const CuckooSearch &rhs) = delete;

    CuckooSearch &operator=(const CuckooSearch &rhs) = delete;

    CuckooSearch(uint eggs, uint nd, double lb, double ub, float pa, const fn_T_2_double<T> &_fn,
                 const fn__2_double &_gen, const fn_T_2_bool<T> &_stop);

    void shuffle();

    virtual ~CuckooSearch();

    virtual const T search();

    virtual const T search(Operators<T> ops);

    const T getBestNest() const;

    virtual void checkBestNest();
};

//---------------------------------------------------------------------
template<typename T>
CuckooSearch<T>::CuckooSearch(uint eggs, uint nd, double lb, double ub, float pa, const fn_T_2_double<T> &_fn,
                              const fn__2_double &_gen, const fn_T_2_bool<T> &_stop) : fn(_fn), gen(_gen), stop(_stop) {
    this->eggs = eggs;
    this->nd = nd;
    this->pa = pa;
    this->ub = ub;
    this->lb = lb;
    this->perm1.resize(eggs);
    this->perm2.resize(eggs);
    IOTA(perm1, 0)
    IOTA(perm2, 0)
    this->nest.resize(eggs);
    this->newNest.resize(eggs);
    GENERATE(nest, [this]() { return std::make_unique<T>(this->fn, this->gen, this->nd, this->lb, this->ub); })
}

//---------------------------------------------------------------------
template<typename T>
CuckooSearch<T>::~CuckooSearch() { }

//---------------------------------------------------------------------
template<typename T>
const T CuckooSearch<T>::search() {
    return search({std::make_unique<GetCuckoos<T>>(),
                   std::make_unique<BestNest<T>>(),
                   std::make_unique<EmptyNest<T>>(),
                   std::make_unique<BestNest<T>>()});
}

//---------------------------------------------------------------------
template<typename T>
const T CuckooSearch<T>::search(Operators<T> ops) {
    checkBestNest();

    while (!stop(getBestNest())) {
        FOR_EACH(ops, [this](const auto &op) { op->apply(*this); })
        niter++;
    }

    return getBestNest();
}

//---------------------------------------------------------------------
template<typename T>
const T CuckooSearch<T>::getBestNest() const {
    return *nest[bestNest];
}

//---------------------------------------------------------------------
template<typename T>
void CuckooSearch<T>::shuffle() {
    _shuffle(perm1, perm2);
}

//---------------------------------------------------------------------
template<typename T>
void CuckooSearch<T>::checkBestNest() {
    auto pos = 0u;
    bestNest = REDUCE(nest, bestNest, ([&pos, this](auto bestnest, const auto &x) {
        return (*this->nest[pos++] < *this->nest[bestnest] ? pos - 1 : bestnest);
    }))
}

//---------------------------------------------------------------------

