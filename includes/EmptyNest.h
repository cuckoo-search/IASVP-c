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
class EmptyNest : public Operator<T> {
public:
    EmptyNest() = default;

    virtual ~EmptyNest() {
    }

    virtual void apply(CuckooSearch<T> &cs) const override;
};

//------------------------------------------------------------
template<typename T>
void EmptyNest<T>::apply(CuckooSearch<T> &cs) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    auto rand = dis(gen);

    cs.shuffle();

    for (auto i = 0u; i < cs.eggs; i++) {
        if (dis(gen) > cs.pa) {
            for (auto j = 0u; j < cs.nd; j++) {
                cs.newNest[i]->solution[j] = cs.nest[i]->solution[j] +
                                            rand * (cs.nest[cs.perm1[i]]->solution[j] -
                                                    cs.nest[cs.perm2[i]]->solution[j]);
            }
            cs.newNest[i]->evaluate();
        } else {
            *cs.newNest[i] = *cs.nest[i];
        }
    }
}

