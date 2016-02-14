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

#include <utility>

template<typename T>
class Operator;

template<typename T>
class BestNest : public Operator<T> {
public:
    BestNest() = default;

    virtual ~BestNest() {
    }

    virtual void apply(CuckooSearch<T> &cs) const override;
};

//-------------------------------------------------------------
template<typename T>
void BestNest<T>::apply(CuckooSearch<T> &cs) const {
    for (auto i = 0u; i < cs.nest.size(); i++) {
        if (*cs.newNest[i] <= *cs.nest[i]) {
            *cs.nest[i] = *cs.newNest[i];
            if (*cs.nest[i] < *cs.nest[cs.bestNest]) {
                cs.bestNest = i;
            }
        }
    }
}



















