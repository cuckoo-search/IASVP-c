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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

#define MAP_2(coll1, coll2, coll3, fn) \
std::transform(std::cbegin(coll1), std::cend(coll1), \
std::cbegin(coll2), std::begin(coll3), fn);

//----------------------------------------------------------------------------------------------
#define MAP(coll1, coll2, fn) \
std::transform(std::cbegin(coll1), std::cend(coll1), \
std::begin(coll2), fn);

//----------------------------------------------------------------------------------------------
#define INNER_MAP_2(coll1, coll2, fn) \
MAP_2(coll1, coll2, coll1, fn)

//----------------------------------------------------------------------------------------------
#define INNER_MAP(coll1, fn) \
MAP(coll1, coll1, fn)

//----------------------------------------------------------------------------------------------
#define REDUCE(coll, init, fn) \
std::accumulate(std::cbegin(coll), std::cend(coll), init, fn);

//----------------------------------------------------------------------------------------------
#define NORM2(coll) \
sqrt(std::accumulate(std::cbegin(coll), std::cend(coll), 0.0, \
        [](auto acc, auto x){return acc + x * x;}));

//----------------------------------------------------------------------------------------------
#define FOR_EACH(coll, fn) \
std::for_each(std::cbegin(coll), std::cend(coll), fn);

//----------------------------------------------------------------------------------------------
#define GENERATE(coll, fn) \
std::generate(std::begin(coll), std::end(coll), fn);

//----------------------------------------------------------------------------------------------
#define PRINT_ROW(coll) \
std::for_each(std::cbegin(coll), std::cend(coll), [](const auto& x){std::cout << x << " ";}); \
std::cout << std::endl;

//----------------------------------------------------------------------------------------------
#define PRINT_COL(coll) \
std::for_each(std::cbegin(coll), std::cend(coll), [](const auto& x){std::cout << x << std::endl;});

//----------------------------------------------------------------------------------------------
#define PRINT_FN(coll, fn) \
std::for_each(std::cbegin(coll), std::cend(coll), fn); \
std::cout << std::endl;

//----------------------------------------------------------------------------------------------
#define IOTA(coll, v) \
std::iota(std::begin(coll), std::end(coll), v);

//----------------------------------------------------------------------------------------------


/**
 * implements Weirstrass's form (infinite product)
 * Mathematical methods for Physicists, 4th ed. page 594.
 * see std::tgamma()
 */
auto Gamma(double z, uint Nterms = 1000) {
    auto g = 0.577216; // Euler-Mascheroni constant
    auto retVal = z * exp(g * z);

    // apply products
    for (auto n = 1u; n <= Nterms; ++n)
        retVal *= (1 + z / n) * exp(-z / n);

    return 1.0 / retVal; // invert
}

template<typename T>
void _shuffle(std::vector<T> &perm1, std::vector<T> &perm2) {
    std::random_device rd;
    srand(rd());

    for (auto i = 1u; i < perm1.size(); i++) {
        std::swap(perm1[i], perm1[rand() % (i + 1)]);
        std::swap(perm2[i], perm2[rand() % (i + 1)]);
    }
}

