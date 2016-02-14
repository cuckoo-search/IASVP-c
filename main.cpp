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

#include <iostream>
#include <chrono>
#include <ratio>

#include <CuckooSearch.h>
#include <Funtions.h>
#include <IASVP.h>

#include <Problem.h>

int main(int argc, char *argv[]) {

    if (argc < 2) {
        std::cout << "./cuckoo-search <pos>" << std::endl;
        return EXIT_SUCCESS;
    }

    int pos = std::numeric_limits<int>::max();
    try { pos = std::atoi(argv[1]); } catch (...) { }

    if (pos > 14) {
        std::cout << "0 <= pos < 15" << std::endl;
        return EXIT_SUCCESS;
    }

    std::string test[] = {"input/c1x10", //0
                          "input/c1x20", //1
                          "input/c1x30", //2
                          "input/c1x40", //3
                          "input/c1x50", //4
                          "input/c2x10", //5
                          "input/c2x20", //6
                          "input/c2x30", //7
                          "input/c2x40", //8
                          "input/c2x50", //9
                          "input/c3x10", //10
                          "input/c3x20", //11
                          "input/c3x30", //12
                          "input/c3x40", //13
                          "input/c3x50"}; //14

    uint nds[] = {10u, 20u, 30u, 40u, 50u,
                  10u, 20u, 30u, 40u, 50u,
                  10u, 20u, 30u, 40u, 50u};

    auto nd = nds[pos];
    const auto lb = -32.0;
    const auto ub = 32.0;
    const double tol = 1.0e-5;
    const uint eggs = 25u;
    const float pa = 0.25f;

    auto seed = load(test[pos], nd, 1);
    IASVP iasvp(seed, makeToeplitz);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lb, ub);

    const auto fn = [&iasvp](const auto &p) { return iasvp.FIASVPToeplitzTriInf(p.solution); };
    const auto fn_gen = [&gen, &dis]() { return dis(gen); };
    const auto stop = [&tol](const auto &p) { return p.fitness < tol; };

    CuckooSearch<Problem> cs(eggs, nd, lb, ub, pa, fn, fn_gen, stop);

    auto start = std::chrono::system_clock::now();

    auto p = cs.search({std::make_unique<GetCuckoos<Problem>>(),
                        std::make_unique<BestNest<Problem>>(),
                        std::make_unique<HybridEmptyNest<Problem>>(iasvp),
                        std::make_unique<BestNest<Problem>>()});

    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration<double>(end - start).count();

    //printf("Elapsed Time,Fitness,R. Error,Iterations,ND\n");
    printf("%lf,%e,%e,%d,%d\n", elapsed, p.fitness, iasvp.RelativeError(p.solution), cs.niter, nd);

    return EXIT_SUCCESS;
}
