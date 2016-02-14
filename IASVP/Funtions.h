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

#include <cblas.h>

#include <lapack.h>

#include <Utils.h>

using vdouble = std::vector<double>;
using vint = std::vector<int>;
using fn_vdouble_2_vdouble = std::function<vdouble(const vdouble &)>;

void print_matrix(double *A, int m, int n) {
    for (auto i = 0; i < m; i++) {
        for (auto j = 0; j < n; j++) {
            std::cout << A[j * m + i] << " ";
        }
        std::cout << std::endl;
    }
}

void print_matrix(const double *A, int m, int n) {
    for (auto i = 0; i < m; i++) {
        for (auto j = 0; j < n; j++) {
            std::cout << A[j * m + i] << " ";
        }
        std::cout << std::endl;
    }
}

void print_matrix(int *A, int m, int n) {
    for (auto i = 0; i < m; i++) {
        for (auto j = 0; j < n; j++) {
            std::cout << A[j * m + i] << " ";
        }
        std::cout << std::endl;
    }
}

vdouble makeToeplitz(const vdouble &seed) {
    auto n = seed.size();
    vdouble Ac(n * n);
    auto it = Ac.begin();
    auto cit = Ac.begin();

    INNER_MAP(Ac, ([&it, &cit, &seed, n](auto x) {
        auto d = std::distance(cit, it++);
        auto c = static_cast<ulong>(d / n);
        auto e = d % n;
        if (e < c) return 0.0;
        else return seed[e - c];
    }))

    return Ac;
}

//----------------------------------------------------------------------------------------------
void newtonBiseccionNLES(const fn_vdouble_2_vdouble &F, vdouble &seed, const fn_vdouble_2_vdouble &Jac, double rel_tol,
                         double abs_tol, int maxIt, int &it) {
    auto trans = 'N';
    auto n = static_cast<int>(seed.size());
    int i, info;
    double r0, n2fx, n2fnewx;
    vint ipiv(n);
    vdouble s(n);
    vdouble newx(n);

    auto fx = F(seed);
    n2fx = r0 = NORM2(fx)

    it = 0;
    while ((n2fx > (rel_tol * r0 + abs_tol)) && (it < maxIt)) {
        auto J = Jac(seed);
        INNER_MAP_2(s, fx, [](auto a, auto b) { return -b; })

        dgetrf_(&n, &n, J.data(), &n, ipiv.data(), &info);
        dgetrs_(&trans, &n, &ione, J.data(), &n, ipiv.data(), s.data(), &n, &info);

        newx = seed;
        INNER_MAP_2(newx, s, [](auto a, auto b) { return a + b; })
        auto fnewx = F(newx);
        n2fx = NORM2(fx)
        n2fnewx = NORM2(fnewx)

        i = 0;
        while ((n2fnewx - n2fx) > 1.0e-10) {
            INNER_MAP_2(newx, seed, [](auto a, auto b) { return a + b; })
            INNER_MAP(newx, [](auto x) { return x * 0.5; })
            fnewx = F(newx);
            n2fnewx = NORM2(fnewx)
            i++;
        }

        seed = newx;
        fx = F(seed);
        n2fx = NORM2(fx)
        it++;
    }
}

//----------------------------------------------------------------------------------------------
vdouble JacIASVPToeplitzTriInf(const vdouble &seed, const fn_vdouble_2_vdouble &matrixMaker) {
    auto n = static_cast<int>(seed.size());
    auto jobu = 'A';
    auto jobvt = 'A';
    auto lwork = 2 * n * n;
    int info;
    auto alpha = 1.0;
    auto beta = 0.0;
    vdouble J(n * n);
    vdouble work(lwork);
    vdouble Ai(n * n);
    vdouble P(n * n);
    vdouble Q(n * n);
    vdouble s(n);
    vdouble pi;
    pi.reserve(n);
    vdouble qi;
    qi.reserve(n);
    vdouble temp(n);
    auto sumA = matrixMaker(seed);

    dgesvd_(&jobu, &jobvt, &n, &n, sumA.data(), &n, s.data(), P.data(), &n, Q.data(), &n,
            work.data(), &lwork, &info);

    for (auto i = 0; i < n; i++) {
        for (auto j = i; j < n; j++) {
            std::swap(Q[j * n + i], Q[i * n + j]);
        }
    }

    for (auto i = 0; i < n; i++) {
        pi.clear();
        qi.clear();
        std::copy_n(std::begin(P) + i * n, n, std::back_inserter(pi));
        std::copy_n(std::begin(Q) + i * n, n, std::back_inserter(qi));

        for (auto j = 0; j < n; j++) {
            for (auto k = 0; k < n - j; k++) {
                Ai[k * n + j + k] = 1.0;
            }

            cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, Ai.data(), n, qi.data(), ione,
                        beta, temp.data(), ione);

            auto it = temp.begin();
            J[j * n + i] = REDUCE(pi, 0.0, [&it](auto acc, auto x) { return acc + x * (*(it++)); })
            for (auto k = 0; k < n - j; k++)
                Ai[k * n + j + k] = 0.0;
        }
    }

    return J;
}

//----------------------------------------------------------------------------------------------
std::vector<double> load(std::string name, int m, int n) {
    std::string line;
    std::ifstream file(name);
    std::shared_ptr<std::ifstream> filePtr(&file, [](auto f) { if (f->is_open()) f->close(); });
    std::vector<double> A(m * n);
    auto i = 0u;
    auto j = 0u;

    if (file) {
        while (getline(file, line)) {
            if (!line.empty()) {
                std::string token;
                std::istringstream tokens(line);
                j = 0u;
                while (tokens >> token) {
                    A[j * m + i] = std::stod(token);
                    j++;
                }
                i++;
            }
        }
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }

    return A;
}

//----------------------------------------------------------------------------------------------
vdouble CalcSV(const vdouble &seed, const fn_vdouble_2_vdouble &matrixMaker) {
    vdouble sigma(seed.size());
    auto n = static_cast<int>(seed.size());
    auto jobu = 'N';
    auto jobvt = 'N';
    double *U = nullptr, *VT = nullptr;
    auto lwork = 2 * n * n;
    int info;
    vdouble work(lwork);
    auto Ac = matrixMaker(seed);

    dgesvd_(&jobu, &jobvt, &n, &n, Ac.data(), &n, sigma.data(), U, &n, VT, &n, work.data(), &lwork,
            &info);

    return sigma;
}
//----------------------------------------------------------------------------------------------

