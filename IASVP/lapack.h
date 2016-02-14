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

/**
 * Common constants
 * */
int iminus_one = -1, izero = 0, ione = 1;
double dminus_one = -1.0, dzero = 0.0, done = 1.0;
char TRANS = 'T', NOTRANS = 'N';
char UPPER = 'U', LOWER = 'L';
char DIAG = 'U', NODIAG = 'N';

extern "C" {

/**
 * LAPACK
 * */
void dgetrf_( int*, int*, double*, int*, int*, int* );
void dgetrs_( char* TRANS, int* n, int* m, double* A, int* lda, int* ipiv, double *b, int* ldb, int* info );

void dgesvd_( char *jobu, char *jobvt, int *m, int *n, double *A, int *lda, 
              double *S, double *U, int *ldu, double *VT, int *ldvt, double *work, int *lwork, 
              int *info );
}

