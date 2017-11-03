/*
 * MMC-Phi-MaxLikelihood-LookUp is used to calculate the expected SFS and expected
 * relative branch length for the Psi-coalescent with exponential growth
 * and its contained subclasses.
 *
 * Copyright (C) 2016 Sebastian Matuszewski & Marcel Hildebrandt
 *
 * This file is part of MMC-Phi-MaxLikelihood-LookUp.
 *
 * MMC-ExpectedSF is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef mExpectedSFS_hpp
#define mExpectedSFS_hpp

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <boost/math/special_functions/binomial.hpp>
#include <vector>
#include <iostream>
#include <gsl/gsl_sf_expint.h>
#include "Eigen/Dense"
#include "Eigen/MPRealSupport"
#include "mpreal.h"

using mpfr::mpreal;

extern mpreal g;

// Declare matrix and vector types with multi-precision scalar type
typedef Eigen::Matrix<mpreal,Eigen::Dynamic,Eigen::Dynamic> MatrixXmp;
typedef Eigen::Matrix<mpreal,Eigen::Dynamic,1> VectorXmp;

void buildBinomCoeffMatrix(unsigned int SampleSize, MatrixXmp &BinomCoeffMatrix);
void buildC(unsigned int SampleSize, MatrixXmp &C);
void buildB(unsigned int SampleSize, MatrixXmp &BinomCoeffMatrix, MatrixXmp &B);
void buildQ(unsigned int SampleSize, double Psi, double Gamma, MatrixXmp &BinomCoeffMatrix, MatrixXmp &Q, VectorXmp &lambda);
void buildU(unsigned int SampleSize, double Psi, double Gamma, MatrixXmp &BinomCoeffMatrix, MatrixXmp &Q, MatrixXmp &U);
void computeInverse1stCol(MatrixXmp &U, MatrixXmp &diagA, VectorXmp &compInvA);

std::vector<double> computePhi(unsigned int SampleSize, long double &totalTreeLength, double Psi, double Gamma, double rho, MatrixXmp &A, MatrixXmp &BinomCoeffMatrix, MatrixXmp &Q, MatrixXmp &U, MatrixXmp &coalD, MatrixXmp &coalV, MatrixXmp &coalL, VectorXmp &compInvA, VectorXmp &c, VectorXmp &a, VectorXmp &tau, bool &success);

#endif /* mExpectedSFS_hpp */