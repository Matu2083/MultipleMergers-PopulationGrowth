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



#include "myExpectedSFS.hpp"


///////////////
// helpers  //
/////////////

mpreal g = mpfr::const_euler(mpreal::get_default_prec(), mpreal::get_default_rnd());

mpreal expint(mpreal x)
{
	mpreal result;
	mpreal one = 1.0;
	mpreal zero = 0.0;
	
	if (x >= 0)
	{
		result = eint(x);
	}
	else
	{
		mpreal eps = mpfr::machine_epsilon();
		if (x > 0) // :: LOOK ::
		{
			mpreal xn = -x;
			mpreal Sn = -x;
			mpreal Sm1 = zero;
			mpreal hsum = one;
			mpreal y = one;
			mpreal factorial = one;
			mpreal epsilon = one;
			if ( x == zero )
			{
				mpfr_free_cache();
				return ( (long double) -DBL_MAX);
			}
			
			while ( abs(Sn - Sm1) > eps * abs(Sm1) )
			{
				Sm1 = Sn;
				++y;
				xn *= (-x);
				factorial *= y;
				hsum += (one / y);
				Sn += hsum * xn / factorial;
			}
			result = (g + log(abs(x)) - exp(x) * Sn);
		}
		else
		{
			
			mpreal Am1 = one;
			mpreal A0 = zero;
			mpreal Bm1 = zero;
			mpreal B0 = one;
			mpreal a = exp(x);
			mpreal b = -x + one;
			mpreal Ap1 = b * A0 + a * Am1;
			mpreal Bp1 = b * B0 + a * Bm1;
			int j = 1;
			
			a = one;
			while ( abs(Ap1 * B0 - A0 * Bp1) > eps * abs(A0 * Bp1) )
			{
				if ( abs(Bp1) > 1.0 )
				{
					Am1 = A0 / Bp1;
					A0 = Ap1 / Bp1;
					Bm1 = B0 / Bp1;
					B0 = one;
				}
				else
				{
					Am1 = A0;
					A0 = Ap1;
					Bm1 = B0;
					B0 = Bp1;
				}
				a = -j * j;
				b += 2.0*one;
				Ap1 = b * A0 + a * Am1;
				Bp1 = b * B0 + a * Bm1;
				++j;
			}
			result = (-Ap1 / Bp1);
		}
	}
	
	mpfr_free_cache();
	return (result);
}


mpreal n_choose_k_boost(int N,  int k)
{
	if(N < k || N < 0 || k < 0)
	{
		return (0);
	}
	else
	{
		return (boost::math::binomial_coefficient<long double>(N, k));
	}
}

mpreal n_choose_k_(unsigned int n, unsigned int k, MatrixXmp &B)
{
	if (B(n,k) >= 0)
	{
		return B(n,k);
	}
	if (k == 0)
	{
		return 1;
	}
	if (n <= 60)
	{
		return n_choose_k_boost(n,k);
	}
	else
	{
		mpreal kk = k;
		mpreal nn = n;
		mpreal value = (nn * n_choose_k_(n - 1, k - 1, B)) / kk;
		B(n,k) = value;
		return value;
	}
}


void computeInverse1stCol(MatrixXmp &U, MatrixXmp &diagA, VectorXmp &compInvA)
{
	long SampleSize = U.rows();
	compInvA(0) = 1.;
	
	for (int i = 1; i < SampleSize; ++i)
	{
		compInvA(i) = 0.;
		for (int j = 0; j < i; ++j)
		{
			compInvA(i) -= U(i,j)*compInvA(j);
		}
	}
	
	diagA = compInvA.asDiagonal();
	mpfr_free_cache();
}

mpreal expIntModified(mpreal lambda, mpreal &rho)
{
	if(rho == 0)
	{
		mpreal one = 1.0;
		return one/lambda;
	}
	
	mpreal x = lambda/rho;
	
	return ( -(exp(x)/rho) * expint(-x) );
}

void buildBinomCoeffMatrix(unsigned int SampleSize, MatrixXmp &BinomCoeffMatrix)
{
	BinomCoeffMatrix.setConstant(-1);
	mpreal value;
	
	for (int j = 0; j < SampleSize + 1; ++j)
	{
		for (int i = 0; i <= j; ++i)
		{
			value = n_choose_k_(j, i, BinomCoeffMatrix); // :: THIS SHOULD NOT BE NECESSARY :: BUT BETTER BE SAVE
			BinomCoeffMatrix(j,i) = value;
		}
	}
	
	mpfr_free_cache();
}


//////////////////
// compute SFS  //
/////////////////


void buildB(unsigned int SampleSize, MatrixXmp &BinomCoeffMatrix, MatrixXmp &B)
{
	mpreal neg_one = -1.0;
	mpreal numerator;
	mpreal two = 2.0;
	mpreal element;
	for (int j = 0; j < SampleSize-1; ++j)
	{
		for (int i = 0; i <= j; ++i)
		{
			numerator = j + two;
			element = (pow(neg_one,i-j))/(numerator) * BinomCoeffMatrix(SampleSize-i-2, j-i) * BinomCoeffMatrix(SampleSize,i+1);
			B(i,j)= element;
		}
	}
	
	mpfr_free_cache();
}


void buildC(unsigned int SampleSize, MatrixXmp &C)
{
	C(0,0) = 2.;
	mpreal two = 2.;
	mpreal element;
	
	for (int k = 1; k < SampleSize-1; ++k)
	{
		element = k + two;
		C(k,k) = element;
		C(k,k-1) = -element;
	}
	
	mpfr_free_cache();
}

void computeCoalescenceRate(double psi, double gamma, int i, MatrixXmp &BinomCoeffMatrix, VectorXmp &lambda)
{
	mpreal Psi = psi;
	mpreal Psi_sq = pow(Psi, 2);
	mpreal xPrime;
	mpreal one = 1.;
	mpreal two = 2.;
	
	if (gamma == 2)
	{
		lambda(1) = exp( log(BinomCoeffMatrix(i, 2)) + log(two/(two+Psi_sq) + Psi_sq/(two+Psi_sq)*pow(one-Psi,i-two)));
		lambda(0) = lambda(1);
		for (int x = 3; x <= i; ++x)
		{
			xPrime = x;
			lambda(x-1) = exp( log(BinomCoeffMatrix(i, x)) - log(two+Psi_sq) + xPrime*log(Psi) + (i-xPrime)*log(one-Psi) );
			lambda(0) += lambda(x-1);
		}
	}
	if (gamma < 2)
	{
		lambda(1) = exp( log(BinomCoeffMatrix(i, 2)) + (i-two)*log(one-Psi) );
		lambda(0) = lambda(1);

		for (int x = 3; x <= i; ++x)
		{
			xPrime = x;
			lambda(x-1) = exp( log(BinomCoeffMatrix(i, x)) + (xPrime-two)*log(Psi) + (i-xPrime)*log(one-Psi) );
			lambda(0) += lambda(x-1);
		}
	}
	else
	{
		lambda(1) = BinomCoeffMatrix(i, 2);
		lambda(0) = lambda(1);
		for (int x = 3; x <= i; ++x)
		{
			lambda(x-1) = 0;
		}
	}
	
	mpfr_free_cache();
}

void buildQ(unsigned int SampleSize, double Psi, double Gamma, MatrixXmp &BinomCoeffMatrix, MatrixXmp &Q, VectorXmp &lambda)
{
	for (int i = 1; i < SampleSize; ++i)
	{
		computeCoalescenceRate(Psi, Gamma, i + 1, BinomCoeffMatrix, lambda);
		Q(i, i) = -lambda(0);
		for (int j = 1; j <= i; ++j)
		{
			Q(i, i-j) = lambda(j);
		}
	}
	
	mpfr_free_cache();
}


void buildU(unsigned int SampleSize, double Psi, double Gamma, MatrixXmp &BinomCoeffMatrix, MatrixXmp &Q, MatrixXmp &U)
{
	mpreal sum;
	mpreal one = 1.;
	mpreal entry;
	
	for (int j = 0; j < SampleSize; ++j)
	{
		U(j, j) =  one;
		
		for (int i = j+1; i < SampleSize; ++i)
		{
			sum = 0.0;
			for (int k = j; k < i; ++k)
			{
				sum += Q(i, k)*U(k, j);
			}
			
			entry = one / (-Q(i, i) + Q(j, j)) * sum;
			U(i, j) =  entry;
		}
	}
	
	mpfr_free_cache();
}


void build_coalescentTime(unsigned int SampleSize, double Psi, double Gamma, double rho, MatrixXmp &BinomCoeffMatrix, MatrixXmp &Q, MatrixXmp &U, MatrixXmp &coalD, MatrixXmp &coalV, MatrixXmp &coalL, VectorXmp &compInvA, VectorXmp &c, VectorXmp &a)
{
	mpreal rho_tilde;
	if (Gamma > 2)
	{
		rho_tilde = rho;
	}
	else
	{
		rho_tilde = rho * Gamma;
	}
	
	mpreal coalescentTime;
	
	for (int l = 0; l < SampleSize - 1; ++l)
	{
		coalescentTime = expIntModified( -Q(l+1, l+1), rho_tilde);
		c(l) = coalescentTime;
	}
	
	a = coalL * c;
	
	mpfr_free_cache();
}


std::vector<double> computePhi(unsigned int SampleSize, long double &totalTreeLength, double Psi, double Gamma, double rho, MatrixXmp &A, MatrixXmp &BinomCoeffMatrix, MatrixXmp &Q, MatrixXmp &U, MatrixXmp &coalD, MatrixXmp &coalV, MatrixXmp &coalL, VectorXmp &compInvA, VectorXmp &c, VectorXmp &a, VectorXmp &tau, bool &success)
{
	std::vector<double> _Phi (SampleSize-1);
	if (Psi != 1)
	{
		build_coalescentTime(SampleSize, Psi, Gamma, rho, BinomCoeffMatrix, Q, U, coalD, coalV, coalL, compInvA, c, a);
		tau = A * a;
		unsigned long N = tau.size();
		long double element;
		
		long double normalize = 0.;
		for (int i = 0; i < N; ++i)
		{
			element = (long double) tau(i);
			if (element < 0)
			{
				success = false;
			}
			normalize += element;
			_Phi[i] = element;
		}
		std::transform(_Phi.begin(), _Phi.end(), _Phi.begin(), std::bind1st(std::multiplies<double>(),1.0/normalize));
		totalTreeLength = normalize;
	}
	else
	{
		mpreal rho_tilde;
		if (Gamma > 2)
		{
			rho_tilde = rho;
		}
		else
		{
			rho_tilde = rho * Gamma;
		}

		totalTreeLength = ((long double) expIntModified( 1, rho_tilde)) * SampleSize;
		 _Phi[0] = 1;
		for(int i = 1; i < SampleSize - 1; i++)
		{
			_Phi[i] = 0;
		}
	}
	
	mpfr_free_cache();
	
	return (_Phi);
}


