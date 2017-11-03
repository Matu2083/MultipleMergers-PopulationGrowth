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

#include <iostream>
#include <fstream>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"


#include "myIOParser.hpp"
#include "myExpectedSFS.hpp"
#include "Eigen/Dense"
#include "mpreal.h"
#include "Eigen/MPRealSupport"

void calculateExpectedSFS();

int main(int argc, const char * argv[])
{
	try
	{
		readcommandLineArguments(argc, argv);
		checkArgvInput();
	}
	catch (std::string e)
	{
		std::cerr << "Error: " << e << std::endl;
		printHelpMessage();
		return EXIT_FAILURE;
	}
	
	set_filename(outfile_name);
	initiateoutput();
	calculateExpectedSFS();
	
	return (0);
}

void calculateExpectedSFS()
{
	MatrixXmp BinomCoeffMatrix(sampleSize + 1, sampleSize + 1);
	MatrixXmp A(sampleSize - 1, sampleSize - 1);
	MatrixXmp coalD(sampleSize - 1, sampleSize - 1); // :: FOR MEMORY REASONS 'MISUSED' AS B MATRIX
	MatrixXmp coalV(sampleSize - 1, sampleSize - 1); // :: FOR MEMORY REASONS 'MISUSED' AS C MATRIX
	MatrixXmp Q(sampleSize, sampleSize);
	MatrixXmp U(sampleSize, sampleSize);
	MatrixXmp coalL(sampleSize - 1, sampleSize - 1);
	
	VectorXmp lambda(sampleSize);
	VectorXmp compInvA(sampleSize);
	VectorXmp a(sampleSize - 1);
	VectorXmp c(sampleSize - 1);
	VectorXmp tau (sampleSize - 1);
	
	buildBinomCoeffMatrix(sampleSize, BinomCoeffMatrix);
	buildB(sampleSize, BinomCoeffMatrix, coalD);
	buildC(sampleSize, coalV);
		
	A = coalD * coalV;
	
	std::vector<double> Phi;
	double Psi, Rho, gamma;
	long double totalTreeLength;
	
	for (int psiIndex = noStepsPsi; psiIndex >= 0; --psiIndex)
	{
		if(noStepsPsi == 0)
		{
			Psi = minPsi;
		}
		else
		{
			Psi = minPsi + psiIndex * (maxPsi - minPsi)/(noStepsPsi);
		}
		
		for (int rhoIndex = noStepsRho; rhoIndex >= 0 ; --rhoIndex)
		{
			if(noStepsRho == 0)
			{
				Rho = minRho;
			}
			else
			{
				Rho = minRho + rhoIndex * (maxRho - minRho)/(noStepsRho);
			}
			
			totalTreeLength = 0;
			
			gamma = 2.5;
			if (Psi > 0)
			{
				gamma = 1.5;
			}
			
			buildQ(sampleSize, Psi, gamma, BinomCoeffMatrix, Q, lambda);
			buildU(sampleSize, Psi, gamma, BinomCoeffMatrix, Q, U);
			computeInverse1stCol(U, coalD, compInvA);
			coalV = U * coalD;
			coalL = -coalV.block(1, 1, sampleSize - 1, sampleSize - 1);
						
			if(gamma > 2)
			{
				Phi = computePhi(sampleSize, totalTreeLength, Psi, gamma, Rho, A, BinomCoeffMatrix, Q, U, coalD, coalV, coalL, compInvA, c, a, tau, success);
			}
			else
			{
				Phi = computePhi(sampleSize, totalTreeLength, Psi, gamma, Rho/gamma, A, BinomCoeffMatrix, Q, U, coalD, coalV, coalL, compInvA, c, a, tau, success);
			}
			
			if (success == false)
			{
				std::cout << "Error: computeRelativeBranchLengths. Negative values encountered. Please increase precision. See ReadMe file for input arguments.\n";
				exit(1);
			}
			printPhi(Phi, totalTreeLength);
		}
	}
	
	_offPhi.close();
	
}
