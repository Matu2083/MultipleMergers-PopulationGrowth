/*
 * MMC-MaxLikelihood-Grid is used to estimate the MMC parameter psi and/or the
 * exponential growth parameter from SFS data.
 *
 * Copyright (C) 2016 Sebastian Matuszewski & Marcel Hildebrandt
 *
 * This file is part of MMC-MaxLikelihood-Grid.
 *
 * MMC-MaxLikelihood-Grid is free software: you can redistribute it and/or modify
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


#include "myEstimationMethod.hpp"

////////////////////
//   grid search  //
////////////////////

myEstimationMethod::myEstimationMethod(std::vector<std::vector<std::vector<double>>> &Phi,std::vector<std::vector<double>> &T_TotalGrid, std::vector<double> SFS, double MinPsi, double MaxPsi, int noStepsPsi, double MinRho, double MaxRho, int noStepsRho, int lumping):
    mLLGrid(noStepsPsi + 1, std::vector<double>(noStepsRho + 1, -std::numeric_limits<double>::infinity())), ml1Grid_abs(noStepsPsi + 1, std::vector<double>(noStepsRho + 1, std::numeric_limits<double>::infinity())), ml2Grid_abs(noStepsPsi + 1, std::vector<double>(noStepsRho + 1, std::numeric_limits<double>::infinity())), mMaxLogLikelihood (- std::numeric_limits<double>::infinity()), mMaxLogLikelihood_logScale(- std::numeric_limits<double>::infinity()), mRhoRange(noStepsRho + 1), mPsiRange(noStepsPsi + 1), mMinl1_abs(std::numeric_limits<double>::infinity()), mMinl2_abs(std::numeric_limits<double>::infinity()),mThetaML(NAN), mThetaML_logScale(NAN), mThetal2_abs(NAN), mThetal1_abs(NAN)
{
	int sampleSize = (int) SFS.size() + 1;
	
	if (lumping != 0)
	{
		while (SFS.size() > lumping)
		{
			SFS[lumping-1] += SFS.back();
			SFS.pop_back();
		}
	}
	
	if (noStepsRho == 0)
	{
		mRhoRange[0] = MinRho;
	}
	else
	{
		for (int i = 0; i < noStepsRho + 1; ++i)
		{

			mRhoRange[i] = MinRho + i*(MaxRho - MinRho)/(noStepsRho);
		}
	}
	
	if (noStepsPsi == 0)
	{
		mPsiRange[0] = MinPsi;
	}
	else
	{
		for (int j = 0; j < noStepsPsi + 1; ++j)
		{
			mPsiRange[j] = MinPsi + j*(MaxPsi - MinPsi)/(noStepsPsi);
		}
	}
	
    int cutoff;
    if (lumping == 0)
    {
        cutoff = sampleSize-1;
    }
    else
    {
        cutoff = lumping;
    }
	
    unsigned int segSites = 0;
    std::vector<double> relBranchLength_observed(cutoff);
	
	for (int j = 0; j < cutoff; ++j)
    {
        segSites += SFS[j];
        relBranchLength_observed[j] = SFS[j];
    }
	
	if (lumping == 0)
	{
		if (folded == true)
		{
			cutoff = (sampleSize) / 2;
		}
		else
		{
			cutoff = sampleSize - 1;
		}
	}
	else
	{
		if (folded == true)
		{
			cutoff = (sampleSize) / 2;
			cutoff = std::min(cutoff, lumping);
		}
		else
		{
			cutoff = lumping;
		}
	}

    std::transform(relBranchLength_observed.begin(), relBranchLength_observed.end(), relBranchLength_observed.begin(), std::bind1st(std::multiplies<double>(),1.0/((double) segSites)));
	
	std::vector<double> _Phi;
	std::vector<double> _passPhi(cutoff, 0.);
	std::vector<double> _SFS;
	std::vector<double> _passSFS(cutoff, 0.);
	double tTotal, theta;
	int index;
	
	for (int k = 0; k < noStepsPsi + 1; ++k)
	{
        int m = 0;
		for (int l = 0; l < noStepsRho + 1; ++l)
		{

            tTotal = T_TotalGrid[k][l];
            theta = segSites / tTotal;
            _Phi = Phi[k][l];
            _SFS = _Phi;
			
			std::fill(_passPhi.begin(), _passPhi.end(), 0.);
			std::fill(_passSFS.begin(), _passSFS.end(), 0.);
            std::transform(_SFS.begin(), _SFS.end(), _SFS.begin(), std::bind1st(std::multiplies<double>(), segSites));
			
			for (int j = 0; j < sampleSize-1; ++j)
			{
				index = std::min(cutoff-1, j);
				_passPhi[index] += _Phi[j];
				_passSFS[index] += _SFS[j];
				
			}
			
            mLLGrid[k][l] = computeLogLikelihood(_passPhi, sampleSize, SFS, segSites, cutoff);
            ml1Grid_abs[k][l] = computeDistance(_passSFS, sampleSize, SFS, segSites, cutoff, 1.0);
            ml2Grid_abs[k][l] = computeDistance(_passSFS, sampleSize, SFS, segSites, cutoff, 2.0);
			
            if (l == ((int) std::pow(2, m)) || l == 0 )
            {
                if (l > 0)
                {
                ++m;
                }
                if (mMaxLogLikelihood_logScale < mLLGrid[k][l])
                {
                    mMaxLogLikelihood_logScale = mLLGrid[k][l];
                    mRhoML_logScale = mRhoRange[l];
                    mPsiML_logScale = mPsiRange[k];
                    mThetaML_logScale = theta;
                }
            }
            
			if (mMaxLogLikelihood < mLLGrid[k][l])
			{
				mRhoML = mRhoRange[l];
				mPsiML = mPsiRange[k];
				mMaxLogLikelihood = mLLGrid[k][l];
				mNoSteps = 1;
				mThetaML = theta;
			}
			
            if (mMinl1_abs > ml1Grid_abs[k][l])
            {
                mRhol1_abs = mRhoRange[l];
                mPsil1_abs = mPsiRange[k];
                mMinl1_abs = ml1Grid_abs[k][l];
                mNoSteps = 1;
				mThetal1_abs = theta;
            }
            if (mMinl2_abs > ml2Grid_abs[k][l])
            {
                mRhol2_abs = mRhoRange[l];
                mPsil2_abs = mPsiRange[k];
                mMinl2_abs = ml2Grid_abs[k][l];
                mNoSteps = 1;
				mThetal2_abs = theta;
            }
			
		}
	}
}


/////////////////////////
//    help functions   //
/////////////////////////

double computeLogLikelihood(std::vector<double> &_passPhi,  int sampleSize, std::vector<double> &SFS, int segSites, int cutoff)
{
		double l = 0;
		for (int k = 0; k < cutoff; ++k)
		{
			l += -_passPhi[k] * segSites + std::log( segSites * _passPhi[k] ) * SFS[k];
		}
		return (l);
}


double computeDistance(std::vector<double> &expected, int sampleSize, std::vector<double> &observed, int segSites, int cutoff, double p)
{
    double d = 0;
    for (int k = 0; k < cutoff; ++k)
    {
        //d += std::pow( (std::abs(observed[k] - expected[k]) ),p) / expected[k];
		d += std::pow( (std::abs(observed[k] - expected[k]) ),p);
    }
    return (std::pow(d, 1.0 / p));
}


/////////////
// getter //
///////////

double myEstimationMethod::getMaxLogLikelihood()
{
	return (this->mMaxLogLikelihood);
}

std::vector<double> myEstimationMethod::getRhoRange()
{
	return (this->mRhoRange);
}

std::vector<double> myEstimationMethod::getPsiRange()
{
	return (this->mPsiRange);
}

std::vector<std::vector<double> > myEstimationMethod::getLLGrid()
{
	return (this->mLLGrid);
}

std::vector<std::vector<double> > myEstimationMethod::getGrid_Psi()
{
	return (this->mGrid_Psi);
}

std::vector<std::vector<bool> > myEstimationMethod::getGrid_Success()
{
	return (this->mSuccessGrid);
}

std::vector<std::vector<double> > myEstimationMethod::getGrid_Rho()
{
	return (this->mGrid_Rho);
}

bool myEstimationMethod::getSuccess()
{
	return (this->mSuccess);
}

double myEstimationMethod::getRhoML()
{
	return (this->mRhoML);
}

double myEstimationMethod::getPsiML()
{
	return (this->mPsiML);
}

double myEstimationMethod::getRhoML_logScale()
{
	return (this->mRhoML_logScale);
}

double myEstimationMethod::getPsiML_logScale()
{
	return (this->mPsiML_logScale);
}

double myEstimationMethod::getMaxLogLikelihood_logScale()
{
	return (this->mMaxLogLikelihood_logScale);
}

unsigned long myEstimationMethod::getNoSteps()
{
	return (this->mNoSteps);
}

double myEstimationMethod::getMinl2_abs()
{
	return (this->mMinl2_abs);
}

double myEstimationMethod::getMinl1_abs()
{
	return (this->mMinl1_abs);
}

double myEstimationMethod::getPsil2_abs()
{
	return (this->mPsil2_abs);
}

double myEstimationMethod::getPsil1_abs()
{
	return (this->mPsil1_abs);
}

double myEstimationMethod::getRhol2_abs()
{
	return (this->mRhol2_abs);
}

double myEstimationMethod::getRhol1_abs()
{
	return (this->mRhol1_abs);
}

std::vector<std::vector <bool > > myEstimationMethod::getSuccessGrid()
{
	return (this->mSuccessGrid);
}

std::vector<std::vector<double> > myEstimationMethod::getml2Grid_abs()
{
	return (this->ml2Grid_abs);
}

std::vector<std::vector<double> > myEstimationMethod::getml1Grid_abs()
{
	return (this->ml1Grid_abs);
}

std::vector<double> myEstimationMethod::getLogLikelihood()
{
	return (this->mLogLikelihood);
}

double myEstimationMethod::getThetaML()
{
	return (this->mThetaML);
}

double myEstimationMethod::getThetaML_logScale()
{
	return (this->mThetaML_logScale);
}

double myEstimationMethod::getThetal1_abs()
{
	return (this->mThetal1_abs);
}
double myEstimationMethod::getThetal2_abs()
{
	return (this->mThetal2_abs);
}
