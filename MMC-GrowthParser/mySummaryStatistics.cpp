/*
 * MMC_Growth is used to simulate gene trees under the psi- and kingman
 * coalescent process with exponential growth.
 *
 * Copyright (C) 2016 Sebastian Matuszewski & Marcel Hildebrandt
 *
 * This file is part of MMC_Growth.
 *
 * MMC_Growth is free software: you can redistribute it and/or modify
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

#include "mySummaryStatistics.hpp"


////////////////////
// help functions //
////////////////////

double compute_pHarmonicNumber(unsigned int n, double p);

////////////////////
//   constructor  //
////////////////////


mySummaryStatistics::mySummaryStatistics(myTree tree)
{
    unsigned int segSites = tree.getNoSegSites();
	if(segSites == 0)
	{
		mTajimasD = 0;
		mFuAndLisD = 0;
		mFayAndWusH	= 0;
	}
	else
	{
		unsigned int sampleSize = tree.getSampleSize();
		std::vector<int> SFS = tree.getSFS();
		double h = compute_pHarmonicNumber(sampleSize-1,1);
		double g = compute_pHarmonicNumber(sampleSize-1,2);
		
		/////////////////////
		//   Tajima's D   //
		///////////////////
		
		mThetaW = segSites / h;
		mThetaT = 0;
		for (int i = 1 ; i < sampleSize; ++i)
		{
			mThetaT += i*(sampleSize-i)*SFS[i-1];
		}
		mThetaT *= 2./(sampleSize*(sampleSize-1));
		double b1 = ((double) sampleSize+1. )/(3.*(sampleSize-1.));
		double b2 = ( 2*(std::pow(sampleSize, 2)+sampleSize +3.) )/(9.*sampleSize*(sampleSize-1.));
		double c1 = b1-1./h;
		double c2 = b2 - (sampleSize+2)/(h*sampleSize)+g/std::pow(h, 2);
		double e1 = c1/h;
		double e2 = c2/(std::pow(h, 2)+ g);
		
		mTajimasD = (mThetaT-mThetaW)/(std::sqrt(e1*segSites + e2*segSites*(segSites-1)));
		
		////////////////////////
		//   Fu and Li's D   //
		//////////////////////
		
		double c = (2.*sampleSize*h-4.*(sampleSize-1.))/((sampleSize-1.)*(sampleSize-2.));
		double v = 1.+std::pow(h, 2)/(g+std::pow(h, 2)) *(c-(sampleSize+1.)/(sampleSize-1.));
		double u = h - 1 - v;
		mFuAndLisD = (segSites-h*SFS[0])/std::sqrt(u*segSites+v*std::pow(segSites, 2));
		
		//////////////////////////
		//     Fay and Wu's D  //
		////////////////////////
		
		mThetaH = 0;
		for (int i = 1; i < sampleSize-1; ++i)
		{
			mThetaH += SFS[i-1]*std::pow(i, 2);
		}
		mThetaH *= 2./(sampleSize*(sampleSize-1.));
		mFayAndWusH = mThetaT - mThetaH;
	}
}


/////////////
// getter //
///////////



double mySummaryStatistics::getTheta_W()
{
    return this->mThetaW;
}

double mySummaryStatistics::getThetaT()
{
    return this->mThetaT;
}
double mySummaryStatistics::getTajimasD()
{
    return this->mTajimasD;
}
double mySummaryStatistics::getFuAndLisD()
{
    return this->mFuAndLisD;
}

double mySummaryStatistics::getThetaH()
{
    return this->mThetaH;
}

double mySummaryStatistics::getFayAndWusH()
{
    return this->mFayAndWusH;
}




////////////////////
// help functions //
////////////////////

double compute_pHarmonicNumber(unsigned int n, double p)
{
    double H = 1;
    for (int i = n; i > 0; --i)
    {
        H *= exp(1/(std::pow(i, p)));
    }
    return log(H);
}