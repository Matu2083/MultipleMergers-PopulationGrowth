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



#ifndef myTree_hpp
#define myTree_hpp

#include <iostream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <cmath>
#include <boost/math/special_functions/binomial.hpp>

#include "myNode.hpp"
#include "myRandomNumbers.hpp"



class myTree{
	int mSampleSize;
	int mNoSegSites;
	double mTMRC;
	double mTotalBranchlength;
	double mTimeScale;
	std::string mNewick;
	std::vector<double> mBranchlengths;
	std::vector<int> mNumberOfDescendents;
	std::vector<double> mBranchlengthsCumulative;
	std::vector<int> mSFS;
	std::vector<std::vector<double> > mMerger;
	std::vector<std::vector<bool> > mSequences;
	
	void mTreeWalk(int node_lable, std::vector<myNode>& AllNodes,long double rho, double theta, double gamma);
public:
	
	// PSI
	myTree(int SampleSize, double Psi, double gamma, double mu, long double rho, unsigned int N);
	// MORAN
	myTree(int SampleSize, double Psi, double gamma, double mu, unsigned int N, long double rho);
	// NEWICK
	myTree(std::string Newick, double Psi, double gamma, double mu, unsigned int N, long double rho);	
	
	/////////////
	// getter //
	///////////
	int getSampleSize();
	int getNoSegSites();
	std::string getNewick();
	std::vector<double> getBranchlengths();
	std::vector<int> getNumberOfDescendents();
	std::vector<int> getSFS();
	std::vector<double> getBranchlengthsCumulative();
	double getTMRCA();
	double getTotalBranchlength();
	std::vector<std::vector<bool> > getSequence();
	std::vector<std::vector<double> > getMerger();
	
};

template <typename _typeName>
std::string convert_NumberToStringTree ( _typeName Number );

#endif /* myTree_hpp */