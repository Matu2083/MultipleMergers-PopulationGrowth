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



#include "myTree.hpp"

////////////////////
// help functions //
////////////////////

int  removeRandomElement(std::vector<int>& v);
std::vector<long double> computeCoalescenceRates(double Psi, double gamma, int i);
std::vector<long double> computeCoalescenceProbabilities(double Psi, double gamma, unsigned int N, unsigned int i);
long double logbinomial(long double N, long double k);
long double n_choose_k(int N, int k);
void ReadNewick(int& current, std::string Newick, std::vector<myNode> &AllNodes, int & SampleSize);

////////////////////
//   constructor  //
////////////////////



///////////////////////
//   Psi-Coalescent  //
///////////////////////


myTree::myTree(int SampleSize, double Psi, double gamma, double mu, long double rho, unsigned int N)
:mBranchlengthsCumulative(SampleSize-1,0.0), mSFS(SampleSize-1)
{
	
	double time = 0;                                //initial time
	mSampleSize = SampleSize;
	mTimeScale = 1.;
	std::vector<myNode> AllNodes (SampleSize);
	std::vector<int> ActiveNodes (SampleSize);      //keeps tracks of the indices of the active nodes (indices correspond to location in AllNodes)
	for (int i = 0; i < SampleSize; ++i)
	{
		AllNodes[i] = myNode(0, i,1);               // create leaves with time 0, label i, 1 descendent
		AllNodes[i].addDescendent(i);               // add itsself as Descendent
		ActiveNodes[i] = i;                         //activate leaves
	}

	double theta;
	if (gamma <= 1)
	{
	 theta = 0;
	}
	else if (gamma < 2)
	{
		theta = 2.*mu*std::pow(N, gamma-1)/std::pow(Psi, 2.);
	}
	else if (gamma == 2)
	{
		theta = 2.*mu*N*(1./(2.+std::pow(Psi, 2)));
	}
	else
	{
		theta = mu*N;
	}

	
	int i = SampleSize;                             //number of ac active nodes
	int lable = i ;                                 //label for the next node
	while (i > 1)
	{
		
		////////////////////////
		// coalescent events  //
		////////////////////////
		
		std::vector< long double> lambda = computeCoalescenceRates(Psi, gamma, i);
		double dt = randexp(lambda[0]);              //sample the waiting time until the next coalescensce event
		lambda.erase(lambda.begin());
		int x = randdiscrete(lambda) + 2; time += dt;                                 //coalescensce time
		AllNodes.push_back(myNode(time, lable,0));
		mMerger.push_back(std::vector<double>(4));
		mMerger.back()[0] = i;
		mMerger.back()[1] = x;
		mMerger.back()[2] = time;
		mMerger.back()[3] = N;
		
		for (int j = 0; j < x; ++j)
		{
			int random_index = removeRandomElement(ActiveNodes); //pick one active node at random and deactivate/ remove from active set
			AllNodes.back().addChild(AllNodes[random_index]);
		}
		ActiveNodes.push_back(lable); //add new node to active set
		++lable;
		i = int(ActiveNodes.size()); //update size of active nodes
	}
	mTotalBranchlength = 0.;
	mNoSegSites = 0;
	mSequences.resize(mSampleSize);
	for (int j = 0; j < AllNodes.size(); ++j)
	{
		AllNodes[j].sortDescendent();
		if (j < mSampleSize)
		{
			mSequences[j].push_back(0);
		}
		
	}
	mTreeWalk(int(AllNodes.size()-1), AllNodes, rho, theta, gamma);
	mNewick += ";";
	
	if (rho > 0)
	{
		if (gamma < 2)
		{
			mTMRC = std::log(1. + time*rho*gamma)/(rho*gamma); //transformation/ time scaling
		}
		else
		{
			mTMRC = std::log(1. + time*rho)/(rho); //transformation/ time scaling
		}
	}
	else
	{
		mTMRC = time;
	}
}


//////////////
//  MORAN   //
//////////////


myTree::myTree(int SampleSize, double Psi, double gamma, double mu, unsigned int N, long double rho)
:mBranchlengthsCumulative(SampleSize-1,0.0), mSFS(SampleSize-1)
{
    
	mNoSegSites = 0;
	long double temp_popsize = (long double) N;
	mTimeScale = (N*(N-1.))/(2.*(1.-std::pow(N, -gamma))+ std::pow(N, -gamma)*N*Psi*(N*Psi-1.));
    long double Growthrate = rho/mTimeScale;

	double time = 0.;                                    //initial time
	mSampleSize = SampleSize;
	mSequences.resize(mSampleSize);
	std::vector<myNode> AllNodes (SampleSize);
	std::vector<int> ActiveNodes (SampleSize);          // keeps tracks of the indices of the active nodes (indices                  correspond to location in AllNodes)
	for (int i = 0; i < SampleSize; ++i)
	{
		AllNodes[i] = myNode(0, i, 1);                   // create leaves with time 0, label i, 1 descendent
		AllNodes[i].addDescendent(i);					// add itsself as Descendent
		ActiveNodes[i] = i;                             // activate leaves
		mSequences[i].push_back(0);
	}
	
	unsigned int i = SampleSize;
	int lable = i ;                                     // label for the next new node
	temp_popsize = N;

	while (i > 1)
	{                                     // iterate as long as more than 1 node active
		++time;
		
		//////////////////////////
		// coalescensce events  //
		//////////////////////////
		
		int prevGenPopsize = std::ceil(temp_popsize*std::exp(-Growthrate));
		int deltaN = N-prevGenPopsize;
		double probSweepstake = std::pow(std::max((int) prevGenPopsize,(int)i), -gamma);  //prob. of sweepstake event
        bool sweepstake =  randbern(probSweepstake);                                  //draw reproduction event
		int NumberOfChildren;
		if (sweepstake)
		{
			NumberOfChildren = std::max(deltaN + 1, (int) std::floor(prevGenPopsize*Psi));
		}
		else
		{
			NumberOfChildren = std::max(2,deltaN + 1) ;
		}
		int x_merger = 0;                                                //counts the number of nodes that merge
		int TotalPopulation = N;
		double NumberOfActiveNodes = ((double)ActiveNodes.size());
		int first_index_active = 0;                         // index of the first child in active nodes
		int first_index_all = 0;                            // index of the first child in AllNodes
		for (int j = 0; j < NumberOfChildren && ActiveNodes.size() != 0; ++j)
		{
			
			double prob = NumberOfActiveNodes/TotalPopulation;      //prob. that child is active
			if (randbern(prob))                                     //child active
			{
				
				--NumberOfActiveNodes;                              //removes element from active set
				++x_merger;
				if (x_merger == 1)
				{
					int n = int(ActiveNodes.size());
					first_index_active = randint(n);
					first_index_all = ActiveNodes[first_index_active];
					if (randbern(mu) && j > 0)                  // true if mutates and not the "first child" = parent
					//if (randbern(mu))							// Alternative mutation model: True if mutates. Parents die such that also the first child can carry a mutation.
					{
						mSFS[AllNodes[first_index_all].getNumberOfDescendents()-1] += 1;
						mNoSegSites ++;
						int counterDescendents = 0;
						for (int k = 0; k< mSampleSize; ++k)
						{
							AllNodes[first_index_all].sortDescendent();
							if (AllNodes[first_index_all].getDescendents()[counterDescendents] == k)
							{
								mSequences[k].push_back(1);         // if descendent at 1
								++counterDescendents;
							}
							else
							{
								mSequences[k].push_back(0);        // else add a 0
							}
						}
						
					}
				}
				
				if (x_merger == 2 )													//coalescence event (at least binary merger)
				{
					ActiveNodes.erase(ActiveNodes.begin() +first_index_active);		// deactivate first child
					AllNodes.push_back(myNode(time, lable,0));						//create new node
					AllNodes.back().addChild(AllNodes[first_index_all]);
					
				}
				if (x_merger > 1)
				{
					int random_index = removeRandomElement(ActiveNodes); //pick one active node at random and deactivate/ remove from active set
					AllNodes.back().addChild(AllNodes[random_index]);
					if (randbern(mu))                  // true if mutates
					{
						mSFS[AllNodes[random_index].getNumberOfDescendents()-1] += 1;
						mNoSegSites++;
						AllNodes[random_index].sortDescendent();
						int counterDescendents = 0;
						for (int k = 0; k< mSampleSize; ++k)
						{
							if (AllNodes[random_index].getDescendents()[counterDescendents] == k)
							{
								mSequences[k].push_back(1);
								++counterDescendents;
							}
							else
							{
								mSequences[k].push_back(0);
							}
						}
						
					}
				}
			}
			--TotalPopulation;           //one element less to choose from
		}
		if (x_merger > 1)                       //coalescence event happend in the current iteration
		{
			ActiveNodes.push_back(lable);       //activate new node
			++lable;
			mMerger.push_back(std::vector<double>(4));
			mMerger.back()[0] = i;
			mMerger.back()[1] = x_merger;
			mMerger.back()[2] = time/mTimeScale;
			mMerger.back()[3] = N;
			i = int(ActiveNodes.size());        //update size of active nodes
			
		}
		temp_popsize *= std::exp(-Growthrate);//update population size
		N =  std::max((int)std::ceil(temp_popsize),(int)i);
	}
	
	mTotalBranchlength = 0;
	for (int j = 0; j < AllNodes.size(); ++j)
	{
		AllNodes[j].sortDescendent();
	}
	
	mTreeWalk(int(AllNodes.size()-1), AllNodes, 0, 0, gamma); // rho = 0, since growth has been accounted for already; theta = 0, since mutations are generated above
	mNewick += ";";
	mTMRC = time/mTimeScale;
    std::cout << " Moran T_MRACA = " << mTMRC << "\n" ;
    
}


//////////////
//  Newick  //
//////////////

myTree::myTree(std::string Newick, double Psi, double gamma, double mu, unsigned int N, long double rho)
{
	mTimeScale = 1.;
	std::vector<myNode> AllNodes ;
	int current = 0;
	int SampleSize = 0;
	ReadNewick(current, Newick, AllNodes, SampleSize);
	mSampleSize = SampleSize;
	mBranchlengthsCumulative.resize(SampleSize-1);
	mSFS.resize(SampleSize-1);
	mSequences.resize(mSampleSize);
	mTotalBranchlength = 0.;
	mNoSegSites = 0;
	for (int j = 0; j < AllNodes.size(); ++j)
	{
		AllNodes[j].sortDescendent();
		if (j < mSampleSize)
		{
			mSequences[j].push_back(0);
		}
	}
	
	double theta;
	if (gamma <= 1)
	{
	 theta = 0;
	}
	else if (gamma < 2)
	{
		theta = 2*mu*std::pow(N, gamma-1)/std::pow(Psi, 2.);
	}
	else if (gamma == 2)
	{
		theta = 2*mu*N*(1./(2+std::pow(Psi, 2)));
	}
	else
	{
		theta = mu*N;
	}

	mTreeWalk(0, AllNodes, rho, theta, gamma);
	mNewick += ";";
	if (rho > 0)
	{
		if (gamma < 2)
		{
			mTMRC = std::log(AllNodes[0].getTime()*rho*gamma+1.)/(rho*gamma); //transformation/ time scaling
		}
		else
		{
			mTMRC = std::log(AllNodes[0].getTime()*rho+1.)/(rho); //transformation/ time scaling
		}
	}
	else
	{
		mTMRC = AllNodes[0].getTime();
	}
    
}


////////////////////
// help functions //
////////////////////

void ReadNewick(int& current, std::string Newick, std::vector<myNode> &AllNodes, int& sampleSize)
{
	if (Newick[current] == 40 ) //opening bracket
	{
		AllNodes.push_back(myNode());
		AllNodes.back().setLable((int) AllNodes.size()-1);
		int i = (int) AllNodes.size()-1;
		++current;
		while (Newick[current] != 41)                              //check if closing bracket
		{
			if  ((Newick[current] >= 48 &&  Newick[current] <= 57) || (Newick[current] >= 65 &&  Newick[current] <= 90) || (Newick[current] >= 97 &&  Newick[current] <= 122)) //check if digit
			{
				sampleSize++;
				AllNodes.push_back(myNode(0, (int) AllNodes.size(),1));
				AllNodes.back().addDescendent((int) AllNodes.size()-1);
				AllNodes[i].addChild(AllNodes.back());
				while(Newick[current] != 58) //as long as next element is not a colon iterate
					++current;
				++current;
				std::string time_string;
				while(Newick[current] != 41 && Newick[current] != 44) //as long as next element is not a colon iterate
				{
					time_string += Newick[current];
					++current;
				}
				--current;
				double time = atof(time_string.c_str());
				AllNodes[i].setTime(time);
			}
			if (Newick[current] == 40) //opening bracket
			{
				int next_index = (int) AllNodes.size();
				ReadNewick(current, Newick, AllNodes, sampleSize);
				AllNodes[i].addChild(AllNodes[next_index]);
				while(Newick[current] != 58)    //as long as next element is not a colon iterate
					++current;
				++current;
				std::string time_string;
				while(Newick[current] != 41 && Newick[current] != 44) //as long as next element is not a colon iterate
				{
					time_string += Newick[current];
					++current;
				}
				--current;
				double time = atof(time_string.c_str())+AllNodes[next_index].getTime();
				AllNodes[i].setTime(time);
			}
			++current;
		}
	}
}


int removeRandomElement(std::vector<int>& v)
{
	int n = int(v.size());
	int i = randint(n);
	int x = v[i];
	v.erase(v.begin() +i);
	return x;
}


std::vector<long double> computeCoalescenceRates(double Psi, double gamma, int i)
{
	double Psi_sq = std::pow(Psi,2);
	std::vector<long double> lambda (i);
	if (gamma == 2)
	{
		
		lambda[1] = std::exp( std::log(n_choose_k(i, 2)) + std::log(2./(2.+Psi_sq) + Psi_sq/(2.+Psi_sq)*std::pow(1.-Psi,i-2)));
		lambda[0] = lambda[1];
		for (int x = 3; x <= i; ++x)
		{
			lambda[x-1] = std::exp( std::log(n_choose_k(i, x)) - std::log(2.0+Psi_sq) + x*std::log(Psi) + (i-x)*std::log(1.-Psi) );
			lambda[0] += lambda[x-1];
		}
	}
	if (gamma < 2)
	{
		lambda[1] = std::exp( std::log(n_choose_k(i, 2)) + (i-2.)*std::log(1.-Psi) );
		lambda[0] = lambda[1];
		for (int x = 3; x <= i; ++x)
		{
			lambda[x-1] = std::exp( std::log(n_choose_k(i, x)) + (x-2.)*std::log(Psi) + (i-x)*std::log(1.-Psi) );
			lambda[0] += lambda[x-1];
		}
	}
	else{
		lambda[1] = n_choose_k(i, 2);
		lambda[0] = lambda[1];
		for (int x = 3; x <= i; ++x)
		{
			lambda[x-1] = 0.;
		}
	}
	return lambda;
}


void myTree::mTreeWalk(int node_lable, std::vector<myNode>& AllNodes, long double rho, double theta, double gamma)
{
	if (AllNodes[node_lable].getChildren().size() == 0) // node is a leave
	{
		mNewick += "A_";
		mNewick+= convert_NumberToStringTree(AllNodes[node_lable].getLable());
	}
	else
	{   mNewick+= "(";
		int i = 0;
		while (i < AllNodes[node_lable].getChildren().size())
		{
			int IndexChild = AllNodes[node_lable].getChildren()[i]; // location of the i-th child in AllNodes
			mTreeWalk(IndexChild, AllNodes, rho, theta, gamma);
			mNewick+= ":";
			double branch_length;
			if (rho > 0)
			{
				if(gamma < 2)
				{
					branch_length= std::log(1+AllNodes[node_lable].getTime()*rho*gamma)/(rho*gamma) - std::log(1+AllNodes[(AllNodes[node_lable].getChildren()[i])].getTime()*rho*gamma)/(rho*gamma);
				}
				else
				{
					branch_length= std::log(1+AllNodes[node_lable].getTime()*rho)/(rho) - std::log(1+AllNodes[(AllNodes[node_lable].getChildren()[i])].getTime()*rho)/(rho);
				}
			}
			else
			{
				branch_length = AllNodes[node_lable].getTime()-AllNodes[(AllNodes[node_lable].getChildren()[i])].getTime();
			}
			
			double MutationRate = branch_length*theta/2.;
			int noMutations = randpoisson(MutationRate);
			mSFS[AllNodes[IndexChild].getNumberOfDescendents()-1] += noMutations;
			mNoSegSites += noMutations;
			if (noMutations > 0)
			{
				std::vector<bool> zeros(noMutations,0);
				std::vector<bool> ones (noMutations,1);
				
				int counterDescendents = 0;
				for (int k = 0; k < mSampleSize; ++k)
				{
					if (AllNodes[IndexChild].getDescendents()[counterDescendents] == k)
					{
						mSequences[k].insert(mSequences[k].end(),ones.begin(),ones.end());
						++counterDescendents;
					}
					else
					{
						mSequences[k].insert(mSequences[k].end(),zeros.begin(),zeros.end());
					}
				}
			}
			branch_length/=mTimeScale;
			mBranchlengths.push_back(branch_length);
			mTotalBranchlength += branch_length;
			mNumberOfDescendents.push_back(AllNodes[IndexChild].getNumberOfDescendents());
			mBranchlengthsCumulative[AllNodes[IndexChild].getNumberOfDescendents()-1] += branch_length;
			mNewick+= convert_NumberToStringTree(branch_length);
			if (i < AllNodes[node_lable].getChildren().size()-1)
			{
				mNewick+= ",";
			}
			++i;
		}
		mNewick+= ")";
	}
}


/////////////
// getter //
///////////


std::vector<double> myTree::getBranchlengths()
{
	return this->mBranchlengths;
}


std::vector<int> myTree::getNumberOfDescendents()
{
	return this->mNumberOfDescendents;
}

std::string myTree::getNewick()
{
	return this->mNewick;
}


int myTree::getSampleSize()
{
	return this->mSampleSize;
}

std::vector<int> myTree::getSFS()
{
	return this->mSFS;
}

std::vector<double> myTree::getBranchlengthsCumulative()
{
	return this->mBranchlengthsCumulative;
}

double myTree::getTMRCA()
{
	return this->mTMRC;
}

double myTree::getTotalBranchlength()
{
	return this->mTotalBranchlength;
}

int myTree::getNoSegSites()
{
	return this->mNoSegSites;
}

std::vector<std::vector<bool> > myTree::getSequence()
{
	return this->mSequences;
}

std::vector<std::vector<double> > myTree::getMerger()
{
	return this->mMerger;
}

long double n_choose_k( int N,  int k)
{
	if(N < k || N < 0 || k < 0)
		return 0;
	else
	{
		return boost::math::binomial_coefficient<long double>(N, k);
	}
}


// Alternative way to calculate the log-binomial (taken from Hybrid-Lambda)
long double logbinomial(long double N, long double k)
{
	assert(N >= 0.0L);
	assert(k >= 0.0L);
	assert(k <= N);
	
	long double x = 0.0L;
	long double bi = 0.0L;
	
	while (k > x)
	{
		bi += logl(N-x);
		bi -= logl(k-x);
		x++;
	}
	return bi;
}


//***convert_NumberToString*****************************
// Converts number to string
// Needed for setting output file name
//
template <typename _typeName>
std::string convert_NumberToStringTree ( _typeName Number )
{
	std::string _temp;
	std::ostringstream _convert;
	_convert << Number;
	_temp = _convert.str();
	
	return (_temp);
}

