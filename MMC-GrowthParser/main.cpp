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


#include <iostream>
#include <stdexcept>
#include <iomanip>    
#include <fstream>
#include "myNode.hpp"
#include "myTree.hpp"
#include "myRandomNumbers.hpp"
#include "myIOParser.hpp"

void simulateCoalescent();
void simulateMoran();

int main(int argc, const char * argv[])
{

	///////////////////
	// parsing input //
	///////////////////
	

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
	
	initializeRNG(seed);
	set_filename(outfile_name);
	
// Initialize random number generator
	
	if (sim_coal)
	{
		if (Gamma > 2)
		{
			initiateoutput(0);
		}
		else
		{
			initiateoutput(1);
		}
		simulateCoalescent();
        closeOffstreams();
	}
	if (sim_moran)
	{
		initiateoutput(2);
		simulateMoran();
        closeOffstreams();

	}
	if (newickFileInput != "-1")
	{
		readNewick();
        closeOffstreams();
	}
	
 return 0;
 }




void simulateCoalescent()
{
	double t_totalAll = 0;
	std::vector<double> cumulBranchLengthsAll;
	cumulBranchLengthsAll.resize(noSamples-1);
	for (int i = 0; i < noSamples-1; ++i)
	{
		cumulBranchLengthsAll[i] = 0.;
	}
	
	for (int i = 0; i < noReplicates; ++i)
	{
		myTree _tree(noSamples, psi, Gamma, mu, rho, initialPopSize);
		
		if (noInf == false)
		{
			if (print_Sequences)
			{
				initiateoutputSequences(1, i+1);
			}
			
			t_totalAll += _tree.getTotalBranchlength();
			for (int j = 0; j < noSamples-1; ++j)
			{
				cumulBranchLengthsAll[j] += _tree.getBranchlengthsCumulative()[j];
			}
			
			printoutput(_tree, i+1);
		}
		else
		{
			if (_tree.getTotalBranchlength() != _tree.getBranchlengthsCumulative()[0])
			{
				if (print_Sequences)
				{
					initiateoutputSequences(1, i+1);
				}
				
				t_totalAll += _tree.getTotalBranchlength();
				for (int j = 0; j < noSamples-1; ++j)
				{
					cumulBranchLengthsAll[j] += _tree.getBranchlengthsCumulative()[j];
				}
				
				printoutput(_tree, i+1);
			}
			else // redo and sample another tree
			{
				i--;
			}
		}
	}
	if(print_relBranchLength) printRelativeBranchLength(cumulBranchLengthsAll, t_totalAll, noSamples);
}


void simulateMoran()
{
	double t_totalAll = 0;
	std::vector<double> cumulBranchLengthsAll;
	cumulBranchLengthsAll.resize(noSamples-1);
	for (int i = 0; i < noSamples-1; ++i)
	{
		cumulBranchLengthsAll[i] = 0.;
	}
	
	for (int i = 0; i < noReplicates; ++i)
	{
		myTree _tree(noSamples, psi, Gamma, mu, initialPopSize, rho);
		if (print_Sequences)
		{
			initiateoutputSequences(2, i+1);
		}
		
		t_totalAll += _tree.getTotalBranchlength();
		for (int j = 0; j < noSamples-1; ++j)
		{
			cumulBranchLengthsAll[j] += _tree.getBranchlengthsCumulative()[j];
		}
		
		printoutput(_tree, i+1);
	}
	if(print_relBranchLength) printRelativeBranchLength(cumulBranchLengthsAll, t_totalAll, noSamples);
}

