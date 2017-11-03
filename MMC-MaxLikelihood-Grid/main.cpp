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

#include "myIOParser.hpp"
#include "myEstimationMethod.hpp"

void estimationRoutine();

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
	checkdataSFS();
	checkdataPhi();
	readSFSFile();
	readPhiFile();
	
	set_filename(outfile_name);
	initiateoutput();
	
	estimationRoutine();
	
	return (0);
}

      

void estimationRoutine()
{
		for (int i = 0; i < noDataSets; i++)
		{
			myEstimationMethod GRID(Phi, T_Total, SFS[i], minPsi, maxPsi, noStepsPsi, minRho, maxRho, noStepsRho, lumping);
			if(print_Grid == true)
			{
				printGrid(GRID, i);
			}
			printMLEstimate(GRID, i);
			printL2AbsEstimate(GRID, i);
			printL1AbsEstimate(GRID, i);
			printMLLogEstimate(GRID, i);
		}
	
	if(print_Grid == true)
	{
		_offGrid.close();
	}

	_offMLEstimate.close();
	_offL2AbsEstimate.close();
	_offL1AbsEstimate.close();
	_offMLLogEstimate.close();
}
