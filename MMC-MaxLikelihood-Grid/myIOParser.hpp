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


#ifndef __MMC_test__myIOParser__
#define __MMC_test__myIOParser__

#include "myEstimationMethod.hpp"

//**********************
//** Output Functions **
//**********************
void initiateoutput();
void printparameters(std::ofstream * outfile); //prints parameters to outfile
void printMLEstimate(myEstimationMethod myEstimate, int dataSet);
void printL2AbsEstimate(myEstimationMethod myEstimate, int dataSet);
void printL1AbsEstimate(myEstimationMethod myEstimate, int dataSet);
void printMLLogEstimate(myEstimationMethod myEstimate, int dataSet);
void printGrid(myEstimationMethod myEstimate, int dataSet);

std::string get_date();
template <typename _typeName>
std::string convert_NumberToString ( _typeName Number );
char * convert_StringToChar (std::string _string);
char * set_filename (char outfile_Name[]);

void readcommandLineArguments(int argc, const char *argv[]);
template < class T > T readNextInput(int argc, int& argc_i, const char *argv[]);
void readNextStringto(std::string &readto , int& argc_i, int argc_, char const *argv[]);
void printHelpMessage();
void checkArgvInput();
void checkdataSFS();
void checkdataPhi();
void readSFSFile();
void readPhiFile();
void initializeValues();

//*********************
//** Output Switches **
//*********************

extern char outfile_name[255];
extern std::string outfile_prefix;
extern std::string dateconst;

//**********************
//*** Output Streams ***
//**********************
extern std::ofstream _offMLEstimate;
extern std::ofstream _offL2AbsEstimate;
extern std::ofstream _offL1AbsEstimate;
extern std::ofstream _offMLLogEstimate;
extern std::ofstream _offGrid;

////////////////
// parameters //
////////////////

extern int sampleSizeSFS, sampleSizePhi, noDataSets;
extern double minPsi, maxPsi, minRho, maxRho;
extern int noStepsPsi, noStepsRho;
extern int lumping;

extern std::vector<std::vector<double> > SFS;							// vector of SFS data sets
extern std::vector<std::vector<std::vector<double> > > Phi;			// vector of Phi data
extern std::vector<std::vector<double> > T_Total;					// vector of T_Total data
extern std::string SFSInputFile;
extern std::string PhiInputFile;
extern bool print_Grid;

#endif 
