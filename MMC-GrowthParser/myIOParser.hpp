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


#ifndef __MMC_test__myIOParser__
#define __MMC_test__myIOParser__


#include <iostream>
#include <iomanip>    
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <stdlib.h>    
#include <sys/stat.h>

#include "myRandomNumbers.hpp"
#include "myTree.hpp"
#include "mySummaryStatistics.hpp"


//**********************
//** Output Functions **
//**********************
void initiateoutput(int whichSim);
void initiateoutputSequences(int whichSim, int rep);
void printparameters(std::ofstream * outfile, int whichSim); //prints parameters to outfile
void printoutput(myTree _tree, int rep);

void printSFS(myTree _tree);																	//prints Site Frequency Spectrum
void printCumulativeBranchLength(myTree _tree);
void printRelativeBranchLength(std::vector<double> cumulBranchLengthsAll, double t_totalAll, int noSamples);
void printNewickTree(myTree _tree);
void printSequences(myTree _tree, int rep);
void printRates(myTree _tree);
void printSumStats(myTree _tree, int rep);
void closeOffstreams();

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
void readNewick();


extern bool sim_coal;
extern bool sim_moran;
extern bool sim_newick;

//*********************
//** Output Switches **
//*********************
extern bool print_SFS;
extern bool noInf;
extern bool print_cumulBranchLength;
extern bool print_relBranchLength;
extern bool print_NewickTree;
extern bool print_Sequences;
extern bool print_Rates;
extern bool print_sumStats;

extern char outfile_name[255];
extern std::string dateconst;
extern std::string outfile_prefix;

//**********************
//*** Output Streams ***
//**********************
extern std::ofstream _offSFS;					// Site frequency spectrum
extern std::ofstream _offCumulBranchLength;
extern std::ofstream _offRelBranchLength;
extern std::ofstream _offNewickTree;
extern std::ofstream _offSequence;
extern std::ofstream _offRates;
extern std::ofstream _offSumStats;

////////////////
// parameters //
////////////////

extern unsigned int initialPopSize, noReplicates;
extern int noSamples;
extern double psi, mu, Gamma;
extern long double rho;
extern std::string newickFileInput;


#endif /* defined(__MMC_test__myIOParser__) */
