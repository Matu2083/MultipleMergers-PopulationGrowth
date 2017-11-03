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



#ifndef __MMC_test__myIOParser__
#define __MMC_test__myIOParser__


#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <string>
#include <stdlib.h>
#include <sys/stat.h>

#include "myExpectedSFS.hpp"

//**********************
//** Output Functions **
//**********************
void initiateoutput();
void printPhi(std::vector<double> &Phi, double totalTreeLength);

template <typename _typeName>
std::string convert_NumberToString ( _typeName Number );
char * convert_StringToChar (std::string _string);
char * set_filename (char outfile_Name[]);

void readcommandLineArguments(int argc, const char *argv[]);
template < class T > T readNextInput(int argc, int& argc_i, const char *argv[]);
void readNextStringto(std::string &readto , int& argc_i, int argc_, char const *argv[]);
void printHelpMessage();
void checkArgvInput();
void checkdata();

//*********************
//** Output Switches **
//*********************

extern char outfile_name[255];
extern std::string OutputPath;
extern std::string outfile_prefix;

//**********************
//*** Output Streams ***
//**********************
extern std::ofstream _offPhi;


////////////////
// parameters //
////////////////

extern double minPsi;
extern double maxPsi;
extern unsigned int noStepsPsi;
extern double minRho;
extern double maxRho;
extern unsigned int noStepsRho;
extern int precisionPhi;
extern unsigned int sampleSize;
extern bool success;


#endif 
