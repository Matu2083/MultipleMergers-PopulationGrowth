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


#include "myIOParser.hpp"

#define MAX_DATE 20

#ifndef VERSION
	#define VERSION "v1.0"
#endif


//*********************
//** Output Switches **
//*********************

char outfile_name[255];
std::string OutputPath = "./";
std::string outfile_prefix;

//**********************
//*** Output Streams ***
//**********************
std::ofstream _offPhi;

////////////////
// parameters //
////////////////

double minPsi = 0;
double maxPsi = 1;
unsigned int noStepsPsi = 100;

double minRho = 0;
double maxRho = 5;
unsigned int noStepsRho = 100;

unsigned int sampleSize = 100;
int precisionPhi = 1;
bool success = true;


void readcommandLineArguments(int argc, const char *argv[])
{
	int argc_i = 1;
	while (argc_i < argc)
	{
		std::string argv_i(argv[argc_i]);
		if ( argv_i == "-h" || argv_i == "-help" ){ printHelpMessage(); exit(0);}
		else if ( argv_i == "-MinPsi" ){ minPsi = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-MaxPsi" ){ maxPsi = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-NoStepsPsi" ){ noStepsPsi = readNextInput<unsigned int>(argc, argc_i, argv); }
		
		else if ( argv_i == "-MinRho" ){ minRho = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-MaxRho" ){ maxRho = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-NoStepsRho" ){ noStepsRho = readNextInput<unsigned int>(argc, argc_i, argv); }
		
		else if ( argv_i == "-sampleSize" ){ sampleSize = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-precPhi" ){ precisionPhi = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-OUT" ){ readNextStringto(OutputPath , argc_i, argc,  argv); }
		
		else{ throw (std::string("Unknown flag: ") + argv_i); }
		argc_i++;
	}
}

template < class T > T readNextInput(int argc, int& argc_i, const char *argv[])
{
	++argc_i;
	if (argc_i >= argc) throw (std::string( "Not enough parameters when parsing options: ") + argv[argc_i-1]);
	
	char c;
	T input;
	std::stringstream ss( argv[argc_i] );
	ss >> input;
	if (ss.fail() || ss.get(c)) throw ( std::string( "Failed to parse option: ") + argv[argc_i]);
	return input;
}

void readNextStringto(std::string &readto , int& argc_i, int argc, char const *argv[])
{
	argc_i++;
	if (argc_i >= argc) throw (std::string("Not enough parameters when parsing options: ") + argv[argc_i-1]);
	readto = std::string(argv[argc_i]);
	if ( readto[0] == '-' ) throw (std::string("Not enough parameters when parsing options: ") + argv[argc_i-1]);
}

void printHelpMessage()
{
	std::cout << std::endl << "MMC-Expected SFS " << VERSION << std::endl << std::endl;
	std::cout << "Usage:" << std::endl;
	std::cout << std::setw(20) << "-h/-help"         << "  --  " << "Help. List the following content." << std::endl;
	std::cout << std::setw(20) << "-MinPsi FLT"           << "  --  " << "Minimal multiple merger parameter psi. By default Psi = 0 is used corresponding to Kingman's Coalescent." << std::endl;
	std::cout << std::setw(20) << "-MaxPsi FLT"           << "  --  " << "Maximal multiple merger parameter psi. By default Psi = 1 is used." << std::endl;
	std::cout << std::setw(20) << "-NoStepsPsi INT"           << "  --  " << "Number of steps between minPsi and maxPsi. By default noStepsPsi = 20 is used." << std::endl;
	std::cout << std::setw(20) << "-MinRho FLT"           << "  --  " << "Minimal exponential growth rate rho. By default Rho = 0 is used corresponding to a constant population." << std::endl;
	std::cout << std::setw(20) << "-MaxRho FLT"           << "  --  " << "Maximal exponential growth rate rho. By default Rho = 2 is used corresponding to a constant population." << std::endl;
	std::cout << std::setw(20) << "-NoStepsRho INT"           << "  --  " << "Number of steps between minRho and maxRho. By default noStepsRho = 20 is used." << std::endl;
	std::cout << std::setw(20) << "-precPhi INT"           << "  --  " << "Number of bits used for the multiple precision calculation. By default 100 is used." << std::endl;
	std::cout << std::setw(20) << "-OUT STR"                << "  --  " << "Specify the FILE PATH where output is written."  << std::endl;

	std::cout << std::endl;
	
	std::cout << "Examples:" << std::endl	<< std::endl;
	std::cout << "MMC-Phi-MaxLikelihood-LookUp.out -OUT /path/To/PhiData " << std::endl;
	std::cout << "MMC-Phi-MaxLikelihood-LookUp.out " << std::endl;
	std::cout << std::endl;
	
}


//******************************************************
//**  INITIALIZING FUNCTIONS  **************************
//******************************************************


void checkArgvInput()
{

	if(minPsi < 0 || minPsi > 1)
	{
		std::cout << "Error: minPsi = " << minPsi << ", but needs to be ]0,1[. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(maxPsi < 0 || maxPsi > 1)
	{
		std::cout << "Error: maxPsi = " << maxPsi << ", but needs to be ]0,1[. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(minPsi > maxPsi)
	{
		std::cout << "Error: minPsi > maxPsi. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(minRho < 0)
	{
		std::cout << "Error: minRho = " << minRho << ", but needs to be positive. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(minRho > maxRho)
	{
		std::cout << "Error: minRho > maxRho. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if (sampleSize <= 0)
	{
		std::cout << "Error: sampleSize = " << sampleSize << " but needs to be positive.\n";
		exit(1);
	}
	
	if(precisionPhi <= 0)
	{
		std::cout << "Error: precisionPhi = " << precisionPhi << ", but needs to be positive. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	checkdata();
}


//***checkdata************************************
//  Pre-scans the data to infer the number of
//	data sets and the number of samples
//
void checkdata()
{
	if(precisionPhi == 1)
	{
		if (sampleSize <= 40)
		{
			precisionPhi = 100;
		}
		else if (sampleSize <= 110)
		{
			precisionPhi = 500;
		}
		else if (sampleSize <= 500)
		{
			precisionPhi = 2000;
		}
		else
		{
			precisionPhi = 5000;
		}
		mpreal::set_default_prec(mpfr::digits2bits(precisionPhi));
	}
	else
	{
		mpreal::set_default_prec(mpfr::digits2bits(precisionPhi));
		std::cout << "WARNING: You have manually set precision to " << precisionPhi << ". The chosen precision might be too low for the given sample size.\n";
	}
	
}


			//////////////////////////////////////////////////
		   /////////            OUTPUT            ///////////
		  //////////////////////////////////////////////////


//**initiateoutput*******************************************
//	Initates all output files
//
void initiateoutput()
{
	
	char _PhiFile[255];
	strcpy(_PhiFile, outfile_name);
	strcat (_PhiFile, "_expectedPhi.txt");
	_offPhi.open(_PhiFile);

}


void printPhi(std::vector<double> &Phi, double totalTreeLength)
{
	_offPhi.precision(20);
	_offPhi << totalTreeLength;
	for (int i = 0; i < sampleSize - 1; i++)
	{
		_offPhi << "\t" << Phi[i];
	}
	_offPhi << std::endl;
}

//***convert_NumberToString*****************************
// Converts number to string
// Needed for setting output file name
//
template <typename _typeName>
std::string convert_NumberToString ( _typeName Number )
{
	std::string _temp;
	std::ostringstream _convert;
	_convert << Number;
	_temp = _convert.str();
	
	return (_temp);
}

//***convert_StringToChar*****************************
// Converts string to char
// Needed for setting output file name
//
char * convert_StringToChar (std::string _string)
{
	return (const_cast<char*> (_string.c_str()));
}

char * set_filename (char outfile_Name[])
{
	outfile_prefix = OutputPath;

	strcpy (outfile_Name, convert_StringToChar(outfile_prefix));
	
	strcat (outfile_Name, "MinPsi-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(minPsi)));

	strcat (outfile_Name, "-MaxPsi-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(maxPsi)));

	strcat (outfile_Name, "-NoStepPsi-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(noStepsPsi)));

	strcat (outfile_Name, "-MinRho-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(minRho)));
	
	strcat (outfile_Name, "-MaxRho-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(maxRho)));
	
	strcat (outfile_Name, "-NoStepRho-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(noStepsRho)));

	strcat (outfile_Name, "-noSamples-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(sampleSize)));
	
	
	return (outfile_Name);
}
