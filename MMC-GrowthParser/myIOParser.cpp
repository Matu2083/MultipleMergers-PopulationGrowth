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

#include "myIOParser.hpp"

#define MAX_DATE 20

#ifndef VERSION
	#define VERSION "v1.0"
#endif

bool sim_coal = true;
bool sim_moran = false;
bool sim_newick = false;

//*********************
//** Output Switches **
//*********************
bool print_SFS = false;
bool noInf = false;
bool print_cumulBranchLength = false;
bool print_relBranchLength = false;
bool print_NewickTree = false;
bool print_Sequences = false;
bool print_Rates = false;
bool print_sumStats = false;

char outfile_name[255];
std::string dateconst;
std::string outfile_prefix = "./";

//**********************
//*** Output Streams ***
//**********************
std::ofstream _offSFS;					// Site frequency spectrum
std::ofstream _offCumulBranchLength;
std::ofstream _offRelBranchLength;
std::ofstream _offNewickTree;
std::ofstream _offSequence;
std::ofstream _offRates;
std::ofstream _offSumStats;

////////////////
// parameters //
////////////////

unsigned int initialPopSize = 10000, noReplicates = 1;
int noSamples = 0;
double psi = 0, mu = 0.00025, Gamma = 3;
long double rho = 0;
std::string newickFileInput = "-1";


void readcommandLineArguments(int argc, const char *argv[])
{
	if (argc == 1) // if input is incomplete stop
	{
		std::cout << "Error: Insufficient number of arguments passed.\n";
		printHelpMessage();
		exit(1);
	}
	
	int argc_i = 1;
	while (argc_i < argc)
	{
		std::string argv_i(argv[argc_i]);
		if ( argv_i == "-h" || argv_i == "-help" ){ printHelpMessage(); exit(0);}
		else if ( argv_i == "-s" || argv_i == "-sampleSize"  ){ noSamples = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-psi" ){ psi = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-mu" ){ mu = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-g" || argv_i == "-gamma"  ){ Gamma = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-num" ){ noReplicates = readNextInput<unsigned int>(argc, argc_i, argv); }
		else if ( argv_i == "-N" || argv_i == "-popSize" ){ initialPopSize = readNextInput<unsigned int>(argc, argc_i, argv); }
		else if ( argv_i == "-r" || argv_i == "-rho" ){ rho = readNextInput<long double>(argc, argc_i, argv); }
		else if ( argv_i == "-seed"){ seed = readNextInput<unsigned long long int>(argc, argc_i, argv); }
		else if ( argv_i == "-SFS"){ print_SFS = true; }
		else if ( argv_i == "-inference"){ noInf = true; }
		else if ( argv_i == "-sumStats"){ print_sumStats = true; }
		else if ( argv_i == "-relBrL"){ print_relBranchLength = true; }
		else if ( argv_i == "-cumulBrL"){ print_cumulBranchLength = true; }
		else if ( argv_i == "-tree"){ print_NewickTree = true; }
		else if ( argv_i == "-seq"){ print_Sequences = true; }
		else if ( argv_i == "-rates"){ print_Rates = true; }
		else if ( argv_i == "-moran"){ sim_coal = false; sim_moran = true; }
		else if ( argv_i == "-newick"){ sim_coal = false; sim_newick = true; readNextStringto(newickFileInput , argc_i, argc,  argv); }
		else if ( argv_i == "-out"){ readNextStringto(outfile_prefix , argc_i, argc,  argv); }
		else{ throw (std::string("Unknown flag: ") + argv_i); }
		argc_i++;
	}
	
	dateconst = get_date();
}

template < class T > T readNextInput(int argc, int &argc_i, const char *argv[])
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
	std::cout << std::endl << "MMC-Growth " << VERSION << std::endl << std::endl;
	std::cout << "Usage:" << std::endl;
	std::cout << std::setw(20) << "-h/-help"         << "  --  " << "Help. List the following content." << std::endl;
	std::cout << std::setw(20) << "-s/-sampleSize INT"           << "  --  " << "Specify the number of samples." << std::endl;
	std::cout << std::setw(20) << "-psi FLT"           << "  --  " << "Multiple merger parameter psi defined by a single numerical constant." << std::endl;
	std::cout << std::setw(20) << "-mu FLT"             << "  --  " << "User defined constant mutation rate per locus. By default mutation rate 0.00025 is used." << std::endl;
	std::cout << std::setw(26) << " "                               << "Assumes infinite sites model." << std::endl;
	std::cout << std::setw(20) << "-g/-gamma FLT"         << "  --  " << "Time-scaling parameter. By default gamma 3 is used corresponding to a Kingman coalescent." << std::endl;
	std::cout << std::setw(26) << " "                               << "Scales the probability of sweepstake reproductive events." << std::endl;
	std::cout << std::setw(20) << "-num INT"            << "  --  " << "The number of gene trees to be simulated." << std::endl;
	std::cout << std::setw(20) << "-N/popSize INT"        << "  --  " << "Population sizes are defined by a single numerical constant." << std::endl;
	std::cout << std::setw(26) << " "                               << "By default, population size 10,000 is used." << std::endl;
	std::cout << std::setw(20) << "-r/-rho FLT"           << "  --  " << "User defined constant (exponential) growth rate. By default rho 0 is used." << std::endl;
	std::cout << std::setw(20) << "-seed INT"           << "  --  " << "User defined random number SEED." << std::endl;
	std::cout << std::setw(20) << "-SFS"                << "  --  " << "Generate site frequency data." << std::endl;
	std::cout << std::setw(20) << "-inference"                << "  --  " << "Generate site frequency data only if SFS[1] != segSites." << std::endl;
	std::cout << std::setw(20) << "-sumStats"                << "  --  " << "Generate summary statistics data." << std::endl;
	std::cout << std::setw(20) << "-cumulBrL"                << "  --  " << "Generate cumulative branch length data." << std::endl;
	std::cout << std::setw(20) << "-relBrL"                << "  --  " << "Generate relative branch length data." << std::endl;
	std::cout << std::setw(20) << "-tree"                << "  --  " << "Generate simulated gene tree data." << std::endl;
	std::cout << std::setw(20) << "-seq"                << "  --  " << "Generate sequence data." << std::endl;
	std::cout << std::setw(20) << "-rates"                << "  --  " << "Generate coalescence rate data." << std::endl;
	std::cout << std::setw(20) << "-moran"                << "  --  " << "Simulate extended Moran model." << std::endl;
	std::cout << std::setw(26) << " "                               << "Note that these are computationally intensive. " << std::endl;
	std::cout << std::setw(20) << "-newick STR"                << "  --  " << "Specify the FILE NAME of trees to generate data." << std::endl;
	std::cout << std::setw(20) << "-out STR"                << "  --  " << "Specify the FILE PATH where generated data is to be put." << std::endl;
	std::cout << std::endl;

	std::cout << "Examples:" << std::endl	<< std::endl;
	std::cout << "MMC-Growth.out -g 1.5 -s 10 -N 50000 -mu 0.0001 -r 0.2 -psi 0.2 -SFS -seq" << std::endl;
	std::cout << "MMC-Growth.out -g 1.5 -s 10 -N 50000 -mu 0.0001 -r 0.2 -psi 0.2 -SFS -seq -moran" << std::endl;
	std::cout << std::endl;
	
}

//******************************************************
//**  INITIALIZING FUNCTIONS  **************************
//******************************************************


void checkArgvInput()
{
	
	if(noSamples <= 1)
	{
		std::cout << "Error: noSamples = " << noSamples << ", but needs to be larger than 1. See ReadMe file for input arguments.\n";
		printHelpMessage();
		exit(1);
	}
	
	if(psi < 0 || psi > 1)
	{
		std::cout << "Error: psi = " << psi << ", but needs to be [0,1]. See ReadMe file for input arguments.\n";
		printHelpMessage();
		exit(1);
	}
	
	if (psi == 0)
	{
		Gamma = 2.5;
	}
	
	if(mu <= 0)
	{
		std::cout << "Error: mu = " << mu << ", but needs to be positive. See ReadMe file for input arguments.\n";
		printHelpMessage();
		exit(1);
	}
	
	if(Gamma <= 0)
	{
		std::cout << "Error: gamma = " << Gamma << ", but needs to be positive. See ReadMe file for input arguments.\n";
		printHelpMessage();
		exit(1);
	}
	
	if(noReplicates <= 0)
	{
		std::cout << "Error: noReplicates = " << noReplicates << ", but needs to be positive. See ReadMe file for input arguments.\n";
		printHelpMessage();
		exit(1);
	}
	
	if(initialPopSize <= 0)
	{
		std::cout << "Error: initialPopSize = " << initialPopSize << ", but needs to be positive. See ReadMe file for input arguments.\n";
		printHelpMessage();
		exit(1);
	}
	
	if(rho < 0)
	{
		std::cout << "Error: rho = " << rho << ", but needs to non-negative. See ReadMe file for input arguments.\n";
		printHelpMessage();
		exit(1);
	}
	
	if(sim_newick == 1 && sim_moran == 1)
	{
		std::cout << "Error: Both Moran model and newick-based simulations requested. Program handles only one simulation request at a time. See ReadMe file for input arguments.\n";
		printHelpMessage();
		exit(1);
	}
	
}


void readNewick()
{
	std::ifstream _infile(newickFileInput.c_str());
	std::string _line;
	initiateoutput(3);
	int i = 0;
	while (std::getline(_infile, _line))
	{
		i++;
		std::istringstream _iss(_line);
		std::string readNewick;

		if (!(_iss >> readNewick))
		{
			break;
		} // error
		readNewick.erase(readNewick.size() - 1);
				
		myTree passedTree(readNewick, psi, Gamma, mu, initialPopSize, rho);
		if (print_Sequences)
		{
			initiateoutputSequences(3, i);
		}
		printoutput(passedTree, i);
	}
	
}


			//////////////////////////////////////////////////
		   /////////            OUTPUT            ///////////
		  //////////////////////////////////////////////////


//**initiateoutput*******************************************
//	Initates all output files
//
void initiateoutput(int whichSim)
{
	if (print_SFS)
	{
		char _SFS_file[255];
		strcpy(_SFS_file, outfile_name);
		if (whichSim == 1)
		{
			strcat (_SFS_file, "Coalescent_");
		}
		else if (whichSim == 2)
		{
			strcat (_SFS_file, "Moran_");
		}
		else if (whichSim == 3)
		{
			strcat (_SFS_file, "Newick_");
		}
		else
		{
			strcat (_SFS_file, "Coalescent_");
		}
		
		strcat (_SFS_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (_SFS_file, "_SFS.txt");
		_offSFS.open(_SFS_file);
		printparameters(& _offSFS, whichSim);
		_offSFS << "#t_total";
		_offSFS << "\tt_MRCA";
		_offSFS << "\tS_total";
		for(int i = 1; i < noSamples; i++)
		{
			_offSFS << "\t" << i;
		}
		_offSFS << std::endl;
	}
	
	if (print_cumulBranchLength)
	{
		char _cumulBL_file[255];
		strcpy(_cumulBL_file, outfile_name);
		if (whichSim == 1)
		{
			strcat (_cumulBL_file, "Coalescent_");
		}
		else if (whichSim == 2)
		{
			strcat (_cumulBL_file, "Moran_");
		}
		else if (whichSim == 3)
		{
			strcat (_cumulBL_file, "Newick_");
		}
		else
		{
			strcat (_cumulBL_file, "Coalescent_");
		}
		strcat (_cumulBL_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (_cumulBL_file, "_cumulativeBranchLength.txt");
		_offCumulBranchLength.open(_cumulBL_file);
		printparameters(& _offCumulBranchLength, whichSim);
		
		_offCumulBranchLength << "#t_total";
		for(int i = 1; i < noSamples; i++)
		{
			_offCumulBranchLength << "\t" << i;
		}
		_offCumulBranchLength << std::endl;
	}
	
	if (print_relBranchLength)
	{
		char _relBL_file[255];
		strcpy(_relBL_file, outfile_name);
		if (whichSim == 1)
		{
			strcat (_relBL_file, "Coalescent_");
		}
		else if (whichSim == 2)
		{
			strcat (_relBL_file, "Moran_");
		}
		else if (whichSim == 3)
		{
			strcat (_relBL_file, "Newick_");
		}
		else
		{
			strcat (_relBL_file, "Coalescent_");
		}
		strcat (_relBL_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (_relBL_file, "_relBranchLength.txt");
		_offRelBranchLength.open(_relBL_file);
		printparameters(& _offRelBranchLength, whichSim);
		
		_offRelBranchLength << "#t_total";
		for(int i = 1; i < noSamples; i++)
		{
			_offRelBranchLength << "\t" << i;
		}
		_offRelBranchLength << std::endl;
	}
	
	if (print_NewickTree)
	{
		char _newick_file[255];
		strcpy(_newick_file, outfile_name);
		if (whichSim == 1)
		{
			strcat (_newick_file, "Coalescent_");
		}
		else if (whichSim == 2)
		{
			strcat (_newick_file, "Moran_");
		}
		else if (whichSim == 3)
		{
			strcat (_newick_file, "Newick_");
		}
		else
		{
			strcat (_newick_file, "Coalescent_");
		}
		strcat (_newick_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (_newick_file, "_NewickTrees.tree");
		_offNewickTree.open(_newick_file);
		printparameters(& _offNewickTree, whichSim);
		
	}
	
	if (print_Rates)
	{
		char _Rates_file[255];
		strcpy(_Rates_file, outfile_name);
		if (whichSim == 1)
		{
			strcat (_Rates_file, "Coalescent_");
		}
		else if (whichSim == 2)
		{
			strcat (_Rates_file, "Moran_");
		}
		else if (whichSim == 3)
		{
			strcat (_Rates_file, "Newick_");
		}
		else
		{
			strcat (_Rates_file, "Coalescent_");
		}
		strcat (_Rates_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (_Rates_file, "_Rates.txt");
		_offRates.open(_Rates_file);
		printparameters(& _offRates, whichSim);
		
		_offRates << "#ActiveLineages" << "\tx-merger" << "\ttime" << "\tPopSize";
		_offRates << std::endl;
	}
	
	if (print_sumStats)
	{
		char _sumStats_file[255];
		strcpy(_sumStats_file, outfile_name);
		if (whichSim == 1)
		{
			strcat (_sumStats_file, "Coalescent_");
		}
		else if (whichSim == 2)
		{
			strcat (_sumStats_file, "Moran_");
		}
		else if (whichSim == 3)
		{
			strcat (_sumStats_file, "Newick_");
		}
		else
		{
			strcat (_sumStats_file, "Coalescent_");
		}
		
		strcat (_sumStats_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (_sumStats_file, "_sumStats.txt");
		_offSumStats.open(_sumStats_file);
		printparameters(& _offSumStats, whichSim);
		
		_offSumStats << "#rep";
		_offSumStats << "\tTheta_W";
		_offSumStats << "\tTheta_T";
		_offSumStats << "\tTheta_H";
		_offSumStats << "\tFayWu_H";
		_offSumStats << "\tTajimasD";
		_offSumStats << "\tFuLiD";
		
		_offSumStats << std::endl;
	}
}

void initiateoutputSequences(int whichSim, int rep)
{
	if (print_Sequences)
	{
		char _Seq_file[255];
		if(rep == 1)
		{
			if (whichSim == 1)
			{
				std::string rm_OldFolder = "rm -rf " + (std::string) outfile_prefix + "SequenceFilesCoalescent";
				int sys = system( rm_OldFolder.c_str() );
				std::string mkdir_NewFolder="mkdir " + (std::string) outfile_prefix + "SequenceFilesCoalescent";
				sys=system( mkdir_NewFolder.c_str() );
			}
			else if (whichSim == 2)
			{
				std::string rm_OldFolder = "rm -rf " + (std::string) outfile_prefix + "SequenceFilesMoran";
				int sys = system( rm_OldFolder.c_str() );
				std::string mkdir_NewFolder="mkdir " + (std::string) outfile_prefix + "SequenceFilesMoran";
				sys=system( mkdir_NewFolder.c_str() );
			}
			else if (whichSim == 3)
			{
				std::string rm_OldFolder = "rm -rf " + (std::string) outfile_prefix + "SequenceFilesNewick";
				int sys = system( rm_OldFolder.c_str() );
				std::string mkdir_NewFolder="mkdir " + (std::string) outfile_prefix + "SequenceFilesNewick";
				sys=system( mkdir_NewFolder.c_str() );
			}
			else
			{
				std::string rm_OldFolder = "rm -rf " + (std::string) outfile_prefix + "SequenceFilesCoalescent";
				int sys = system( rm_OldFolder.c_str() );
				std::string mkdir_NewFolder="mkdir " + (std::string) outfile_prefix + "SequenceFilesCoalescent";
				sys=system( mkdir_NewFolder.c_str() );
			}

		}
		
		strcpy(_Seq_file, convert_StringToChar(outfile_prefix));
		if (whichSim == 1)
		{
			strcat(_Seq_file, "SequenceFilesCoalescent/");
		}
		else if (whichSim == 2)
		{
			strcat(_Seq_file, "SequenceFilesMoran/");
		}
		else if (whichSim == 3)
		{
			strcat(_Seq_file, "SequenceFilesNewick/");
		}
		else
		{
			strcat(_Seq_file, "SequenceFilesCoalescent/");
		}

		
		strcat (_Seq_file, "Psi-");
		strcat (_Seq_file, convert_StringToChar (convert_NumberToString(psi)));
		
		strcat (_Seq_file, "-N-");
		strcat (_Seq_file, convert_StringToChar (convert_NumberToString(initialPopSize)));
		
		strcat (_Seq_file, "-rho-");
		strcat (_Seq_file, convert_StringToChar (convert_NumberToString(rho)));
		
		strcat (_Seq_file, "-MutRate-");
		strcat (_Seq_file, convert_StringToChar (convert_NumberToString(mu)));
		
		strcat (_Seq_file, "-gamma-");
		strcat (_Seq_file, convert_StringToChar (convert_NumberToString(Gamma)));
		
		strcat (_Seq_file, "-noSamples-");
		strcat (_Seq_file, convert_StringToChar (convert_NumberToString(noSamples)));
		strcat (_Seq_file, "-");

		if (whichSim == 1)
		{
			strcat (_Seq_file, "Coalescent_");
		}
		else if (whichSim == 2)
		{
			strcat (_Seq_file, "Moran_");
		}
		else if (whichSim == 3)
		{
			strcat (_Seq_file, "Newick_");
		}
		else
		{
			strcat (_Seq_file, "Coalescent_");
		}
		strcat (_Seq_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (_Seq_file, "_SeqFile_");
		strcat (_Seq_file, convert_StringToChar(convert_NumberToString(rep)));
		strcat (_Seq_file, ".txt");
		_offSequence.open(_Seq_file);
		printparameters(& _offSequence, whichSim);
		_offSequence << "#ID";
		_offSequence << "\t Sequence";
		_offSequence << std::endl;
	}
}

//***printparameters************************************
//  Prints parameters to output file
//
void printparameters(std::ofstream * outfile, int whichSim)
{
	if (whichSim == 0)
	{
		*outfile << "# Kingman Coalescent with exponential growth";
		*outfile << "\n# Psi: " << psi;
		*outfile << "\n# Population Size: " << initialPopSize;
		*outfile << "\n# Growth rate: " << rho;
		*outfile << "\n# Gamma: " << Gamma;
		*outfile << "\n# Mutation rate: " << mu;
		*outfile << "\n# noSamples: " << noSamples;
		*outfile << "\n# Random seed: " << seed;
		*outfile << std::endl << std::endl;
	}
	else
	{
		*outfile << "# Psi-Coalescent with exponential growth";
		*outfile << "\n# Psi: " << psi;
		*outfile << "\n# Population Size: " << initialPopSize;
		*outfile << "\n# Growth rate: " << rho;
		*outfile << "\n# Gamma: " << Gamma;
		*outfile << "\n# Mutation rate: " << mu;
		*outfile << "\n# noSamples: " << noSamples;
		*outfile << "\n# Random seed: " << seed;
		*outfile << std::endl << std::endl;
	}

}

void printoutput(myTree _tree, int rep)
{
	if (print_SFS) printSFS(_tree);
	if (print_cumulBranchLength) printCumulativeBranchLength(_tree);
	if (print_NewickTree) printNewickTree(_tree);
	if (print_Sequences) printSequences(_tree, rep);
	if (print_Rates) printRates(_tree);
	if (print_sumStats) printSumStats(_tree, rep);
}

void printSFS(myTree _tree)
{
	_offSFS << _tree.getTotalBranchlength();
	_offSFS << "\t" << _tree.getTMRCA();
	_offSFS << "\t" << _tree.getNoSegSites();
	for (int i = 0; i < _tree.getSampleSize()-1; i++)
	{
		_offSFS << "\t" << _tree.getSFS()[i];
	}
	_offSFS << "\n";
	
}

void printCumulativeBranchLength(myTree _tree)
{
	_offCumulBranchLength << _tree.getTotalBranchlength();
	for (int i = 0; i < _tree.getSampleSize()-1; i++)
	{
		_offCumulBranchLength << "\t" << _tree.getBranchlengthsCumulative()[i];
	}
	_offCumulBranchLength << "\n";
}


void printRelativeBranchLength(std::vector<double> cumulBranchLengthsAll, double t_totalAll, int sampleSize)
{
	_offRelBranchLength.precision(20);
	_offRelBranchLength << t_totalAll/ (double) noReplicates;
	for (int i = 0; i < sampleSize-1; i++)
	{
		_offRelBranchLength << "\t" << cumulBranchLengthsAll[i] / t_totalAll;
	}
	_offRelBranchLength << "\n";
}


void printSequences(myTree _tree, int rep)
{
	for(int _sample = 0; _sample < _tree.getSequence().size(); _sample++)
	{
		_offSequence << "A_" << _sample + 1 << "\t";
		for (int _pos = 1; _pos < _tree.getSequence()[0].size(); _pos++)
		{
			_offSequence << _tree.getSequence()[_sample][_pos];
		}
		_offSequence << "\n";
	}
	_offSequence << "\n";
	_offSequence.close();
}


void printNewickTree(myTree _tree)
{
	_offNewickTree << _tree.getNewick();
	_offNewickTree << "\n";
}


void printRates(myTree _tree)
{
	for (int i = 0; i < _tree.getMerger().size(); i++)
	{
		_offRates << _tree.getMerger()[i][0] << "\t" << _tree.getMerger()[i][1] << "\t" << _tree.getMerger()[i][2] << "\t" <<  _tree.getMerger()[i][3] << "\n";
	}
}

void printSumStats(myTree _tree, int rep)
{
	mySummaryStatistics _sumStats(_tree);
	
	_offSumStats << rep << "\t";
	_offSumStats << _sumStats.getTheta_W() << "\t";
	_offSumStats << _sumStats.getThetaT() << "\t";
	_offSumStats << _sumStats.getThetaH() << "\t";
	_offSumStats << _sumStats.getFayAndWusH() << "\t";
	_offSumStats << _sumStats.getTajimasD() << "\t";
	_offSumStats << _sumStats.getFuAndLisD();
	
	_offSumStats << std::endl;
}


//***get_date********************************************
// Gets current date and time for naming the outputfile uniquely
// Format is dd_mm_yyyy_hh_minmin_secsec
//
std::string get_date()
{
	time_t now;
	char the_date[MAX_DATE];
	
	now = time(NULL);
	
	if (now != -1)
	{
		strftime(the_date, MAX_DATE, "%d_%m_%Y_%H%M%S", gmtime(&now));
	}
	
	return std::string (the_date);
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
	strcpy (outfile_Name, convert_StringToChar(outfile_prefix));
	
	strcat (outfile_Name, "Psi-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(psi)));
	
	strcat (outfile_Name, "-N-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(initialPopSize)));
	
	strcat (outfile_Name, "-rho-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(rho)));
	
	strcat (outfile_Name, "-MutRate-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(mu)));
	
	strcat (outfile_Name, "-gamma-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(Gamma)));
	
	strcat (outfile_Name, "-noSamples-");
	strcat (outfile_Name, convert_StringToChar (convert_NumberToString(noSamples)));
	strcat (outfile_Name, "-");
	
	return (outfile_Name);
}

void closeOffstreams()
{
	if (print_SFS) _offSFS.close();
	if (print_cumulBranchLength) _offCumulBranchLength.close();
	if (print_relBranchLength) _offRelBranchLength.close();
	if (print_NewickTree) _offNewickTree.close();
	if (print_Rates) _offRates.close();
	if (print_sumStats) _offSumStats.close();
}
