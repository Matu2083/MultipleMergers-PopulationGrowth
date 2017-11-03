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


#ifndef __MMC_test__myRandomNumbers__
#define __MMC_test__myRandomNumbers__

#include <stdio.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//*****************************
//** Random Number Functions **
//*****************************
unsigned long long rdtsc();								//generates random seed based on rdtsc() function
int randint(int n);										//generates an integer random number from a uniform distribution between 0 and n-1
double randnum();										//generates a randon number from a uniform distribution between 0 and 1
double randexp(double lambda);							//generates a randon number from an exponential distribution with parameter 1/lambda
int randpoisson(double lambda);							//generates a randon number from a Possion distribution with mean lambda distribution
int randbern(double pi);								//generates a randon number from a Bernoulli distribution with parameter pi
int randdiscrete(std::vector<long double> lambda);


extern unsigned long long int seed;
void initializeRNG(unsigned long long int &seed);



#endif /* defined(__MMC_test__myRandomNumbers__) */


