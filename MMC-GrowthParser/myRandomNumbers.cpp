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



#include "myRandomNumbers.hpp"

unsigned long long int seed = 0;
gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);

//*****************************
//** Random Number Functions **
//*****************************

void initializeRNG(unsigned long long int &seed)
{
	if(seed <= 0)
	{
		seed = rdtsc();
	}
	gsl_rng_set (rng, seed);                        //this seeds the RNG
}

//***randnum********************************************
//  Generates a random number from a uniform distribution between 0 and 1

double randnum()
{
	return (double) gsl_rng_uniform (rng);
}

//***randint********************************************
//  Generates a discrete random number from a uniform distribution between 0 and n-1

int randint(int n)
{
	return (int) gsl_rng_uniform_int (rng, (unsigned long int) n);
}


//***randexp******************************************
//  Generates a random number from an exponential distribution with
//  parameter (1/lambda)

double randexp(double lambda)
{
	return (double) gsl_ran_exponential (rng, 1/lambda);
}



//   Generates a random seed
unsigned long long rdtsc()
{
	unsigned int lo,hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return ((unsigned long long) hi << 32) | lo;
}



//***randpoisson******************************************
//  Generates a random number from a Poisson distribution with mean
//  lambda

int randpoisson(double lambda)
{
	return (int) gsl_ran_poisson (rng, lambda);
}

//***randbern********************************************
//  Generates a random number from a Bernoulli distribution with  probability pi

int randbern(double pi)
{
	return (int) gsl_ran_bernoulli (rng, pi);
}



//***randmultinomial********************************************
//  This function computes a random sample n[lambda.size()] from the multinomial distribution formed by a single trial from an underlying distribution lambda.
//
int randdiscrete(std::vector<long double> lambda)
{
	int _tempPos = 0;
	double _p[lambda.size()];
	std::copy(lambda.begin(), lambda.end(), _p);

	unsigned int _n[lambda.size()];
	gsl_ran_multinomial (rng, lambda.size(), 1, _p, _n);

	while(_n[_tempPos] == 0)
	{
		_tempPos++;
	}
	
	return (_tempPos);
}
