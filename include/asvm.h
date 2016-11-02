/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 *
 *  asvm.h
 *
 *  Created on : Aug 14, 2012
 *  Author     : Ashwini Shukla and Saurav Aryan
 *  Email      : ashwini.shukla@epfl.ch
 *  Website    : lasa.epfl.ch
 *
 * Permission is granted to copy, distribute, and/or modify this program
 * under the terms of the GNU General Public License, version 2 or any
 * later version published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details
 */
#ifndef _asvm_H_
#define _asvm_H_


#include "util.h"
#include "asvmdata.h"

#define MYMAX(a,b) ((a)>=(b) ? (a):(b))
#define MYMIN(a,b) ((a)<=(b) ? (a):(b))

class asvm
{
	friend class ASVM_SMO_Solver;
	friend class ASVM_NLopt_Solver;
	friend class ASVM_ALGLIB_Solver;
public:

	double* alpha;
	double* beta;
	double* gamma;
	int *y;
	unsigned int numAlpha;
	unsigned int numBeta;
	unsigned int dim;
	double** svalpha;
	double** svbeta;
	double* lambda;
    char type[1024];
	double b0;
    double* target;

	void setSVMFromOptimizedSol(double *x, asvmdata* dat);
public:

    asvm()
        : numAlpha(0), numBeta(0), dim(0), b0(0), lambda(0),
          alpha(NULL), beta(NULL), gamma(NULL), y(NULL),
          svalpha(NULL), svbeta(NULL), target(NULL),
          temp(NULL), temp1(NULL), temp2(NULL), tempMat(NULL)
    {}
    asvm(const asvm& o);
    asvm& operator=(const asvm& o);
    ~asvm();

	double getclassifiervalue(double *pt);
	void getclassifierderivative(double *point, double* derivative);
	void printinfo();
	void saveToFile(const char* filename);
	void calcb0();

private:
	double *temp;
	double *temp1;
	double *temp2;
	double ** tempMat;
};

#endif
