/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * 
 *  mxSMOSolver.cpp
 *
 *  Created on : Aug 14, 2012
 *  Author     : Ashwini Shukla
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


#include "mex.h"
#include "mxUtil.h"
#include "../../include/asvm_smo_solver.h"



#ifdef MX_API_VER
#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif
#endif


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    
	if(nrhs == 0)
	{
		printf("Usage: mxSMOSolver ( data_cell_object, target_class, kernel_width, [parameter_file_name]\n");
		return;
	}

	size_t tclass = (size_t)mxGetScalar(prhs[1])-1;

	double sigma = mxGetScalar(prhs[2]);
	char paramfilename[1025];
	if(nrhs < 4)
			strcpy(paramfilename, " ");
	else
		mxGetString(prhs[3], paramfilename, 1025);


	 asvmdata* dataobj = mxUtil_Data_mat2c(prhs[0]);


	dataobj->setParams("rbf", sigma, 1.0);
	asvm svmobj;
	ASVM_SMO_Solver smo_solver;
	smo_solver.configure(paramfilename);

	smo_solver.learn(*dataobj, tclass, &svmobj);

	//Converting asvm class instance to matlab struct
	mxArray* learned_svm = mxUtil_SVM_c2mat(&svmobj);

	nlhs=1;
	plhs[0] = learned_svm;

	delete dataobj;
//	delete &svmobj;

     
}
