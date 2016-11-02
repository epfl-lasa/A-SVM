/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 *
 *  asvm_nlopt_solver.h
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

#ifndef _ASVM_NLOPT_SOLVER_H_
#define _ASVM_NLOPT_SOLVER_H_

#include <nlopt.hpp>
#include "asvm.h"


typedef struct {
	unsigned int num_alpha;
	double** NLModulationKernel;
	bool bVerbose;
}obj_func_data;

typedef struct {
	unsigned int num_alpha;
	int *data_labels;
	bool bVerbose;
} constr_func_data;



class ASVM_NLopt_Solver
{
private:
	nlopt_opt 			main_opt;
	nlopt_algorithm 	main_algo;
	asvmdata* copy_data;
	nlopt_algorithm 	sub_algo;



	double 	constr_tol;
	double 	fun_rel_tol;
	double 	x_rel_tol;
	long	max_iter;
	long	max_time_sec;
	double	paramC;
	bool bVerbose;

public:

	ASVM_NLopt_Solver()
	{
		main_algo = NLOPT_LD_AUGLAG;
			sub_algo = NLOPT_LD_LBFGS;
			paramC = 1e6;
			constr_tol = 1e-5;
			fun_rel_tol = 1e-10;
			x_rel_tol = 1e-10;
			max_iter = 1e6;
			max_time_sec = 120;
			bVerbose = false;
	}

	int learn(asvmdata& input1, unsigned int tclass, asvm* svmobj);
	void configure(const char* paramFile);
	void force_stop() {  nlopt_set_force_stop(main_opt,1);}

private:
	double forward(int index, double* x_smo, asvmdata* dat);
};


#endif
