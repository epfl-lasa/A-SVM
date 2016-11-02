/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * 
 *  mxUtil.cpp
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

#include "mxUtil.h"
#ifdef MX_API_VER
#if MX_API_VER < 0x07030000
typedef int mwIndex;
#endif
#endif

asvmdata* mxUtil_Data_mat2c(const mxArray *mat_data)
{
	size_t num_targets = mxGetNumberOfElements(mat_data);
	size_t num_traj, dim, traj_len;
	double* trajdata;
	asvmdata *cpp_data = new asvmdata();
	mxArray* curr_target = mxGetCell(mat_data, 0);
	mxArray* curr_traj = mxGetCell(curr_target, 0);
	dim = mxGetM(curr_traj);
	cpp_data->dim = dim;

	unsigned int i,j,k,l;
	for(i=0;i<num_targets;i++)
	{
		target t;
		t.targ = new double[dim];
		t.dim = dim;
		for(j=0;j<dim;j++)
			t.targ[j] = 0;
		curr_target = mxGetCell(mat_data, i);
		num_traj = mxGetM(curr_target);

		for(j=0;j<num_traj;j++)
		{
			trajectory tr;
			curr_traj = mxGetCell(curr_target, j);
			trajdata = (double*)mxGetPr(curr_traj);
			traj_len = mxGetN(curr_traj);
			dim = mxGetM(curr_traj);

			tr.dim = dim;
			tr.nPoints = traj_len;
			tr.y = new int[traj_len];
			tr.coords = new double*[traj_len];
			tr.vel = new double*[traj_len];
			for(k=0;k<traj_len;k++)
			{
				tr.coords[k] = new double[dim];
				tr.vel[k] = new double[dim];
			}

			for(k=0;k<dim;k++)
			{
				for(l=0;l<traj_len;l++)
					tr.coords[l][k] = trajdata[k+l*dim];

			}
			for(l=0;l<traj_len;l++)
				tr.y[l] = i;

			for(l=0;l<dim;l++)
				t.targ[l] += tr.coords[traj_len-1][l];

			t.traj.push_back(tr);
		}

		for(l=0;l<dim;l++)
			t.targ[l] /= (double)num_traj;

		cpp_data->addTarget(t);
	}

	cpp_data->isOkay = true;
	return cpp_data;
}

asvm* mxUtil_SVM_mat2c(const mxArray *mat_svm)
{
	mxArray* tmp;


	asvm* svm_obj = new asvm;

	tmp = mxGetField(mat_svm, 0, "alpha");
	svm_obj->numAlpha = mxGetNumberOfElements(tmp);


	svm_obj->alpha = (double*)mxGetPr(tmp);
	tmp = mxGetField(mat_svm, 0, "Sva");
	double* tmpdat = mxGetPr(tmp);
	svm_obj->dim = mxGetM(tmp);
	svm_obj->svalpha = new double*[svm_obj->numAlpha];
	for(unsigned int i =0; i<svm_obj->numAlpha; i++)
		svm_obj->svalpha[i] = new double[svm_obj->dim];
	unsigned int cnt=0;
	for(unsigned int i=0;i<svm_obj->numAlpha;i++)
		for(unsigned int j=0;j<svm_obj->dim;j++)
			svm_obj->svalpha[i][j] = tmpdat[cnt++];

	tmp = mxGetField(mat_svm, 0, "y");
	double* tmplabel = (double*)mxGetPr(tmp);
	svm_obj->y = new int[svm_obj->numAlpha];
	for(unsigned int i=0;i<svm_obj->numAlpha;i++)
		svm_obj->y[i] = (int)tmplabel[i];

	tmp = mxGetField(mat_svm, 0, "beta");
	svm_obj->numBeta = mxGetNumberOfElements(tmp);
	svm_obj->beta = (double*)mxGetPr(tmp);
	tmp = mxGetField(mat_svm, 0, "Svb");
	svm_obj->svbeta = new double*[svm_obj->numBeta];
    for(unsigned int i =0; i<svm_obj->numBeta; i++)
    	svm_obj->svbeta[i] = new double[svm_obj->dim*2];
	tmpdat = mxGetPr(tmp);
	cnt=0;
	for(unsigned int i=0;i<svm_obj->numBeta;i++)
		for(unsigned int j=0;j<2*svm_obj->dim;j++)
			svm_obj->svbeta[i][j] = tmpdat[cnt++];

	tmp = mxGetField(mat_svm, 0, "gamma");
	svm_obj->gamma = (double*)mxGetPr(tmp);

	tmp = mxGetField(mat_svm, 0, "lambda");
	if(mxGetNumberOfElements(tmp) == 1)
	{
		svm_obj->lambda = new double[svm_obj->dim];
		double l=mxGetScalar(tmp);
		for(unsigned int i=0;i<svm_obj->dim;i++)
			svm_obj->lambda[i] = l;
	}
	else if(mxGetNumberOfElements(tmp) == svm_obj->dim)
		svm_obj->lambda = (double*)mxGetPr(tmp);
	else
	{
		printf("ERROR: lambda can have length = 1 or equal to the number of dimensions \n");
		return NULL;
	}

	tmp = mxGetField(mat_svm, 0, "type");
	mxGetString(tmp, svm_obj->type, 1025);

	tmp = mxGetField(mat_svm, 0, "b0");
	svm_obj->b0 = mxGetScalar(tmp);

	tmp = mxGetField(mat_svm, 0, "target");
	svm_obj->target = (double*)mxGetPr(tmp);

	return svm_obj;

}

mxArray* mxUtil_SVM_c2mat(asvm* cpp_svm)
{
	const char** fnames = (const char**)mxCalloc((size_t)10, sizeof(*fnames));
	fnames[0] = "type";
	fnames[1] = "lambda";
	fnames[2] = "b0";
	fnames[3] = "target";
	fnames[4] = "alpha";
	fnames[5] = "y";
	fnames[6] = "beta";
	fnames[7] = "gamma";
	fnames[8] = "Sva";
	fnames[9] = "Svb";

	mxArray* mat_svm = mxCreateStructMatrix(1,1,10, fnames);

	mxSetField(mat_svm, 0, "type", mxCreateString("rbf"));
	mxSetField(mat_svm, 0, "b0", mxCreateDoubleScalar(cpp_svm->b0));



	mxArray *mxptr= mxCreateDoubleMatrix(cpp_svm->dim, 1, mxREAL);
	double *cptr = mxGetPr(mxptr);
	for(unsigned int i=0;i<cpp_svm->dim;i++)
		cptr[i] = cpp_svm->target[i];
	mxSetField(mat_svm, 0, "target", mxptr);

	mxptr = mxCreateDoubleMatrix(cpp_svm->dim,1,mxREAL);
	cptr = mxGetPr(mxptr);
	for(unsigned int i=0;i<cpp_svm->dim;i++)
		cptr[i] = cpp_svm->lambda[i];
	mxSetField(mat_svm, 0, "lambda", mxptr);

	mxptr = mxCreateDoubleMatrix(cpp_svm->numAlpha, 1 ,mxREAL);
	cptr = mxGetPr(mxptr);
	for(unsigned int i=0;i<cpp_svm->numAlpha;i++)
		cptr[i] = cpp_svm->alpha[i];
	mxSetField(mat_svm, 0, "alpha", mxptr);

	mxptr = mxCreateDoubleMatrix(cpp_svm->numAlpha, 1 ,mxREAL);
	cptr = mxGetPr(mxptr);
	for(unsigned int i=0;i<cpp_svm->numAlpha;i++)
		cptr[i] = (double)cpp_svm->y[i];
	mxSetField(mat_svm, 0, "y", mxptr);

	mxptr = mxCreateDoubleMatrix(cpp_svm->numBeta, 1 ,mxREAL);
	cptr = mxGetPr(mxptr);
	for(unsigned int i=0;i<cpp_svm->numBeta;i++)
		cptr[i] = cpp_svm->beta[i];
	mxSetField(mat_svm, 0, "beta", mxptr);

	mxptr = mxCreateDoubleMatrix(cpp_svm->dim, 1 ,mxREAL);
	cptr = mxGetPr(mxptr);
	for(unsigned int i=0;i<cpp_svm->dim;i++)
		cptr[i] = cpp_svm->gamma[i];
	mxSetField(mat_svm, 0, "gamma", mxptr);

	mxptr = mxCreateDoubleMatrix(cpp_svm->dim, cpp_svm->numAlpha ,mxREAL);
	cptr = mxGetPr(mxptr);
	int cnt=0;
	for(unsigned int i=0;i<cpp_svm->numAlpha;i++)
		for(unsigned int j=0;j<cpp_svm->dim;j++)
			cptr[cnt++] = cpp_svm->svalpha[i][j];
	mxSetField(mat_svm, 0, "Sva", mxptr);

	mxptr = mxCreateDoubleMatrix(2*cpp_svm->dim, cpp_svm->numBeta ,mxREAL);
	cptr = mxGetPr(mxptr);
	cnt=0;
	for(unsigned int i=0;i<cpp_svm->numBeta;i++)
		for(unsigned int j=0;j<2*cpp_svm->dim;j++)
			cptr[cnt++] = cpp_svm->svbeta[i][j];
	mxSetField(mat_svm, 0, "Svb", mxptr);

	return mat_svm;

}

