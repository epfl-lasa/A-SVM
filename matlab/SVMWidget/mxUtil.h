/*
 * Copyright (C) 2012 Learning Algorithms and Systems Laboratory, EPFL, Switzerland
 * 
 *  mxUtil.h
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


#ifndef MXUTIL_H_
#define MXUTIL_H_

#include "mex.h"
#include "../../include/asvm.h"

asvmdata* mxUtil_Data_mat2c(const mxArray *mat_data);
mxArray* mxUtil_SVM_c2mat(asvm* cpp_svm);
asvm* mxUtil_SVM_mat2c(const mxArray *mat_svm);

#endif /* MXUTIL_H_ */
