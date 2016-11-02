function Q = get_bsb_value_osb_derivative_kernel( input_data, d, type )
% This function calculates the kernel matrix for surface reconstruction 
% using A-SVM
%
%   Inputs ----------------------------------------------------------------
%   o input_data :  D x N matrix of N data points of dimension D
%   o target     :  D x 1 vector representing the target location
%   o labels     :  N x 1 label vector of +1/-1
%   o d          :  Scalar kernel parameter
%   o type       :  String kernel type - 'rbf', 'poly' etc.
%
%   Outputs ---------------------------------------------------------------
%   o K, G, H   :  Block matrices which form the full kernel matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Copyright (c) 2012 Ashwini SHUKLA, LASA, EPFL,          %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The program is free for non-commercial academic use. Please contact the
% author if you are interested in using the software for commercial purposes.
% The software must not be modified or distributed without prior permission
% of the authors. Please acknowledge the authors in any academic publications
% that have made use of this code or part of it. Please use this BibTex
% reference:
%
%
% To get latest upadate of the software please visit
%                          http://asvm.epfl.ch
%
% Please send your feedbacks or questions to:
%                           ashwini.shukla_at_epfl_dot_ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = size(input_data,2);
N = size(input_data,1)/2;
lambda=d;
ktype=type;

% Set up the kernel matrix
disp('Calculating non-linear kernel...');
K=zeros(M,M);
G=zeros(M,M);
H = zeros(M, M);

for i=1:M
    xi = input_data(1:N,i);
    for j=1:M
        xj = input_data(1:N,j);
        K(i,j) = getKernel(xi, xj, lambda, ktype);
        G(i,j) = getKernelFirstDerivative(xi, xj, lambda, ktype, 2)'*input_data(N+1:end,j);
        H(i,j) = input_data(N+1:end,i)'*getKernelSecondDerivative(xi, xj, lambda, ktype)*input_data(N+1:end,j);
    end
end

Q = [K   -K   G;
    -K    K  -G;
     G'  -G'  H];




