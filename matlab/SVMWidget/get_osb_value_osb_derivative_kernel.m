function Q = get_osb_value_osb_derivative_kernel( input_data, target, labels, d, type )
% This function calculates the kernel matrix for A-SVM
%
%   Inputs ----------------------------------------------------------------
%   o input_data :  D x N matrix of N data points of dimension D
%   o target     :  D x 1 vector representing the target location
%   o labels     :  N x 1 label vector of +1/-1
%   o d          :  Scalar kernel parameter
%   o type       :  String kernel type - 'rbf', 'poly' etc.
%
%   Outputs ---------------------------------------------------------------
%   o Q          :  kernel matrix
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

pl = find(labels==1);

M = size(input_data,2);
N = size(input_data,1)/2;
P = length(pl);

K = zeros(M,M);
G = zeros(M,P);
G_s = zeros(M,N);
H = zeros(P,P);
H_s = zeros(P,N);

for i=1:M
    yi = labels(i);
    xi = input_data(1:N,i);
    for j=1:M
        K(i,j) = yi*labels(j)*getKernel(xi, input_data(1:N,j), d, type);
    end
end

for i=1:M
    yi = labels(i);
    xi = input_data(1:N,i);
    for j=1:P
        G(i,j) = yi*getKernelFirstDerivative(xi,input_data(1:N,pl(j)), d, type, 2)'*input_data(N+1:end,pl(j));
    end
end


for i=1:M
    G_s(i,1:N) = labels(i)*getKernelFirstDerivative(input_data(1:N,i), target, d, type, 2)';
end

for i=1:P
    xi = input_data(1:N,pl(i));
    xidot = input_data(N+1:end,pl(i));
    for j=1:P
        H(i,j) = xidot'*getKernelSecondDerivative(xi, input_data(1:N,pl(j)), d, type)*input_data(N+1:end,pl(j));
    end
end

for i=1:P
    H_s(i,1:N) = input_data(N+1:end,pl(i))'*getKernelSecondDerivative(input_data(1:N, pl(i)), target, d, type);
end


H_ss = getKernelSecondDerivative(target, target, d, type);

Q=[K G -G_s;G' H -H_s;-G_s' -H_s' H_ss];
