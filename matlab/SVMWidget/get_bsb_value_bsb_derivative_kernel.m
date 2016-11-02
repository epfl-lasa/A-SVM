function Q = get_bsb_value_bsb_derivative_kernel( input_data, d, type )
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


K=zeros(M,M);
G=zeros(M,N*M);
H = zeros(N*M, N*M);
for i=1:M
    xi = input_data(1:N,i);
    for j=1:M
        xj = input_data(1:N,j);
        K(i,j) = getKernel(xi, xj, lambda, ktype);
    end
end


    for i=1:M
        xi = input_data(1:N,i);
        for j=1:M
            xj = input_data(1:N,j);
            G(i,j+(0:N-1)*M) = getKernelFirstDerivative(xi, xj, lambda, ktype, 2)';
        end
    end



        for i=1:M
            xi = input_data(1:N,i);
            for j=1:M
                xj = input_data(1:N,j);
                H(i+(0:N-1)*M,j+(0:N-1)*M) = getKernelSecondDerivative(xi, xj, lambda, ktype);
            end
        end

Qtmp = [K -G; -G' H];
Q=[];
for i=1:N+1
    tmp=[];
    for j=1:N+1
        ri=(i-1)*M+1:i*M;
        rj=(j-1)*M+1:j*M;
        tmp = [tmp,Qtmp(ri,rj), -Qtmp(ri,rj)];
    end
    Q=[Q;tmp;-tmp];
end



