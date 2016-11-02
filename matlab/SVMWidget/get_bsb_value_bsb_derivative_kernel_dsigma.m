function Q = get_bsb_value_bsb_derivative_kernel_dsigma( input_data, d, type)
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

ktype=type;

% Set up the kernel matrix

K=zeros(M,M);
Kp=zeros(M,M);
Gp=zeros(M,N*M);
Hp = zeros(N*M, N*M);
D=diag(d);

for i=1:M
    x1 = input_data(1:N,i);
    for j=1:M
        K(i,j) = getKernel(x1, input_data(1:N,j), d, ktype);
    end
end

Q=cell(N,1);
for l=1:N
    
    for i=1:M
        Kp(i,:) = -K(i,:).*(input_data(l,:)-input_data(l,i)).^2;
    end
    
    for i=1:M
        x1 = input_data(1:N,i);
        for j=1:M
            x2 = input_data(1:N,j);
            tmpvec = -2*(x1(l)-x2(l))^2*K(i,j)*((x1-x2).*d)';
            tmpvec(l) = tmpvec(l) + 2*K(i,j)*(x1(l)-x2(l));
            Gp(i,j+(0:N-1)*M) = tmpvec;
        end
    end
    
    
    
    for i=1:M
        x1 = input_data(1:N,i);
        for j=1:M
            x2 = input_data(1:N,j);
            tmpval = (x1(l)-x2(l))^2;
            v=(x1-x2).*d;
            tmpmat = 2*K(i,j)*tmpval*(2*(v*v') - D);
            tmpmat(l,:) = -4*(x1(l)-x2(l))*K(i,j)*(1-d(l)*tmpval)*v';
            tmpmat(:,l) = tmpmat(l,:);
            tmpmat(l,l) = -4*tmpval*K(i,j)*d(l)*(2-d(l)*tmpval) + 2*K(i,j)*(1-d(l)*tmpval);
            Hp(i+(0:N-1)*M,j+(0:N-1)*M) = tmpmat;
        end
    end
    
    Qtmp = [Kp -Gp; -Gp' Hp];
    Qtmp2=[];
    for i=1:N+1
        tmp=[];
        for j=1:N+1
            ri=(i-1)*M+1:i*M;
            rj=(j-1)*M+1:j*M;
            tmp = [tmp,Qtmp(ri,rj), -Qtmp(ri,rj)];
        end
        Qtmp2=[Qtmp2;tmp;-tmp];
    end
    
    Q{l} = Qtmp2;
end

