function Q = get_osb_value_osb_derivative_kernel_dsigma( input_data, target, labels, d, type )
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
Kp = zeros(M,M);
Kt=zeros(M,1);
G = zeros(M,P);
G_s = zeros(M,N);
H = zeros(P,P);
H_s = zeros(P,N);

for i=1:M
    yi = labels(i);
    x1 = input_data(1:N,i);
    Kt(i) = getKernel(x1, target,d, type);
    for j=1:M
        K(i,j) = getKernel(x1, input_data(1:N,j), d, type);
    end
end

Q=cell(N,1);
D=diag(d);
for l=1:N
    for i=1:M
        Kp(i,:) = -labels(i)*labels'.*K(i,:).*(input_data(l,i) - input_data(l,:)).^2;
    end
    
    for i=1:M
        yi = labels(i);
        x1 = input_data(1:N,i);
        %     for j=1:P
        %         kij=K(i,pl(j));
        %         v = x1-input_data(1:N,pl(j));
        %         tmpvec = -2*(v(l))^2*kij*(v.*d);
        %         tmpvec(l) = tmpvec(l) + 2*kij*(v(l));
        %         G(i,j) = yi*tmpvec'*input_data(N+1:end,pl(j));
        %     end
        ind = pl(1:P);
        x2=input_data(1:N,ind);
        kij = K(i,ind);
        tmpvec = x1(l)-x2(l,:);
        tmpmat= -2*repmat((tmpvec).^2,N,1).*repmat(kij,N,1).*((repmat(x1,1,P)-x2).*repmat(d,1,P));
        tmpmat(l,:) = tmpmat(l,:)+2*kij.*(tmpvec);
        G(i,:) = yi*dot(tmpmat,input_data(N+1:end,ind));
    end
    
    
    for i=1:M
        x1=input_data(1:N,i);
        kij=Kt(i);
        tmpvec = -2*(x1(l)-target(l))^2*kij*((x1-target).*d)';
        tmpvec(l) = tmpvec(l) + 2*kij*(x1(l)-target(l));
        G_s(i,1:N) = labels(i)*tmpvec';
    end
    
    
    for i=1:P
        x1 = input_data(1:N,pl(i));
        x1dot = input_data(N+1:end,pl(i));
        for j=1:P
            kij=K(pl(i),pl(j));
            x2 = input_data(1:N,pl(j));
            tmpval = (x1(l)-x2(l))^2;
            v=(x1-x2).*d;
            tmpmat = 2*kij*tmpval*(2*(v*v') - D);
            tmpmat(l,:) = -4*(x1(l)-x2(l))*kij*(1-d(l)*tmpval)*v';
            tmpmat(:,l) = tmpmat(l,:);
            tmpmat(l,l) = -4*tmpval*kij*d(l)*(2-d(l)*tmpval) + 2*kij*(1-d(l)*tmpval);
            
            H(i,j) = x1dot'*tmpmat*input_data(N+1:end,pl(j));
        end
    end
    
    for i=1:P
        x1=input_data(1:N,pl(i));
        kij = Kt(pl(i));
        tmpval = (x1(l)-target(l))^2;
        v=(x1-target).*d;
        tmpmat = 2*kij*tmpval*(2*(v*v') - D);
        tmpmat(l,:) = -4*(x1(l)-target(l))*kij*(1-d(l)*tmpval)*v';
        tmpmat(:,l) = tmpmat(l,:);
        tmpmat(l,l) = -4*tmpval*kij*d(l)*(2-d(l)*tmpval) + 2*kij*(1-d(l)*tmpval);
        H_s(i,1:N) = input_data(N+1:end,pl(i))'*tmpmat;
    end
    
    H_ss = zeros(N);
    H_ss(l,l) =  2*getKernel(target, target, d, type);
    
    Q{l}=[Kp G -G_s;G' H -H_s;-G_s' -H_s' H_ss];
end