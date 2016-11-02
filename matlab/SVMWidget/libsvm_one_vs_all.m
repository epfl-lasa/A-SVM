function [ learned_svm ] = libsvm_one_vs_all( all_data, lambda, target_class )
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

libsvm_data = [];
data_labels=[];
dim = size(all_data{1}{1},1);

target = zeros(dim,1);
for i = 1:length(all_data)
    tmp = all_data{i};
    num_points = 0;
    for j=1:length(tmp)
        libsvm_data = [ libsvm_data ; tmp{j}(:,1:end-1)' ];
        num_points = num_points + size(tmp{j},2)-1;
        if(i == target_class)
            target = target + tmp{j}(:,end);
        end
    end
    
    if(i==target_class)
        one_of_k_target = 1;
        target = target/length(tmp);
    else
        one_of_k_target = -1;
    end
    
    data_labels = [data_labels ; repmat(one_of_k_target, num_points,1)];
    
end

%% For libsvm
stroptions = '-s 0 -t 4  -c 1e6 -e 1e-4 coef0 0 ';
M = length(data_labels);
K = zeros(M,M);
for i=1:M
    for j=i:M
        K(i,j) = data_labels(i)*data_labels(j)*getKernel(libsvm_data(i,:)', libsvm_data(j,:)', lambda, 'rbf');
        K(j,i)=K(i,j);
    end
end
model = svmtrain(data_labels, [[1:M]' K], stroptions);
structSVM=[];
structSVM.y = data_labels(model.SV_indices+1);
structSVM.Sva =  libsvm_data(model.SV_indices+1,:)';
structSVM.alpha = abs(model.sv_coef);
structSVM.Svb=[];
structSVM.beta=[];
structSVM.gamma=[];
structSVM.b0 = -model.rho;
structSVM.lambda = lambda;
structSVM.type='rbf';
structSVM.target = target;
if(model.Label(1) == -1)
   structSVM.b0 = -structSVM.b0;
end
learned_svm = structSVM;

