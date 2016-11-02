function [ x, learned_svm ] = ipopt_one_vs_all( dynamics_data , data_labels, lambda, C, tol, max_iter, target, initial_sol )
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
ktype = 'rbf';
use_gamma = 1;

M = size(dynamics_data,2);
N = size(dynamics_data,1)/2;

% Set up the kernel matrix
disp('Calculating non-linear kernel...');
[K, G, G_s, H, H_s, H_ss] = getNLModulationKernel(dynamics_data, target, data_labels, lambda, ktype);


P = length(find(data_labels == 1));

if(use_gamma == 0)
    
    Q = [       K     G;
                G'    H            ];
    
else
    
    Q = [       K                      G                      -G_s;
        G'                    H                        -H_s
        -G_s'                       -H_s'                   H_ss            ];
    
end

[ee,vv] =eig(Q+Q');
disp(['Minimum eigen value of Q : ' num2str(min(vv(find(vv<0))))]);
pause(2);
% Q = Q+1e-11*eye(size(Q));
if(~isscalar(initial_sol))

    x0 = initial_sol;
else
    
    if(initial_sol == 0)
        x0 = zeros(M+P+N,1);
    end
    if(initial_sol == 1)
        x0 = -2000*rand(M+P+N,1) + 1000;
    end
    if(initial_sol == 2)
            x0 = zeros(M,1);
            ker = kernel(ktype,lambda);
            algo =  svm(ker);
            [tr_data, algo] = train(algo, data(dynamics_data(1:N,:)', data_labels));
            invalid_alpha_index = find((algo.alpha));
            disp('Calculating initial solution...');
            disp(['Got ' num2str(length(invalid_alpha_index)) ' support vectors']);
            x0(invalid_alpha_index) = abs(algo.alpha(invalid_alpha_index) );
            x0 = [x0;zeros(P+N,1)];
    end
end

% IPOPT options
options = [];
if(use_gamma == 0)
    options.lb = zeros(M+P,1);
    options.ub = [inf*ones(M,1);C*ones(P,1)];
else
    options.lb = [zeros(M+P,1);-inf*ones(N,1)];
    options.ub = [inf*ones(M,1);C*ones(P,1);inf*ones(N,1)];
end
options.cl = 0;
options.cu = 0;

if(use_gamma == 0)
    options.auxdata = {Q (Q+Q') data_labels M P N [ones(M,1);zeros(P,1)]            use_gamma};
else
    options.auxdata = {Q (Q+Q') data_labels M P N [ones(M,1);zeros(P,1);zeros(N,1)] use_gamma};
end

options.ipopt.hessian_approximation         = 'limited-memory';
options.ipopt.mu_strategy                   = 'adaptive';
% options.ipopt.derivative_test             = 'first-order';
options.ipopt.jac_c_constant                = 'yes';
options.ipopt.jac_d_constant                = 'yes';
options.ipopt.max_iter                      = max_iter;
options.ipopt.tol                           = tol;
options.ipopt.print_level                   = 5;
options.ipopt.max_iter         
% The callback functions.
funcs = [];
funcs.objective         = @ipopt_objective;
funcs.constraints       = @ipopt_constraints;
funcs.gradient          = @ipopt_gradient;
funcs.jacobian          = @ipopt_jacobian;
funcs.jacobianstructure = @ipopt_jacobianstructure;
funcs.iterfunc          = @ipopt_callback;

% Run IPOPT.
if(options.ipopt.max_iter)
    [x info] = ipopt(x0,funcs,options);
else
    x=x0;
end

%%

disp('Parsing ipopt solution to SVM object');

sv_tol_alpha = max(x(1:M))*1e-4;
sv_tol_beta = max(x(M+1:M+P))*1e-4;

positive_labels = find(data_labels == 1);
alpha_index = find(abs(x(1:M))>sv_tol_alpha);       % find SVs
beta_index = positive_labels(find(x(M+1:M+P) > sv_tol_beta));

learned_svm = [];
learned_svm.alpha = x(alpha_index);  
learned_svm.beta = x(M+find(x(M+1:M+P) > sv_tol_beta));
learned_svm.gamma = x(M+P+1:end);
learned_svm.y = data_labels(alpha_index);
learned_svm.Sva = dynamics_data(1:N,alpha_index);
learned_svm.Svb = dynamics_data(:,beta_index);
learned_svm.lambda = lambda;
learned_svm.type = ktype;
learned_svm.target = target;



% Calculating b0
learned_svm.b0 = 0;
f_svp = [];
f_svn = [];
n_p = 0;
n_n = 0;
for i=1:size(dynamics_data,2)
    tmp = 0;
    pt = dynamics_data(1:N,i);
   tmp = calculateClassifier(learned_svm, pt);
    
    if(data_labels(i) == 1)
       n_p = n_p + 1;
       f_svp = [f_svp;tmp];
    else
        n_n = n_n + 1;
        f_svn = [f_svn;tmp];
    end
end

learned_svm.b0 = -(max(f_svn) + min(f_svp))/2;


