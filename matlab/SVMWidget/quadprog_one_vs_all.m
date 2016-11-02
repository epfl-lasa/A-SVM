function [ x, learned_svm, exitflag ] = quadprog_one_vs_all( dynamics_data , data_labels, lambda, C, nu, tol, R, max_iter, target, initial_sol, nu_ASVM )
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

M = size(dynamics_data,2);
N = size(dynamics_data,1)/2;

% Set up the kernel matrix
disp('Calculating non-linear kernel...');
[K, G, G_s, H, H_s, H_ss] = getNLModulationKernel(dynamics_data, target, data_labels, lambda, ktype);


P = length(find(data_labels == 1));

    Q = [       K                      G                      -G_s;
        G'                    H                        -H_s
        -G_s'                       -H_s'                   H_ss            ];

[ee,vv] =eig(Q+Q');
disp(['Minimum eigen value of Q : ' num2str(min(vv(find(vv<0))))]);
% [ee,vv]= eig((G+G')/2);
% disp(['G : ' num2str(min(vv))]);

Q = Q+1e-11*eye(size(Q));
if(~isscalar(initial_sol))

    x0 = initial_sol;
else
    
    if(initial_sol == 0)
        x0 = zeros(M+P+N,1);
    end
    if(initial_sol == 1)
        x0 = C*rand(M+P+N,1);
    end
    if(initial_sol == 2)
            x0 = zeros(M,1);
            stroptions = ['-s 0 -t 2 -g ' num2str(lambda) ' -c 1e6 -e 1e-4 coef0 0 '];
            model = svmtrain(data_labels, dynamics_data(1:N,:)', stroptions);
            x0(model.SV_indices+1,1) = abs(model.sv_coef);
            x0 = [x0;zeros(P+N,1)];
    end
end

auxdata_1 = Q;
auxdata_3 = data_labels;
auxdata_4 = M;

% FMINCON options
lb = [zeros(M+P,1);-inf*ones(N,1)];
if(nu_ASVM)
    ub = [inf*ones(M,1);C/P*ones(P,1);inf*ones(N,1)];
else
    ub = [inf*ones(M,1);C*ones(P,1);inf*ones(N,1)];
end

obj = @(x)fmincon_objective(x, auxdata_4, auxdata_1);
options=[];


% pause;
% Run FMINCON.
tic


if(nu_ASVM)

    options=optimset('TolX',tol,'TolCon',tol,'Algorithm','interior-point-convex','Display','iter','TolFun',tol,'MaxIter',max_iter,'OutputFcn',@fmincon_callback,'MaxFunEvals',inf);
    [x, fval, exitflag] = quadprog((auxdata_1+auxdata_1')/2, [-ones(M,1);-R*ones(P,1);zeros(N,1)], -[zeros(1,M),ones(1,P), zeros(1,N)], -nu*C, [auxdata_3;zeros(P+N,1)]', 0, lb, ub, x0, options);
    
else

    %  C-ASVM with quadprog
    options=optimset('TolX',tol,'TolCon',tol,'Algorithm','interior-point-convex','Display','iter','TolFun',tol,'MaxIter',max_iter,'OutputFcn',@fmincon_callback,'MaxFunEvals',inf);
    [x, fval, exitflag] = quadprog((auxdata_1+auxdata_1')/2, [-ones(M,1);-R*ones(P,1);zeros(N,1)],[], [], [auxdata_3;zeros(P+N,1)]', 0, lb, ub, x0, options);

end
if(nu_ASVM)
    disp(['Expected number of SV    = ' num2str(nu*P)]);
    for i=1:M
       for j=1:M
          Q(i,j) = Q(i,j)*data_labels(i)*data_labels(j); 
       end
       for j=M+1:M+P+N
          Q(i,j) = Q(i,j)*data_labels(i);
          Q(j,i) = Q(i,j);
       end
    end
    er = Q(M+1:M+P,:)*(x.*[data_labels(1:M);ones(P+N,1)]);
    er_ind = find(er < 0);
    disp(['Number of violating points = ' int2str(length(er_ind))]);
    disp(['sum(beta) , C*nv           = ' num2str(sum(x(M+1:M+P))) ' , ' num2str(nu*C)]);
    disp(' ');
end
toc
%%

disp('Parsing fmincon solution to SVM object');

sv_tol_alpha = max(x(1:M+P))*1e-5;
sv_tol_beta = sv_tol_alpha;

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


    plt = [];
for i=1:size(dynamics_data,2)
    if(data_labels(i) == 1)
       dr = calculateClassifierDerivative(learned_svm, dynamics_data(1:2,i));
       tmp = dr'*dynamics_data(3:4,i);
       plt(i,:) = [tmp , tmp/norm(dr)];
    end
end
% figure
% plot(plt,'*-');
% hold on
% 
% legend('Un-Normalized dot product', 'Normalized dot product');
% if(~nu_ASVM)
% lm = get(gca,'XLim');
% plot(lm',[eps;eps],'k-','Linewidth',2); 
% end
