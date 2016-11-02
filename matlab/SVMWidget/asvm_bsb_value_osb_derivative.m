function learned_svm= asvm_bsb_value_osb_derivative( dynamics_data, y, Q, lambda, ktype,  C, tol, eps_nu, use_nu )
%LEARN_SURFACE_ASVM Summary of this function goes here
%   Detailed explanation goes here

max_iter = 1000;
use_derivatives = 1;
if(isscalar(C))
   C=[C;C]; 
end

if(isscalar(eps_nu))
   eps_nu=[eps_nu;eps_nu]; 
end


M = size(dynamics_data,2);
N = size(dynamics_data,1)/2;

y = reshape(y,M,1);

if(use_nu)
    f= [y;-y; zeros(M,1)] ;
else
    f= [y;-y; zeros(M,1)] + eps_nu(1)*[ones(M,1);ones(M,1); zeros(M,1)] - eps_nu(2)*[zeros(2*M,1);ones(M,1)];
end
% FMINCON options
lb = zeros(3*M,1);
if(use_derivatives)
    if(use_nu)
        ub = [C(1)/M*ones(2*M,1);C(2)/M*ones(M,1)];
    else
        ub = [C(1)*ones(2*M,1);C(2)*ones(M,1)];
    end
else
    if(use_nu)
        ub = [C(1)/M*ones(2*M,1);zeros(M,1)];
    else
        ub = [C(1)*ones(2*M,1);zeros(M,1)];
    end
end
options=optimset('TolX',tol,'TolCon',tol,'Algorithm','interior-point-convex','Display','iter','TolFun',tol,'MaxIter',max_iter,'OutputFcn',@fmincon_callback,'MaxFunEvals',inf);

tic
if(use_nu)
    [x, fval, exitflag] = quadprog((Q+Q')/2, f, [ones(1,2*M),zeros(1,M);zeros(1,2*M), -ones(1,M)],[eps_nu(1)*C(1);-eps_nu(2)*C(2)],[],[],lb,ub, [], options);
else
    [x, fval, exitflag] = quadprog((Q+Q')/2, f, [],[],[],[],lb,ub, [], options);
end
toc
if(use_nu)
x'*[ones(1,2*M),zeros(1,M);zeros(1,2*M), ones(1,M)]'
[eps_nu(1)*C(1), eps_nu(2)*C(2)]
end
disp('Parsing fmincon solution to SVM object');

learned_svm = [];
e=eye(N);

alphasub = x(1:M)-x(M+1:2*M);
sv_tol_alpha = max(abs(alphasub))*1e-5;
alpha_index = find(abs(alphasub)>sv_tol_alpha & abs(alphasub) > tol);       % find SVs
learned_svm.alpha = alphasub(alpha_index);  
learned_svm.y = ones(length(learned_svm.alpha),1);
learned_svm.Sva = dynamics_data(1:N,alpha_index);


learned_svm.beta = [];
learned_svm.Svb = [];
sv_tol_beta = max(abs(x(2*M+1:end)))*1e-5;
beta_index = find((x(2*M+1:end))>sv_tol_beta);
learned_svm.beta = x(2*M+beta_index);
learned_svm.Svb = dynamics_data(:,beta_index);

learned_svm.gamma = zeros(N,1);
learned_svm.lambda = lambda;
learned_svm.type = ktype;
learned_svm.target = zeros(N,1);

% Calculating b0
learned_svm.b0 = 0;
tmp = 0;
count=0;
for i=1:length(alpha_index)
  
    if(x(alpha_index(i)) < sv_tol_alpha && x(alpha_index(i)+M) < C(1))
        tmp = tmp + (y(alpha_index(i))+eps_nu(1) - calculateClassifier(learned_svm, dynamics_data(1:N, alpha_index(i))));
        count=count+1;
    else if(x(alpha_index(i)+M) < sv_tol_alpha && x(alpha_index(i)) < C(1))
       tmp = tmp + (y(alpha_index(i))-eps_nu(1) - calculateClassifier(learned_svm, dynamics_data(1:N, alpha_index(i))));
       count=count+1;
        end
    end
end
if(count ~=0)
learned_svm.b0 = tmp/count;
end


end

