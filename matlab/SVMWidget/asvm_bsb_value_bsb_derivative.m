function learned_svm= asvm_bsb_value_bsb_derivative( dynamics_data, y, Q, lambda, ktype,  C, tol, eps_nu, use_nu )
%LEARN_SURFACE_ASVM Summary of this function goes here
%   Detailed explanation goes here

use_derivatives = 1;
max_iter = 1000;

if(isscalar(C))
   C=[C;C]; 
end

if(isscalar(eps_nu))
   eps_nu=[eps_nu;eps_nu]; 
end
M = size(dynamics_data,2);
N = size(dynamics_data,1)/2;

y = reshape(y,M,1);
tmp = [];
for i=1:N
   tmp = [tmp,dynamics_data(N+i,:),-dynamics_data(N+i,:)]; 
end
if(use_nu)
    f=[y;-y; zeros(2*N*M,1)] + [zeros(2*M,1);tmp'];
else
f=eps_nu(2)*[zeros(2*M,1);ones(2*N*M,1)] + [y;-y; zeros(2*N*M,1)] + eps_nu(1)*[ones(2*M,1);...
    zeros(2*N*M,1)] + [zeros(2*M,1);tmp'];
end


% x0 = zeros(2*M+2*N*M,1);


% FMINCON options
lb = zeros(2*M+2*N*M,1);
if(use_derivatives)
    if(use_nu)
        ub = [C(1)*ones(2*M,1);C(2)*ones(2*N*M,1)];
    else
        ub = [C(1)*ones(2*M,1);C(2)*ones(2*N*M,1)];
    end
else
    if(use_nu)
        ub = [C(1)*ones(2*M,1);zeros(2*N*M,1)];
    else
        ub = [C(1)*ones(2*M,1);zeros(2*N*M,1)];
    end
end

options=optimset('TolX',tol,'TolCon',tol,'Algorithm','interior-point-convex','Display','iter','TolFun',tol,'MaxIter',max_iter,'MaxFunEvals',inf);

tic
if(use_nu)
     [x, fval, exitflag] = quadprog((Q+Q')/2, f, [ones(1,2*M),zeros(1,2*N*M);zeros(1,2*M), ones(1,2*N*M)],...
         [C(1)*eps_nu(1);C(2)*eps_nu(2)],[],[],lb,ub, [], options);
else
    [x, fval, exitflag] = quadprog((Q+Q')/2, f, [],[],[],[],lb,ub, [], options);
%         [x, fval, exitflag] = mx_quadprog(Q+1e-4*eye(2*M+2*N*M), f, lb, ub, tol);

end
toc

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


for i=1:N
    betasub = x(2*M*i+1:2*M*i+M)-x(2*M*i+M+1:2*M*(i+1));
    sv_tol_beta = max(abs(betasub))*1e-5; 
    beta_index = find(abs(betasub)  > sv_tol_beta & abs(betasub) > tol);
    beta = betasub(beta_index);
    Svb = dynamics_data(1:N,beta_index);
    
    learned_svm.beta = [learned_svm.beta;-beta(:)];
    learned_svm.Svb = [learned_svm.Svb, [Svb;repmat(e(:,i),1,size(Svb,2))]];

end

learned_svm.gamma = zeros(N,1);
learned_svm.lambda = lambda;
learned_svm.type = ktype;
learned_svm.target = zeros(N,1);

% Calculating b0
learned_svm.b0 = 0;
tmp = 0;
count=0;
indp=[];
indn=[];
for i=1:length(alpha_index)
  
    if(x(alpha_index(i)) < sv_tol_alpha && x(alpha_index(i)+M) < C(1))
        indp = [indp,i];
%         tmp = tmp + (y(alpha_index(i))+eps_nu(1) - calculateClassifier(learned_svm, dynamics_data(1:N, alpha_index(i))));
        count=count+1;
    else if(x(alpha_index(i)+M) < sv_tol_alpha && x(alpha_index(i)) < C(1))
%        tmp = tmp + (y(alpha_index(i))-eps_nu(1) - calculateClassifier(learned_svm, dynamics_data(1:N, alpha_index(i))));
        indn=[indn,i];
       count=count+1;
        end
    end
end
tmp = eps_nu(1) + sum(y(alpha_index(indp)) - mx_calculate_classifier(learned_svm, dynamics_data(1:N, alpha_index(indp))));
tmp = tmp - eps_nu(1) + sum(y(alpha_index(indn)) - mx_calculate_classifier(learned_svm, dynamics_data(1:N, alpha_index(indn))));

if(count ~=0)
learned_svm.b0 = tmp/count;
end

if(use_nu)
disp(num2str([length(learned_svm.alpha), M*eps_nu(1)]));
disp(num2str([length(learned_svm.beta), M*eps_nu(2)]));
end

end

