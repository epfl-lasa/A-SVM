function learned_svm  = smo_one_vs_all( all_data, target_class, lambda, C, tol )
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
global x_smo dlabels err_cache_alpha err_cache_beta err_cache_gamma tl Cparam Bparam ...
    ker_matrix M minimum_alpha maximum_alpha P N maximum_beta maximum_gamma H beta_tol d beta_relaxation

 [ dynamics_data, data_labels, indices, target_list ] = ipopt_preprocess_data( all_data,...
                                                        target_class );
target = target_list(:,target_class);
dlabels = data_labels;
tl = tol;
Cparam = C;
Bparam = 0;
beta_tol = .1;
d = lambda;
beta_relaxation = 0*pi/180;

positive_labels = find(data_labels == 1);
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
    

Q = Q+1e-08*eye(size(Q));
[ee,vv] =eig(Q+Q');
v=[];for i=1:size(vv,1)
v=[v;vv(i,i)];
end
disp(['minimum eigen value of Q : ' num2str(min(v))]);



ker_matrix = Q;
for i=1:M
    for j=1:M
        ker_matrix(i,j) = ker_matrix(i,j)*dlabels(i)*dlabels(j);
    end
     for j=M+1:M+P+N
        ker_matrix(i,j) = ker_matrix(i,j)*dlabels(i);
        ker_matrix(j,i)=ker_matrix(i,j);
    end
end


disp('Warm starting using libsvm (classifier-only) solution...');
stroptions = ['-s 0 -t 2 -g ' num2str(lambda) ' -c 1e6 -e 1e-4 '];
model = svmtrain(data_labels, dynamics_data(1:N,:)', stroptions);
x_smo = [zeros(M,1); zeros(P,1); zeros(N,1)];   % remove this when solving full problem
x_smo(model.SV_indices+1) = abs(model.sv_coef);               % check if this even is a valid classifier
Bparam = model.rho;
disp(' ');

% disp('Cold starting ...');
% x_smo = [ones(M,1); ones(P,1); ones(N,1)];
% ind = find(data_labels > 0);
% x_smo(ind)  = 1/P;
% ind = find(data_labels < 0);
% x_smo(ind)  = 1/(M-P);

err_cache_alpha = zeros(M,1);
for i=1:M
    err_cache_alpha(i) =  forward_alpha(i) - dlabels(i);
end
err_cache_beta = zeros(P,1);
for i=1:P
    err_cache_beta(i) =  forward_beta(i+M);
end
err_cache_gamma = zeros(N,1);
for i=1:N
    err_cache_gamma(i) = forward_gamma(i+M+P);
end

[tmp, minimum_alpha] = min(err_cache_alpha);
[tmp, maximum_alpha] = max(err_cache_alpha);
[tmp, maximum_beta] = max(err_cache_beta);
[tmp, maximum_gamma] = max(err_cache_gamma);

% updateB0();
%% Solve SMO here starting with only x_smo. Optimal sol as x
disp('Starting SMO');
tic
doSMO();
toc

disp('Parsing SMO solution to SVM object');


alpha_index = find(x_smo(1:M)>tl);       % find SVs
beta_index = positive_labels(find(x_smo(M+1:M+P) > beta_tol));

learned_svm = [];
learned_svm.alpha = x_smo(alpha_index);
learned_svm.beta = x_smo(M+find(x_smo(M+1:M+P) > beta_tol));
learned_svm.gamma = x_smo(M+P+1:end);
learned_svm.y = data_labels(alpha_index);
learned_svm.Sva = dynamics_data(1:N,alpha_index);
learned_svm.Svb = dynamics_data(:,beta_index);
learned_svm.lambda = lambda;
learned_svm.type = ktype;
learned_svm.target = target;
learned_svm.b0 = -Bparam;

disp(['Alpha err     : ' num2str(1-abs(calculateClassifier(learned_svm, learned_svm.Sva(:,1))))]);
disp(['Beta err      : ' num2str(abs(sum(min(forward_beta(M+1:M+P),0))))]);
disp(['Gamma err     : ' num2str(forward_gamma(M+P+1:M+P+N)')]);
disp(['sum(alpha.*y) : ' num2str(sum(x_smo(1:M).*data_labels))]);


function isChanged = examineExampleForAlpha(i2)
global x_smo dlabels err_cache_alpha tl Cparam minimum_alpha maximum_alpha M

isChanged = 0;
y2 = dlabels(i2);
alph2 = x_smo(i2);
if alph2 > 0.0 && alph2 < Cparam
    E2 = err_cache_alpha(i2);
else
    E2 = forward_alpha(i2) - y2;
end
r2 = E2*y2;

if( (r2 < -tl && alph2 < Cparam) || (r2 > tl && alph2 > 0))
    if minimum_alpha > -1
        if(abs(E2-err_cache_alpha(minimum_alpha)) > abs(E2 - err_cache_alpha(maximum_alpha)))
            i1 = minimum_alpha;
        else
            i1 = maximum_alpha;
        end

        if takeStepForAlpha(i1, i2, E2)
            isChanged = 1;
            return;
        end
    end

    active_alpha_index = find(x_smo(1:M) > 0 & x_smo(1:M) < Cparam);
    for j=1:length(active_alpha_index)
        if takeStepForAlpha(active_alpha_index(j), i2, E2)
            isChanged = 1;
            return;
        end
    end

    for j = 1:M
        if x_smo(j) == 0 || x_smo(j) == Cparam
            if takeStepForAlpha(j, i2, E2)
                isChanged = 1;
                return;
            end
        end
    end
    
end


function isStepTaken = takeStepForAlpha(i1, i2, E2)
global x_smo dlabels err_cache_alpha tl P Cparam ai ker_matrix M ...
    minimum_alpha maximum_alpha err_cache_beta

isStepTaken = 1;
ai = i2;
if(i1 == i2)
    isStepTaken = 0;
    return;
end

alph1 = x_smo(i1);
alph2 = x_smo(i2);
y1 = dlabels(i1);
y2 = dlabels(i2);
if alph1 > 0.0 && alph1 < Cparam
    E1 = err_cache_alpha(i1);
else
     E1 = forward_alpha(i1) - y1;
end


s = y1*y2;
if(y1 ~= y2)
    L = max(0, alph2-alph1);
    H = min(Cparam, Cparam+alph2-alph1);
else
    L = max(0, alph2+alph1-Cparam);
    H = min(Cparam, alph2+alph1);
end

if(abs(L-H) < tl)
    isStepTaken = 0;
    return;
end

k11 = ker_matrix(i1, i1);
k12 = ker_matrix(i1, i2);
k22 = ker_matrix(i2, i2);
eta = k11 + k22 - 2*k12;

if(eta > 0)
    a2 = alph2 + y2*(E1-E2)/eta;
    if(a2 < L)
        a2 = L ;
    else if(a2 > H)
            a2 = H;
        end
    end
else
    disp('Eta < 0 !!');
    isStepTaken = 0;
    return;
end


if(abs(a2-alph2) < (a2+alph2+tl)*tl)      % check this
    isStepTaken = 0;
    return;
end

a1 = alph1 + s*(alph2-a2);
if a1 < tl
    a1=0;
end

x_smo (i1) = a1;x_smo (i2) = a2;

w1 = y1*(a1-alph1);
w2 = y2*(a2-alph2);

if(a1 > 0 && a1 < Cparam)
    err_cache_alpha(i1) = forward_alpha(i1) - y1;
end
if(a2 > 0 && a2 < Cparam)
    err_cache_alpha(i2) = forward_alpha(i2) - y2;
end

if (err_cache_alpha(i1) > err_cache_alpha(i2))
    minimum_alpha = i2;
    maximum_alpha = i1;
end
if (err_cache_alpha(i1) < err_cache_alpha(i2))
    minimum_alpha = i1;
    maximum_alpha = i2;
end

for i=1:M
    if i~=i1 && i~=i2 && x_smo(i) >0 && x_smo(i) < Cparam
            err_cache_alpha(i) =  err_cache_alpha(i) + w1*ker_matrix(i1, i) + w2*ker_matrix(i2, i);
        if (err_cache_alpha(i)  > err_cache_alpha(maximum_alpha))
            maximum_alpha = i;
        end
        
        if (err_cache_alpha(i) < err_cache_alpha(minimum_alpha))
            minimum_alpha = i;
        end

    end
    
end

for i=M+1:M+P
    if x_smo(i) > 0 && x_smo(i) < Cparam
        err_cache_beta(i-M) = err_cache_beta(i-M) + w1*ker_matrix(i, i1) + w2*ker_matrix(i, i2);
    end
end

function isChanged = examineExampleForBeta(i1)
global x_smo err_cache_beta M Cparam beta_tol
isChanged = 0;


beta1 = x_smo(i1);
if beta1 > 0 && beta1 < Cparam
   E1 = err_cache_beta(i1-M);
else
   E1 = forward_beta(i1);
end

if (E1 < -beta_tol  && x_smo(i1) < Cparam) || (E1 > beta_tol && x_smo(i1) > 0)  % if this beta is in violation
   if takeStepForBeta(i1, E1)
      isChanged = 1;
      return;
   end        
end


function isStepTaken = takeStepForBeta(i1, E1)
global beta_tol x_smo Cparam H err_cache_beta M P ker_matrix err_cache_alpha ...
     maximum_alpha minimum_alpha  beta_relaxation

isStepTaken = 1;
beta1 = x_smo(i1);

% if E1 > -beta_relaxation && E1 < 0
%     beta_new = 0;
% %         isStepTaken = 0;
% %     return;
% else
quad_term = H(i1-M, i1-M);
if(quad_term > 0)
    beta_new = beta1 - E1/quad_term;
else
    disp(['H(i,i) = ' num2str(quad_term) ' for i=' int2str(i1-M) '!! Expected positive!!']);
    isStepTaken = 0;
    return;
end
if beta_new < 0
    beta_new = 0;
else if beta_new > Cparam
        beta_new = Cparam;
    end
end
% end

bdiff = beta_new - beta1;
if(abs(bdiff) < (beta_new+beta1+beta_tol)*beta_tol)
    isStepTaken = 0;
    return;
end
x_smo(i1) = beta_new;

if beta_new > 0 && beta_new <Cparam
    err_cache_beta(i1-M) = forward_beta(i1);
end


for i=M+1:M+P
    if i ~= i1 && x_smo(i) > 0 && x_smo(i) < Cparam
        err_cache_beta(i-M) = err_cache_beta(i-M) +  ker_matrix(i1, i)*bdiff;
    end
end

for i=1:M
    if x_smo(i) > 0 && x_smo(i) < Cparam
        err_cache_alpha(i) = err_cache_alpha(i) +  ker_matrix(i, i1)*bdiff;
        
        if (err_cache_alpha(i)  > err_cache_alpha(maximum_alpha))
            maximum_alpha = i;
        end

        if (err_cache_alpha(i) < err_cache_alpha(minimum_alpha))
            minimum_alpha = i;
        end
    end
    
end


function isChanged = examineExampleForGamma(i1)
global  tl
isChanged = 0;

E1 = forward_gamma(i1);
if  abs(E1) > tl
    if takeStepForGamma(i1, E1)
        isChanged = 1;
    end
end


function isStepTaken = takeStepForGamma(i1, E1)
global x_smo maximum_alpha minimum_alpha err_cache_alpha err_cache_beta...
    ker_matrix Cparam d M P 

isStepTaken = 1;
g1 = x_smo(i1);
gamma_new = g1-E1/(2*d);

gdiff = gamma_new - g1;
if abs(gdiff) < 0.0001 || abs(gamma_new) < 0.001
    isStepTaken = 0;
    return;
end

x_smo(i1) = gamma_new;

% update err cache for alpha 
for i=1:M
    if x_smo(i) > 0 && x_smo(i) < Cparam
        err_cache_alpha(i) = err_cache_alpha(i) +  ker_matrix(i, i1)*gdiff;
        
        if (err_cache_alpha(i)  > err_cache_alpha(maximum_alpha))
            maximum_alpha = i;
        end

        if (err_cache_alpha(i) < err_cache_alpha(minimum_alpha))
            minimum_alpha = i;
        end
    end
    
end

% update err cache for beta
for i=M+1:M+P
    if x_smo(i) > 0 && x_smo(i) < Cparam
        err_cache_beta(i-M) = err_cache_beta(i-M) +  ker_matrix(i, i1)*gdiff;
    end
end


function fval = forward_alpha(ind)
global ker_matrix M dlabels x_smo Bparam P N 

fval=-Bparam;
for i=1:M
    if x_smo(i) > 0
        fval = fval + ker_matrix(ind, i)*dlabels(i)*x_smo(i);
    end
end
for i=M+1:M+P
    if x_smo(i) > 0
    fval = fval + ker_matrix(ind, i)*x_smo(i);
    end
end
for i=M+P+1:M+P+N
    fval = fval + ker_matrix(ind, i)*x_smo(i);
end

function fval = forward_beta(ind)
global ker_matrix M dlabels x_smo  P N 

fval=0;
for i=1:M
    if x_smo(i) > 0
    fval = fval + ker_matrix(ind, i)*dlabels(i)*x_smo(i);
    end
end
for i=M+1:M+P
    if x_smo(i) > 0
    fval = fval + ker_matrix(ind, i)*x_smo(i);
    end
end
for i=M+P+1:M+P+N
    fval = fval + ker_matrix(ind, i)*x_smo(i);
end

function fval = forward_gamma(ind)
global ker_matrix M dlabels x_smo  P N
fval=0;
for i=1:M
    if x_smo(i) > 0
    fval = fval + ker_matrix(ind, i)*dlabels(i)*x_smo(i);
    end
end

for i=M+1:M+P
    if x_smo(i) > 0
    fval = fval + ker_matrix(ind, i)*x_smo(i);
    end
end


for i=M+P+1:M+P+N
    fval = fval + ker_matrix(ind, i)*x_smo(i);
end





function updateB0()
global Bparam Cparam maximum_alpha minimum_alpha err_cache_alpha M  dlabels x_smo
b_old = Bparam;

ind = find(x_smo(1:M) >0  & x_smo(1:M) < Cparam);
fnc=forward_alpha(ind) + Bparam;
Bparam = mean(-dlabels(ind) + fnc);

% update err cache for alpha 
for i=1:M
    if x_smo(i) > 0 && x_smo(i) < Cparam
        err_cache_alpha(i) = err_cache_alpha(i)  + b_old - Bparam;
        
        if (err_cache_alpha(i)  > err_cache_alpha(maximum_alpha))
            maximum_alpha = i;
        end

        if (err_cache_alpha(i) < err_cache_alpha(minimum_alpha))
            minimum_alpha = i;
        end
    end
    
end

function doSMO()
global M P N x_smo Cparam
numChanged = 0;
examineAll = 1;

while (numChanged > 0 || examineAll)
    
    numChanged = 0;
    if (examineAll)
        for i = 1:M
                numChanged = numChanged + examineExampleForAlpha(i);
        end
        
        for i = M+1:M+P
                numChanged = numChanged + examineExampleForBeta(i);
        end
        for i = M+P+1:M+P+N
                numChanged = numChanged + examineExampleForGamma(i);
        end

        
        examineAll = 0;
    else
        for i = 1:M
            if x_smo(i) > 0 && x_smo(i) < Cparam
                numChanged = numChanged + examineExampleForAlpha(i);
            end
        end

        for i = M+1:M+P
            if x_smo(i) > 0 && x_smo(i) < Cparam
                numChanged = numChanged + examineExampleForBeta(i);
            end
        end

                for i = M+P+1:M+P+N
                numChanged = numChanged + examineExampleForGamma(i);
                end
        
        if numChanged == 0
            examineAll = 1;
        end
        

    end
    updateB0();

end
