function learned_svm = GMKL_ASVM_bsb_OptimizeKernel(data_set, eps_value, eps_derivative, gamma_init)

data_set = [zeros(size(data_set,1),1),  data_set];
labels =  [0;ones(size(data_set,2)-1,1)];
tol = 1e-5;
TYPEreg=0;
NUMkernels=size(data_set,1)/2;
TYPEker=2;

parms=GMKL_initparms(NUMkernels,TYPEker,TYPEreg);
parms.initd = gamma_init(:);
parms.TYPEsolver='ASVM-bsb';
parms.outer_tol=1e-3;
parms.EPSv = eps_value;
parms.EPSd = eps_derivative;
parms.Cv = 1;
parms.Cd = 1;

tic
svm=GMKL_optimize(data_set,labels,parms);
toc
x=svm.x;
N = size(data_set,1)/2;
M = size(data_set,2);

learned_svm = [];
e=eye(N);

alphasub = x(1:M)-x(M+1:2*M);
sv_tol_alpha = max(abs(alphasub))*tol;
alpha_index = find(abs(alphasub)>sv_tol_alpha & abs(alphasub) > tol);       % find SVs
learned_svm.alpha = alphasub(alpha_index);  
learned_svm.y = ones(length(learned_svm.alpha),1);
learned_svm.Sva = data_set(1:N,alpha_index);
learned_svm.alpha_index=alpha_index;

learned_svm.beta = [];
learned_svm.Svb = [];
beta_index=[];
for i=1:N
    betasub = x(2*M*i+1:2*M*i+M)-x(2*M*i+M+1:2*M*(i+1));
    sv_tol_beta = max(abs(betasub))*tol; 
    tmp_beta_index = find(abs(betasub)  > sv_tol_beta & abs(betasub) > tol);
    beta_index = [beta_index;tmp_beta_index];
    beta = betasub(tmp_beta_index);
    Svb = data_set(1:N,tmp_beta_index);
    
    learned_svm.beta = [learned_svm.beta;-beta(:)];
    learned_svm.Svb = [learned_svm.Svb, [Svb;repmat(e(:,i),1,size(Svb,2))]];

end
learned_svm.beta_index=beta_index;
learned_svm.gamma = zeros(N,1);
learned_svm.lambda = svm.d;
learned_svm.type = 'rbf';
learned_svm.target = zeros(N,1);

% Calculating b0
learned_svm.b0 = 0;
tmp = 0;
count=0;
for i=1:length(alpha_index)
  
    if(x(alpha_index(i)) < sv_tol_alpha && x(alpha_index(i)+M) < parms.Cv)
        tmp = tmp + (labels(alpha_index(i))+parms.EPSv - calculateClassifier(learned_svm, data_set(1:N, alpha_index(i))));
        count=count+1;
    else if(x(alpha_index(i)+M) < sv_tol_alpha && x(alpha_index(i)) < parms.Cv)
       tmp = tmp + (labels(alpha_index(i))-parms.EPSv - calculateClassifier(learned_svm, data_set(1:N, alpha_index(i))));
       count=count+1;
        end
    end
end
if(count ~=0)
learned_svm.b0 = tmp/count;
end
learned_svm.b0=learned_svm.b0-1;
