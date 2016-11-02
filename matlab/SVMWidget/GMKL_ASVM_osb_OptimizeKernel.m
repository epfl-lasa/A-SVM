function learned_svm = GMKL_ASVM_osb_OptimizeKernel(data_set, labels, target, R, gamma_init)


M=size(data_set,2);
N=size(data_set,1)/2;
P=length(find(labels==1));
tol = 1e-5;
TYPEreg=0;
NUMkernels = N;
TYPEker = 3;

parms = GMKL_initparms(NUMkernels, TYPEker, TYPEreg);
parms.initd=gamma_init;
parms.TYPEsolver='ASVM-osb';
parms.outer_tol=1e-3;
parms.R = R;
parms.C = 1;
parms.target = target;
tic
svm=GMKL_optimize(data_set,labels,parms);
toc
gamma_opt = svm.d;
alpha_opt = svm.x;
x=svm.x;

sv_tol_alpha = max(x(1:M+P))*tol;
sv_tol_beta = sv_tol_alpha;

positive_labels = find(labels == 1);
alpha_index = find(abs(x(1:M))>sv_tol_alpha);       % find SVs
beta_index = positive_labels(find(x(M+1:M+P) > sv_tol_beta));

learned_svm = [];
learned_svm.alpha = x(alpha_index);  
learned_svm.beta = x(M+find(x(M+1:M+P) > sv_tol_beta));
learned_svm.gamma = x(M+P+1:end);
learned_svm.y = labels(alpha_index);
learned_svm.Sva = data_set(1:N,alpha_index);
learned_svm.Svb = data_set(:,beta_index);
learned_svm.lambda = svm.d;
learned_svm.type = 'rbf';
learned_svm.target = parms.target;
% Calculating b0
learned_svm.b0 = 0;
f_svp = [];
f_svn = [];
n_p = 0;
n_n = 0;
for i=1:M
    tmp = 0;
    pt = data_set(1:N,i);
   tmp = calculateClassifier(learned_svm, pt);
    
    if(labels(i) == 1)
       n_p = n_p + 1;
       f_svp = [f_svp;tmp];
    else
        n_n = n_n + 1;
        f_svn = [f_svn;tmp];
    end
end

learned_svm.b0 = -(max(f_svn) + min(f_svp))/2;