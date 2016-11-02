function learned_svm = GMKL_SVC_OptimizeKernel(data_set, labels, gamma_init)

TYPEreg=0;
NUMkernels=size(data_set,1);
% The following code shows how to call COMPGDoptimize with
% feature vectors as input
TYPEker=1;
parms=GMKL_initparms(NUMkernels,TYPEker,TYPEreg);
% parms.MAXITER = 500;
% parms.C=10;
% parms.gamma = [1/(1)];
parms.TYPEsolver='LIB';
parms.outer_tol=1e-3;
parms.initd = gamma_init;
parms.C=1;
% parms.sigma=parms.C*ones(NUMkernels,1);
tic
svm=GMKL_optimize(data_set,labels,parms);
toc
tol=1e-5;
alpha_index = find(abs(svm.x)>tol);
learned_svm=[];
learned_svm.alpha = abs(svm.x(alpha_index));
learned_svm.y=sign(labels(alpha_index));
learned_svm.Sva = data_set(:,alpha_index);
learned_svm.beta=[];
learned_svm.Svb=[];
learned_svm.gamma=zeros(size(data_set,1),1);
learned_svm.target=zeros(size(data_set,1),1);
learned_svm.type='rbf';
learned_svm.b0 = svm.b;
learned_svm.lambda = svm.d;
