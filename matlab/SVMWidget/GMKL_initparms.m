function parms=GMKL_initparms(NUMk,TYPEker,TYPEreg)
  % Kernel parameters
  switch (TYPEker),
  
    case 1, parms.KERname='Product of RBF kernels across features';
            parms.gamma=1;
	    parms.fncK=@KProdRBF;
	    parms.fncKdash=@KdashProdRBF;
	    parms.fncKdashint=@KdashProdRBFintermediate;
      case 2, parms.KERname='ASVM bsb-value bsb-gradient';
          parms.gamma=1;
          parms.fncK=@KASVM_bsb_value_bsb_gradient;
          parms.fncKdash=@KdashASVM_bsb_value_bsb_gradient;
          parms.fncKdashint=@KdashASVM_bsb_value_bsb_gradient_intermediate;
       case 3, parms.KERname='ASVM osb-value osb-gradient';
          parms.gamma=1;
          parms.fncK=@KASVM_osb_value_osb_gradient;
          parms.fncKdash=@KdashASVM_osb_value_osb_gradient;
          parms.fncKdashint=@KdashASVM_osb_value_osb_gradient_intermediate;
    otherwise, fprintf('Unknown kernel type.\n'); keyboard; 
  end;
      
  % Regularization parameters
  switch (TYPEreg),
      case 0, parms.REGname='l1';          % L1 Regularization
          parms.sigma=ones(NUMk,1);
          parms.fncR=@Rl1;
          parms.fncRdash=@Rdashl1;
      case 1, parms.REGname='l2';          % L2 Regularization
          parms.mud=ones(NUMk,1);
          parms.covd=1e-1*eye(NUMk);
          parms.invcovd=inv(parms.covd);
          parms.fncR=@Rl2;
          parms.fncRdash=@Rdashl2;
      otherwise, fprintf('Unknown regularisation type.\n'); keyboard;
  end;
    
  % Standard SVM parameters
 
  parms.C=1;                 % Misclassification penalty for SVC/SVR ('-c' in libSVM)
  
  % Gradient descent parameters
  parms.initd=rand(NUMk,1);   % Starting point for gradient descent
  parms.TYPEsolver='MAT';     % Use Matlab's quadprog as an SVM solver.
  parms.TYPEstep=1;           % 0 = Armijo, 1 = Variant Armijo, 2 = Hessian (not yet implemented)
  parms.MAXITER=40;           % Maximum number of gradient descent iterations
  parms.MAXEVAL=200;          % Maximum number of SVM evaluations
  parms.MAXSUBEVAL=20;        % Maximum number of SVM evaluations in any line search
  parms.SIGMATOL=0.3;         % Needed by Armijo, variant Armijo for line search
  parms.BETAUP=2.1;           % Needed by Armijo, variant Armijo for line search
  parms.BETADN=0.3;           % Needed by variant Armijo for line search
  parms.BOOLverbose=true;     % Print debug information at each iteration
  parms.SQRBETAUP=parms.BETAUP*parms.BETAUP; 
  parms.outer_tol = 1e-5;
end



function r=Rl1(parms,d), r=parms.sigma'*d; end
function rdash=Rdashl1(parms,d), rdash=parms.sigma; end

function r=Rl2(parms,d), delta=parms.mud-d;r=0.5*delta'*parms.invcovd*delta; end
function rdash=Rdashl2(parms,d), rdash=parms.invcovd*(d-parms.mud); end

function K=KASVM_bsb_value_bsb_gradient(TRNfeatures,TSTfeatures, parms,d)
NUMtrn=size(TRNfeatures,2);
NUMtst=size(TSTfeatures,2);


if (~isempty(find(d>1e-4,1))),
    K = mx_get_bsb_value_bsb_derivative_kernel(TRNfeatures, parms.gamma*d, 'rbf');
else
    K=zeros(NUMtrn,NUMtst);
end;

end


function Kdashint=KdashASVM_bsb_value_bsb_gradient_intermediate(TRNfeatures,TSTfeatures,parms,d)
Kdashint = mx_get_bsb_value_bsb_derivative_kernel_dsigma(TRNfeatures, parms.gamma*d, 'rbf');
end


function Kdash=KdashASVM_bsb_value_bsb_gradient(TRNfeatures,TSTfeatures,parms,d, l, Kdashint)
    Kdash = Kdashint{l};
end

function K=KASVM_osb_value_osb_gradient(TRNfeatures,TSTfeatures, parms,d)
NUMtrn=size(TRNfeatures,2);
NUMtst=size(TSTfeatures,2);

  
if (~isempty(find(d>1e-4,1)))
    K = mx_get_osb_value_osb_derivative_kernel(TRNfeatures, parms.target, parms.labels, parms.gamma*d, 'rbf');
else
    K=zeros(NUMtrn,NUMtst);
end
  
end

function Kdashint=KdashASVM_osb_value_osb_gradient_intermediate(TRNfeatures,TSTfeatures,parms,d)
Kdashint = mx_get_osb_value_osb_derivative_kernel_dsigma(TRNfeatures, parms.target, parms.labels, parms.gamma*d,'rbf');
end

function Kdash=KdashASVM_osb_value_osb_gradient(TRNfeatures,TSTfeatures,parms,d, l, Kdashint)
Kdash = Kdashint{l};
end


function K=KProdRBF(TRNfeatures,TSTfeatures,parms,d)
  [NUMdims,NUMtrn]=size(TRNfeatures);
  [NUMdims,NUMtst]=size(TSTfeatures);
  if (NUMdims~=length(d)), fprintf('NUMdims not equal to NUMk\n'); keyboard; end;
  
  nzind=find(d>1e-4);
  K=zeros(NUMtrn,NUMtst);
  if (~isempty(nzind)),
    for i=1:length(nzind),
      k=nzind(i);
      K=K+d(k)*(repmat(TRNfeatures(k,:)',1,NUMtst)-repmat(TSTfeatures(k,:),NUMtrn,1)).^2;
    end;
  end;
  K=exp(-parms.gamma*K);
end

function Kdashint=KdashProdRBFintermediate(TRNfeatures,TSTfeatures,parms,d)
  Kdashint=-parms.gamma*parms.fncK(TRNfeatures,TSTfeatures,parms,d);
end

function Kdash=KdashProdRBF(TRNfeatures,TSTfeatures,parms,d,k,Kdashint)
  [NUMdims,NUMtrn]=size(TRNfeatures);
  [NUMdims,NUMtst]=size(TSTfeatures);

  Kdash=1/(parms.C)*Kdashint.*(repmat(TRNfeatures(k,:)',1,NUMtst)-repmat(TSTfeatures(k,:),NUMtrn,1)).^2;
end

