% function [ x, learned_svm, exitflag ] = classification_dynamics_asvm( all_data, target_class, lambda, C, tol, max_iter )
%%

clear
use_alpha=1;
load data/examples/mat/ex5.mat 
sig = 0.5;
all_data = saved_data;
target_class = 1;
lambda = 1/(2*sig*sig);
C=1000;
tol=1e-6;
max_iter = 1000;

ktype = 'rbf';
[dynamics_data, data_labels, ind, tl] = ipopt_preprocess_data_unnormalized(all_data, target_class);
dynamics_data = [dynamics_data, [tl(:,target_class);zeros(size(dynamics_data,1)/2,1)]];
data_labels = [data_labels;1];
M = size(dynamics_data,2);
N = size(dynamics_data,1)/2;
P = length(find(data_labels == 1));
positive_labels = find(data_labels == 1);

dynamics_data(N+1:end,:) = dynamics_data(N+1:end,:)/0.005;
% Set up the kernel matrix
disp('Calculating non-linear kernel...');
K=zeros(M,M);
G=zeros(M,N*P);
H = zeros(N*P, N*P);
e=eye(N);
for i=1:M
    yi = data_labels(i);
    xi = dynamics_data(1:N,i);
    for j=1:M
        yj = data_labels(j);
        xj = dynamics_data(1:N,j);
        K(i,j) = yi*yj*getKernel(xi, xj, lambda, ktype);
    end
end

for b1=1:N
    for i=1:M
        yi = data_labels(i);
        xi = dynamics_data(1:N,i);
        for j=1:P
            xj = dynamics_data(1:N,positive_labels(j));
            G(i,j+(b1-1)*P) = yi*getKernelFirstDerivative(xi, xj, lambda, ktype, 2)'*e(:,b1);
        end
    end
end

for b1=1:N
    for b2=1:N
        for i=1:P
            xi = dynamics_data(1:N,positive_labels(i));
            for j=1:P
                xj = dynamics_data(1:N,positive_labels(j));
                H(i+(b1-1)*P,j+(b2-1)*P) = e(:,b1)'*getKernelSecondDerivative(xi, xj, lambda, ktype)*e(:,b2);
            end
        end
    end
end
if(use_alpha)
Q = [K -G; -G' H];
else
Q=H;
end

[ee,vv] =eig(Q+Q');
disp(['Minimum eigen value of Q : ' num2str(min(vv(find(vv<0))))]);


Q = Q+1e-11*eye(size(Q));

if(use_alpha)
x0 = zeros(M+2*N*P,1);
else
x0 = zeros(2*N*P,1);
end

auxdata_1 = Q;
auxdata_2 = lambda;
auxdata_3 = dynamics_data;
auxdata_4 = M;
auxdata_5 = P;
auxdata_6 = N;
auxdata_7 = find(data_labels == 1);
auxdata_8=use_alpha;

% FMINCON options
if(use_alpha)
lb = zeros(M+2*N*P,1);
ub = C*ones(M+2*N*P,1);
else
lb = zeros(2*N*P,1);
ub = C*ones(2*N*P,1);
% for i=1:N
%    ub(i*P)=inf; 
% end
end

obj = @(x)classification_dynamics_objective(x, auxdata_4, auxdata_1, auxdata_2, auxdata_5, auxdata_6, auxdata_3, auxdata_7, auxdata_8);
options=[];
options=optimset('TolX',tol,'TolCon',tol,'Algorithm','sqp','Display','iter','TolFun',tol,'MaxIter',max_iter,'OutputFcn',@fmincon_callback,'MaxFunEvals',inf);

% Run FMINCON.
tic

% Only C-SV
if(use_alpha)
[x, fval, exitflag] = fmincon(obj, x0, [], [], [data_labels;zeros(2*N*P,1)]', 0, lb, ub, [], options);
else
[x, fval, exitflag] = fmincon(obj, x0, [], [], [], [], lb, ub, [], options);
end

if~(use_alpha)
x=[zeros(M,1);x];
end

toc


disp('Parsing fmincon solution to SVM object');

learned_svm = [];
e=eye(N);

sv_tol_alpha = max(x(1:M))*1e-5;
alpha_index = find(abs(x(1:M))>sv_tol_alpha);       % find SVs
learned_svm.alpha = x(alpha_index);  
learned_svm.y = data_labels(alpha_index);
learned_svm.Sva = dynamics_data(1:N,alpha_index);


learned_svm.beta = [];
learned_svm.Svb = [];


for i=1:N
    betasub = x(M+(i-1)*P + 1:M+i*P)-x(M+N*P+(i-1)*P+1:M+N*P+i*P);
    sv_tol_beta = max(abs(betasub))*1e-5; 
    beta_index = find(abs(betasub)  > sv_tol_beta & abs(betasub) > tol);
    beta = betasub(beta_index);
    Svb = dynamics_data(1:N,positive_labels(beta_index));
    
    learned_svm.beta = [learned_svm.beta;-beta(:)];
    learned_svm.Svb = [learned_svm.Svb, [Svb;repmat(e(:,i),1,size(Svb,2))]];
end

learned_svm.gamma = [];
learned_svm.lambda = lambda;
learned_svm.type = ktype;
learned_svm.target = tl(:,target_class);

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
learned_svm.b0 = -(max(f_svn) + min(f_svp))/2

%%
errs=[];
for i=1:length(dynamics_data)
if(data_labels(i) > 0)
errs=[errs,calculateClassifierDerivative(learned_svm,dynamics_data(1:N,i)) - dynamics_data(N+1:2*N,i)];
end
end
disp(['Velocity Error: [' num2str(mean(errs')) '] +- [' num2str(std(errs')) ']']);

figure(1)
cla; hold on;
hh=plot(dynamics_data(1,:)',dynamics_data(2,:)','k.');
axis tight
lims=[get(gca,'XLim'), get(gca,'YLim')];
lims = lims - [(lims(2)-lims(1))*0.2 -(lims(2)-lims(1))*0.2 (lims(4)-lims(3))*0.2 -(lims(4)-lims(3))*0.2];
[tx,ty]=meshgrid(linspace(lims(1), lims(2), 20),linspace(lims(3), lims(4), 20));
vx=[];vy=[];
hx=[];
for i=1:size(tx,1)
    for j=1:size(tx,2)
        pt=[tx(i,j);ty(i,j)];
        tmp = calculateClassifierDerivative(learned_svm,pt);
        hx(i,j) = calculateClassifier(learned_svm, pt);
        vx(i,j) = tmp(1);
        vy(i,j) = tmp(2);
    end
end
delete(hh);

displaySVM(learned_svm, 'range',lims, 'sv',1);

plot(dynamics_data(1,:)',dynamics_data(2,:)','k.');
streamslice(tx, ty, vx, vy, 1.5);
axis tight
colormap autumn
shading interp



figure(2)
subplot(1,2,1); cla;hold on; grid on;
plot3(dynamics_data(1,positive_labels)', dynamics_data(2,positive_labels)',dynamics_data(3,positive_labels)','k.');

lims=[get(gca,'XLim'), get(gca,'YLim')];
[tx,ty]=meshgrid(linspace(lims(1), lims(2), 20),linspace(lims(3), lims(4), 20));
vx=[];
for i=1:size(tx,1)
    for j=1:size(tx,2)
        pt=[tx(i,j);ty(i,j)];
        tmp = calculateClassifierDerivative(learned_svm,pt);
        vx(i,j) = tmp(1);
    end
end
surf(tx, ty, vx);
xlabel('X');ylabel('Y');zlabel('Xdot');
subplot(1,2,2); cla;hold on; grid on;
plot3(dynamics_data(1,positive_labels)', dynamics_data(2,positive_labels)',dynamics_data(4,positive_labels)','k.');
lims=[get(gca,'XLim'), get(gca,'YLim')];
[tx,ty]=meshgrid(linspace(lims(1), lims(2), 20),linspace(lims(3), lims(4), 20));
vx=[];
for i=1:size(tx,1)
    for j=1:size(tx,2)
        pt=[tx(i,j);ty(i,j)];
        tmp = calculateClassifierDerivative(learned_svm,pt);
        vx(i,j) = tmp(2);
    end
end
surf(tx, ty, vx);
xlabel('X');ylabel('Y');zlabel('Ydot');



% figure(3)
% cla
% plot(dynamics_data(1,positive_labels)', dynamics_data(2,positive_labels)','k.');








