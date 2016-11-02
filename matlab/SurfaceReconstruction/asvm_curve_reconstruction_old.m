% function [ x, learned_svm, exitflag ] = classification_dynamics_asvm( all_data, target_class, lambda, C, tol, max_iter )
%%

clear
close all
figure
axis([-1 1 -1 1]);
rawdata = grabDataFromCursor(1, 0.01);
%%
figure
axis([-1 1 -1 1]);
cla
hold on
posdata=[];
all_dynamics_data=[];
for i=1:length(rawdata)
    veldata =  diff(rawdata{i},1,2);
    vdata=[];
    for j=1:length(veldata)
        if(norm(veldata(:,j)) > 1e-10)
        normv = veldata(:,j)/norm(veldata(:,j));
         all_dynamics_data = [all_dynamics_data,[rawdata{i}(:,j);-normv(2);normv(1)]];
        end
    end
    plot(rawdata{i}(1,:),rawdata{i}(2,:),'k-*');
end


for i=1:size(all_dynamics_data,2)
   line([all_dynamics_data(1,i);all_dynamics_data(1,i) + 0.1*all_dynamics_data(3,i)] ,[all_dynamics_data(2,i);all_dynamics_data(2,i) + 0.1*all_dynamics_data(4,i)]);
end

%%
dynamics_data = all_dynamics_data(:,1:5:end);
use_derivatives =1;
sig = 0.5;
lambda = 1/(2*sig*sig);
C=10000;
tol=1e-9;
max_iter = 1000;
eps_p = 0.01;
eps_d=0.2;

ktype = 'rbf';

M = size(dynamics_data,2);
N = size(dynamics_data,1)/2;

% Set up the kernel matrix
disp('Calculating non-linear kernel...');
K=zeros(M,M);
G=zeros(M,N*M);
H = zeros(N*M, N*M);
e=eye(N);
for i=1:M
    xi = dynamics_data(1:N,i);
    for j=1:M
        xj = dynamics_data(1:N,j);
        K(i,j) = getKernel(xi, xj, lambda, ktype);
    end
end

for b1=1:N
    for i=1:M
        xi = dynamics_data(1:N,i);
        for j=1:M
            xj = dynamics_data(1:N,j);
            G(i,j+(b1-1)*M) = getKernelFirstDerivative(xi, xj, lambda, ktype, 2)'*e(:,b1);
        end
    end
end

for b1=1:N
    for b2=1:N
        for i=1:M
            xi = dynamics_data(1:N,i);
            for j=1:M
                xj = dynamics_data(1:N,j);
                H(i+(b1-1)*M,j+(b2-1)*M) = e(:,b1)'*getKernelSecondDerivative(xi, xj, lambda, ktype)*e(:,b2);
            end
        end
    end
end

Q = [K -G; -G' H];


[ee,vv] =eig(Q+Q');
disp(['Minimum eigen value of Q : ' num2str(min(vv(find(vv<0))))]);
% Q = Q+1e-11*eye(size(Q));

x0 = zeros(2*M+2*N*M,1);

auxdata_1 = Q;
auxdata_2 = lambda;
auxdata_3 = dynamics_data;
auxdata_4 = M;
auxdata_5 = N;
auxdata_6 = eps_p;
auxdata_7 = eps_d;

% FMINCON options
lb = zeros(2*M+2*N*M,1);
if(use_derivatives)
ub = [inf*ones(2*M,1);C*ones(2*N*M,1)];
else
ub = [inf*ones(2*M,1);zeros(2*N*M,1)];
end
obj = @(x)asvm_curve_objective(x, auxdata_1, auxdata_2, auxdata_3, auxdata_4, auxdata_5, auxdata_6, auxdata_7);
options=[];
options=optimset('TolX',tol,'TolCon',tol,'Algorithm','sqp','Display','iter','TolFun',tol,'MaxIter',max_iter,'OutputFcn',@fmincon_callback,'MaxFunEvals',inf);

tic
% [x, fval, exitflag] = fmincon(obj, x0, [], [], [ones(1,M),-ones(1,M),zeros(1,2*N*M)], 0, lb, ub, [], options);
[x, fval, exitflag] = fmincon(obj, x0, [], [], [], [], lb, ub, [], options);
% asvm_curve_objective(x, auxdata_1, auxdata_2, auxdata_3, auxdata_4, auxdata_5, auxdata_6, auxdata_7);
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
    betasub = x(M+i*M + 1:2*M+i*M)-x(M+N*M+i*M+1:2*M+N*M+i*M);
    sv_tol_beta = max(abs(betasub))*1e-5; 
    beta_index = find(abs(betasub)  > sv_tol_beta & abs(betasub) > tol);
    beta = betasub(beta_index);
    Svb = dynamics_data(1:N,beta_index);
    
    learned_svm.beta = [learned_svm.beta;-beta(:)];
    learned_svm.Svb = [learned_svm.Svb, [Svb;repmat(e(:,i),1,size(Svb,2))]];
end

learned_svm.gamma = [];
learned_svm.lambda = lambda;
learned_svm.type = ktype;
learned_svm.target = zeros(N,1);

% Calculating b0
learned_svm.b0 = 0;
tmp = 0;
count=0;
for i=1:length(alpha_index)
  
    if(x(alpha_index(i)) < tol && x(alpha_index(i)+M) < C)
        tmp = tmp + (1+eps_p - calculateClassifier(learned_svm, dynamics_data(1:N, alpha_index(i))));
        count=count+1;
    else if(x(alpha_index(i)+M) < tol && x(alpha_index(i)) < C)
       tmp = tmp + (1-eps_p - calculateClassifier(learned_svm, dynamics_data(1:N, alpha_index(i))));
       count=count+1;
        end
    end
end
if(count ~=0)
learned_svm.b0 = tmp/count-1;
end
learned_svm


errs=[];
for i=1:length(dynamics_data)
errs=[errs,calculateClassifier(learned_svm,dynamics_data(1:N,i)) ];
end
disp(['Value Error: [' num2str(mean(errs')) '] +- [' num2str(std(errs')) ']']);
errs=[];
for i=1:length(dynamics_data)
errs=[errs,calculateClassifierDerivative(learned_svm,dynamics_data(1:N,i)) - dynamics_data(N+1:2*N,i)];
end
disp(['Gradient Error: [' num2str(mean(errs')) '] +- [' num2str(std(errs')) ']']);

figure(1)
cla; hold on;
hh=plot(dynamics_data(1,:)',dynamics_data(2,:)','k.');

axis tight

lims=[get(gca,'XLim'), get(gca,'YLim')];
lims = lims - [(lims(2)-lims(1))*0.2 -(lims(2)-lims(1))*0.2 (lims(4)-lims(3))*0.2 -(lims(4)-lims(3))*0.2];
[tx,ty]=meshgrid(linspace(lims(1), lims(2), 50),linspace(lims(3), lims(4), 50));
vx=[];vy=[];
hx=[];
for i=1:size(tx,1)
    for j=1:size(tx,2)
        pt=[tx(i,j);ty(i,j)];
        hx(i,j) = calculateClassifier(learned_svm, pt);

    end
end
delete(hh);
displaySVM(learned_svm, 'range',lims, 'sv',1, 'target',0, 'num_contours',0);
contour(tx, ty, hx, [eps_p, -eps_p],'k--');
plot(dynamics_data(1,:)',dynamics_data(2,:)','k.');
for i=1:size(dynamics_data,2)
   line([dynamics_data(1,i);dynamics_data(1,i) + 0.1*dynamics_data(3,i)] ,[dynamics_data(2,i);dynamics_data(2,i) + 0.1*dynamics_data(4,i)]);
end
axis tight
colormap autumn
shading interp


figure(2)

subplot(1,2,1); cla;hold on; grid on;
plot3(dynamics_data(1,:)', dynamics_data(2,:)',dynamics_data(3,:)','k.');
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
plot3(dynamics_data(1,:)', dynamics_data(2,:)',dynamics_data(4,:)','k.');
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

