
%%not finish
clc
clear
close all
%%creat a dataset represent the surface.
ns=120;
r=1;
thetar=2*pi*rand(ns,1);
xs=r*cos(thetar);
ys=r*sin(thetar);
zs=1.25*rand(ns,1);
nu=25;
r=1*rand(nu,1);
thetar=2*pi*rand(nu,1);
xu=r.*cos(thetar);
yu=r.*sin(thetar);
zu=1.25*ones(nu,1);
nd=25;
r=1*rand(nd,1);
thetar=2*pi*rand(nd,1);
xd=r.*cos(thetar);
yd=r.*sin(thetar);
zd=zeros(nd,1);
cyld_set=[xs,ys,zs;xu,yu,zu;xd,yd,zd];
%%

ns = [xs,ys,zeros(length(xs),1)];
for i=1:length(ns)
   ns(i,:) = ns(i,:)/norm(ns(i,:)); 
end
nu = [zeros(length(xu),2),ones(length(xu),1) ];
nd = [zeros(length(xd),2),-ones(length(xd),1) ];

all_dynamics_data = [xs',xu',xd';ys',yu',yd';zs',zu',zd';ns',nu',nd'];
%%
% plot3(cyld_set(:,1),cyld_set(:,2),cyld_set(:,3),'r.');
% cc
 xon=cyld_set;
yon=zeros(size(xon,1),1);
xin=[0,0,0.065;];
yin=[-1;];

n=20;
thetar=pi*rand(1,n);
phi=2*pi*rand(1,n);
[x,y,z]=sph2cart(thetar,phi,0.2*ones(1,n));
xout=[x;y;z]';
yout=ones(n,1);
x=[xin;xon;xout];
y=[yin;yon;yout];
% size(x);
% size(y);
%%
[xt,yt,zt]=meshgrid(-0.1:0.01:0.1,-0.1:0.01:0.1,-0.02:0.02:0.16);
size_testdata=size(xt);
x_test=[xt(:),yt(:),zt(:)];
%%
[theta0,fval]=fminunc(@(theta)est_lik(theta,x,y),[0.1,0.1,0.1,0.1])
mu=theta0(1);
kernel_width=theta0(2:end);
[ver,f_var,f_mean,invk,alp]=Est_gps(x,y,mu,kernel_width,x_test,size_testdata);
xlabel('X');
ylabel('Y');
zlabel('Z');
qq=reshape(f_mean,size(xt));
psurf=patch(isosurface(xt,yt,zt,qq,0));
hold on;
set(psurf,'FaceColor','green','EdgeColor','none',...
    'FaceAlpha',0.8);
daspect([1,1,1]),view(150,30);
camlight 
lighting phong
plot3(xon(:,1),xon(:,2),xon(:,3),'r.','MarkerSize',16),hold on;
xtrain_cyl=x;
ytrain_cyl=y;
save  xtrain_cyl;
save ytrain_cyl;
save alp;
save theta0;
max(f_var)