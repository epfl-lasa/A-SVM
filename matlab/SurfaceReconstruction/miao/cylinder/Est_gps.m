function [ver,f_var,f_mean,invk,alp]=Est_gps(x,y,mu,kernel_width,x_test,size_testdata);



% plot3(x(:,1),x(:,2),x(:,3),'g.')
var_noise=0.01;
datasize=size(x,1);
K_mat=Kernel_func(x,x,mu,kernel_width)+var_noise*eye(datasize);
invk=inv(K_mat);
[vecmat,valmat]=eig(K_mat);
L=chol(K_mat,'lower');

K_test_train=Kernel_func(x_test,x,mu,kernel_width);
K_test=Kernel_func(x_test,x_test,mu,kernel_width);
f_mean=1+K_test_train*(L'\(L\(y-1)));
alp=L'\(L\(y-1));
qq=reshape(f_mean,size_testdata);
xt=reshape(x_test(:,1),size_testdata);
yt=reshape(x_test(:,2),size_testdata);
zt=reshape(x_test(:,3),size_testdata);
 [face,ver]=isosurface(xt,yt,zt,qq,0);
%  size(ver);
 x_testsurf=ver;
% K_test_train=Kernel_func(x_testsurf,x,mu,kernel_width);
% K_test=Kernel_func(x_testsurf,x_testsurf,mu,kernel_width);
% fver_mean=1+K_test_train*(L'\(L\(y-1)));

if size(ver,1)==0
error('no verticle')
end
 for i=1:size(x_testsurf,1)    
K_ss=Kernel_func(x_testsurf(i,:),x_testsurf(i,:),mu,kernel_width);
K_s=Kernel_func(x_testsurf(i,:),x,mu,kernel_width);
v=L\K_s';
  f_var1(i)=K_ss-v'*v;
end
f_var=f_var1';
end
