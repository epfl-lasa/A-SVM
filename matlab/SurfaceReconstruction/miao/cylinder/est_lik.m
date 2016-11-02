function lik=est_lik(theta,x,y)
var_noise=0.01;
datasize=size(x,1);
k=Kernel_func(x,x,theta(1),theta(2:end))+var_noise*eye(datasize);
KL=chol(k,'lower');
v=KL\y;
% lik=0.5*log(det(k))+0.5*y'*inv(k)*y;
lik=0.5*log(det(k)+realmin)+0.5*v'*v;
end

