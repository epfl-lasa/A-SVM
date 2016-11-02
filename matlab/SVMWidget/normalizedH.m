function f = normalizedH(lsvm, point)

f=lsvm.b0;
d=lsvm.lambda;
type = lsvm.type;
dim = length(point);
e = eye(dim);
for i=1:length(lsvm.alpha)
    f = f + lsvm.alpha(i)*lsvm.y(i)*getKernel(point, lsvm.Sva(:,i),d, type);
end

for i=1:length(lsvm.beta1)
    f = f - lsvm.beta1(i)*getKernelFirstDerivative(point, lsvm.Svb1(1:dim,i), d, type, 2)'*e(:,1);
end

for i=1:length(lsvm.beta2)
    f = f - lsvm.beta2(i)*getKernelFirstDerivative(point, lsvm.Svb2(1:dim,i), d, type, 2)'*e(:,2);
end