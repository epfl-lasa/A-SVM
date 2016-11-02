function gradient = normalizedHgrad(lsvm, point)


gradient =zeros(length(point),1);
d=lsvm.lambda;
type = lsvm.type;
dim = length(point);
e = eye(dim);
for i=1:length(lsvm.alpha)
    gradient = gradient + lsvm.alpha(i)*lsvm.y(i)*getKernelFirstDerivative(point, lsvm.Sva(:,i),d, type, 1);
end

for i=1:length(lsvm.beta1)
    gradient = gradient - lsvm.beta1(i)*getKernelSecondDerivative(point, lsvm.Svb1(1:dim,i), d, type)*e(:,1);
end

for i=1:length(lsvm.beta2)
    gradient = gradient - lsvm.beta2(i)*getKernelSecondDerivative(point, lsvm.Svb2(1:dim,i), d, type)*e(:,2);
end