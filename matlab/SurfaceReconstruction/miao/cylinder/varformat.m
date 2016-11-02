function varformat

load theta0.mat;
 A=1:size(theta0,1);
x=[A',theta0];
dlmwrite('theta0_cyld.txt',x,'delimiter', '\t')
% load Rh.mat;
% 
% Rh=[A',Rh'];
% dlmwrite('rh.txt',Rh,'delimiter', '\t')