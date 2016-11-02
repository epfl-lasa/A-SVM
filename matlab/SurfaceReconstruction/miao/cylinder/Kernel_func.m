function K_mat=Kernel_func(x,y,mu,lamda)
xdatasize=size(x,1);
ydatasize=size(y,1);

% for i=1:xdatasize
%         for j=1:ydatasize
%        K_mat(i,j)=mu^2*exp(-(norm(x(i,:)-y(j,:)))^2/(2*lamda^2));
%         end
%     end            
% end

    var_mat=inv(diag(lamda.^2));
for i=1:xdatasize
        for j=1:ydatasize
       K_mat(i,j)=mu^2*exp(-(x(i,:)-y(j,:))*var_mat*(x(i,:)-y(j,:))');
        end
    end            
end
