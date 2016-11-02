%% Create data by hand
clear
close all
figure
axis([-1 1 -1 1]);
rawdata = grabDataFromCursor(1, 0.01);

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

%% Reading from obj files
[V,F]=read_vertices_and_faces_from_obj_file('matlab/SurfaceReconstruction/miao/cup.obj');
dd=convert_vf_to_data(V,F);
all_dynamics_data=center_and_normalize_data(dd);

%% Resampling and making testing/training sets
resample=30;
figure
subplot(1,2,1);
N = size(all_dynamics_data,1)/2;
plot3(all_dynamics_data(1,:), all_dynamics_data(2,:), all_dynamics_data(3,:), 'k*');
hold on
grid on
for i=1:size(all_dynamics_data,2)
   line([all_dynamics_data(1,i);all_dynamics_data(1,i) + 0.2*all_dynamics_data(4,i)] ,...
       [all_dynamics_data(2,i);all_dynamics_data(2,i) + 0.2*all_dynamics_data(5,i)], ...
       [all_dynamics_data(3,i);all_dynamics_data(3,i) + 0.2*all_dynamics_data(6,i)],'Color','r','Linewidth',1 );
end
axis equal
rand_ind =  randperm(size(all_dynamics_data,2));
test_size = floor(0.1*size(all_dynamics_data,2));
test_data = all_dynamics_data(:,rand_ind(1:test_size));
dynamics_data = all_dynamics_data(:,rand_ind(test_size+1:resample:end));
subplot(1,2,2);
plot3(dynamics_data(1,:), dynamics_data(2,:), dynamics_data(3,:), 'k*');
hold on
grid on
for i=1:size(dynamics_data,2)
   line([dynamics_data(1,i);dynamics_data(1,i) + 0.2*dynamics_data(4,i)] ,...
       [dynamics_data(2,i);dynamics_data(2,i) + 0.2*dynamics_data(5,i)], ...
       [dynamics_data(3,i);dynamics_data(3,i) + 0.2*dynamics_data(6,i)],'Color','r','Linewidth',1 );
end
axis equal
%% Creating the kernel matrix 
y = ones(size(dynamics_data,2),1);

% Creating a fake data point at the origin
dynamics_data = [zeros(size(dynamics_data,1),1),dynamics_data];
y=[0;y];

sig = 0.3;  
ktype = 'rbf';

Q=mx_get_bsb_value_bsb_derivative_kernel(dynamics_data, 1/(2*sig*sig), ktype);

%% Optimizing 
C=[1;1];
tol=1e-12;
eps = [0.005;0.1];
learned_svm = asvm_bsb_value_bsb_derivative(dynamics_data, y, Q,  1/(2*sig*sig), ktype, C, tol,eps, 0)

% Since we learned a value of h(x)=1 on the data points. Shift it back to
% zero
learned_svm.b0 = learned_svm.b0-1;

N=size(dynamics_data,1)/2;
errs = mx_calculate_classifier(learned_svm, dynamics_data(1:N,:));
disp(['Value Training Error: [' num2str(mean(abs(errs)')) '] +- [' num2str(std(abs(errs)')) ']']);
errs_tst = mx_calculate_classifier(learned_svm, test_data(1:N,:));
disp(['Value Testing Error: [' num2str(mean(abs(errs_tst)')) '] +- [' num2str(std(abs(errs_tst)')) ']']);

errs = mx_calculate_classifier_derivative(learned_svm, dynamics_data(1:N,:))-dynamics_data(N+1:end,:);
disp(['Derivative Training Error: [' num2str(mean(abs(errs)')) '] +- [' num2str(std(abs(errs)')) ']']);
errs_tst = mx_calculate_classifier_derivative(learned_svm, test_data(1:N,:));
disp(['Derivative Testing Error: [' num2str(mean(abs(errs_tst)')) '] +- [' num2str(std(abs(errs_tst)')) ']']);

%% Plot result
plotASVMSurface(learned_svm, dynamics_data, 0.7, 1, 1);
