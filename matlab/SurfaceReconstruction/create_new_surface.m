clear
close all
figure
axis([-1 1 -1 1]);
unprocessed_rawdata = grabDataFromCursor(1, 0.01);
%%
figure
axis([-1 1 -1 1]);
cla
hold on
posdata=[];
rawdata = unprocessed_rawdata;
all_dynamics_data=[];
for i=1:length(rawdata)
    rawdata{i}(1,:) = smooth(rawdata{i}(1,:)',5)';
    rawdata{i}(2,:) = smooth(rawdata{i}(2,:)',5)';
    veldata =  diff(rawdata{i},1,2);
    vdata=[];
    for j=1:length(veldata)
        if(norm(veldata(:,j)) > 1e-10)
        normv = veldata(:,j)/norm(veldata(:,j));
        grd = [-normv(2);normv(1)];
         all_dynamics_data = [all_dynamics_data,[rawdata{i}(:,j);grd(1);grd(2)]];
        end
    end
    plot(rawdata{i}(1,:),rawdata{i}(2,:),'k-*');
end


for i=1:size(all_dynamics_data,2)
   line([all_dynamics_data(1,i);all_dynamics_data(1,i) + 0.1*all_dynamics_data(3,i)] ,[all_dynamics_data(2,i);all_dynamics_data(2,i) + 0.1*all_dynamics_data(4,i)]);
end