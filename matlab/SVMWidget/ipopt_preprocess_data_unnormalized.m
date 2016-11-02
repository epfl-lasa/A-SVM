function [ dynamics_data, data_labels, indices, target_list ] = ipopt_preprocess_data_unnormalized( all_data, target_class )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Copyright (c) 2012 Ashwini SHUKLA, LASA, EPFL,          %%%
%%%          CH-1015 Lausanne, Switzerland, http://lasa.epfl.ch         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The program is free for non-commercial academic use. Please contact the
% author if you are interested in using the software for commercial purposes.
% The software must not be modified or distributed without prior permission
% of the authors. Please acknowledge the authors in any academic publications
% that have made use of this code or part of it. Please use this BibTex
% reference:
%
%
% To get latest upadate of the software please visit
%                          http://asvm.epfl.ch
%
% Please send your feedbacks or questions to:
%                           ashwini.shukla_at_epfl_dot_ch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_labels = [];
velocity_data = [];
dynamics_data = [];

dim = size(all_data{1}{1},1);
target_list = zeros(dim, length(all_data));
indices = [];
for i = 1:length(all_data)
    data = all_data{i};
    num_points = 0;
    velocity_data=[];
    pos_data = [];
    
    target = zeros(dim,1);
    for j=1:length(data)
        target = target + data{j}(:,end);
    end
    target = target/length(data);
    
    target_list(:,i) = target;
    
    for j=1:length(data)

        read_data = data{j};

        tmp = diff(read_data,1,2)';
        
%         for k=1:size(tmp,1)
%             if(norm(tmp(k,:)) > 1e-20)
%                 tmp(k,:) = tmp(k,:)/norm(tmp(k,:));%*(1 + (k-1)*(500-1)/(size(tmp,1)-1));
%             else
%                 tmp(k,:) = [0,0];
%             end
%         end
        if(i==target_class)
            indices = [indices;length(dynamics_data)+length(pos_data)+1];
        end
        
        velocity_data = [velocity_data,tmp'];
        pos_data = [pos_data,read_data(:,1:end-1)];
        num_points = num_points + size(read_data,2)-1;
        
        
        
    end
    
    if(i==target_class)
        one_of_k_target = 1;
        indices = [indices;length(dynamics_data)+length(pos_data)+1];
    else
        one_of_k_target = -1;
    end
    
    data_labels = [data_labels ; repmat(one_of_k_target, num_points,1)];
    dynamics_data = [dynamics_data , [pos_data;velocity_data]];
end

end

