function output_data = center_and_normalize_data( dynamics_data )
%CENTER_AND_NORMALIZE_DATA Summary of this function goes here
%   Detailed explanation goes here

N = size(dynamics_data,1)/2;
M = size(dynamics_data,2);
output_data = dynamics_data;
output_data(1:N,:) = output_data(1:N,:) - repmat(mean(output_data(1:N,:)')',1,M);


mx = max(max(output_data(1:N,:)'));
mn = min(min(output_data(1:N,:)'));
output_data(1:N,:) = output_data(1:N,:)./max(abs(mx), abs(mn));
end

