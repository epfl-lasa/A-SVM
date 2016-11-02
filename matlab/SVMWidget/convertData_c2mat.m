function convertData_c2mat(varargin)
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

txtfilepath='';
savepath='';

for i=1:length(varargin)
    if( strcmp(varargin{i},'input'))
        txtfilepath=varargin{i+1};
    end
    if( strcmp(varargin{i},'output'))
        savepath=varargin{i+1};
    end
end

if(isempty(txtfilepath))
    [FileName,PathName,FilterIndex] =   uigetfile('*.txt','Select data file');
    txtfilepath = [PathName FileName];
end
if(isempty(savepath))
    [FileName,PathName,FilterIndex] =   uiputfile('*.mat','Select output file');
    savepath = [PathName FileName];
end

disp([txtfilepath ' ----> ' savepath]);

[file msg] = fopen(txtfilepath,'r');
if(file<0)
    disp('Error while opening data file');
    return;
end


num_att = fscanf(file, '%d', 1);
toSave=cell(num_att,1);

dim = fscanf(file, '%d', 1);

for i=1:num_att
    num_traj = fscanf(file, '%d', 1);
toSave{i} = cell(num_traj,1);
    for j=1:num_traj
        num_pts = fscanf(file, '%d', 1);
        toSave{i}{j} = zeros(dim,num_pts);
        for k = 1:num_pts
           for l=1:dim
                toSave{i}{j}(l,k) = fscanf(file, '%f', 1);
           end
        end

    end
end

save(savepath,'toSave');
fclose(file);
