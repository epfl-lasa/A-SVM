function convertData_mat2c(varargin)
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
matfilepath='';
savepath='';

for i=1:length(varargin)
    if( strcmp(varargin{i},'input'))
        matfilepath=varargin{i+1};
    end
    if( strcmp(varargin{i},'output'))
        savepath=varargin{i+1};
    end
end

if(isempty(matfilepath))
    [FileName,PathName,FilterIndex] =   uigetfile('*.mat','Select data mat file');
    matfilepath = [PathName FileName];
end
if(isempty(savepath))
    [FileName,PathName,FilterIndex] =   uiputfile('*.txt','Select output file');
    savepath = [PathName FileName];
end


if(ischar(matfilepath))
    tmp = load(matfilepath);
    ce = fieldnames(tmp);
    if(isempty(ce))
        disp('Bad input file!');
        return;
    end
    toSave = tmp.(ce{1});
else if(iscell(matfilepath))
    toSave = matfilepath;
    else
        disp('ERROR: Bad input!!');
        return;
    end
end

[file msg] = fopen(savepath,'w');
if(file<0)
    disp('Error while opening output file');
    return;
end

% disp([matfilepath ' ----> ' savepath]);
tmp = length(toSave);
fprintf(file, '%d\n', tmp);     %no. of attractors (two)
dim = size(toSave{1}{1},1);
fprintf(file, '%d\n', dim);

for i = 1:length(toSave)
    tmp = length(toSave{i});
    fprintf(file, '%d\n', tmp);     %no. of trajectories
    
    for j= 1:length(toSave{i})
        
        [dim, pointsno] = size(toSave{i}{j});   %dimension and no. of points
        fprintf(file, '%d\n', pointsno);
        for k=1:pointsno
            for l = 1:dim
                fprintf(file, '%4.10f\t', toSave{i}{j}(l,k));
            end
            fprintf(file, ' \n');
        end
    end
end


fclose(file);
