function [retval] = convertSVM_mat2c(varargin)
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
retval = -1;
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

tmp = load(matfilepath);
ce = fieldnames(tmp);
if(isempty(ce))
    disp('Bad input file!');
    return;
end
toSave = tmp.(ce{1});

[file msg] = fopen(savepath,'w');
if(file<0)
    disp('Error while opening output file');
    return;
end

disp([matfilepath ' ----> ' savepath]);

fprintf(file, '%s\n', toSave.type);
dim = size(toSave.Sva,1);
fprintf(file, '%d\n', dim);
fprintf(file, '%f\n', toSave.lambda);
fprintf(file, '%f\n', toSave.b0);
fprintf(file, '%d\n', length(toSave.alpha));
fprintf(file, '%d\n\n', length(toSave.beta));
for i=1:dim
    fprintf(file, '%f\t',toSave.target(i));
end
fprintf(file, '\n');

for i=1:length(toSave.alpha)
    fprintf(file, '%f\t',toSave.alpha(i));
end
fprintf(file, '\n\n');

for i=1:length(toSave.y)
    fprintf(file, '%d\t',toSave.y(i));
end
fprintf(file, '\n\n');

for i=1:length(toSave.beta)
    fprintf(file, '%f\t',toSave.beta(i));
end
fprintf(file, '\n\n');

for i=1:length(toSave.gamma)
    fprintf(file, '%f\t',toSave.gamma(i));
end
fprintf(file, '\n\n');

for i=1:dim
    for j=1:length(toSave.alpha)
        fprintf(file, '%f\t',toSave.Sva(i,j));
    end
    fprintf(file, '\n');
end
fprintf(file, '\n\n');

for i=1:2*dim
    for j=1:length(toSave.beta)
        fprintf(file, '%f\t',toSave.Svb(i,j));
    end
    fprintf(file, '\n');
end
fprintf(file, '\n\n');

fclose(file);
retval = 0;
