function [retval] = convertSVM_c2mat(varargin)
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
    [FileName,PathName,FilterIndex] =   uigetfile('*.txt','Select SVM result file');
    txtfilepath = [PathName FileName];
end
if(isempty(savepath))
    [FileName,PathName,FilterIndex] =   uiputfile('*.mat','Select output file');
    savepath = [PathName FileName];
end

[file msg] = fopen(savepath,'w');
if(file<0)
    disp('Error while opening output file');
    return;
end

    
file = fopen(txtfilepath,'r');
 if(file<0)
    disp('Error while opening output file');
   
    return;
 end

disp([txtfilepath ' ----> ' savepath]);

    learned_svm = [];
    learned_svm.type = fscanf(file, '%s',1);
    dim = fscanf(file, '%d', 1);
    learned_svm.lambda = fscanf(file, '%f', 1);
     learned_svm.b0 = fscanf(file, '%f', 1);     
    numalpha = fscanf(file, '%d', 1);
    numbeta = fscanf(file, '%d', 1);
    learned_svm.target = fscanf(file, '%f', dim);
    learned_svm.alpha =fscanf(file, '%f', numalpha);


    learned_svm.y = fscanf(file, '%d',[1 numalpha])';

    learned_svm.beta = fscanf(file, '%f', numbeta);
    learned_svm.gamma = fscanf(file, '%f', dim);
    for i=1:numalpha
        for j=1:dim
            learned_svm.Sva(j,i) = fscanf(file, '%f', 1);
        end
    end
            
   for i=1:numbeta
        for j=1:2*dim
            learned_svm.Svb(j,i) = fscanf(file, '%f', 1);
        end
    end
            
    
    save(savepath, 'learned_svm');
  fclose(file);
  retval = 0;
