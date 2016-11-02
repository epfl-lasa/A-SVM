function func = classification_dynamics_objective( x, auxdata_4, auxdata_1, auxdata_2, auxdata_5, auxdata_6,auxdata_3, auxdata_7, auxdata_8 )
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


M = auxdata_4;
Q = auxdata_1;
lambda = auxdata_2;
dd = auxdata_3;
P = auxdata_5;
N = auxdata_6;
eps = [5*ones(P,1);1*ones(P,1)];
%     for i=1:N
%         eps(i*P)=0;
%     end

if(auxdata_8)

betasub = x(M+1:M+N*P)-x(M+N*P+1:end);
vec = [x(1:M);betasub];
func = 0.5*vec'*Q*vec - sum(x(1:M)) + (x(M+1:M+N*P)+x(M+N*P+1:end))'*eps + betasub'*reshape(dd(N+1:end,auxdata_7)',N*P,1);
else

    
betasub = x(1:N*P)-x(N*P+1:end);
func = 0.5*betasub'*Q*betasub  + (x(1:N*P)+x(N*P+1:end))'*eps + betasub'*reshape(dd(N+1:end,auxdata_7)',N*P,1);
end

end

