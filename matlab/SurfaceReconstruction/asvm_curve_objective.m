function func = asvm_curve_objective( x, auxdata_1, auxdata_2, auxdata_3, auxdata_4, auxdata_5, auxdata_6, auxdata_7)
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


Q = auxdata_1;
lambda = auxdata_2;
dd = auxdata_3;
M = auxdata_4;
N = auxdata_5;
eps_p = auxdata_6;
eps_d = auxdata_7;
yi=ones(1,M);

betasub = x(2*M+1:2*M+N*M)-x(2*M+N*M+1:end);
alphasub=x(1:M)-x(M+1:2*M);
vec = [alphasub;betasub];

func = 0.5*vec'*Q*vec + sum(x(2*M+1:2*M+N*M)+x(2*M+N*M+1:end))*eps_d + betasub'*reshape(dd(N+1:end,:)',N*M,1) + yi*alphasub + eps_p*sum(x(1:M)+x(M+1:2*M));


end

