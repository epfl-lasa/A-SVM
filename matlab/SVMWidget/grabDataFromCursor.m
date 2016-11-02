function data = grabDataFromCursor(numDemos,dt_timer,varargin)
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
global nbNewData D dt timer_dt fig
nbSamplesMax = numDemos;
if ~exist('nbSamplesMax', 'var'), 
    nbSamplesMax = 3;
end
% D=[-1 1 -1 1];
if isempty(varargin)
    nbNewData=100;
    fig=[];
elseif length(varargin)==1
    nbNewData=varargin{1};
    fig=[];
else
    nbNewData=varargin{1};
    fig=varargin{2};
    D=varargin{3};
end

nbSamples = 0;
dt=dt_timer;

fig = figure(gcf);

hold on
grid on
setappdata(gcf,'motion',[]);
setappdata(gcf,'data',[]);
setappdata(gcf,'nbSamples',nbSamples);
setappdata(gcf,'nbSamplesMax',nbSamplesMax);

% axis(D);
set(fig,'WindowButtonDownFcn',{@wbd});

while(nbSamples~=nbSamplesMax)
  figure(fig);
  nbSamples = getappdata(gcf,'nbSamples');
  if nbSamples==nbSamplesMax
    disp('Finished.');
    data = getappdata(gcf,'data');
    delete(timer_dt);
    clear timer_dt;
    set(fig,'WindowButtonDownFcn','');
    break;
  end
end

% -----------------------------------------------------------------------
function wbd(h,evd) % executes when the mouse button is pressed
global D dt timer_dt flag cur_point fig
disp(['Recording demo ' num2str(getappdata(fig,'nbSamples')+1) ' of ' num2str(getappdata(fig,'nbSamplesMax')) '...']);
cur_point = get(gca,'Currentpoint');
motion = cur_point(1,1:2)';
setappdata(gcf,'motion',motion);
%hold off;
% plot(-1,-1);
% axis(D);
hold on;

%get the values and store them in the figure's appdata
props.WindowButtonMotionFcn = get(fig,'WindowButtonMotionFcn');
props.WindowButtonUpFcn = get(fig,'WindowButtonUpFcn');
setappdata(fig,'TestGuiCallbacks',props);
set(fig,'WindowButtonMotionFcn',{@wbm});
set(fig,'WindowButtonUpFcn',{@wbu});
timer_dt=timer;
set(timer_dt,'ExecutionMode','fixedSpacing');
set(timer_dt,'Period',dt);
set(timer_dt,'TimerFcn',@timerfnc);
cur_point = get(gca,'Currentpoint');
start(timer_dt);
flag=false;

% -----------------------------------------------------------------------
function wbm(h,evd) % executes while the mouse moves
global flag timer_dt cur_point 
cur_point = get(gca,'Currentpoint');
% motion = getappdata(gcf,'motion');
% motion = [motion cur_point(1,1:2)'];
% setappdata(gcf,'motion',motion);

% plot(motion(1,end-1:end),motion(2,end-1:end), 'k', 'lineWidth', 2);
% set(timer_dt,'TimerFcn',@timerfnc);
% start(timer_dt);



% -----------------------------------------------------------------------
function timerfnc(obj, event, string_arg)
global flag timer_dt cur_point fig
% stop(timer_dt);
motion = getappdata(fig,'motion');
motion = [motion cur_point(1,1:2)'];
setappdata(fig,'motion',motion);
plot(motion(1,end-1:end),motion(2,end-1:end), 'k', 'lineWidth', 2);

function wbu(h,evd) % executes when the mouse button is released

global nbNewData timer_dt fig

stop(timer_dt);
motion = getappdata(fig,'motion');
if size(motion,2)==1
  motion = [motion motion];
end

disp('Resampling data...');
%Resampling
nbData = size(motion,2);
%xx = 1:(nbData-1)/(nbNewData-1):nbData;
motion = pchip(1:nbData, motion, 1:nbData);
%motion = spline(1:nbData, motion, xx);

motion_smooth(1,:) = smooth(motion(1,:),1);
motion_smooth(2,:) = smooth(motion(2,:),1);
plot(motion_smooth(1,:),motion_smooth(2,:), 'r', 'lineWidth', 2);

data = getappdata(fig,'data');
nbSamples = getappdata(fig,'nbSamples');
nbSamples = nbSamples + 1;
%data(:,:,nbSamples) = motion_smooth;
%data(:,:,nbSamples) = [1:nbNewData; motion_smooth];
data{nbSamples} = motion_smooth;
setappdata(fig,'nbSamples',nbSamples);
setappdata(fig,'data',data);

%get the properties and restore them
props = getappdata(fig,'TestGuiCallbacks');
set(fig,props);



