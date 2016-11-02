function varargout = svm_widget(varargin)
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


% SVM_WIDGET MATLAB code for svm_widget.fig
%      SVM_WIDGET, by itself, creates a new SVM_WIDGET or raises the existing
%      singleton*.
%
%      H = SVM_WIDGET returns the handle to a new SVM_WIDGET or the handle to
%      the existing singleton*.
%
%      SVM_WIDGET('CALLBACK',hObject,eventData,handles,...) calls the local
%      functioni named CALLBACK in SVM_WIDGET.M with the given input arguments.
%
%      SVM_WIDGET('Property','Value',...) creates a new SVM_WIDGET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before svm_widget_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to svm_widget_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help svm_widget

% Last Modified by GUIDE v2.5 04-Mar-2013 14:01:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @svm_widget_OpeningFcn, ...
    'gui_OutputFcn',  @svm_widget_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before svm_widget is made visible.
function svm_widget_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to svm_widget (see VARARGIN)


global numAtt all_data  learned_svm real_data_ranges all_data_ranges scroll_count  ipopt_force_quit data_path model_path

ipopt_force_quit = false;


numAtt = 0;
all_data={};
classifier_val = [];
ch = get(handles.attPanel,'Children');

for i=1:length(ch)
    set(ch(i),'Enable','off');
end
updateTargetGroup(handles);

% Choose default command line output for svm_widget
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% colorbar ('location','WestOutside');
axes(handles.axes1);
real_data_ranges = [-1 1 -1 1];
axis(real_data_ranges);
scroll_count=0;
all_data_ranges = real_data_ranges;
% axis(all_data_ranges);
hold on
grid on

learned_svm=[];
data_path = '.';
model_path='.';

% UIWAIT makes svm_widget wait for user response (see UIRESUME)
% uiwait(handles.main);


% --- Outputs from this function are returned to the command line.
function varargout = svm_widget_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global all_data numAtt real_data prev_ipopt_sol all_data_ranges real_data_ranges scroll_count data_path 

[FileName,PathName,FilterIndex] = uigetfile({'*.mat';'*.txt'},'Select data file',data_path);

if(FilterIndex == 0)
    return;
else if(FilterIndex == 1)
        data_path = PathName;
        s = importdata([PathName FileName]);
        if(isempty(s))
            return;
        end
        tmp = s;
    else if(FilterIndex == 2)
            data_path = PathName;
            convertData_c2mat('input',[PathName FileName],'output','.tmprappdata.mat');
            s = importdata('.tmprappdata.mat');
            if(isempty(s))
                return;
            end
            tmp = s;
        end
    end
    
end

ch = get(handles.attPanel,'Children');

if(length(tmp)>length(ch))
    disp(['Cannot support more than' num2str(length(ch)) ' attractors']);
    %addOutputString(['Cannot support more than' num2str(length(ch)) ' attractors'],handles.cmdWnd);
    return;
end

if(length(tmp) == 0)
    disp('Need atleast one attractor');
    %  addOutputString('Need atleast one attractor' , handles.cmdWnd);
    return;
end


all_data = tmp;
numAtt = length(all_data);
for j=1:length(ch)
   set(ch(j),'Enable','off'); 
end
for i=1:length(all_data)
    for j=1:length(ch)
        if(~isempty(strfind(get(ch(j), 'String'), int2str(i))))
            set(ch(j),'Enable','on');
            set(ch(j),'Value',1);
        end
    end
end

axes(handles.axes1);
cla
hold on
grid on
plot_all_data([1:length(all_data)]');
updateTargetGroup(handles);


for i=1:length(all_data)
    dim = size(all_data{i}{1},1);
    tar = zeros(dim,1);
    for j=1:length(all_data{i})
        tar = tar + all_data{i}{j}(:,end);
    end
    tar = tar/length(all_data{i});
    
    for j=1:length(all_data{i})
        all_data{i}{j} = all_data{i}{j} + repmat(tar-all_data{i}{j}(:,end),1,size(all_data{i}{j},2));
    end
end

real_data = all_data;
real_data_ranges = getDisplayRange();
all_data_ranges = real_data_ranges;
scroll_count=0;
prev_ipopt_sol = [];
cent = zeros(size(all_data{1}{1},1),1);
for i=1:length(all_data)
   cent = cent+all_data{i}{1}(:,end); 
end
cent  = cent/length(all_data);

if(abs(real_data_ranges(1)-real_data_ranges(2)) >= abs(real_data_ranges(3)-real_data_ranges(4)))
    len = abs(real_data_ranges(1)-real_data_ranges(2));
else
    len = abs(real_data_ranges(3)-real_data_ranges(4));
end
axis([cent(1)-len/2.0 cent(1)+len/2.0 cent(2)-len/2.0 cent(2)+len/2.0]);

updateTargetValBox(handles);
refreshDisplay(handles);
% axis(all_data_ranges);
% axis equal


% --------------------------------------------------------------------
function save_data_Callback(hObject, eventdata, handles)
% hObject    handle to save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global all_data toSave data_path

toSave = all_data(getActiveIndices(handles.attPanel));
saved_data = toSave;


[FileName,PathName,FilterIndex] = uiputfile({'*.mat';'*.txt'},'Save data as...', data_path);
if(FilterIndex == 0)
    return;
else if(FilterIndex == 1)
        data_path = PathName;
        save([PathName FileName],'saved_data');
        
    else if(FilterIndex == 2)
            data_path = PathName;
            save('.tmprappdata.mat','saved_data');
            convertData_mat2c('input','.tmprappdata.mat','output',[PathName FileName]);
        end
    end
    
end



% --- Executes on button press in Att1.
function Att1_Callback(hObject, eventdata, handles)
% hObject    handle to Att1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Att1
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
updateTargetGroup(handles);

% --- Executes on button press in Att2.
function Att2_Callback(hObject, eventdata, handles)
% hObject    handle to Att2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Att2
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
updateTargetGroup(handles)

% --- Executes on button press in Att3.
function Att3_Callback(hObject, eventdata, handles)
% hObject    handle to Att3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Att3
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
updateTargetGroup(handles);

% --- Executes on button press in Att4.
function Att4_Callback(hObject, eventdata, handles)
% hObject    handle to Att4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Att4
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
updateTargetGroup(handles);

% --- Executes on button press in Att5.
function Att5_Callback(hObject, eventdata, handles)
% hObject    handle to Att5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Att5
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
updateTargetGroup(handles);

% --- Executes on button press in Att6.
function Att6_Callback(hObject, eventdata, handles)
% hObject    handle to Att6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Att6
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
updateTargetGroup(handles);


function plot_all_data(indices)
global all_data 

if(~isempty(all_data))
    hl = [];
    if(size(all_data{1}{1},1) == 2)
        
        for i=indices
            tmp = all_data{i};
            for j=1:length(tmp)
                hl = [hl plot(tmp{j}(1,:)',tmp{j}(2,:)','k.')];
            end
        end
    else if(size(all_data{1}{1},1) == 3)
            for i=indices
                tmp = all_data{i};
                for j=1:length(tmp)
                    hl = [hl plot3(tmp{j}(1,:)',tmp{j}(2,:)',tmp{j}(3,:),'k.')];
                end
            end
        end
        
    end
    hg = hggroup;
    set(hl,'Parent',hg);
    
    set(get(get(hg,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on');
    %     axis equal
    box on
end



function indices = getActiveIndices(panelHandle)

ch = get(panelHandle,'Children');
indices = [];
for i=1:length(ch)
    if(get(ch(i),'Value') == 1)
        indices = [indices,eval(get(ch(i), 'String'))];
    end
end
indices = sort(indices);

function updateTargetGroup(handles)
ch1 = get(handles.attPanel,'Children');

ch3 = handles.tclass_list;

strdata={};
cnt=1;
for i=1:length(ch1)
    if(strcmp(get(ch1(i),'Enable'), 'on') && get(ch1(i),'Value') == 1)
       strdata{cnt} = get(ch1(i), 'String');
       cnt=cnt+1;
    end
end

if(~isempty(strdata))
    set(ch3, 'Enable','on');
    set(ch3,'String',strdata); 
else
    set(ch3, 'Enable','off');
end


function class_id = getTargetClass(panelHandle)

algo_selection = get(panelHandle,'Value');
all_str = get(panelHandle, 'String');
select_str = all_str{algo_selection};

if(strcmp(get(panelHandle, 'Enable'),'off'))
    class_id = 0;
else
    class_id = eval(select_str);
end

% --- Executes on button press in classify.
function classify_Callback(hObject, eventdata, handles)
% hObject    handle to classify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global all_data learned_svm ipopt_force_quit app_handle prev_ipopt_sol


indices = getActiveIndices(handles.attPanel);

if(getTargetClass(handles.tclass_list) == 0)
    disp('Please choose a target class before classifying.');
    return
end

algo_selection = get(handles.select_algo,'Value');
all_str = get(handles.select_algo, 'String');
select_str = all_str{algo_selection};
app_handle = handles.main;
target_class = getTargetClass(handles.tclass_list);
sigma = eval(get(handles.box_lambda,'String'));
lambda = 1/(2*sigma*sigma);

C = eval(get(handles.box_C,'String'));
tol = eval(get(handles.box_tol,'String'));
btol = eval(get(handles.box_btol,'String'));
max_iter = eval(get(handles.box_maxiter,'String'));
R = eval(get(handles.box_eps, 'String'));

set(hObject, 'Enable','off');
exepflag=false;


if(~isempty(strfind(select_str, 'LibSVM')))
    disp('Using LibSVM to classify (only)');
    try
        if( get(handles.chk_Auto, 'Value') )
            [dd, dl, ind, tl] = ipopt_preprocess_data(all_data(indices),find(indices == getTargetClass(handles.tclass_list)) );
            N = size(dd,1)/2;
            perm = [find(dl == 1);find(dl == -1)];
           lsvm = GMKL_SVC_OptimizeKernel(dd(1:N,perm), dl(perm), ones(N,1));
           lambda = lsvm.lambda;
        end
        learned_svm = libsvm_one_vs_all( all_data(indices) , lambda, find(indices == getTargetClass(handles.tclass_list)) );
        
    catch exception
        set(hObject, 'Enable','on');
        disp(['ERROR: ' exception.message]);
        exepflag=true;
    end
    
else if(~isempty(strfind(select_str, 'Ipopt')))
        ipopt_force_quit = false;
        disp('Using Ipopt');
        
        % Preparing data
        [dynamics_data, data_labels, indices, target_list] = ipopt_preprocess_data(all_data(indices), find(indices == getTargetClass(handles.tclass_list)));
        
        tim = timer('TimerFcn',@mycallback, 'Period', 1,'ExecutionMode','fixedSpacing');
        start(tim);
        
        % Learning
        
        try
            if(getTargetClass(handles.ipopt_init) == 1)
                [x,learned_svm] = ipopt_one_vs_all( dynamics_data, data_labels, lambda, C, tol, max_iter, target_list(:,target_class), 0 );
            end
            
            if(getTargetClass(handles.ipopt_init) == 4)
                [x,learned_svm] = ipopt_one_vs_all( dynamics_data, data_labels, lambda, C, tol, max_iter, target_list(:,target_class), 1 );
            end
            
            if(getTargetClass(handles.ipopt_init) == 3)
                [x,learned_svm] = ipopt_one_vs_all( dynamics_data, data_labels, lambda, C, tol, max_iter, target_list(:,target_class), 2 );
            end
            
            if(getTargetClass(handles.ipopt_init) == 2)
                if(isempty(prev_ipopt_sol))
                    [x,learned_svm] = ipopt_one_vs_all( dynamics_data, data_labels, lambda, C, tol, max_iter, target_list(:,target_class), 0 );
                else
                    [x,learned_svm] = ipopt_one_vs_all( dynamics_data, data_labels, lambda,  C, tol, max_iter, target_list(:,target_class), prev_ipopt_sol );
                end
            end
            prev_ipopt_sol = x;
            
        catch exception
            exepflag=true;
            set(hObject, 'Enable','on');
            disp(['ERROR: ' exception.message]);
        end
        
        stop(tim);
        delete(tim);
        
    else if(~isempty(strfind(select_str, 'Quadprog')))
            ipopt_force_quit = false;
            disp('Using Quadprog');
            
            % Preparing data
            [dynamics_data, data_labels, indices, target_list] = ipopt_preprocess_data(all_data(indices), find(indices == getTargetClass(handles.tclass_list)));
            
            tim = timer('TimerFcn',@mycallback, 'Period', 1,'ExecutionMode','fixedSpacing');
            start(tim);
            
            try
                ch = get(handles.ipopt_init,'Children');
                for i=1:length(ch)
                    if get(ch(i),'Value') == 1
                        initstr = get(ch(i), 'String');
                        break;
                    end
                end
                
                if(~isempty(strfind(select_str, 'nu-')))
                    disp('Learning nu-ASVM');
                    use_nu = true;
                    nu=eval(get(handles.box_nu, 'String'));
                else
                    use_nu=false;
                    nu=0;
                end
                
                if(~isempty(strfind(initstr, 'Zeros')))
                    
                    init=0;
                else if (~isempty(strfind(initstr, 'Random')))
                        init=1;
                    else if (~isempty(strfind(initstr, 'SVM')))
                            init=2;
                        else  if (~isempty(strfind(initstr, 'Last')))
                                if(isempty(prev_ipopt_sol))
                                    init=0;
                                else
                                    init=prev_ipopt_sol;
                                end
                             end
                        end
                    end
                end
                
                 if( get(handles.chk_Auto, 'Value') )
                     disp('Tuning kernel width...');
                     lsvm = GMKL_ASVM_osb_OptimizeKernel(dynamics_data, data_labels, target_list(:,target_class), R, ones(size(dynamics_data,1)/2,1));
                     lambda = lsvm.lambda;
                 end
                [x,learned_svm, exitflag] = quadprog_one_vs_all( dynamics_data, data_labels, lambda, C, nu, tol, R, max_iter, target_list(:,target_class), init, use_nu );
                prev_ipopt_sol = x;
                
            catch exception
                set(hObject, 'Enable','on');
                disp(['ERROR: ' exception.message]);
                exepflag=true;
            end
            
            stop(tim);
            delete(tim);
            
        else if(~isempty(strfind(select_str, 'NLOPT')))
                
                toSave = all_data(getActiveIndices(handles.attPanel));
                
                save('.tmprappdata.mat','toSave');
                convertData_mat2c('input', '.tmprappdata.mat','output','.matlabdata.txt');
                file = fopen('.optparams_matlab.ini','w');
                fprintf(file, 'algo\t\tLD_AUGLAG\n');
                if(~isempty(strfind(select_str, 'LBFGS')))
                    fprintf(file, 'sub_algo\t\tLD_LBFGS\n');
                else if(~isempty(strfind(select_str, 'SQP')))
                        fprintf(file, 'sub_algo\t\tLD_SLSQP\n');
                    else
                        disp('Invalid selection in algorithm!!');
                        set(hObject, 'Enable','on');
                        return;
                    end
                end
                fprintf(file, 'C\t\t%f\n', C);
                fprintf(file, 'ftol_rel\t\t%s\n', num2str(tol));
                fprintf(file, 'xtol_rel\t\t%s\n', num2str(tol));
                fprintf(file, 'constr_tol\t\t%f\n', 1e-4);
                fprintf(file, 'max_eval\t\t%d\n', max_iter);
                fprintf(file, 'max_time_sec\t\t%d\n', 300);
                fprintf(file, 'verbose\t\t%s\n', 'off');
                fclose(file);

                disp('Solving...');
                drawnow;
            learned_svm = mxNLOPTSolver(all_data(indices), find(indices == getTargetClass(handles.tclass_list))...
                , sigma,'.optparams_matlab.ini' );
            else if(~isempty(strfind(select_str, 'SMO')))
                    if( ~isempty(strfind(select_str, 'C++')))
                        disp('Calling mxSMOSolver...');
                         file = fopen('.optparams_matlab.ini','w');
                         fprintf(file, 'C\t\t%f\n',C);
                         fprintf(file, 'classification_tol\t\t%f\n',tol);
                         fprintf(file, 'lyapunov_tol\t\t%f\n',btol);
                         fprintf(file, 'lyapunov_relaxation\t\t0\n');
                         fprintf(file, 'max_eval\t\t%f\n',max_iter);
                         fprintf(file, 'verbose\t\toff\n');
                         fclose(file);
                         convertData_mat2c('input',all_data(indices), 'output','.matlabdata.txt');
                        learned_svm = mxSMOSolver(all_data(indices), ...
                            find(indices == getTargetClass(handles.tclass_list)), sigma,'.optparams_matlab.ini');
                   
                    else if(~isempty(strfind(select_str, 'MAT')))
                      
                            learned_svm = smo_one_vs_all( all_data(indices),  find(indices == getTargetClass(handles.tclass_list)), lambda, C, tol);
                            
                        else
                            disp('Invalid selection in algorithm!!');
                            exepflag = true;
                        end
                    end
                else
                    disp('ERROR: Algorithm not implemented!!');
                    exepflag = true;
                end
            end
        end
    end
end

set(hObject, 'Enable','on');

if(exepflag)
    
    return;
end
disp(' ');
disp(['Alpha*Y (' num2str(length(learned_svm.alpha)) ') : ' num2str((learned_svm.alpha.*learned_svm.y)')]);
disp(['Beta    (' num2str(length(learned_svm.beta))  ') : ' num2str(learned_svm.beta')]);
disp(['Gamma       : ' num2str(learned_svm.gamma(:)')]);
disp(['Lambda      : ' num2str(learned_svm.lambda(:)')]);
disp(['B0          : ' num2str(learned_svm.b0)]);
disp(' ');disp(' ');
axes(handles.axes1);
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
plot_curr_svm(handles);
% axis tight




function mycallback(object, event, string_arg)
drawnow;



% --- Executes on selection change in select_algo.
function select_algo_Callback(hObject, eventdata, handles)
% hObject    handle to select_algo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_algo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_algo
contents = cellstr(get(hObject,'String'));
selected = contents{get(hObject,'Value')};

ch = get(handles.ipopt_init,'Children');

if(~isempty(strfind(selected, 'nu-')))
%      for i=1:length(ch)
%          set(ch(i) , 'Enable','off');
%      end
     set(handles.box_nu, 'Enable','on');
else
    
     set(handles.box_nu, 'Enable','off');
end

if(~isempty(strfind(selected, 'LibSVM')))
    for i=1:length(ch)
        if(~isempty(strfind(get(ch(i), 'String'), 'Zeros')))
            set(ch(i) , 'Enable','off');
        else if(~isempty(strfind(get(ch(i), 'String'), 'Random')))
                set(ch(i) , 'Enable','off');
            else if (~isempty(strfind(get(ch(i), 'String'), 'Classifier')))
                    set(ch(i) , 'Enable','off');
                else if (~isempty(strfind(get(ch(i), 'String'), 'Last')))
                        set(ch(i) , 'Enable','off');
                    end
                end
            end
        end
    end
    set(handles.box_btol, 'Enable','off');
else if (~isempty(strfind(selected, 'SMO')))
        
        for i=1:length(ch)
            if(~isempty(strfind(get(ch(i), 'String'), 'Zeros')))
                set(ch(i) , 'Enable','off');
            else if(~isempty(strfind(get(ch(i), 'String'), 'Random')))
                    set(ch(i) , 'Enable','off');
                else if (~isempty(strfind(get(ch(i), 'String'), 'Classifier')))
                        set(ch(i) , 'Enable','off');
                    else if (~isempty(strfind(get(ch(i), 'String'), 'Last')))
                            set(ch(i) , 'Enable','off');
                        end
                    end
                end
            end
        end
        set(handles.box_btol, 'Enable','on');
    else if  (~isempty(strfind(selected, 'NLOPT')) || ~isempty(strfind(selected, 'IPOPT')) || ~isempty(strfind(selected, 'Quadprog')))
            
            for i=1:length(ch)
                if(~isempty(strfind(get(ch(i), 'String'), 'Zeros')))
                    set(ch(i) , 'Enable','on');
                else if(~isempty(strfind(get(ch(i), 'String'), 'Random')))
                        set(ch(i) , 'Enable','on');
                    else if (~isempty(strfind(get(ch(i), 'String'), 'Classifier')))
                            set(ch(i) , 'Enable','on');
                        else if (~isempty(strfind(get(ch(i), 'String'), 'Last')))
                                set(ch(i) , 'Enable','on');
                            end
                        end
                    end
                end
            end
            set(handles.box_btol, 'Enable','off');
        end
    end
end

% --- Executes during object creation, after setting all properties.
function select_algo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_algo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  ipopt_force_quit
disp('Force Quitting..');
ipopt_force_quit = true;
drawnow;

set(handles.classify, 'Enable','on');




% --- Executes on slider movement.
function cont_density_Callback(hObject, eventdata, handles)
% hObject    handle to cont_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global  learned_svm

if(isempty(learned_svm))
    disp('Learn a SVM obejct first!');
    return;
end

cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
plot_curr_svm(handles);
% axis tight






% --- Executes during object creation, after setting all properties.
function cont_density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cont_density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Att7.
function Att7_Callback(hObject, eventdata, handles)
% hObject    handle to Att7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Att7
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
updateTargetGroup(handles);

% --- Executes on button press in Att8.
function Att8_Callback(hObject, eventdata, handles)
% hObject    handle to Att8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Att8
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));
updateTargetGroup(handles);



function data_rng = getDisplayRange()
global all_data

N = size(all_data{1}{1},1);
margin = 0.5;
if(N == 3)
    
    min_x = inf; min_y = inf; min_z = inf;
    max_x = -inf; max_y = -inf; max_z = -inf;
    for i=1:length(all_data)
        for j=1:length(all_data{i})
            if(min_x > min(all_data{i}{j}(1,:)))
                min_x = min(all_data{i}{j}(1,:));
            end
            if(max_x < max(all_data{i}{j}(1,:)))
                max_x = max(all_data{i}{j}(1,:));
            end
            
            if(min_y > min(all_data{i}{j}(2,:)))
                min_y = min(all_data{i}{j}(2,:));
            end
            if(max_y < max(all_data{i}{j}(2,:)))
                max_y = max(all_data{i}{j}(2,:));
            end
            
            if(min_z > min(all_data{i}{j}(3,:)))
                min_z = min(all_data{i}{j}(3,:));
            end
            if(max_z < max(all_data{i}{j}(3,:)))
                max_z = max(all_data{i}{j}(3,:));
            end
        end
    end
    
    min_x = min_x - (max_x-min_x)*margin;
    min_y = min_y - (max_y-min_y)*margin;
    min_z = min_z - (max_z-min_z)*margin;
    
    max_x = max_x + (max_x-min_x)*margin;
    max_y = max_y + (max_y-min_y)*margin;
    max_z = max_z + (max_z-min_z)*margin;
    
    data_rng = [min_x max_x min_y max_y min_z max_z];
else if(N ==2)
        
        min_x = inf; min_y = inf;
        max_x = -inf; max_y = -inf;
        for i=1:length(all_data)
            for j=1:length(all_data{i})
                if(min_x > min(all_data{i}{j}(1,:)))
                    min_x = min(all_data{i}{j}(1,:));
                end
                if(max_x < max(all_data{i}{j}(1,:)))
                    max_x = max(all_data{i}{j}(1,:));
                end
                
                if(min_y > min(all_data{i}{j}(2,:)))
                    min_y = min(all_data{i}{j}(2,:));
                end
                if(max_y < max(all_data{i}{j}(2,:)))
                    max_y = max(all_data{i}{j}(2,:));
                end
            end
        end
        
        min_x = min_x - (max_x-min_x)*margin;
        min_y = min_y - (max_y-min_y)*margin;
        
        max_x = max_x + (max_x-min_x)*margin;
        max_y = max_y + (max_y-min_y)*margin;
        
        data_rng = [min_x max_x min_y max_y];
    end
end


% --- Executes on slider movement.
function smoothness_slider_Callback(hObject, eventdata, handles)
% hObject    handle to smoothness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if( getTargetClass(handles.tclass_list) == 0)
    disp('Choose a target class to smooth first!');
    return;
end
updateData(handles);
cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));

% --- Executes during object creation, after setting all properties.
function smoothness_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function resample_slider_Callback(hObject, eventdata, handles)
% hObject    handle to resample_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if( getTargetClass(handles.tclass_list) == 0)
    disp('Choose a target class to smooth first!');
    return;
end

updateData(handles);

cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));


% --- Executes during object creation, after setting all properties.
function resample_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resample_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes during object creation, after setting all properties.
function smooth_lowlim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smooth_lowlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function resample_lowlim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resample_lowlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function smooth_uplim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smooth_uplim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function resample_uplim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resample_uplim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function box_lambda_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function box_maxiter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_maxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function box_tol_Callback(hObject, eventdata, handles)
% hObject    handle to box_tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_tol as text
%        str2double(get(hObject,'String')) returns contents of box_tol as a double
if ( strcmp(get(handles.box_btol, 'Enable'), 'off'))
   set(handles.box_btol, 'String',get(handles.box_tol, 'String')); 
end

% --- Executes during object creation, after setting all properties.
function box_tol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_tol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function load_svm_Callback(hObject, eventdata, handles)
% hObject    handle to load_svm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global learned_svm model_path

[FileName,PathName,FilterIndex] = uigetfile({'*.mat';'*.txt'},'Select A-SVM model file', model_path);
if(FilterIndex == 0)
    return;
else if(FilterIndex == 1)
        model_path = PathName;
        s = importdata([PathName FileName]);
        if(isempty(s))
            return;
        end
        
        learned_svm=s;
    else if(FilterIndex == 2)
            model_path = PathName;
            convertSVM_c2mat('input',[PathName FileName],'output','.tmprappdata.mat');
            s = importdata('.tmprappdata.mat');
            if(isempty(s))
                return;
            end
            
            learned_svm=s;
        end
    end
    
end

learned_svm
refreshDisplay(handles);

% --------------------------------------------------------------------
function save_svm_Callback(hObject, eventdata, handles)
% hObject    handle to save_svm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global learned_svm model_path
saved_svm = learned_svm;

[FileName,PathName,FilterIndex] = uiputfile({'*.mat';'*.txt'},'Save A-SVM as...', model_path);
if(FilterIndex == 0)
    return;
else if(FilterIndex == 1)
        model_path = PathName;
        save([PathName FileName],'saved_svm');
        
    else if(FilterIndex == 2)
            model_path = PathName;
            save('.tmprappdata.mat','saved_svm');
            convertSVM_mat2c('input','.tmprappdata.mat','output',[PathName FileName]);
        end
    end
    
end


% --- Executes during object creation, after setting all properties.
function box_C_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function testModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global learned_svm

% cla
% hold on
% grid on
% plot_all_data(getActiveIndices(handles.attPanel));
% plot_curr_svm(handles);

[traj_list, spur_att, h_val] = testModel(learned_svm, [get(handles.axes1, 'XLim'),get(handles.axes1, 'YLim'),get(handles.axes1, 'ZLim')], eval(get(handles.box_step,'String')),eval(get(handles.box_mesh,'String')), 2);


N = size(learned_svm.Sva,1);

if(N==2)
    %     for i=1:length(traj_list)
    %         if(~isempty(traj_list{i}))
    %             plot(traj_list{i}(:,1), traj_list{i}(:,2),'b-');
    %         end
    %     end
    %
    %     if(~isempty(spur_att))
    %         plot(spur_att(1,:)',spur_att(2,:)','r*','Linewidth',3,'Markersize',10);
    %     end
end
if(N==3)
    %     for i=1:length(traj_list)
    %         if(~isempty(traj_list{i}))
    %             plot3(traj_list{i}(:,1), traj_list{i}(:,2),traj_list{i}(:,3),'b-');
    %             plot3(traj_list{i}(1,1),traj_list{i}(1,2),traj_list{i}(1,3),'ko','Linewidth',2);
    %         end
    %     end
    %
    %     if(~isempty(spur_att))
    %         plot3(spur_att(1,:)',spur_att(2,:)',spur_att(3,:)','r*','Linewidth',3,'Markersize',10);
    %     end
end

% --- Executes on scroll wheel click while the figure is in focus.
function main_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to main (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
% hObject

% axis(all_data_ranges);
curr_lim = [get(handles.axes1, 'XLim') get(handles.axes1, 'YLim') get(handles.axes1, 'ZLim')];
xlen = curr_lim(2) - curr_lim(1);
ylen = curr_lim(4) - curr_lim(3);
zlen = curr_lim(6) - curr_lim(5);

xc = (curr_lim(2) + curr_lim(1))/2.0;
yc = (curr_lim(4) + curr_lim(3))/2.0;
zc = (curr_lim(6) + curr_lim(5))/2.0;

xlen =  xlen + eventdata.VerticalScrollCount*0.1*xlen;
ylen =  ylen + eventdata.VerticalScrollCount*0.1*ylen;
zlen =  zlen + eventdata.VerticalScrollCount*0.1*zlen;
set(handles.axes1,'XLim',[xc - xlen/2.0, xc + xlen/2.0],'YLim',[yc - ylen/2.0, yc + ylen/2.0],'ZLim',[zc - zlen/2.0, zc + zlen/2.0]);
% axis([xc - xlen/2.0, xc + xlen/2.0 , yc - ylen/2.0, yc + ylen/2.0 , zc - zlen/2.0, zc + zlen/2.0])
% axis equal

% --- Executes during object creation, after setting all properties.
function cd_lowlim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cd_lowlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cd_uplim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cd_uplim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function plot_curr_svm(handles)
global learned_svm

if(~isempty(learned_svm))
    displ=[];
    N = size(learned_svm.Sva,1);
    
    v = get(handles.cont_density,'Value');

    disp('Plotting SVM...');
    if(N==2)
        ranges = [get(handles.axes1,'XLim'),get(handles.axes1,'YLim')];
        num_contours = eval(get(handles.cd_lowlim,'String')) + v*(eval(get(handles.cd_uplim,'String'))-eval(get(handles.cd_lowlim,'String')));
        displaySVM(learned_svm,'range', ranges, 'step', eval(get(handles.box_step,'String')),...
                    'num_contours', num_contours, 'sv', get(handles.SV_check,'value'), 'boundary',...
                    eval(get(handles.box_mesh,'String')),  'sv_contour',0);
    else if(N==3)
            %             view(3);
            displ.ranges = [get(handles.axes1,'XLim'),get(handles.axes1,'YLim'),get(handles.axes1,'ZLim')];
            display3DSVM(learned_svm, 'targets',true,'contour_list',eval(get(handles.box_mesh,'String')),'range',...
                displ.ranges,'step', eval(get(handles.box_step,'String')), 'sv',get(handles.SV_check,'value'));
            camlight
            lighting phong
        end
    end
    disp('done.');
%     axis tight
end


% --- Executes on button press in SV_check.
function SV_check_Callback(hObject, eventdata, handles)
% hObject    handle to SV_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SV_check
refreshDisplay(handles);

% --- Executes during object creation, after setting all properties.
function box_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function box_mesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function targetValBox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to targetValBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function refreshDisplay(handles)
global all_data  
% if(size(all_data{1}{1},1) == 2)
%     tmp_range = [get(handles.axes1,'XLim'),get(handles.axes1,'YLim')];
% else
%     tmp_range = [get(handles.axes1,'XLim'),get(handles.axes1,'YLim'),get(handles.axes1,'ZLim')];
% end
% tmp_view=get(handles.axes1,'view');
% axis(handles.axes1);
currlims= [get(handles.axes1,'XLim') get(handles.axes1,'YLim')];
cla
hold on
grid on

if(~isempty(all_data))
    plot_all_data(getActiveIndices(handles.attPanel));
    %     axis(all_data_ranges);
end

plot_curr_svm(handles);
set(handles.axes1, 'XLim',currlims(1:2));
set(handles.axes1, 'YLim',currlims(3:4));

box on

% axis(tmp_range);


% --------------------------------------------------------------------
function export_eps_Callback(hObject, eventdata, handles)
% hObject    handle to export_eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile('*.eps','Export Display') ;
if(FileName == 0)
    return
end
newfig=figure;
h=copyobj(handles.axes1, newfig);
box on
print(newfig,'-dpsc',[PathName '/' FileName]);


% --- Executes on slider movement.
function trim_slider_Callback(hObject, eventdata, handles)
% hObject    handle to trim_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
if( getTargetClass(handles.tclass_list) == 0)
    disp('Choose a target class to smooth first!');
    return;
end

updateData(handles);

cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));

% --- Executes during object creation, after setting all properties.
function trim_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trim_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes during object creation, after setting all properties.
function trim_lowlim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trim_lowlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function trim_uplim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trim_uplim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateData(handles)
global all_data real_data


if(~isempty(real_data))
    
    if(getTargetClass(handles.tclass_list) ~=0)
        
        tar2 = eval(get(handles.targetValBox,'String'))';
        tclass = getTargetClass(handles.tclass_list);
        
        
        dim = size(real_data{tclass}{1},1);
        tar = zeros(dim,1);
        for i = 1:length(real_data{tclass})
            tar = tar + real_data{tclass}{i}(:,end);
        end
        tar = tar/length(real_data{tclass});
        for i = 1:length(real_data{tclass})
            all_data{tclass}{i} = real_data{tclass}{i} + repmat(tar2-tar,1,size(real_data{tclass}{i},2));
        end
        
        if(get(handles.trim_slider,'value') ~=0)
            tr = ceil(get(handles.trim_slider,'value')*eval(get(handles.trim_uplim, 'String')));
            for j=1:length(all_data{tclass})
                all_data{tclass}{j} =  all_data{tclass}{j}(:,tr:end);
            end
        end
        
        sf = eval(get(handles.smooth_lowlim, 'String')) + get(handles.smoothness_slider, 'value')*(eval(get(handles.smooth_uplim, 'String'))-eval(get(handles.smooth_lowlim, 'String')));
        if(get(handles.resample_slider,'value') ~=0)
            nd = eval(get(handles.resample_lowlim, 'String')) + get(handles.resample_slider,'value')*eval(get(handles.resample_uplim, 'String'));
            nd = max(nd,5);
        else
            nd=-1;
        end
        all_data{tclass} = resampleAndSmoothData(all_data{tclass}, sf, nd);
        
        
    end
end


% --- Executes on button press in kboard.
function kboard_Callback(hObject, eventdata, handles)
% hObject    handle to kboard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global learned_svm all_data  

keyboard;


% --- Executes during object creation, after setting all properties.
function box_btol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_btol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function box_maxtiter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_maxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in refresh.
function refresh_Callback(hObject, eventdata, handles)
% hObject    handle to refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global all_data

if(~isempty(all_data))
    updateData(handles);
end
refreshDisplay(handles);

% --- Executes during object creation, after setting all properties.
function tclass_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tclass_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in add_class.
function add_class_Callback(hObject, eventdata, handles)
% hObject    handle to add_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global numAtt all_data real_data all_data_ranges
axes(handles.axes1);

if(numAtt == 8)
   disp('Only upto 8 attractors are supported!!');
   return;
end

numAtt = numAtt + 1;
all_data{numAtt} = {};


ch = get(handles.attPanel,'Children');
for j=1:length(ch)
    if(~isempty(strfind(get(ch(j),'String'), int2str(numAtt))))
        set(ch(j),'Enable','on');
        set(ch(j),'Value',1);
        break;
    end
end



all_data{numAtt} = grabDataFromCursor(3,0.001);



dim = size(all_data{numAtt}{1},1);
tar = zeros(dim,1);
for j=1:length(all_data{numAtt})
    tar = tar + all_data{numAtt}{j}(:,end);
end
tar = tar/length(all_data{numAtt});

for j=1:length(all_data{numAtt})
    all_data{numAtt}{j} = all_data{numAtt}{j} + repmat(tar-all_data{numAtt}{j}(:,end),1,size(all_data{numAtt}{j},2));
end


real_data{numAtt} = all_data{numAtt};
indices = [];
for i=1:length(ch)
    if(get(ch(i),'Value') == 1)
        indices = [indices;eval(get(ch(i),'String'))];
    end
end
plot_all_data(indices);

updateTargetGroup(handles);
updateTargetValBox(handles);




% --- Executes on button press in clr_model.
function clr_model_Callback(hObject, eventdata, handles)
% hObject    handle to clr_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global learned_svm
learned_svm=[];

refreshDisplay(handles);

% --- Executes on button press in clr_data.
function clr_data_Callback(hObject, eventdata, handles)
% hObject    handle to clr_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global all_data real_data numAtt
all_data = [];
real_data=[];
ch1 = get(handles.attPanel,'Children');
for i=1:length(ch1)
    set(ch1(i),'Enable','off');
    set(ch1(i),'Value',0);
end
numAtt = 0;
cla
hold on
grid on
updateTargetGroup(handles);
refreshDisplay(handles);


% --- Executes on selection change in tclass_list.
function tclass_list_Callback(hObject, eventdata, handles)
% hObject    handle to tclass_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tclass_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tclass_list
updateTargetValBox(handles);

function updateTargetValBox(handles)
global all_data

allstr = get(handles.tclass_list,'String');
selectstr = allstr{get(handles.tclass_list,'Value')};

if(~isempty(all_data))
    tclass = eval(selectstr);
    dim = size(all_data{tclass}{1},1);
    tar = zeros(dim,1);
    for i = 1:length(all_data{tclass})
        tar = tar + all_data{tclass}{i}(:,end);
    end
    tar = tar/length(all_data{tclass});
    str='[ ';
    for i=1:dim
        str = [str sprintf('%1.2f ', tar(i))];
    end
    str = [str ']'];
    set(handles.targetValBox,'String', str);
end



function box_btol_Callback(hObject, eventdata, handles)
% hObject    handle to box_btol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_btol as text
%        str2double(get(hObject,'String')) returns contents of box_btol as a double



function box_lambda_Callback(hObject, eventdata, handles)
% hObject    handle to box_lambda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_lambda as text
%        str2double(get(hObject,'String')) returns contents of box_lambda as a double



function box_C_Callback(hObject, eventdata, handles)
% hObject    handle to box_C (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_C as text
%        str2double(get(hObject,'String')) returns contents of box_C as a double



function box_maxiter_Callback(hObject, eventdata, handles)
% hObject    handle to box_maxiter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_maxiter as text
%        str2double(get(hObject,'String')) returns contents of box_maxiter as a double



function box_step_Callback(hObject, eventdata, handles)
% hObject    handle to box_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_step as text
%        str2double(get(hObject,'String')) returns contents of box_step as a double



function box_mesh_Callback(hObject, eventdata, handles)
% hObject    handle to box_mesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_mesh as text
%        str2double(get(hObject,'String')) returns contents of box_mesh as a double



function targetValBox_Callback(hObject, eventdata, handles)
% hObject    handle to targetValBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of targetValBox as text
%        str2double(get(hObject,'String')) returns contents of targetValBox as a double


% --- Executes on button press in simds.
function simds_Callback(hObject, eventdata, handles)
% hObject    handle to simds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global learned_svm all_data

if(isempty(learned_svm))
    disp('Learn a SVM object first!!');
    return;
end

disp('Preparing nominal regression model');
sigma = eval(get(handles.box_lambda,'String'));
dim = size(all_data{1}{1},1);
tc = getTargetClass(handles.tclass_list);
rgs = [get(handles.axes1, 'XLim'),get(handles.axes1, 'YLim')];
stp = eval(get(handles.box_step,'String'));
%% Using SVR
% lambda = 1/(2*sigma*sigma)
% tc = getTargetClass(handles.tclass_list);
% lsvmdata=[];
% for i=1:length(all_data{tc})
%     lsvmdata = [lsvmdata ; [all_data{tc}{i}(:,1:end-1)', 20*[diff(all_data{tc}{i},1,2)']]];
% end
% stroptions = ['-s 3 -t 2 -g ' num2str(15) ' -c 1000 -e 1e-4'];
% 
% nominal_models = cell(dim,1);
% for i=1:dim
% nominal_models{i} = svmtrain(lsvmdata(:,i+dim), lsvmdata(:,1:dim) , stroptions);
% end
% disp('Done');

% [tempx, tempy] = meshgrid(rgs(1):stp:rgs(2),rgs(3):stp:rgs(4));
% tempz=zeros(size(tempx));
% for i=1:size(tempx,1)
%     for j=1:size(tempx,2)
%         pt = [tempx(i,j);tempy(i,j)];
%         tempz(i,j) = 
%     end
% end

%% Using EM
%initializing with EM

gmrdata=[];
for i=1:length(all_data{tc})
    gmrdata = [gmrdata , [all_data{tc}{i}(:,1:end-1); [diff(all_data{tc}{i},1,2)]]];
end
[Priors_0, Mu_0, Sigma_0] = EM_init_kmeans(gmrdata, 6);
[Priors_0, Mu_0, Sigma_0] = EM(gmrdata, Priors_0, Mu_0, Sigma_0);

[tempx, tempy] = meshgrid(rgs(1):stp:rgs(2),rgs(3):stp:rgs(4));
u=zeros(size(tempx));
v=zeros(size(tempx));
for i=1:size(tempx,1)
    for j=1:size(tempx,2)
        pt = [tempx(i,j);tempy(i,j)];
        if(calculateClassifier(learned_svm, pt) >= 0)
            tmp = GMR(Priors_0, Mu_0, Sigma_0, pt, 1:2,3:4);
            tmp = getModulatedVelocity(learned_svm, tmp, pt);
            u(i,j) = tmp(1);
            v(i,j) = tmp(2);
        end
    end
end
disp('Done');
refreshDisplay(handles);
h=streamslice(tempx, tempy, u,v,2);
set(h, 'Linewidth',2);



function box_nu_Callback(hObject, eventdata, handles)
% hObject    handle to box_nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_nu as text
%        str2double(get(hObject,'String')) returns contents of box_nu as a double


% --- Executes during object creation, after setting all properties.
function box_nu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateData(handles);

cla
hold on
grid on
plot_all_data(getActiveIndices(handles.attPanel));


% --- Executes on button press in chk_Auto.
function chk_Auto_Callback(hObject, eventdata, handles)
% hObject    handle to chk_Auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chk_Auto
if(get(hObject,'Value'))
    set(handles.box_lambda, 'Enable','off');
else
    set(handles.box_lambda, 'Enable','on');
end



function box_eps_Callback(hObject, eventdata, handles)
% hObject    handle to box_eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of box_eps as text
%        str2double(get(hObject,'String')) returns contents of box_eps as a double


% --- Executes during object creation, after setting all properties.
function box_eps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to box_eps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
