function varargout = neural_constraint_gui(varargin)
% NEURAL_CONSTRAINT_GUI M-file for neural_constraint_gui.fig
%      NEURAL_CONSTRAINT_GUI, by itself, creates a new NEURAL_CONSTRAINT_GUI or raises the existing
%      singleton*.
%
%      H = NEURAL_CONSTRAINT_GUI returns the handle to a new NEURAL_CONSTRAINT_GUI or the handle to
%      the existing singleton*.
%
%      NEURAL_CONSTRAINT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURAL_CONSTRAINT_GUI.M with the given input arguments.
%
%      NEURAL_CONSTRAINT_GUI('Property','Value',...) creates a new NEURAL_CONSTRAINT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before neural_constraint_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to neural_constraint_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help neural_constraint_gui

% Last Modified by GUIDE v2.5 18-Aug-2011 13:05:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @neural_constraint_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @neural_constraint_gui_OutputFcn, ...
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


% --- Executes just before neural_constraint_gui is made visible.
function neural_constraint_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to neural_constraint_gui (see VARARGIN)

% Choose default command line output for neural_constraint_gui
set(handles.OK_button,'UserData',2);
handles.output = hObject;
set(handles.algorithm_panel,'SelectionChangeFcn',@algorithm_panel_SelectionChangeFcn);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes neural_constraint_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = neural_constraint_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes on button press in OK_button.
function OK_button_Callback(hObject, eventdata, handles)
% hObject    handle to OK_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
if (get(handles.OK_button,'UserData') ==2.0);                   % If it is the "start" button that is clicked...
    
    %------------------Get Input Values from Input Boxes-------------------
    N = str2num(get(handles.No_pattern_neurons,'String'));      % N is the number of pattern nodes.
    N_const = str2num(get(handles.No_constraints,'String'));    % N_const is the number of constraint nodes.
    deg_row_H = str2num(get(handles.deg_row,'String'));         % deg_row_H is the row degree of the constraint matrix.
    deg_column_H = str2num(get(handles.deg_column,'String'));   % deg_column_H is the column degree of the constraint matrix.
    K = N-N_const;                                              % K is the number of the message nodes.
    
    
    convergence_flag = (get(handles.conv_flag,'Value'));        % This flag determines if the algorithm waits for convergence to stop or performs a pre-determined number of iterations.
    max_itr = str2num(get(handles.maximum_itr,'String'));       % max_itr is the pre-determined number of iterations in case this option was chosen.
    
    max_simulated_instances=str2num(get(handles.max_instances,'String'));   % Determines how many instances should be simulated for each network.    
    no_ensemble_generated = str2num(get(handles.no_ensembles,'String'));    % Determines how many different networks should be simulated.
    
    algorithm_option = (get(handles.bit_flip,'Value'));         % Determines if winner-take-all or bit-flipping should be used.
    gamma = str2num(get(handles.gamma_box,'String'));           % gamma is the fraction of neighbors that trigger an update in neurons' states in the bit-flipping algorithm.
        
    y_max = str2num(get(handles.max_fire,'String'));            % y_max is the maximum firing rate of a neuron in the model.
    y_min = str2num(get(handles.min_fire,'String'));            % y_min is the minimum firing rate of a neuron in the model.

    error_bits = str2num(get(handles.no_err_bits,'String'));    % Stores the number of initial erroneous nodes.
    max_noise_amp = str2num(get(handles.noise_amp,'String'));   % Determines the maximum absolute value of noise for each erroneous node.
    %----------------------------------------------------------------------    
    
    
    %-------------Check the Validity of the Input Arguments----------------
    error_flag = 0;                 
    if (N*deg_column_H ~= N_const*deg_row_H)
        h = msgbox('Invalid input arguments. Number of pattern nodes times column degree must be equal to the number of constraint nodes times row degree.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (N <= 0)
        h = msgbox('The number of pattern nodes must be positive.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (N_const <= 0)
        h = msgbox('The number of constraint nodes must be positive.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (N <= N_const)
        h = msgbox('The number of pattern nodes must be bigger than the number of constraint nodes.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (deg_row_H <= 0)
        h = msgbox('The row degree must be positive.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (deg_column_H <= 0)
        h = msgbox('The column degree must be positive.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (max_simulated_instances <= 0)
        h = msgbox('The maximum number of simulated instances must be positive.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (no_ensemble_generated <= 0)
        h = msgbox('The number of generated ensembles must be positive.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (gamma <= 0)
        h = msgbox('gamma must be positive.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (y_max <=y_min)
        h = msgbox('The Maximum neural firing rate must be greater than the minimum neural firing rate.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (max_noise_amp <= 0)
        h = msgbox('The maximum noise amplitude must be positive.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    
    if (min(error_bits)<0)
        h = msgbox('The initial number of erroneous nodes must be non-negative.','Error','error');                            
        uiwait(h);
        error_flag = 1;
    end
    %----------------------------------------------------------------------    
    
    if (error_flag == 0)
        set(handles.OK_button,'UserData',1);
        initialization_done = 1;
        pause_flag = 0;
        
        %-----------Make Invisible Buttons and Text Boxes Visible----------
        set(handles.Pause_button,'Visible','on');
        set(handles.text29,'Visible','on');
        set(handles.current_ensemble,'Visible','on');
        set(handles.text31,'Visible','on');
        set(handles.cur_simul_itr,'Visible','on');
        set(handles.curr_error_rate,'Visible','on');
        set(handles.text30,'Visible','on');
        set(handles.initial_errors,'Visible','on');
        set(handles.text32,'Visible','on');
        set(handles.uipanel11,'Visible','on');
        set(handles.OK_button,'String','Stop');
        %------------------------------------------------------------------
        
        %------------------------Initialize Displays-----------------------
        set(handles.current_ensemble,'String',num2str(1));                                           
        set(handles.initial_errors,'String',num2str(error_bits(1)));                                   
        set(handles.cur_simul_itr,'String',num2str(1));                                    
        set(handles.curr_error_rate,'String',num2str(0));            
        %------------------------------------------------------------------
    
        pause(1)
        run Neural_constraint;              % Execute the neural network code.
    end
else    
    
    set(handles.OK_button,'UserData',2);        %Reset the UserData value
    
    %-----------------------Reset the Displays-----------------------------
    set(handles.Pause_button,'Visible','off');
    set(handles.text29,'Visible','on');
    set(handles.current_ensemble,'Visible','off');
    set(handles.text31,'Visible','off');
    set(handles.cur_simul_itr,'Visible','off');
    set(handles.curr_error_rate,'Visible','off');
    set(handles.text30,'Visible','off');
    set(handles.initial_errors,'Visible','off');
    set(handles.uipanel11,'Visible','off');
    set(handles.text32,'Visible','off');    
    set(handles.OK_button,'String','Start');    
    set(handles.current_ensemble,'String','');                                          
    set(handles.initial_errors,'String','');                    
    set(handles.cur_simul_itr,'String','');                    
    set(handles.curr_error_rate,'String','');
     set(handles.Pause_button,'String','Pause');
    %----------------------------------------------------------------------
end

function No_pattern_neurons_Callback(hObject, eventdata, handles)
% hObject    handle to No_pattern_neurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of No_pattern_neurons as text
%        str2double(get(hObject,'String')) returns contents of No_pattern_neurons as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function No_pattern_neurons_CreateFcn(hObject, eventdata, handles)
% hObject    handle to No_pattern_neurons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function No_constraints_Callback(hObject, eventdata, handles)
% hObject    handle to No_constraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of No_constraints as text
%        str2double(get(hObject,'String')) returns contents of No_constraints as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function No_constraints_CreateFcn(hObject, eventdata, handles)
% hObject    handle to No_constraints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function deg_row_Callback(hObject, eventdata, handles)
% hObject    handle to deg_row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deg_row as text
%        str2double(get(hObject,'String')) returns contents of deg_row as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function deg_row_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deg_row (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function deg_column_Callback(hObject, eventdata, handles)
% hObject    handle to deg_column (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of deg_column as text
%        str2double(get(hObject,'String')) returns contents of deg_column as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function deg_column_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deg_column (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_instances_Callback(hObject, eventdata, handles)
% hObject    handle to max_instances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_instances as text
%        str2double(get(hObject,'String')) returns contents of max_instances as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function max_instances_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_instances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function no_ensembles_Callback(hObject, eventdata, handles)
% hObject    handle to no_ensembles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of no_ensembles as text
%        str2double(get(hObject,'String')) returns contents of no_ensembles as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function no_ensembles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to no_ensembles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in winner.

function winner_Callback(hObject, eventdata, handles)
% hObject    handle to winner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of winner
111

% --- Executes on button press in conv_flag.
function conv_flag_Callback(hObject, eventdata, handles)
% hObject    handle to conv_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = get(handles.conv_flag,'Value');
if (val == 0)
    set(handles.text20,'Visible','on');
    set(handles.maximum_itr,'Visible','on');
else
    set(handles.text20,'Visible','off');
    set(handles.maximum_itr,'Visible','off');
end

% Hint: get(hObject,'Value') returns toggle state of conv_flag



function maximum_itr_Callback(hObject, eventdata, handles)
% hObject    handle to maximum_itr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maximum_itr as text
%        str2double(get(hObject,'String')) returns contents of maximum_itr as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function maximum_itr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maximum_itr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over conv_flag.
function conv_flag_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to conv_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text20,'Visible','on');



function min_fire_Callback(hObject, eventdata, handles)
% hObject    handle to min_fire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_fire as text
%        str2double(get(hObject,'String')) returns contents of min_fire as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function min_fire_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_fire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_fire_Callback(hObject, eventdata, handles)
% hObject    handle to max_fire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_fire as text
%        str2double(get(hObject,'String')) returns contents of max_fire as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function max_fire_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_fire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gamma_box_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma_box as text
%        str2double(get(hObject,'String')) returns contents of gamma_box as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gamma_box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function algorithm_panel_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to algorithm_panel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function noise_amp_Callback(hObject, eventdata, handles)
% hObject    handle to noise_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noise_amp as text
%        str2double(get(hObject,'String')) returns contents of noise_amp as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function noise_amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noise_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function no_err_bits_Callback(hObject, eventdata, handles)
% hObject    handle to no_err_bits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of no_err_bits as text
%        str2double(get(hObject,'String')) returns contents of no_err_bits as a double
input = str2num(get(hObject,'String'));
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function no_err_bits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to no_err_bits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function algorithm_panel_SelectionChangeFcn(hObject, eventdata)
 
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'winner'
      %execute this code when fontsize08_radiobutton is selected
    set(handles.text23,'Visible','off');
    set(handles.gamma_box,'Visible','off');
 
    case 'bit_flip'
      %execute this code when fontsize12_radiobutton is selected
    set(handles.text23,'Visible','on');
    set(handles.gamma_box,'Visible','on');

    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);



% --- Executes on button press in Pause_button.
function Pause_button_Callback(hObject, eventdata, handles)
% hObject    handle to Pause_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (get(handles.OK_button,'UserData') ==1)
    set(handles.OK_button,'UserData',0);
    guidata(hObject, handles);    
else
    set(handles.OK_button,'UserData',1);
    guidata(hObject, handles);       
end
% close all;
