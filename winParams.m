function varargout = winParams(varargin)
% WINPARAMS MATLAB code for winParams.fig
%      WINPARAMS, by itself, creates a new WINPARAMS or raises the existing
%      singleton*.
%
%      H = WINPARAMS returns the handle to a new WINPARAMS or the handle to
%      the existing singleton*.
%
%      WINPARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINPARAMS.M with the given input arguments.
%
%      WINPARAMS('Property','Value',...) creates a new WINPARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before winParams_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to winParams_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help winParams

% Last Modified by GUIDE v2.5 22-Jan-2018 16:28:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @winParams_OpeningFcn, ...
                   'gui_OutputFcn',  @winParams_OutputFcn, ...
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

% --- All my functions are here --- %

function handles = set_fields(hObject, handles, winParams)
% update fields according to given winParams
handles.winParams = winParams;
handles.fieldWinSize.String = num2str(handles.winParams.size);
handles.fieldWinCenter.String = num2str(handles.winParams.center);
handles.checkStart.Value = handles.winParams.by_start;
handles.checkStop.Value = handles.winParams.by_stop;
handles.checkBlockEnd.Value = handles.winParams.truncate_seq;
guidata(hObject, handles);


function handles = update_figure(hObject, handles)
% update winParams and axPreview according to fields

% parse params
handles.winParams.by_start = handles.checkStart.Value;
handles.winParams.by_stop = handles.checkStop.Value;
handles.winParams.truncate_seq = handles.checkBlockEnd.Value;
tmp = str2double(handles.fieldWinSize.String);
if isnan(tmp)
    errordlg('invalid window size value', 'window size error', 'modal');
    handles.fieldWinSize.String = num2str(handles.winParams.size);
else
    handles.winParams.size = round(tmp);
end
tmp = str2double(handles.fieldWinCenter.String);
if isnan(tmp)
    errordlg('invalid window center value', 'window center error', 'modal');
    handles.fieldWinCenter.String = num2str(handles.winParams.center);
else
    handles.winParams.center = round(tmp);
end
guidata(hObject, handles);

% preview
draw_win(handles.axPreview, ...
         handles.winParams.size, ...
         handles.winParams.center, ...
         handles.winParams.by_start, handles.winParams.by_stop, ...
         handles.winParams.truncate_seq);

% --- End of my functions --- %

% --- Executes just before winParams is made visible.
function winParams_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to winParams (see VARARGIN)

if isempty(varargin)
    handles.winParams = struct;
else
    handles = set_fields(hObject, handles, varargin{1});
end
handles = update_figure(hObject, handles);

% Choose default command line output for winParams
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes winParams wait for user response (see UIRESUME)
uiwait(handles.figParams);


% --- Outputs from this function are returned to the command line.
function varargout = winParams_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.winParams;

delete(handles.figParams);


% --- Executes on button press in checkBlockStart.
function checkBlockStart_Callback(hObject, eventdata, handles)
% hObject    handle to checkBlockStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkBlockStart


% --- Executes on button press in checkBlockEnd.
function checkBlockEnd_Callback(hObject, eventdata, handles)
% hObject    handle to checkBlockEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkBlockEnd
update_figure(hObject, handles);

% --- Executes on button press in checkStart.
function checkStart_Callback(hObject, eventdata, handles)
% hObject    handle to checkStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkStart
if ~handles.winParams.by_stop
    handles.winParams.by_stop = 1;
    handles.checkStop.Value = 1;
end
update_figure(hObject, handles);


% --- Executes on button press in checkStop.
function checkStop_Callback(hObject, eventdata, handles)
% hObject    handle to checkStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkStop
if ~handles.winParams.by_start
    handles.winParams.by_start = 1;
    handles.checkStart.Value = 1;
end
update_figure(hObject, handles);


function fieldWinSize_Callback(hObject, eventdata, handles)
% hObject    handle to fieldWinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldWinSize as text
%        str2double(get(hObject,'String')) returns contents of fieldWinSize as a doubl
update_figure(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fieldWinSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldWinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fieldWinCenter_Callback(hObject, eventdata, handles)
% hObject    handle to fieldWinCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldWinCenter as text
%        str2double(get(hObject,'String')) returns contents of fieldWinCenter as a double
update_figure(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fieldWinCenter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldWinCenter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butOK.
function butOK_Callback(hObject, eventdata, handles)
% hObject    handle to butOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figParams);


% --- Executes when user attempts to close figParams.
function figParams_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end
