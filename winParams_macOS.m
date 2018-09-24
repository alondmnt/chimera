function varargout = winParams_macOS(varargin)
% WINPARAMS_MACOS MATLAB code for winParams_macOS.fig
%      WINPARAMS_MACOS, by itself, creates a new WINPARAMS_MACOS or raises the existing
%      singleton*.
%
%      H = WINPARAMS_MACOS returns the handle to a new WINPARAMS_MACOS or the handle to
%      the existing singleton*.
%
%      WINPARAMS_MACOS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WINPARAMS_MACOS.M with the given input arguments.
%
%      WINPARAMS_MACOS('Property','Value',...) creates a new WINPARAMS_MACOS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before winParams_macOS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to winParams_macOS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help winParams_macOS

% Last Modified by GUIDE v2.5 25-Aug-2018 19:48:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @winParams_macOS_OpeningFcn, ...
                   'gui_OutputFcn',  @winParams_macOS_OutputFcn, ...
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
% update fields according to given winParams_macOS
handles.winParams = winParams;
handles.fieldWinSize.String = num2str(handles.winParams.size);
handles.fieldWinCenter.String = num2str(handles.winParams.center);
handles.checkStart.Value = handles.winParams.by_start;
handles.checkStop.Value = handles.winParams.by_stop;
handles.checkBlockEnd.Value = handles.winParams.truncate_seq;
if isfinite(handles.winParams.max_len) && (handles.winParams.max_len > 0)
    handles.checkHomolog1.Value = 1;
    handles.fieldMaxLen.String = num2str(handles.winParams.max_len);
else
    handles.checkHomolog1.Value = 0;
end
if isfinite(handles.winParams.max_pos) && (handles.winParams.max_pos > 0)
    handles.checkHomolog2.Value = 1;
    handles.fieldMaxPos.String = num2str(handles.winParams.max_pos);
else
    handles.checkHomolog2.Value = 0;
end

guidata(hObject, handles);


function handles = update_figure(hObject, handles)
% update winParams_macOS and axPreview according to fields

% parse params
handles.winParams.by_start = handles.checkStart.Value;
handles.winParams.by_stop = handles.checkStop.Value;
handles.winParams.truncate_seq = handles.checkBlockEnd.Value;
tmp = str2double(handles.fieldWinSize.String);
if ~isfinite(tmp) || (tmp < 0)
    errordlg(sprintf('invalid window size value: %g', tmp), 'window size error', 'modal');
else
    handles.winParams.size = round(tmp);
end
handles.fieldWinSize.String = num2str(handles.winParams.size);

tmp = str2double(handles.fieldWinCenter.String);
if ~isfinite(tmp)
    errordlg(sprintf('invalid window center value: %g', tmp), 'window center error', 'modal');
else
    handles.winParams.center = round(tmp);
end
handles.fieldWinCenter.String = num2str(handles.winParams.center);

tmp = str2double(handles.fieldMaxLen.String);
if handles.checkHomolog1.Value && (~isfinite(tmp) || tmp <= 0)
    errordlg(sprintf('invalid max len value: %g', tmp), 'max block length error', 'modal');
elseif ~handles.checkHomolog1.Value
    handles.winParams.max_len = NaN;
else
    handles.winParams.max_len = ceil(tmp);
end
handles.fieldMaxLen.String = num2str(handles.winParams.max_len);

tmp = str2double(handles.fieldMaxPos.String);
if handles.checkHomolog2.Value && (~isfinite(tmp) || tmp > 1 || tmp <= 0)
    errordlg(sprintf('invalid max fraction value: %g', tmp), 'max fraction error', 'modal');
elseif ~handles.checkHomolog2.Value
    handles.winParams.max_pos = NaN;
else
    handles.winParams.max_pos = tmp;
end
handles.fieldMaxPos.String = num2str(handles.winParams.max_pos);

guidata(hObject, handles);

% preview
draw_win(handles.axPreview, ...
         handles.winParams.size, ...
         handles.winParams.center, ...
         handles.winParams.by_start, handles.winParams.by_stop, ...
         handles.winParams.truncate_seq);

% --- End of my functions --- %

% --- Executes just before winParams_macOS is made visible.
function winParams_macOS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to winParams_macOS (see VARARGIN)

if isempty(varargin)
    handles.winParams = struct;
else
    handles = set_fields(hObject, handles, varargin{1});
end
handles = update_figure(hObject, handles);

% Choose default command line output for winParams_macOS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes winParams_macOS wait for user response (see UIRESUME)
uiwait(handles.figParams);


% --- Outputs from this function are returned to the command line.
function varargout = winParams_macOS_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in checkHomolog1.
function checkHomolog1_Callback(hObject, eventdata, handles)
% hObject    handle to checkHomolog1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkHomolog1
if handles.checkHomolog1.Value
    handles.fieldMaxLen.String = '40';
end
update_figure(hObject, handles);


function fieldMaxLen_Callback(hObject, eventdata, handles)
% hObject    handle to fieldMaxLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldMaxLen as text
%        str2double(get(hObject,'String')) returns contents of fieldMaxLen as a double
update_figure(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fieldMaxLen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldMaxLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkHomolog2.
function checkHomolog2_Callback(hObject, eventdata, handles)
% hObject    handle to checkHomolog2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkHomolog2
if handles.checkHomolog2.Value
    handles.fieldMaxPos.String = '0.5';
end
update_figure(hObject, handles);


function fieldMaxPos_Callback(hObject, eventdata, handles)
% hObject    handle to fieldMaxPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldMaxPos as text
%        str2double(get(hObject,'String')) returns contents of fieldMaxPos as a double
update_figure(hObject, handles);


% --- Executes during object creation, after setting all properties.
function fieldMaxPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldMaxPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
