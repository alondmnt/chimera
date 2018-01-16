function varargout = chimeraGUI(varargin)
% CHIMERAGUI MATLAB code for chimeraGUI.fig
%      CHIMERAGUI, by itself, creates a new CHIMERAGUI or raises the existing
%      singleton*.
%
%      H = CHIMERAGUI returns the handle to a new CHIMERAGUI or the handle to
%      the existing singleton*.
%
%      CHIMERAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHIMERAGUI.M with the given input arguments.
%
%      CHIMERAGUI('Property','Value',...) creates a new CHIMERAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before chimeraGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to chimeraGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help chimeraGUI

% Last Modified by GUIDE v2.5 16-Jan-2018 17:12:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chimeraGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @chimeraGUI_OutputFcn, ...
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
function [region, err] = get_region(fieldFrom, fieldTo, L)
err = 0;
region = [int64(str2double(fieldFrom.String)), ...
          int64(str2double(fieldTo.String))];
if any(~isfinite(region) | region == 0)
    errordlg(sprintf(['expecting [from] and [to] to be signed integers. \n', ...
                      '>= 1 coordinates measured from START. \n', ...
                      '<= -1 coordinates measured from STOP.']));
    err = 1;
    return
end

if any(region < 0) && ~isfinite(L)
    errordlg('negative coordinates are undefined when target protein is unknown.')
    err = 3;
    return
end
region(region < 0) = L + 1 + region(region < 0);

if region(1) > region(2)
    errordlg('expecting [to] >= [from].');
    err = 2;
    return
end

if ~isfinite(L) && any(region > L)
    errordlg(sprint('new region exceeds target protein length (%d > %d).', ...
                    region(find(region > L, 1)), L))
    err = 4;
    return
end


function [err, merge, trunc_new, trunc_old] = check_region(new_reg, regions)
new_reg = region2set(new_reg);
regions = region2set(regions);

merge = set2region(union(new_reg, regions));
ovlp = intersect(new_reg, regions);

if isempty(ovlp)
    err = 0;
    trunc_new = new_reg;
    trunc_old = regions;
    return
elseif length(ovlp) == length(new_reg)
    % new region is a subset of regions
    err = 2;
else
    err = 1;
end
trunc_new = set2region(setdiff(new_reg, regions));
trunc_old = set2region(setdiff(regions, new_reg));


function x = cellcat(x, dim)
if nargin < 1
    dim = 1;
end
if ~isempty(x)
    x = cat(dim, x{:});
else
    x = [];
end


function set = region2set(region)
% [region] is a 2-column matrix with start/end coordinates
% [set] is a sorted array of a union of coordinates
if isempty(region)
    set = [];
    return
end
set = cellcat(arrayfun(@(x, y) {x:y}, region(:, 1), region(:, 2)), 2);
n = length(set);
set = unique(set);
assert(length(set) == n, 'overlapping regions')


function region = set2region(set)
% [set] is an array of coordinates
% [region] is a 2-column matrix with start/end coordinates
if isempty(set)
    region = zeros(0, 2);
    return
elseif length(set) == 1
    region = [set, set];
    return
end
set = sort(set);
ireg = [0, find(diff(set) > 1), length(set)];
region = cellcat(arrayfun(@(x, y) {[set(x+1), set(y)]}, ireg(1:end-1), ireg(2:end)), 1);


function update_regions(hObject, handles)
% update all list boxes and barplot.
chim_region = region2set(handles.chimera_regions);
codon_region = region2set(handles.codon_regions);
assert(isempty(intersect(chim_region, codon_region)), 'codon/chimera overlap');

if isfinite(handles.target_len)
    L = handles.target_len;
else
    L = max([chim_region, codon_region]);
end
handles.default_regions = set2region(setdiff(setdiff(1:L, chim_region), codon_region));
if ~handles.default_exist && ~isempty(handles.default_regions)
    handles.statDefault.String = 'required';
elseif ~handles.default_exist
    handles.statDefault.String = 'optional';
end

[handles.listChimera.String, handles.listChimera.Value] = update_listbox(...
	handles.chimera_regions, handles.target_len, handles.fieldChimeraFrom, handles.fieldChimeraTo);
[handles.listCodon.String, handles.listCodon.Value] = update_listbox(...
    handles.codon_regions, handles.target_len, handles.fieldCodonFrom, handles.fieldCodonTo);
handles.listDefault.Value = min(handles.listDefault.Value, size(handles.default_regions, 1));
handles.listDefault.String = update_listbox(handles.default_regions, handles.target_len);

update_bar(handles);

guidata(hObject, handles);


function [String, Value] = update_listbox(regions, L, fieldFrom, fieldTo)
if isempty(regions)
    String = '';
    Value = [];
    return
end

String = sprintf('from %d to %d\n', regions');
String = String(1:end-1);

if nargout == 1
    return
end

try
    reg = get_region(fieldFrom, fieldTo, L);
catch
    Value = [];
    return
end
Value = [find(regions(:, 1) == reg(:, 1)); find(regions(:, 2) == reg(:, 2))];


function update_bar(handles)
regions = {handles.chimera_regions, handles.codon_regions, handles.default_regions};
lens = cumsum(cellfun(@(x) size(x, 1), regions));

[regions, categ] = sortrows(cellcat(regions, 1));
regions = diff(regions, 1, 2) + 1;

categ(categ <= lens(1)) = 1;
categ(lens(1) < categ & categ <= lens(2)) = 2;
categ(lens(2) < categ & categ <= lens(3)) = 3;

cmap = colormap;
b = barh([regions, regions]', 1, 'stacked', 'EdgeColor', 'none', ...
         'ButtonDownFcn', @select_region_from_figure);
[b(categ == 1).FaceColor] = deal([1, 0, 0]);
[b(categ == 2).FaceColor] = deal([0, 1, 0]);
[b(categ == 3).FaceColor] = deal([0, 0, 1]);

handles.axPreview.XLim = [0, sum(regions)+eps];
if lens(end) == 0
    handles.axPreview.XTick = [];
end
handles.axPreview.YTick = [];
handles.axPreview.YLim = [0.5, 1.5];


function select_region_from_figure(hObject, eventdata)

handles = guidata(hObject.Parent.Parent);
pos = ceil(eventdata.IntersectionPoint(1));

ichim = handles.chimera_regions(:, 1) <= pos & ...
        pos <= handles.chimera_regions(:, 2);
if any(ichim)
    handles.listChimera.Value = find(ichim);
    handles.fieldChimeraFrom.String = sprintf('%d', handles.chimera_regions(ichim, 1));
    handles.fieldChimeraTo.String = sprintf('%d', handles.chimera_regions(ichim, 2));
end

icod = handles.codon_regions(:, 1) <= pos & ...
       pos <= handles.codon_regions(:, 2);
if any(icod)
    handles.listCodon.Value = find(icod);
    handles.fieldCodonFrom.String = sprintf('%d', handles.codon_regions(icod, 1));
    handles.fieldCodonTo.String = sprintf('%d', handles.codon_regions(icod, 2));
end

idef = handles.default_regions(:, 1) <= pos & ...
       pos <= handles.default_regions(:, 2);
if any(idef)
    handles.listDefault.Value = find(idef);
    handles.fieldDefaultFrom.String = sprintf('%d', handles.default_regions(idef, 1));
    handles.fieldDefaultTo.String = sprintf('%d', handles.default_regions(idef, 2));
end

guidata(hObject, handles);
% --- End of my functions --- %


% --- Executes just before chimeraGUI is made visible.
function chimeraGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chimeraGUI (see VARARGIN)

handles.axPreview.Box = 'on';

handles.target_len = NaN;
handles.target_exist = false;
handles.default_exist = false;
handles.reference_exist = false;

handles.chimera_regions = zeros(0, 2);
handles.codon_regions = zeros(0, 2);
update_regions(hObject, handles);

% Choose default command line output for chimeraGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chimeraGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = chimeraGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listChimera.
function listChimera_Callback(hObject, eventdata, handles)
% hObject    handle to listChimera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listChimera contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listChimera
for reg = handles.listChimera.Value
    handles.fieldChimeraFrom.String = sprintf('%d', handles.chimera_regions(reg, 1));
    handles.fieldChimeraTo.String = sprintf('%d', handles.chimera_regions(reg, 2));
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listChimera_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listChimera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butChimeraAddRegion.
function butChimeraAddRegion_Callback(hObject, eventdata, handles)
% hObject    handle to butChimeraAddRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_reg, err] = get_region(handles.fieldChimeraFrom, handles.fieldChimeraTo, ...
                            handles.target_len);
if err
    return
end

err = check_region(new_reg, handles.chimera_regions);
if err == 2
    errordlg('new region already exists');
    return
elseif err == 1
    decision = questdlg('new region overlaps with an existing Chimera region' , ...
                        'region overlap', 'merge', 'cancel', 'cancel');
    if strcmp(decision, 'cancel')
        return
    end
end

[err, ~, trunc_new, trunc_old] = check_region(new_reg, handles.codon_regions);
if err
    decision = questdlg('new region overlaps with an existing codon region', ...
                        'region overlap', 'keep chimera', 'keep codon', 'keep chimera');
    switch decision
        case 'keep codon'
            new_reg = trunc_new;
        case 'keep chimera'
            handles.codon_regions = trunc_old;
        otherwise
            error('too many cooks!');
    end
end

[~, handles.chimera_regions] = check_region(new_reg, handles.chimera_regions);

update_regions(hObject, handles);


% --- Executes on button press in butRun.
function butRun_Callback(hObject, eventdata, handles)
% hObject    handle to butRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function fileTarget_Callback(hObject, eventdata, handles)
% hObject    handle to fileTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileTarget as text
%        str2double(get(hObject,'String')) returns contents of fileTarget as a double


% --- Executes during object creation, after setting all properties.
function fileTarget_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butTarget.
function butTarget_Callback(hObject, eventdata, handles)
% hObject    handle to butTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function fileReference_Callback(hObject, eventdata, handles)
% hObject    handle to fileReference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileReference as text
%        str2double(get(hObject,'String')) returns contents of fileReference as a double


% --- Executes during object creation, after setting all properties.
function fileReference_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileReference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butReference.
function butReference_Callback(hObject, eventdata, handles)
% hObject    handle to butReference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listCodon.
function listCodon_Callback(hObject, eventdata, handles)
% hObject    handle to listCodon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listCodon contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listCodon
for reg = handles.listCodon.Value
    handles.fieldCodonFrom.String = sprintf('%d', handles.codon_regions(reg, 1));
    handles.fieldCodonTo.String = sprintf('%d', handles.codon_regions(reg, 2));
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listCodon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listCodon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in butChimeraRemRegion.
function butChimeraRemRegion_Callback(hObject, eventdata, handles)
% hObject    handle to butChimeraRemRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_reg, err] = get_region(handles.fieldChimeraFrom, handles.fieldChimeraTo, ...
                            handles.target_len);
if err
    return
end

[err, ~, ~, handles.chimera_regions] = check_region(new_reg, handles.chimera_regions);
if err == 0
    errordlg('region does not exist');
    return
end

update_regions(hObject, handles);


function fieldChimeraFrom_Callback(hObject, eventdata, handles)
% hObject    handle to fieldChimeraFrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldChimeraFrom as text
%        str2double(get(hObject,'String')) returns contents of fieldChimeraFrom as a double

% --- Executes during object creation, after setting all properties.
function fieldChimeraFrom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldChimeraFrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fieldChimeraTo_Callback(hObject, eventdata, handles)
% hObject    handle to fieldChimeraTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldChimeraTo as text
%        str2double(get(hObject,'String')) returns contents of fieldChimeraTo as a double


% --- Executes during object creation, after setting all properties.
function fieldChimeraTo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldChimeraTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listDefault.
function listDefault_Callback(hObject, eventdata, handles)
% hObject    handle to listDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listDefault contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listDefault


% --- Executes during object creation, after setting all properties.
function listDefault_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fileDefault_Callback(hObject, eventdata, handles)
% hObject    handle to fileDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileDefault as text
%        str2double(get(hObject,'String')) returns contents of fileDefault as a double


% --- Executes during object creation, after setting all properties.
function fileDefault_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butDefault.
function butDefault_Callback(hObject, eventdata, handles)
% hObject    handle to butDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in optSourceRef.
function optSourceRef_Callback(hObject, eventdata, handles)
% hObject    handle to optSourceRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optSourceRef


% --- Executes on button press in optSource.
function optSource_Callback(hObject, eventdata, handles)
% hObject    handle to optSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optSource


% --- Executes on button press in butCodonAddRegion.
function butCodonAddRegion_Callback(hObject, eventdata, handles)
% hObject    handle to butCodonAddRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_reg, err] = get_region(handles.fieldCodonFrom, handles.fieldCodonTo, ...
                            handles.target_len);
if err
    return
end

err = check_region(new_reg, handles.codon_regions);
if err == 2
    errordlg('new region already exists');
    return
elseif err == 1
    decision = questdlg('new region overlaps with an existing codon region' , ...
                        'region overlap', 'merge', 'cancel', 'cancel');
    if strcmp(decision, 'cancel')
        return
    end
end

[err, ~, trunc_new, trunc_old] = check_region(new_reg, handles.chimera_regions);
if err
    decision = questdlg('new region overlaps with an existing codon region', ...
                        'region overlap', 'keep chimera', 'keep codon', 'keep codon');
    switch decision
        case 'keep chimera'
            new_reg = trunc_new;
        case 'keep codon'
            handles.chimera_regions = trunc_old;
        otherwise
            error('too many cooks!');
    end
end

[~, handles.codon_regions] = check_region(new_reg, handles.codon_regions);

update_regions(hObject, handles);


% --- Executes on button press in butCodonRemRegion.
function butCodonRemRegion_Callback(hObject, eventdata, handles)
% hObject    handle to butCodonRemRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_reg, err] = get_region(handles.fieldCodonFrom, handles.fieldCodonTo, ...
                            handles.target_len);
if err
    return
end

[err, ~, ~, handles.codon_regions] = check_region(new_reg, handles.codon_regions);
if err == 0
    errordlg('region does not exist');
    return
end

update_regions(hObject, handles);


function fieldCodonFrom_Callback(hObject, eventdata, handles)
% hObject    handle to fieldCodonFrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldCodonFrom as text
%        str2double(get(hObject,'String')) returns contents of fieldCodonFrom as a double


% --- Executes during object creation, after setting all properties.
function fieldCodonFrom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldCodonFrom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fieldCodonTo_Callback(hObject, eventdata, handles)
% hObject    handle to fieldCodonTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fieldCodonTo as text
%        str2double(get(hObject,'String')) returns contents of fieldCodonTo as a double


% --- Executes during object creation, after setting all properties.
function fieldCodonTo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fieldCodonTo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butChimeraOpts.
function butChimeraOpts_Callback(hObject, eventdata, handles)
% hObject    handle to butChimeraOpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function textStat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to textStat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axPreview_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axPreview
