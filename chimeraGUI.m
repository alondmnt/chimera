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

% Last Modified by GUIDE v2.5 21-Jan-2018 21:13:49

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
function update_figure(hObject, handles)
handles = update_regions(handles);
handles = update_bar(handles);
handles = update_chimera(handles);
handles = update_status(handles);
guidata(hObject, handles);


function [region, err] = get_region(fieldFrom, fieldTo, target)
err = 0;
region = [int64(str2double(fieldFrom.String)), ...
          int64(str2double(fieldTo.String))];
if any(~isfinite(region) | region == 0)
    errordlg(sprintf(['expecting [from] and [to] to be signed integers. \n', ...
                      '>= 1 coordinates measured from START. \n', ...
                      '<= -1 coordinates measured from STOP.']), 'regions', 'modal');
    err = 1;
    return
end
if isempty(target)
    L = NaN;
else
    L = length(target);
end

if any(region < 0) && ~isfinite(L)
    errordlg('negative coordinates are undefined when target protein is unknown.', 'regions', 'modal')
    err = 3;
    return
end
region(region < 0) = L + 1 + region(region < 0);

if region(1) > region(2)
    errordlg('expecting [to] >= [from].', 'regions', 'modal');
    err = 2;
    return
end

if isfinite(L) && any(region > L)
    errordlg(sprintf('new region exceeds target protein length (%d > %d).', ...
                     region(find(region > L, 1)), L), 'regions', 'modal');
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


function handles = update_regions(handles)
% update all list boxes and barplot.
chim_region = region2set(handles.chimera_regions);
codon_region = region2set(handles.codon_regions);
assert(isempty(intersect(chim_region, codon_region)), 'codon/chimera overlap');

if handles.target_exist
    L = length(handles.target_seq);
else
    L = max([chim_region, codon_region]);
end
chim_region = chim_region(1 <= chim_region & chim_region <= L);
codon_region = codon_region(1 <= codon_region & codon_region <= L);

handles.default_regions = set2region(setdiff(setdiff(1:L, chim_region), codon_region));
handles.chimera_regions = set2region(chim_region);
handles.codon_regions = set2region(codon_region);

[handles.listChimera.String, handles.listChimera.Value] = update_listbox(...
	handles.chimera_regions, handles.target_seq, handles.fieldChimeraFrom, handles.fieldChimeraTo);
[handles.listCodon.String, handles.listCodon.Value] = update_listbox(...
    handles.codon_regions, handles.target_seq, handles.fieldCodonFrom, handles.fieldCodonTo);
handles.listDefault.Value = min(max(1, handles.listDefault.Value), size(handles.default_regions, 1));
handles.listDefault.String = update_listbox(handles.default_regions, handles.target_seq);


function [String, Value] = update_listbox(regions, target, fieldFrom, fieldTo)
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
    reg = get_region(fieldFrom, fieldTo, target);
catch
    Value = [];
    return
end
Value = [find(regions(:, 1) == reg(:, 1)); find(regions(:, 2) == reg(:, 2))];


function handles = update_bar(handles)
regions = {handles.chimera_regions, handles.codon_regions, handles.default_regions};
lens = cumsum(cellfun(@(x) size(x, 1), regions));

[regions, categ] = sortrows(cellcat(regions, 1));
regions = diff(regions, 1, 2) + 1;

ind = [categ <= lens(1), lens(1) < categ & categ <= lens(2), lens(2) < categ & categ <= lens(3)];
categ(ind(:, 1)) = 1;
categ(ind(:, 2)) = 2;
categ(ind(:, 3)) = 3;

cmap = colormap;
b = barh(handles.axPreview, [regions, regions]', 1, 'stacked', 'EdgeColor', 'none', ...
         'ButtonDownFcn', @select_region_from_figure);
[b(categ == 1).FaceColor] = deal([160, 197, 95]/255);  % [205, 120, 35]/255);
[b(categ == 2).FaceColor] = deal([123, 59, 59]/255);  % [41, 131, 20]/255);
[b(categ == 3).FaceColor] = deal([28, 135, 162]/255);  % [31, 149, 179]/255);

if handles.target_exist
    handles.axPreview.XLim = [0, length(handles.target_seq)+eps];
else
    handles.axPreview.XLim = [0, sum(regions)+eps];
end
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


function handles = update_status(handles)
is_ready = true;

if (handles.target_exist || handles.default_exist) && ~isempty(handles.target_seq)
    handles.statTarget.String = 'OK';
else
    handles.statRun.String = 'data still missing';
    is_ready = false;
    handles.statTarget.String = 'required';
end

is_default_req = ~isempty(handles.default_regions);
if handles.default_exist
    if strcmp(handles.target_seq, nt2aa(handles.default_seq, 'AlternativeStartCodons', false)) 
        handles.statDefault.String = 'OK';
    else
        handles.statDefault.String = 'aa != nt';
        if is_ready
            handles.statRun.String = 'invalid default';
        end
        is_ready = false;
    end
elseif is_default_req
    if is_ready
        handles.statRun.String = 'default missing';
    end
    is_ready = false;
    handles.statDefault.String = 'required';
else
    handles.statDefault.String = 'optional';
end

if handles.reference_exist
    handles.statReference.String = 'OK';  % sprintf('OK: %d', length(handles.reference_seq));
else
    if is_ready
        handles.statRun.String = 'reference missing';
    end
    is_ready = false;
    handles.statReference.String = 'required';
end

if ~any([length(handles.chimera_regions), length(handles.codon_regions)])
    if is_ready
        handles.statRun.String = 'regions missing';
    end
    is_ready = false;
end

if is_ready
    handles.statRun.String = 'ready';
    handles.butRun.Enable = 'on';
else
    handles.butRun.Enable = 'off';
end


function handles = update_chimera(handles)
tmp = sprintf('win size: %d codons\ncenter: %d\n', handles.winParams.size, handles.winParams.center);
if handles.winParams.by_start
    tmp = [tmp, 'start'];
end
if handles.winParams.by_start && handles.winParams.by_stop
    tmp = [tmp, ' + '];
end
if handles.winParams.by_stop
    tmp = [tmp, 'stop'];
end
handles.paramChimera.String = tmp;

% --- End of my functions --- %


% --- Executes just before chimeraGUI is made visible.
function chimeraGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chimeraGUI (see VARARGIN)
addpath('../');
addpath('../pos-spec/');

handles.axPreview.Box = 'on';

handles.target_seq = '';
handles.target_name = '';
handles.target_exist = false;

handles.default_seq = '';
handles.default_name = '';
handles.default_exist = false;

handles.reference_seq = {};
handles.reference_name = {};
handles.reference_exist = false;

handles.SA = zeros(0, 3);
handles.SA_exist = false;

handles.CUB = [];

handles.chimera_regions = zeros(0, 2);
handles.codon_regions = zeros(0, 2);
handles.default_regions = zeros(0, 2);

handles.winParams = struct('size', 70, 'center', 35, 'by_start', 1, 'by_stop', 1, 'truncate_seq', 1);

tuller_logo(handles.axLogo);

update_figure(hObject, handles);

% Choose default command line output for chimeraGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chimeraGUI wait for user response (see UIRESUME)
% uiwait(handles.figChimera);


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
                            handles.target_seq);
if err
    return
end

err = check_region(new_reg, handles.chimera_regions);
if err == 2
    errordlg('new region already exists', 'chimera regions', 'modal');
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

update_figure(hObject, handles);


% --- Executes on button press in butRun.
function butRun_Callback(hObject, eventdata, handles)
% hObject    handle to butRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 0. select output file
% TODO

% 1. codons optimization
handles.statRun.String = 'codons optim'; guidata(hObject, handles);

if handles.optSourceRef.Value
    codon_seq = maximize_CUB(handles.target_seq, handles.reference_seq);
elseif handles.optSourceTable.Value && isempty(handles.CUB)
    errordlg('missing a codon table', 'codon optimization', 'modal');
    return
else
    codon_seq = maximize_CUB(handles.target_seq, handles.CUB);
end

% 2. chimera optimization
if ~handles.SA_exist
    handles.statRun.String = 'building index'; guidata(hObject, handles);

    handles.reference_aa = nt2aa(handles.reference_seq, 'AlternativeStartCodon', false);  % here false is good
    lens = cellfun(@length, handles.reference_aa);
    handles.SA = build_suffix_array(handles.reference_aa, false);
    handles.SA(:, 3) = handles.SA(:, 1) - lens(handles.SA(:, 2)) - 1;  % equals -1 at end of seq
    handles.SA_exist = true;
end
handles.statRun.String = 'chimera optim'; guidata(hObject, handles);

if handles.winParams.size >= length(handles.target_seq)
    % not position specific
    warndlg('window size is larger than target sequence. that is, optimization is not position specific.');
    [chimera_seq, blocks] = calc_map(handles.target_seq, handles.SA, ...
                                     handles.reference_aa, handles.reference_seq);
else
    [chimera_seq, blocks] = calc_cmap_posspec1(handles.target_seq, handles.SA, ...
                                               handles.reference_aa, handles.reference_seq, ...
                                               handles.winParams);
end
blocks = table(blocks(:, 1), blocks(:, 2), blocks(:, 3), 'VariableNames', {'gene', 'pos', 'block'});
writetable(blocks, '../../output.csv')


% 3. complete construct from regions
handles.statRun.String = 'combining regions'; guidata(hObject, handles);

codon_region = region2set(handles.codon_regions);
chim_region = region2set(handles.chimera_regions);
def_region = region2set(handles.default_regions);
assert(isempty(intersect(codon_region, chim_region)) && ...
       isempty(intersect(codon_region, def_region)) && ...
       isempty(intersect(chim_region, def_region)), 'overlapping regions')

seq_bank = {codon_seq, chimera_seq, handles.default_seq};
seq_source(codon_region) = 1;
seq_source(chim_region) = 2;
seq_source(def_region) = 3;
final_seq = cellcat(arrayfun(@(x, y) {seq_bank{x}(3*(y-1)+1:3*y)}, seq_source, 1:length(seq_source)), 2);
assert(strcmp(nt2aa(final_seq, 'AlternativeStartCodons', false), handles.target_seq), 'final seq error');

outfile = '../../output.fasta';
if exist(outfile, 'file')
    delete(outfile);
end
fastawrite(outfile, sprintf('%s optimized by chimeraDesigner', handles.target_name), final_seq);

handles.statRun.String = 'done';
guidata(hObject, handles);


function fileTarget_Callback(hObject, eventdata, handles)
% hObject    handle to fileTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileTarget as text
%        str2double(get(hObject,'String')) returns contents of fileTarget as a double
aa_alphabet = '*ACDEFGHIKLMNPQRSTVWY';
tmp = upper(hObject.String);
valid = ismember(tmp, aa_alphabet);
tmp = tmp(valid);
if ~isempty(tmp) && tmp(end) ~= '*'
    tmp(end+1) = '*';
end
hObject.String = tmp;
if strcmp(tmp, handles.target_seq)
    return
end
handles.target_seq = tmp;
handles.target_name = 'user_input';

if ~isempty(handles.target_seq)
    handles.target_exist = true;
else
    handles.target_exist = false;
end

update_figure(hObject, handles);


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
[fname, dirname] = uigetfile({'*.fa;*.fasta', 'fasta file'}, 'select a protein sequence file');
if fname == 0
    errordlg('missing file', 'target seq', 'modal');
    return
end
seqs = fastaread(fullfile(dirname, fname));

if length(seqs) > 1
    [ind, ok] = listdlg('ListString', {seqs.Header}, 'Name', 'select a protein', ...
                        'ListSize', [450, 300], 'SelectionMode', 'single');
    if ok == 0
        return
    end
    seqs = seqs(ind);
end

aa_alphabet = '*ACDEFGHIKLMNPQRSTVWY';
valid = ismember(upper(seqs.Sequence), aa_alphabet);
if ~all(valid)
    errordlg(sprintf('protein sequence contains illegal chars ("%s"). \nexpecting an amino acid sequence.', unique(seqs.Sequence(~valid))), 'target seq', 'modal');
    return
end
nt_alphabet = 'ACGTU';
if all(ismember(upper(seqs.Sequence), nt_alphabet))
    warndlg(sprintf('protein sequence may comprise of nucleotides. \nexpecting an amino acid sequence. \nyou may use the default sequence input to select both AA seq and NT seq at once.'));
end

if seqs.Sequence(end) ~= '*'
    seqs.Sequence(end+1) = '*';
end

handles.fileTarget.String = upper(seqs.Sequence);
handles.target_seq = upper(seqs.Sequence);
handles.target_name = seqs.Header;
handles.target_exist = true;

update_figure(hObject, handles);


function check_text(eventdata, alphabet)
if ~ismember(upper(eventdata.Character), alphabet)
    if ismember(eventdata.Key, {'backspace', 'rightarrow', 'leftarrow', 'return', 'control'})
        return
    end
    if any(ismember(eventdata.Modifier, {'control'}));
        return
    end
    errordlg(sprintf('expecting a sequence with alphabet: %s', alphabet), 'sequence error', 'modal');
    return
end


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
[fname, dirname] = uigetfile({'*.fa;*.fasta', 'fasta file'}, 'select a reference sequence file');
if fname == 0
    errordlg('missing file', 'reference seq', 'modal');
    return
end
seqs = fastaread(fullfile(dirname, fname));

nt_alphabet = 'ACGTU';
valid_seq = true(length(seqs), 1);
valid_len = true(length(seqs), 1);
for i = 1:length(seqs)
    valid_seq(i) = all(ismember(upper(seqs(i).Sequence), nt_alphabet));
    valid_len(i) = mod(length(seqs(i).Sequence), 3) == 0;
    seqs(i).Sequence = strrep(upper(seqs(i).Sequence), 'U', 'T');
end

if ~all(valid_seq)
    listdlg('ListString', {seqs(~valid_seq).Header}, 'Name', 'invalid seqs', ...
            'ListSize', [450, 300], ...
            'PromptString', 'the following seqs contain invalid chars and will be discarded.');
end
seqs = seqs(valid_seq);
valid_len = valid_len(valid_seq);

if ~all(valid_len)
    [sel, decision] = listdlg('ListString', {seqs(~valid_len).Header}, ...
                              'Name', 'invalid length', 'ListSize', [450, 300], ...
                              'InitialValue', 1:sum(~valid_len), ...
                              'OKString', 'discard selected', 'CancelString', 'ignore warning', ...
                              'PromptString', 'the following seqs have lengths not divisible by 3.');
    if decision == 1  % 'discard selected'
        seqs(sel) = [];
    end
end

handles.fileReference.String = sprintf('%s: %d seqs', fname, length(seqs));
handles.reference_seq = {seqs.Sequence}';
handles.SA = zeros(0, 3);

handles.reference_exist = true;
handles.SA_exist = false;

update_figure(hObject, handles);


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
                            handles.target_seq);
if err
    return
end

[err, ~, ~, handles.chimera_regions] = check_region(new_reg, handles.chimera_regions);
if err == 0
    errordlg('region does not exist', 'chimera regions', 'modal');
    return
end

update_figure(hObject, handles);


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
nt_alphabet = 'ACGT';
tmp = strrep(upper(hObject.String), 'U', 'T');
valid = ismember(tmp, nt_alphabet);
tmp = tmp(valid);
hObject.String = tmp;
if strcmp(tmp, handles.target_seq)
    return
end
handles.default_seq = tmp;

if isempty(tmp)
    handles.default_exist = false;
else
    % default seq overrides target AA seq
    handles.target_seq = nt2aa(tmp, 'AlternativeStartCodons', false);  % TODO: reconsider alternative
    handles.fileTarget.String = handles.target_seq;
    handles.target_name = 'user_input';
    
    handles.default_exist = true;
    handles.target_exist = true;
end

update_figure(hObject, handles);


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
[fname, dirname] = uigetfile({'*.fa;*.fasta', 'fasta file'}, 'select a default sequence file');
if fname == 0
    errordlg('missing file', 'default seq', 'modal');
    return
end
seqs = fastaread(fullfile(dirname, fname));

if length(seqs) > 1
    [ind, ok] = listdlg('ListString', {seqs.Header}, 'Name', 'select a sequence', ...
                        'ListSize', [450, 300], 'SelectionMode', 'single');
    if ok == 0
        return
    end
    seqs = seqs(ind);
end

nt_alphabet = 'ACGTU';
valid = ismember(upper(seqs.Sequence), nt_alphabet);
if ~all(valid)
    errordlg(sprintf('default sequence contains illegal chars ("%s"). \nexpecting a nucleotide sequence.', unique(seqs.Sequence(~valid))), 'default seq', 'modal');
    return
end

tmp = strrep(upper(seqs.Sequence), 'U', 'T');
tmp_aa = nt2aa(tmp, 'AlternativeStartCodons', false);
if handles.target_exist && ~strcmp(tmp_aa, handles.target_seq)
    decision = questdlg('default sequence does not match the current target sequence.', ...
                        'default sequence error', 'update target', 'cancel', 'cancel');
    switch decision
        case 'update target'
            handles.fileTarget.String = tmp_aa;
            handles.target_seq = tmp_aa;
            handles.target_name = seqs.Header;
        otherwise
            return
    end
elseif ~handles.target_exist
    handles.fileTarget.String = tmp_aa;
    handles.target_seq = tmp_aa;
    handles.target_name = seqs.Header;
    handles.target_exist = true;
end
handles.fileDefault.String = tmp;
handles.default_seq = tmp;

handles.default_exist = true;

update_figure(hObject, handles);


% --- Executes on button press in optSourceRef.
function optSourceRef_Callback(hObject, eventdata, handles)
% hObject    handle to optSourceRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optSourceRef
handles.CUB = [];  % ensure that we don't keep old tables
guidata(hObject, handles);


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
                            handles.target_seq);
if err
    return
end

err = check_region(new_reg, handles.codon_regions);
if err == 2
    errordlg('new region already exists', 'codon regions', 'modal');
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
    decision = questdlg('new region overlaps with an existing chimera region', ...
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

update_figure(hObject, handles);


% --- Executes on button press in butCodonRemRegion.
function butCodonRemRegion_Callback(hObject, eventdata, handles)
% hObject    handle to butCodonRemRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_reg, err] = get_region(handles.fieldCodonFrom, handles.fieldCodonTo, ...
                            handles.target_seq);
if err
    return
end

[err, ~, ~, handles.codon_regions] = check_region(new_reg, handles.codon_regions);
if err == 0
    errordlg('region does not exist', 'codon regions', 'modal');
    return
end

update_figure(hObject, handles);


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
handles.winParams = winParams(handles.winParams);
update_figure(hObject, handles);


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


% --- Executes on key press with focus on fileTarget and none of its controls.
function fileTarget_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to fileTarget (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
aa_alphabet = '*ACDEFGHIKLMNPQRSTVWY';
check_text(eventdata, aa_alphabet);


function fileRegions_Callback(hObject, eventdata, handles)
% hObject    handle to fileRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileRegions as text
%        str2double(get(hObject,'String')) returns contents of fileRegions as a double


% --- Executes during object creation, after setting all properties.
function fileRegions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in butRegions.
function butRegions_Callback(hObject, eventdata, handles)
% hObject    handle to butRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on fileDefault and none of its controls.
function fileDefault_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to fileDefault (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
nt_alphabet = 'ACGTU';
check_text(eventdata, nt_alphabet);


% --- Executes during object creation, after setting all properties.
function axLogo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axLogo


% --- Executes on button press in optSourceTable.
function optSourceTable_Callback(hObject, eventdata, handles)
% hObject    handle to optSourceTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optSourceTable
[fname, dirname] = uigetfile({'*.tsv;*.tab', 'tab separated values'}, 'select a codon score table');
if fname == 0
    errordlg('missing file', 'codon table', 'modal');
    handles.optSourceRef.Value = 1;
    handles.optSourceTable.Value = 0;
    handles.CUB = [];
    guidata(hObject, handles);
    return
end

file_OK = true;
aa_list = fieldnames(aacount(''));
fid = fopen(fullfile(dirname, fname));
try
    T = textscan(fid, '%s\t%f');
    T{3} = nt2aa(T{1}, 'AlternativeStartCodon', false);
    aa_OK = ismember(aa_list, T{3});
catch
    file_OK = false;
end
fclose(fid);

if ~file_OK
    errordlg('codon table format error. expecting a tab-separated file with 2 columns (and no header): [codon, score]', 'codon table', 'modal');
    handles.optSourceRef.Value = 1;
    handles.optSourceTable.Value = 0;
    handles.CUB = [];
    guidata(hObject, handles);
    return
end

if ~all(aa_OK)
    errordlg(sprintf('the following AA are missing from table: %s', char(aa_list(~aa_OK))'), 'codon table', 'modal');
    handles.optSourceRef.Value = 1;
    handles.optSourceTable.Value = 0;
    handles.CUB = [];
    guidata(hObject, handles);
    return
end

handles.CUB = codonbias('');
aa_list = fieldnames(handles.CUB)';
for aa = aa_list
    iaa = find(strcmpi(T{3}, aminolookup(aa{1})))';
    for i = iaa
        icod = strcmpi(T{1}(i), handles.CUB.(aa{1}).Codon);
        handles.CUB.(aa{1}).Freq(icod) = T{2}(i);
    end
end
guidata(hObject, handles);
