function varargout = chimeraGUI_macOS(varargin)
% CHIMERAGUI_MACOS MATLAB code for chimeraGUI_macOS.fig
%
%      ChimeraUGEM: Unsupervised Gene Expression Modeling
%      Diament et al., 2018
%
%      Version 1.0
%      Alon Diament / Tuller Lab, September 2018.
%
%      CHIMERAGUI_MACOS, by itself, creates a new CHIMERAGUI_MACOS or raises the existing
%      singleton*.
%
%      H = CHIMERAGUI_MACOS returns the handle to a new CHIMERAGUI_MACOS or the handle to
%      the existing singleton*.
%
%      CHIMERAGUI_MACOS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHIMERAGUI_MACOS.M with the given input arguments.
%
%      CHIMERAGUI_MACOS('Property','Value',...) creates a new CHIMERAGUI_MACOS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before chimeraGUI_macOS_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to chimeraGUI_macOS_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help chimeraGUI_macOS

% Last Modified by GUIDE v2.5 12-Aug-2018 13:04:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @chimeraGUI_macOS_OpeningFcn, ...
                   'gui_OutputFcn',  @chimeraGUI_macOS_OutputFcn, ...
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
function x = get_i0(x, i)
% ~immitates python 0-indexing

x = x(mod(i, length(x)) + 1);


function handles = update_figure(hObject, handles)
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
% n = length(set);
% set = unique(set);
% assert(length(set) == n, 'overlapping regions')


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
% update all list boxes and barplot (regions map). showing the first set of
% regions when a multiple-region file is known.
if isempty(handles.chimera_regions)
    chim_region = [];
else
    chim_region = region2set(handles.chimera_regions{1});
end
if isempty(handles.codon_regions)
    codon_region = [];
else
    codon_region = region2set(handles.codon_regions{1});
end
assert(isempty(intersect(chim_region, codon_region)), 'codon/chimera overlap');

if strcmp(handles.fieldChimeraTo.Enable, 'off')
    % otherwise this may raise warnings
    handles.fieldChimeraTo.String = '-1';
    handles.fieldCodonTo.String = '-1';
    handles.fieldChimeraFrom.String = '1';
    handles.fieldCodonFrom.String = '1';
end
if isempty(handles.target_seq)
    handles.target_exist = false;
    handles.default_exist = false;
    handles.target_seq = {''};
    handles.target_name = {''};
    handles.default_seq = {''};
    handles.fileTarget.String = '';
    handles.fileRegions.String = '';
end
if handles.target_exist
    L = length(handles.target_seq{1});
else
    L = max([chim_region, codon_region]);
end
chim_region = chim_region(1 <= chim_region & chim_region <= L);
codon_region = codon_region(1 <= codon_region & codon_region <= L);

handles.default_regions{1} = set2region(setdiff(setdiff(1:L, chim_region), codon_region));
handles.chimera_regions{1} = set2region(chim_region);
handles.codon_regions{1} = set2region(codon_region);

[handles.listChimera.String, handles.listChimera.Value] = update_listbox(...
	handles.chimera_regions{1}, handles.target_seq{1}, handles.fieldChimeraFrom, handles.fieldChimeraTo);
[handles.listCodon.String, handles.listCodon.Value] = update_listbox(...
    handles.codon_regions{1}, handles.target_seq{1}, handles.fieldCodonFrom, handles.fieldCodonTo);
handles.listDefault.Value = min(max(1, handles.listDefault.Value), size(handles.default_regions{1}, 1));
handles.listDefault.String = update_listbox(handles.default_regions{1}, handles.target_seq{1});


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
% plots a map of all region types in the currently designed gene.
regions = {handles.chimera_regions{1}, handles.codon_regions{1}, handles.default_regions{1}};
lens = cumsum(cellfun(@(x) size(x, 1), regions));

[regions, categ] = sortrows(cellcat(regions, 1));
regions = diff(regions, 1, 2) + 1;

ind = [categ <= lens(1), lens(1) < categ & categ <= lens(2), lens(2) < categ & categ <= lens(3)];
categ(ind(:, 1)) = 1;
categ(ind(:, 2)) = 2;
categ(ind(:, 3)) = 3;

if ~isempty(regions)
    cmap = colormap;
    b = barh(handles.axPreview, [regions, regions]', 1, 'stacked', 'EdgeColor', 'none', ...
        'ButtonDownFcn', @select_region_from_figure);
    [b(categ == 1).FaceColor] = deal([160, 197, 95]/255);  % [205, 120, 35]/255);
    [b(categ == 2).FaceColor] = deal([28, 135, 162]/255);  % [31, 149, 179]/255);
    [b(categ == 3).FaceColor] = deal([123, 59, 59]/255);  % [41, 131, 20]/255);
else
    cla(handles.axPreview);
end
if handles.target_exist
    handles.axPreview.XLim = [0, length(handles.target_seq{1})+eps];
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

ichim = handles.chimera_regions{1}(:, 1) <= pos & ...
        pos <= handles.chimera_regions{1}(:, 2);
if any(ichim) && strcmp(handles.fieldChimeraFrom.Enable, 'on')
    handles.listChimera.Value = find(ichim);
    handles.fieldChimeraFrom.String = sprintf('%d', handles.chimera_regions{1}(ichim, 1));
    handles.fieldChimeraTo.String = sprintf('%d', handles.chimera_regions{1}(ichim, 2));
end

icod = handles.codon_regions{1}(:, 1) <= pos & ...
       pos <= handles.codon_regions{1}(:, 2);
if any(icod) && strcmp(handles.fieldCodonFrom.Enable, 'on')
    handles.listCodon.Value = find(icod);
    handles.fieldCodonFrom.String = sprintf('%d', handles.codon_regions{1}(icod, 1));
    handles.fieldCodonTo.String = sprintf('%d', handles.codon_regions{1}(icod, 2));
end

idef = handles.default_regions{1}(:, 1) <= pos & ...
       pos <= handles.default_regions{1}(:, 2);
if any(idef) && strcmp(handles.fieldChimeraFrom.Enable, 'on')
    handles.listDefault.Value = find(idef);
    handles.fieldDefaultFrom.String = sprintf('%d', handles.default_regions{1}(idef, 1));
    handles.fieldDefaultTo.String = sprintf('%d', handles.default_regions{1}(idef, 2));
end

guidata(hObject, handles);


function handles = update_status(handles)
% this functions determines whether all prerequisites for the optimization
% are met (and signals the user whether data is missing).
is_ready_2map = true;
is_ready_2ars = true;
handles.statTarget.Visible = 'off';
handles.statReference.Visible = 'off';
handles.statRegionFile.Visible = 'off';
handles.statRegionSel.Visible = 'off';
handles.optNT.ForegroundColor = [0, 0, 0];
handles.optNT.FontWeight = 'normal';
handles.optRefNT.ForegroundColor = [0, 0, 0];
handles.optRefNT.FontWeight = 'normal';
handles.statMap.String = '';
handles.statARS.String = '';
handles.butRegions.Enable = 'off';


if (handles.target_exist || handles.default_exist) && ~isempty(handles.target_seq)
    handles.butRegions.Enable = 'on';
%     handles.statTarget.String = 'OK';
else
    handles.statMap.String = 'target missing';
    is_ready_2map = false;
    is_ready_2ars = false;
    handles.statTarget.Visible = 'on';
end

is_default_req = ~all(cellfun(@isempty, handles.default_regions));
if is_default_req && ~handles.default_exist
    if is_ready_2map
        handles.statMap.String = 'default NT seq missing';
        handles.statTarget.Visible = 'on';
        handles.optNT.ForegroundColor = handles.statTarget.ForegroundColor;
        handles.optNT.FontWeight = 'bold';
    end
    is_ready_2map = false;
end

if handles.reference_exist
    if ~strcmpi(handles.reference_type, 'NT') && is_ready_2map
        handles.statMap.String = 'NT reference needed';
        handles.statReference.Visible = 'on';
        handles.optRefNT.ForegroundColor = handles.statTarget.ForegroundColor;
        handles.optRefNT.FontWeight = 'bold';
        is_ready_2map = false;
        % can still calc cARS
    end
else
    if is_ready_2map || is_ready_2ars
        handles.statMap.String = 'reference missing';
        handles.statReference.Visible = 'on';
    end
    is_ready_2map = false;
    is_ready_2ars = false;
end

if length(handles.target_seq) > 1
    region_editing = 'off';
else
    region_editing = 'on';
end
handles.fieldChimeraFrom.Enable = region_editing;
handles.fieldChimeraTo.Enable = region_editing;
handles.butChimeraAddRegion.Enable = region_editing;
handles.butChimeraRemRegion.Enable = region_editing;
handles.fieldCodonFrom.Enable = region_editing;
handles.fieldCodonTo.Enable = region_editing;
handles.butCodonAddRegion.Enable = region_editing;
handles.butCodonRemRegion.Enable = region_editing;

if ~all(cellfun(@length, {handles.chimera_regions, ...
                          handles.codon_regions, ...
                          handles.default_regions}) == length(handles.target_seq))
    if is_ready_2map
%         handles.statMap.String = 'regions file missing';
        handles.statRegionFile.Visible = 'on';
    end
%     is_ready = false;
end

if any(cellfun(@(x, y) isempty(x) & isempty(y), ...
               handles.chimera_regions, handles.codon_regions))
    if is_ready_2map
        handles.statMap.String = 'regions missing';
        handles.statRegionSel.Visible = 'on';
    end
    is_ready_2map = false;
end

if any(cellfun(@isempty, handles.chimera_regions))
    if is_ready_2ars
        handles.statARS.String = 'no chimera region';
        handles.statRegionSel.Visible = 'on';
    end
    is_ready_2ars = false;
end

if handles.during_ars && is_ready_2map
    handles.statMap.String = 'cARS is on';
    is_ready_2map = false;
end
if handles.during_map && is_ready_2ars
    handles.statARS.String = 'cMap is on';
    is_ready_2ars = false;
end

if is_ready_2map
    handles.statMap.String = 'ready';
    handles.butMap.Enable = 'on';
else
    handles.butMap.Enable = 'off';
end

if is_ready_2ars
    handles.statARS.String = 'ready';
    handles.butARS.Enable = 'on';
else
    handles.butARS.Enable = 'off';
end

if handles.target_exist
    n_seq = length(handles.target_seq);
    if n_seq > 1
        handles.fileTarget.String = sprintf('%s: [%d %s seqs]', handles.target_file, ...
            n_seq, handles.target_type);
    else
        handles.fileTarget.String = sprintf('%s: [%s seq] %s', handles.target_file, ...
            handles.target_type, handles.target_name{1});
    end
else
    handles.fileTarget.String = '';
end


function handles = update_chimera(handles)
% algorithm parameters displayed
str_block = {};
if isfinite(handles.winParams.max_len)
    str_block = [str_block, {sprintf('%d cod', handles.winParams.max_len)}];
end
if isfinite(handles.winParams.max_pos)
    str_block = [str_block, {sprintf('%.0f%%', 100 * handles.winParams.max_pos)}];
end
if isempty(str_block)
    str_block = 'none';
else
    str_block = strjoin(str_block, ', ');
end
tmp = sprintf('win size: %d codons\nlimit: %s\ncenter: %d\n', ...
              handles.winParams.size, str_block, handles.winParams.center);
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


function [seq_res, len_res, iref] = is_target_in_ref(handles, alphabet)
% DEPRECATED. see: filter_target_from_ref()
assert(ischar(alphabet), 'alphabet must be a string');
alphabet = upper(alphabet);

len_res = false;
seq_res = false;
iref = find(strcmp(handles.target_name{1}, handles.reference_name));

for i = 1:length(iref)
    switch alphabet
        case 'AA'
            len_res(i) = length(handles.target_seq{1}) == length(handles.reference_aa{iref(i)});
            seq_res(i) = strcmpi(handles.target_seq{1}, handles.reference_aa{iref(i)});

        case 'NT'
            if strcmpi(handles.target_type, 'NT')
                tlen = length(handles.default_seq{1});
                tseq = handles.default_seq{1};
            else
                tlen = 3*length(handles.target_seq{1});
                tseq = 'z';
            end
            
            if strcmpi(handles.reference_type, 'NT')
                rlen = length(handles.reference_seq{iref(i)});
                rseq = handles.reference_seq{iref(i)};
            else
                rlen = 3*length(handles.reference_aa{iref(i)});
                rseq = 'x';  % different from tseq
            end

            len_res(i) = rlen == tlen;
            seq_res(i) = strcmpi(rseq, tseq);

        otherwise
            error('unknown alphabet type');
    end
end


function [filtered, handles] = filter_target_from_ref(handles, ctarg, cref)
% using the suffix array to search the reference for sequences identical
% to the target.
% accepts handles containing the suffix array (containing redundant,
% non-unique suffixes in field [SA]) and optionally an action to take
% (field [filt_action]), the target sequence, and the refrence set.

handles.SA(:, 3) = 1;  % reset all suffix frequencies to 1 (assuming non-unique SA)
filtered = 0;
prompt_user = ~isfield(handles, 'filt_action');
if ~prompt_user && strcmpi(handles.filt_action, 'ignore')
    return
end

while strcmp(longest_prefix(ctarg, handles.SA, cref), ctarg)
    filtered = filtered + 1;
    if prompt_user
        handles.filt_action = questdlg('some target sequences appear in the reference' , ...
                          'target in reference', 'filter from ref', ...
                          'skip target', 'ignore', 'filter from ref');
        prompt_user = false;
    end

    switch handles.filt_action
        case 'filter from ref'
            % find which ref to remove
            [~, iSA] = longest_prefix(ctarg, handles.SA, cref);
            iref = handles.SA(iSA, 2);  % seq index
            handles.SA(handles.SA(:, 2) == iref, 3) = 0; 
            % set freq to 0 to remove from pos-spec cARS/cMap reference.
            % this is more efficient than copy/slicing and is required to
            % keep the sorted indices in columns 5-6 stable.
            handles.filt_msg = 'targets filtered from reference:';
        case 'skip target'
            handles.target_seq(1) = [];
            handles.default_seq(1) = [];
            handles.target_name(1) = [];
            handles.chimera_regions(1) = [];
            handles.codon_regions(1) = [];
            handles.default_regions(1) = [];
            handles.filt_msg = 'targets skipped:';
            return
        case 'ignore'
            handles.filt_msg = 'targets ignored:';
            return
        otherwise
            error('too many cooks!');
    end
end


function w = CUB2weights(CUB)
    aa_list = fieldnames(CUB);
    n = length(aa_list);
    codon_list = cell(n, 1);
    w = cell(n, 1);
    for aa = 1:n
        codon_list{aa} = CUB.(aa_list{aa}).Codon;
        w{aa} = CUB.(aa_list{aa}).Freq;
    end
    codon_list = cellcat(codon_list, 2);
    w = cellcat(w, 2);
    [~, isort] = sort(codon_list);
    w = w(isort)';

% --- End of my functions --- %


% --- Executes just before chimeraGUI_macOS is made visible.
function chimeraGUI_macOS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to chimeraGUI_macOS (see VARARGIN)

handles.axPreview.Box = 'on';

handles.target_seq = {''};
handles.target_name = {''};
handles.target_exist = false;

handles.default_seq = {''};
handles.default_name = {''};
handles.default_exist = false;

handles.reference_seq = {};
handles.reference_name = {};
handles.reference_exist = false;

handles.SA = zeros(0, 3);
handles.SA_exist = false;

handles.CUB = [];

handles.chimera_regions = {zeros(0, 2)};
handles.codon_regions = {zeros(0, 2)};
handles.default_regions = {zeros(0, 2)};

handles.during_ars = false;
handles.during_map = false;

handles.winParams = struct('size', 40, 'center', 0, ...
                           'by_start', 1, 'by_stop', 1, ...
                           'truncate_seq', 0, ...
                           'max_len', 40, 'max_pos', 0.5);

logo(handles.axLogo);

update_figure(hObject, handles);

% Choose default command line output for chimeraGUI_macOS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes chimeraGUI_macOS wait for user response (see UIRESUME)
% uiwait(handles.figChimera);


% --- Outputs from this function are returned to the command line.
function varargout = chimeraGUI_macOS_OutputFcn(hObject, eventdata, handles) 
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
    handles.fieldChimeraFrom.String = sprintf('%d', handles.chimera_regions{1}(reg, 1));
    handles.fieldChimeraTo.String = sprintf('%d', handles.chimera_regions{1}(reg, 2));
end
guidata(hObject, handles);


% --- Executes on button press in butChimeraAddRegion.
function butChimeraAddRegion_Callback(hObject, eventdata, handles)
% hObject    handle to butChimeraAddRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_reg, err] = get_region(handles.fieldChimeraFrom, handles.fieldChimeraTo, ...
                            handles.target_seq{1});
if err
    return
end

err = check_region(new_reg, handles.chimera_regions{1});
if err == 2
    errordlg('new region already exists', 'chimera regions', 'modal');
    return
elseif err == 1
    decision = questdlg('new region overlaps with an existing Chimera region' , ...
                        'region overlap', 'merge', 'cancel', 'cancel');
    if ~strcmp(decision, 'merge')
        return
    end
end

[err, ~, trunc_new, trunc_old] = check_region(new_reg, handles.codon_regions{1});
if err
    decision = questdlg('new region overlaps with an existing codon region', ...
                        'region overlap', 'keep chimera', 'keep codon', 'keep chimera');
    switch decision
        case 'keep codon'
            new_reg = trunc_new;
        case 'keep chimera'
            handles.codon_regions{1} = trunc_old;
        otherwise
            return
    end
end

[~, handles.chimera_regions{1}] = check_region(new_reg, handles.chimera_regions{1});

update_figure(hObject, handles);


% --- Executes on button press in butMap.
function butMap_Callback(hObject, eventdata, handles)
% hObject    handle to butMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 0. select output file
profiling_run = false;
[fname, dirname] = uiputfile({'*.fa;*.fasta', 'fasta file'}, 'select output file');
if fname == 0
    return
elseif strcmp(fname, 'profiling.fa')
    profiling_run = true;
end
[~, fname, extent] = fileparts(fname);
if isempty(extent)
    extent = '.fasta';
end
outfile = fullfile(dirname, fname);
outfasta = strcat(outfile, extent);
if exist(outfasta, 'file')
    delete(outfasta);
end
if handles.optSourceRef.Value
    handles.CUB = [];  % we will compute CUB later, once for all targets
end
tic;

handles.during_map = true;
handles = update_figure(hObject, handles);

ntarg = length(handles.target_seq);
target_filter = {};
homolog_filter = {};
target_error = {};
cmap_names = cellfun(@(x) get_i0(strsplit(strrep(x, ',', ';'), ' '), 0), handles.target_name);
for t = 1:ntarg
    % 1. codons optimization
    do_codon = ~isempty(handles.codon_regions{1});
    handles.statMap.String = sprintf('%d: codons optim', t); drawnow;

    if do_codon && handles.optSourceRef.Value
        if isempty(handles.CUB)
            [codon_seq, handles.CUB] = maximize_CUB(handles.target_seq{1}, handles.reference_seq);
        else
            codon_seq = maximize_CUB(handles.target_seq{1}, handles.CUB);
        end
    elseif do_codon && handles.optSourceTable.Value && isempty(handles.CUB)
        errordlg('missing a codon table', 'codon optimization', 'modal');
        handles.during_map = false;
        update_figure(hObject, handles);
        return
    elseif do_codon
        codon_seq = maximize_CUB(handles.target_seq{1}, handles.CUB);
    else
        codon_seq = '';
    end

    % 2. chimera optimization
    do_chimera = ~isempty(handles.chimera_regions{1});
    if do_chimera && ~(handles.SA_exist && strcmpi(handles.SA_type, 'AA'))
        handles.statMap.String = sprintf('%d: building index', t); drawnow;
        handles.SA = build_suffix_array(handles.reference_aa, false);
        handles.SA_type = 'AA';
        handles.SA_exist = true;
    end
    handles.statMap.String = sprintf('%d: chimera optim', t); drawnow;

    % check if target is already in reference
    if handles.default_exist
        % re-translating with alternative=False to match the reference
        [filtered, handles] = filter_target_from_ref(handles, ...
                                                nt2aa(handles.default_seq{1}, 'AlternativeStartCodons', false), ...
                                                handles.reference_aa);
    else
        [filtered, handles] = filter_target_from_ref(handles, handles.target_seq{1}, handles.reference_aa);
    end
    if filtered
        switch handles.filt_action
            case 'filter from ref'
                target_filter{end+1} = sprintf('%05d: %s (filtered %d)', t, cmap_names{t}, filtered);
                handles.statMap.String = sprintf('%d: chimera (filtered %d)', t, filtered); drawnow;
            case 'skip target'
                target_filter{end+1} = sprintf('%05d: %s', t, cmap_names{t});
                handles = update_figure(hObject, handles);
                continue
        end
    end

    if do_chimera && (handles.winParams.size == 0)
        % not position specific
        [chimera_seq, mblocks, err, homologs] = calc_cmap(handles.target_seq{1}, handles.SA, ...
            handles.reference_aa, handles.reference_seq, ...
            handles.winParams.max_len, handles.winParams.max_pos);
    elseif do_chimera
        [chimera_seq, mblocks, err, homologs] = calc_cmap_posspec(handles.target_seq{1}, handles.SA, ...
            handles.reference_aa, handles.reference_seq, ...
            handles.winParams, ...
            handles.winParams.max_len, handles.winParams.max_pos);
    else
        chimera_seq = '';
    end

    if homologs
        homolog_filter{end+1} = sprintf('%05d: %s (filtered %d)', t, cmap_names{t}, homologs);
    end
    if any(err)
        err_string = sprintf('%05d: %s', t, cmap_names{t});
        if err(1)
            err_string = sprintf('%s [empty block]', err_string);
        end
        if err(2)
            err_string = sprintf('%s [empty window]', err_string);
        end
        if err(3)
            err_string = sprintf('%s [peptide error]', err_string);
        end
        target_error{end+1} = err_string;
    end

    % 3. complete construct from regions
    handles.statMap.String = sprintf('%d: combining regions', t); drawnow;

    codon_region = region2set(handles.codon_regions{1});
    chim_region = region2set(handles.chimera_regions{1});
    def_region = region2set(handles.default_regions{1});
    assert(isempty(intersect(codon_region, chim_region)) && ...
        isempty(intersect(codon_region, def_region)) && ...
        isempty(intersect(chim_region, def_region)), 'overlapping regions')

    seq_bank = {codon_seq, chimera_seq, handles.default_seq{1}};
    seq_source = [];
    seq_source(codon_region) = 1;
    seq_source(chim_region) = 2;
    seq_source(def_region) = 3;
    final_seq = cellcat(arrayfun(@(x, y) {seq_bank{x}(3*(y-1)+1:3*y)}, seq_source, 1:length(seq_source)), 2);
    if ~do_chimera || ~any(err)
        assert(strcmp(nt2aa(final_seq, 'AlternativeStartCodons', true), handles.target_seq{1}), 'final seq error');
    end
    % NOTE: as long as we allow alternative starts in [default_seq], so do
    %       we need to allow it here

    % generate a block table
    if ~profiling_run
        handles.statMap.String = sprintf('%d: saving', t); drawnow;
        blocks = table(0, 0, {'init'}, {''}, 0, {''}, 'VariableNames', ...
            {'pos_s', 'pos_e', 'type', 'gene', 'gene_loc', 'block'});
        if do_chimera
            mblocks = table(cumsum([1; cellfun(@length, mblocks(1:end-1, 3))/3]), ...
                cumsum(cellfun(@length, mblocks(:, 3))/3), ...
                mblocks(:, 1), cell2mat(mblocks(:, 2)), mblocks(:, 3), 'VariableNames', ...
                {'pos_s', 'pos_e', 'gene', 'gene_loc', 'block'});
            erase = [];
            for b = 1:height(mblocks)
                % keep entire block if it intersects with a chimeric region
                in_region = intersect(mblocks.pos_s(b):mblocks.pos_e(b), chim_region);
                if isempty(in_region)
                    erase = [erase; b];
                end
            end
            mblocks(erase, :) = [];
            valid = ~cellfun(@isnan, mblocks.gene);
            mblocks.gene(valid) = cellfun(@(x) get_i0(strsplit(strrep(handles.reference_name{x}, ',', ';'), ' '), 0), mblocks.gene(valid));
            mblocks.gene(~valid) = {'NA'};
            mblocks.type = repelem({'cMap'}, height(mblocks), 1);
            blocks = outerjoin(blocks, mblocks, 'MergeKeys', true);
        end
        if do_codon
            cblocks = table(handles.codon_regions{1}(:, 1), handles.codon_regions{1}(:, 2), ...
                'VariableNames', {'pos_s', 'pos_e'});
            if handles.optSourceRef.Value
                str_codon = sprintf('codon: %s', handles.reference_file);
            else
                str_codon = sprintf('codon: %s', handles.CUB_file);
            end
            cblocks.type = repelem({str_codon}, height(cblocks), 1);
            blocks = outerjoin(blocks, cblocks, 'MergeKeys', true);
        end
        if ~isempty(def_region)
            dblocks = table(handles.default_regions{1}(:, 1), handles.default_regions{1}(:, 2), ...
                'VariableNames', {'pos_s', 'pos_e'});
            dblocks.type = repelem({'default'}, height(dblocks), 1);
            blocks = outerjoin(blocks, dblocks, 'MergeKeys', true);
        end
        blocks(1, :) = [];

        writetable(blocks, sprintf('%s_%s.csv', outfile, cmap_names{t}))
        fastawrite(outfasta, ...
            sprintf('%s optimized by ChimeraMap', handles.target_name{1}), final_seq);
    end

    handles.target_seq(1) = [];
    handles.default_seq(1) = [];
    handles.target_name(1) = [];
    handles.chimera_regions(1) = [];
    handles.codon_regions(1) = [];
    handles.default_regions(1) = [];
    handles = update_figure(hObject, handles);
end
if handles.optSourceRef.Value
    handles.CUB = [];  % make sure we compute it next time as well
end
toc;
% msgbox(sprintf('%d sequences optimized.', ntarg), 'optimization', 'modal');

if isfield(handles, 'filt_action')
    handles = rmfield(handles, 'filt_action');
end
handles.during_map = false;
update_figure(hObject, handles);

if ~isempty(target_filter)
    listdlg('Name', 'target in reference', 'PromptString', handles.filt_msg, ...
            'ListString', target_filter, 'SelectionMode', 'single', ...
            'ListSize', [240, 300], 'CancelString', 'Cool');
end
if ~isempty(homolog_filter)
    listdlg('Name', 'homologs in reference', 'PromptString', ...
            'homologs suspected for the following targets:', ...
            'ListString', homolog_filter, 'SelectionMode', 'single', ...
            'ListSize', [240, 300], 'CancelString', 'Cool');
end
if ~isempty(target_error)
    listdlg('Name', 'cMap errors', 'PromptString', ...
            'the following errors occured during optimization:', ...
            'ListString', target_error, 'SelectionMode', 'single', ...
            'ListSize', [240, 300], 'CancelString', 'Bummer');
end


% --- Executes on button press in butTarget.
function butTarget_Callback(hObject, eventdata, handles)
% hObject    handle to butTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, dirname] = uigetfile({'*.fa;*.fasta', 'fasta file'}, 'select a sequence file');
if fname == 0
    return
end
try
    seqs = fastaread(fullfile(dirname, fname));
catch
    errordlg('invalid file: %s', fullfile(dirname, fname));
end

if length(seqs) > 1
    ls = arrayfun(@(x) {sprintf('%05d: %s', x, seqs(x).Header)}, 1:length(seqs));
    [ind, ok] = listdlg('ListString', ls, 'Name', 'select a target sequence', ...
                        'ListSize', [450, 300], 'SelectionMode', 'multiple');
    clear('ls');
    if ok == 0
        return
    end
    seqs = seqs(ind);
end
n_seq = length(seqs);

aa_alphabet = '*ACDEFGHIKLMNPQRSTVWY';
nt_alphabet = 'ACGTU';
for s = 1:n_seq
    valid_aa = ismember(upper(seqs(s).Sequence), aa_alphabet);
    valid_nt = ismember(upper(seqs(s).Sequence), nt_alphabet);

    if handles.optAA.Value && ~all(valid_aa)
        errordlg(sprintf('protein %s sequence contains illegal AA chars ("%s"). \nexpecting an amino acid sequence.', ...
            seqs(s).Header, unique(seqs(s).Sequence(~valid_aa))), 'target seq', 'modal');
        return
    end
    if handles.optNT.Value && ~all(valid_nt)
        errordlg(sprintf('protein %s sequence contains illegal NT chars ("%s"). \nexpecting a nucleotide sequence.', ...
            seqs(s).Header, unique(seqs(s).Sequence(~valid_nt))), 'target seq', 'modal');
        return
    end
    if ~all(valid_aa) && ~all(valid_nt)
        errordlg(sprintf('protein %s sequence contains illegal AA chars ("%s") \nand illegal NT chars ("%s"). \nexpecting an amino acid / nucleotide sequence.', ...
            seqs(s).Header, unique(seqs(s).Sequence(~valid_aa)), unique(seqs(s).Sequence(~valid_nt))), 'target seq', 'modal');
        return
    end
end

% sequence type logic
if handles.optAA.Value
    seq_type = 'AA';
elseif handles.optNT.Value
    seq_type = 'NT';
elseif all(valid_aa)
    seq_type = 'AA';
else
    seq_type = 'NT';
end

handles.target_seq = {};
handles.default_seq = {};
handles.target_name = {};
handles.target_NT_alpha = '';  % for alphabet assertion
handles.target_AA_alpha = '';
valid_len = [];
discard_stop = {};
added_stop = {};
for s = 1:n_seq
    seqs(s).Sequence = upper(seqs(s).Sequence);
    switch seq_type
        case 'AA'
            handles.optAA.Value = 1;
            handles.optNT.Value = 0;
            if seqs(s).Sequence(end) ~= '*'
                if ~exist('stop_decis', 'var')
                    stop_decis = questdlg('AA sequence is missing a STOP codon.', ...
                                          sprintf('%s STOP codon', seqs(s).Header), ...
                                          'add STOP', 'discard', 'ignore', 'add STOP');
                end
                if strcmp(stop_decis, 'add STOP')
                    seqs(s).Sequence(end+1) = '*';
                    added_stop{end+1} = seqs(s).Header;
                elseif strcmp(stop_decis, 'discard')
                    discard_stop{end+1} = seqs(s).Header;
                    continue
                end
            end
            handles.target_seq{end+1} = seqs(s).Sequence;
            handles.default_seq{end+1} = '';
            handles.target_name{end+1} = seqs(s).Header;
            handles.default_exist = false;
            valid_len = [valid_len; true];
        case 'NT'
            seqs(s).Sequence = strrep(seqs(s).Sequence, 'U', 'T');
            if nt2aa(seqs(s).Sequence(end-2:end)) ~= '*'
                if ~exist('stop_decis', 'var')
                    stop_decis = questdlg('NT sequence is missing a STOP codon.', ...
                                          sprintf('%s stop codon', seqs(s).Header), ...
                                          'add STOP', 'discard', 'ignore', 'add STOP');
                end
                if strcmp(stop_decis, 'add STOP')
                    seqs(s).Sequence(end+1:end+3) = aa2nt('*');  % random STOP
                    added_stop{end+1} = seqs(s).Header;
                elseif strcmp(stop_decis, 'discard')
                    discard_stop{end+1} = seqs(s).Header;
                    continue
                end
            end
            handles.optNT.Value = 1;
            handles.optAA.Value = 0;
            handles.default_seq{end+1} = seqs(s).Sequence;
            handles.target_seq{end+1} = nt2aa(handles.default_seq{end}, ...
                                              'AlternativeStartCodons', true);
            % NOTE: here we allow alternative starts so that the given
            %       target is translated correctly.
            handles.target_name{end+1} = seqs(s).Header;
            handles.default_exist = true;
            valid_len = [valid_len; mod(length(seqs(s).Sequence), 3) == 0];
        otherwise
            error('too many cooks!');
    end
    handles.target_AA_alpha = union(handles.target_AA_alpha, handles.target_seq{end});
    handles.target_NT_alpha = union(handles.target_NT_alpha, handles.default_seq{end});
end

if ~all(valid_len)
    [sel, decision] = listdlg('ListString', handles.target_name(~valid_len), ...
                              'Name', 'invalid length', 'ListSize', [450, 300], ...
                              'InitialValue', 1:sum(~valid_len), ...
                              'OKString', 'discard selected', 'CancelString', 'ignore warning', ...
                              'PromptString', 'the following seqs have lengths not divisible by 3.');
    if decision == 1  % 'discard selected'
        invalid_len = find(~valid_len);
        handles.default_seq(invalid_len(sel)) = [];
        handles.target_seq(invalid_len(sel)) = [];
        handles.target_name(invalid_len(sel)) = [];
    end
end
if ~isempty(discard_stop)
    listdlg('PromptString', sprintf('%d targets missing a STOP codon (discarded):', ...
                                    length(discard_stop)), ...
            'ListString', discard_stop, 'Name', 'target protein STOP', ...
            'ListSize', [240, 300], 'SelectionMode', 'single', 'CancelString', 'Cool');
end
if ~isempty(added_stop)
    listdlg('PromptString', sprintf('%d targets missing a STOP codon (added): \n%s', ...
                                    length(added_stop)), ...
            'ListString', added_stop ,'Name', 'target protein STOP', ...
            'ListSize', [240, 300], 'SelectionMode', 'single', 'CancelString', 'Cool');
end
if handles.reference_exist
    repres = ismember(handles.target_AA_alpha, handles.reference_AA_alpha);
    if ~all(repres)
        warndlg(sprintf('%d letters in the target AA alphabet are missing from the reference alphabet: \n%s', ...
                        sum(~repres), handles.target_AA_alpha(~repres)), 'AA alphabet assertion');
    end
    if strcmpi(seq_type, 'NT') && strcmpi(handles.reference_type, 'NT')
        repres = ismember(handles.target_NT_alpha, handles.reference_NT_alpha);
        if ~all(repres)
            warndlg(sprintf('%d letters in the target NT alphabet are missing from the reference alphabet: \n%s', ...
                            sum(~repres), handles.target_NT_alpha(~repres)), 'NT alphabet assertion');
        end
    end
end

handles.target_type = seq_type;
handles.target_file = fname;
handles.target_exist = true;
if n_seq > 1
    handles.fieldChimeraFrom.String = '1';
    handles.fieldChimeraTo.String = '-1';
    handles.fieldCodonFrom.String = '1';
    handles.fieldCodonTo.String = '-1';
end
handles.chimera_regions = cellfun(@(x) {[1, length(x)]}, handles.target_seq);
handles.codon_regions = cellfun(@(x) {[]}, handles.target_seq);

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


% --- Executes on button press in butReference.
function butReference_Callback(hObject, eventdata, handles)
% hObject    handle to butReference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, dirname] = uigetfile({'*.fa;*.fasta', 'fasta file'}, 'select a reference sequence file');
if fname == 0
    return
end

handles.statMap.String = 'loading reference'; drawnow;
try
    seqs = fastaread(fullfile(dirname, fname));
catch
    errordlg('invalid file: %s', fullfile(dirname, fname));
end

nt_alphabet = 'ACGTU';
aa_alphabet = '*ACDEFGHIKLMNPQRSTVWY';
if handles.optRefNT.Value
    alphabet = nt_alphabet;
elseif handles.optRefAA.Value
    alphabet = aa_alphabet;
else
    error('too many cooks!');
end
valid_seq = true(length(seqs), 1);
valid_len = true(length(seqs), 1);
for i = 1:length(seqs)
    seqs(i).Sequence = upper(seqs(i).Sequence);
    valid_seq(i) = all(ismember(seqs(i).Sequence, alphabet));
    valid_len(i) = mod(length(seqs(i).Sequence), 3) == 0;
    if handles.optRefNT.Value && valid_seq(i)
        seqs(i).Sequence = strrep(seqs(i).Sequence, 'U', 'T');
        seqs(i).aa_seq = nt2aa(seqs(i).Sequence, 'AlternativeStartCodons', false);  % here false is good
        % we will need the AA seq soon for alphabet assertion
    end
end

if ~all(valid_seq)
    listdlg('ListString', {seqs(~valid_seq).Header}, 'Name', 'invalid seqs', ...
            'ListSize', [450, 300], 'CancelString', 'Cool', ...
            'PromptString', 'the following seqs contain invalid letters and will be discarded.');
end
seqs = seqs(valid_seq);
valid_len = valid_len(valid_seq);

if ~all(valid_len) && handles.optRefNT.Value
    [sel, decision] = listdlg('ListString', {seqs(~valid_len).Header}, ...
                              'Name', 'invalid length', 'ListSize', [450, 300], ...
                              'InitialValue', 1:sum(~valid_len), ...
                              'OKString', 'discard selected', 'CancelString', 'ignore warning', ...
                              'PromptString', 'the following seqs have lengths not divisible by 3.');
    if decision == 1  % 'discard selected'
        invalid_len = find(~valid_len);
        seqs(invalid_len(sel)) = [];
    end
end

handles.reference_NT_alpha = '';
handles.reference_AA_alpha = '';
for i = 1:length(seqs)
    if handles.optRefNT.Value
        handles.reference_NT_alpha = union(handles.reference_NT_alpha, seqs(i).Sequence);
        handles.reference_AA_alpha = union(handles.reference_AA_alpha, seqs(i).aa_seq);
    else
        handles.reference_AA_alpha = union(handles.reference_AA_alpha, seqs(i).Sequence);
    end
end
if handles.target_exist
    repres = ismember(handles.target_AA_alpha, handles.reference_AA_alpha);
    if ~all(repres)
        errordlg(sprintf('%d letters in the target AA alphabet are missing from the reference alphabet: \n%s', ...
                         sum(~repres), handles.target_AA_alpha(~repres)), 'AA alphabet assertion');
        return
    end
    if strcmpi(handles.target_type, 'NT') && handles.optRefNT.Value
        repres = ismember(handles.target_NT_alpha, handles.reference_NT_alpha);
        if ~all(repres)
            errordlg(sprintf('%d letters in the target NT alphabet are missing from the reference alphabet: \n%s', ...
                             sum(~repres), handles.target_NT_alpha(~repres)), 'NT alphabet assertion');
            return
        end
    end
end

if handles.optRefNT.Value
    handles.reference_type = 'NT';
    handles.reference_seq = {seqs.Sequence}';
    handles.reference_aa = {seqs.aa_seq}';
    n_seq = length(handles.reference_seq);
else
    handles.reference_type = 'AA';
    handles.reference_seq = {};
    handles.reference_aa = {seqs.Sequence}';
    n_seq = length(handles.reference_aa);
end
handles.reference_name = {seqs.Header}';
handles.reference_file = fname;
handles.fileReference.String = sprintf('%s: [%d %s seqs]', handles.reference_file, n_seq, handles.reference_type);
handles.SA = zeros(0, 3);

if ~isempty(handles.reference_seq) || ~isempty(handles.reference_aa)
    handles.reference_exist = true;
else
    handles.reference_exist = false;
end
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
    handles.fieldCodonFrom.String = sprintf('%d', handles.codon_regions{1}(reg, 1));
    handles.fieldCodonTo.String = sprintf('%d', handles.codon_regions{1}(reg, 2));
end
guidata(hObject, handles);


% --- Executes on button press in butChimeraRemRegion.
function butChimeraRemRegion_Callback(hObject, eventdata, handles)
% hObject    handle to butChimeraRemRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_reg, err] = get_region(handles.fieldChimeraFrom, handles.fieldChimeraTo, ...
                            handles.target_seq{1});
if err
    return
end

[err, ~, ~, handles.chimera_regions{1}] = check_region(new_reg, handles.chimera_regions{1});
if err == 0
    errordlg('region does not exist', 'chimera regions', 'modal');
    return
end

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
                            handles.target_seq{1});
if err
    return
end

err = check_region(new_reg, handles.codon_regions{1});
if err == 2
    errordlg('new region already exists', 'codon regions', 'modal');
    return
elseif err == 1
    decision = questdlg('new region overlaps with an existing codon region' , ...
                        'region overlap', 'merge', 'cancel', 'cancel');
    if ~strcmp(decision, 'merge')
        return
    end
end

[err, ~, trunc_new, trunc_old] = check_region(new_reg, handles.chimera_regions{1});
if err
    decision = questdlg('new region overlaps with an existing chimera region', ...
                        'region overlap', 'keep chimera', 'keep codon', 'keep codon');
    switch decision
        case 'keep chimera'
            new_reg = trunc_new;
        case 'keep codon'
            handles.chimera_regions{1} = trunc_old;
        otherwise
            return
    end
end

[~, handles.codon_regions{1}] = check_region(new_reg, handles.codon_regions{1});

update_figure(hObject, handles);


% --- Executes on button press in butCodonRemRegion.
function butCodonRemRegion_Callback(hObject, eventdata, handles)
% hObject    handle to butCodonRemRegion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_reg, err] = get_region(handles.fieldCodonFrom, handles.fieldCodonTo, ...
                            handles.target_seq{1});
if err
    return
end

[err, ~, ~, handles.codon_regions{1}] = check_region(new_reg, handles.codon_regions{1});
if err == 0
    errordlg('region does not exist', 'codon regions', 'modal');
    return
end

update_figure(hObject, handles);


% --- Executes on button press in butChimeraOpts.
function butChimeraOpts_Callback(hObject, eventdata, handles)
% hObject    handle to butChimeraOpts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.winParams = winParams_macOS(handles.winParams);
update_figure(hObject, handles);


% --- Executes on button press in butRegions.
function butRegions_Callback(hObject, eventdata, handles)
% hObject    handle to butRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, dirname] = uigetfile({'*.reg', 'region file'}, 'select a regions file');
if fname == 0
    return
end
try
    regions = fastaread(fullfile(dirname, fname));
catch
    errordlg('invalid file: %s', fullfile(dirname, fname));
end

[target_in_file, itarg] = ismember(handles.target_name, {regions.Header});
if ~all(target_in_file)
    str_targets = cellfun(@(x) get_i0(strsplit(strrep(x, ',', ';'), ' '), 0), ...
                          handles.target_name(~target_in_file));
    listdlg('Name', 'regions file error', 'PromptString', ...
            sprintf('%d targets missing from file:', sum(~target_in_file)), ...
            'ListString', str_targets, 'SelectionMode', 'single', ...
            'ListSize', [240, 300], 'CancelString', 'Shame');
    return
end

lens_OK = cellfun(@length, {regions(itarg(itarg>0)).Sequence}) == ...
          cellfun(@length, handles.target_seq);
if ~all(lens_OK)
    str_targets = cellfun(@(x) get_i0(strsplit(strrep(x, ',', ';'), ' '), 0), ...
                          handles.target_name(~lens_OK));
    listdlg('Name', 'regions file error', 'PromptString', ...
            sprintf('%d targets with wrong length', sum(~lens_OK)), ...
            'ListString', str_targets, 'SelectionMode', 'single', ...
            'ListSize', [240, 300], 'CancelString', 'Shame');
    return
end

chars_OK = cellfun(@(x) all(ismember(upper(x), 'BCD')), {regions(itarg(itarg>0)).Sequence});
if ~all(chars_OK)
    str_targets = cellfun(@(x) get_i0(strsplit(strrep(x, ',', ';'), ' '), 0), ...
                          handles.target_name(~chars_OK));
    listdlg('Name', 'regions file error', 'PromptString', ...
            sprintf('%d bad region definitions. consult the user guide.', sum(~chars_OK)), ...
            'ListString', str_targets, 'SelectionMode', 'single', ...
            'ListSize', [240, 300], 'CancelString', 'Shame');
    return
end

handles.chimera_regions = {};
handles.codon_regions = {};
handles.default_regions = {};
for i = 1:length(handles.target_seq)
    if ~target_in_file(i)
        continue
    end
    handles.chimera_regions{i} = set2region(find(upper(regions(itarg(i)).Sequence) == 'B'));
    handles.codon_regions{i} = set2region(find(upper(regions(itarg(i)).Sequence) == 'C'));
    handles.default_regions{i} = set2region(find(upper(regions(itarg(i)).Sequence) == 'D'));
end
handles.fileRegions.String = fname;

update_figure(hObject, handles);


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


% --- Executes on button press in optSourceTable.
function optSourceTable_Callback(hObject, eventdata, handles)
% hObject    handle to optSourceTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optSourceTable
[fname, dirname] = uigetfile({'*.csv', 'comma separated values'}, 'select a codon score table');
if fname == 0
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
    T = textscan(fid, '%s %f', 'Delimiter', ',');
    T{3} = nt2aa(T{1}, 'AlternativeStartCodons', false);
    aa_OK = ismember(aa_list, T{3});
catch
    file_OK = false;
end
fclose(fid);

if ~file_OK
    errordlg('codon table format error. expecting a comma-separated file with 2 columns (and no header): [codon, score]', 'codon table', 'modal');
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

handles.CUB_file = fname;
handles.CUB = codonbias('');
aa_list = fieldnames(handles.CUB)';
for aa = aa_list
    iaa = find(strcmpi(T{3}, aminolookup(aa{1})))';
    for i = iaa
        icod = strcmpi(T{1}(i), handles.CUB.(aa{1}).Codon);
        handles.CUB.(aa{1}).Freq(icod) = T{2}(i);
    end
    % TODO: if we intend to add any algorithm which requires sum over aa to
    %       to be 1, need to add this scaling here.
end
guidata(hObject, handles);


% --- Executes on button press in optNT.
function optNT_Callback(hObject, eventdata, handles)
% hObject    handle to optNT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optNT
handles.optNT.Value = 1;
handles.optAA.Value = 0;


% --- Executes on button press in optNT.
function optAA_Callback(hObject, eventdata, handles)
% hObject    handle to optNT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of optNT
handles.optAA.Value = 1;
handles.optNT.Value = 0;


% --- Executes on button press in butARS.
function butARS_Callback(hObject, eventdata, handles)
% hObject    handle to butARS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 0. select output file
profiling_run = false;
[fname, dirname] = uiputfile({'*.csv', 'comma separated values'}, 'select output file');
if fname == 0
    return
elseif strcmp(fname, 'profiling.csv')
    profiling_run = true;
end
[~, fname, extent] = fileparts(fname);
if isempty(extent)
    extent = '.csv';
end
outfile = fullfile(dirname, fname);
outmain = strcat(outfile, extent);
fmain = fopen(outmain, 'w');
outreport = strcat(outfile, '_summary', extent);
frep = fopen(outreport, 'w');
tic;

handles.during_ars = true;
handles = update_figure(hObject, handles);
handles.statARS.String = 'building index'; drawnow;
win_params = handles.winParams;

% select ARS alphabet
if strcmpi(handles.reference_type, 'AA') || strcmpi(handles.target_type, 'AA')
    ars_type = 'AA';
else
    % NT may be converted to AA
    ars_type = {'NT', 'codon', 'AA'};
end
if iscell(ars_type)
    ars_type = questdlg('select an alphabet for Chimera ARS', 'cARS', ...
                        ars_type{:}, ars_type{1});
end
switch ars_type
    case 'AA'
        cars_ref = handles.reference_aa;
        fprintf(frep, 'target,cARS_total_%s,cARS_region\n', ars_type);
    case 'codon'
        cars_ref = nt2codon(handles.reference_seq);
        fprintf(frep, 'target,cARS_total_%s,cARS_region,codon_score\n', ars_type);
    case 'NT'
        cars_ref = handles.reference_seq;
        win_params.size = 3 * win_params.size;  % aa2nt
        win_params.center = 3 * win_params.center;
        win_params.max_len = 3 * win_params.max_len;
        fprintf(frep, 'target,cARS_total_%s,cARS_region,codon_score\n', ars_type);
    otherwise
        return
end
if handles.SA_exist && strcmpi(ars_type, handles.SA_type)
    % pass
else
    handles.SA = build_suffix_array(cars_ref, false);
    handles.SA_type = ars_type;
    handles.SA_exist = true;
end

if isempty(handles.CUB)
    [~, handles.CUB] = maximize_CUB('', handles.reference_seq);
end
w = CUB2weights(handles.CUB);

ntarg = length(handles.target_seq);
cars_score = nan(ntarg, 1);
cars_score_reg = nan(ntarg, 1);
codon_score = nan(ntarg, 1);
cars_names = cellfun(@(x) get_i0(strsplit(strrep(x, ',', ';'), ' '), 0), handles.target_name);
target_filter = {};
homolog_filter = {};
target_error = {};
for t = 1:ntarg
    % chimera ARS
    do_chimera = ~isempty(handles.chimera_regions{1});
    if ~do_chimera
        error('no chimera region defined.');
    end

    handles.statARS.String = sprintf('target %d', t); drawnow;

    chim_region = double(region2set(handles.chimera_regions{1}));
    switch ars_type
        case 'AA'
            cars_targ = handles.target_seq{1};
        case 'codon'
            cars_targ = nt2codon(handles.default_seq{1});
        case 'NT'
            cars_targ = handles.default_seq{1};
            chim_region = bsxfun(@minus, 3*chim_region, [2; 1; 0]);
            chim_region = chim_region(:);
        otherwise
            error('too many cooks!');
    end

    % check if target is already in reference
    [filtered, handles] = filter_target_from_ref(handles, cars_targ, cars_ref);
    if filtered
        switch handles.filt_action
            case 'filter from ref'
                target_filter{end+1} = sprintf('%05d: %s (filtered %d)', t, cars_names{t}, filtered);
                handles.statARS.String = sprintf('target %d (filtered %d)', t, filtered); drawnow;
            case 'skip target'
                target_filter{end+1} = sprintf('%05d: %s', t, cars_names{t});
                handles = update_figure(hObject, handles);
                continue
        end
    end

    if handles.winParams.size == 0
        [cars_score(t), cars_vec, err, homologs] = calc_cars(cars_targ, ...
            handles.SA, cars_ref, win_params.max_len, win_params.max_pos);
    else
        [cars_score(t), cars_vec, err, homologs] = calc_cars_posspec(cars_targ, ...
            handles.SA, cars_ref, win_params, ...
            win_params.max_len, win_params.max_pos);
    end

    if homologs
        homolog_filter{end+1} = sprintf('%05d: %s (filtered %d)', t, cars_names{t}, homologs);
    end
    if any(err)
        err_string = sprintf('%05d: %s', t, cars_names{t});
        if err(1)
            err_string = sprintf('%s [empty substring]', err_string);
        end
        if err(2)
            err_string = sprintf('%s [empty window]', err_string);
        end
        target_error{end+1} = err_string;
    end

    cars_vec = cars_vec(chim_region);
    cars_score_reg(t) = mean(cars_vec);
    if handles.default_exist
        codon_score(t) = calc_score_from_weights(w, handles.default_seq{1});
    end

    if ~profiling_run
        cars_str = strjoin(arrayfun(@(x) {int2str(x)}, cars_vec), ',');
        fprintf(fmain, sprintf('%s,%s\n', cars_names{t}, cars_str));
        fprintf(frep, sprintf('%s,%.4f,%.4f,%.4f\n', cars_names{t}, cars_score(t), cars_score_reg(t), codon_score(t)));
    end

    handles.target_seq(1) = [];
    handles.default_seq(1) = [];
    handles.target_name(1) = [];
    handles.chimera_regions(1) = [];
    handles.codon_regions(1) = [];
    handles.default_regions(1) = [];
    handles = update_figure(hObject, handles);
end
fclose(fmain);
fclose(frep);
toc;

if isfield(handles, 'filt_action')
    handles = rmfield(handles, 'filt_action');
end
handles.during_ars = false;
update_figure(hObject, handles);

if ~isempty(target_filter)
    listdlg('Name', 'target in reference', 'PromptString', handles.filt_msg, ...
            'ListString', target_filter, 'SelectionMode', 'single', ...
            'ListSize', [240, 300], 'CancelString', 'Cool');
end
if ~isempty(homolog_filter)
    listdlg('Name', 'homologs in reference', 'PromptString', ...
            'homologs suspected for the following targets:', ...
            'ListString', homolog_filter, 'SelectionMode', 'single', ...
            'ListSize', [240, 300], 'CancelString', 'Cool');
end
if ~isempty(target_error)
    listdlg('Name', 'cARS errors', 'PromptString', ...
            'the following errors occured during calculation:', ...
            'ListString', target_error, 'SelectionMode', 'single', ...
            'ListSize', [240, 300], 'CancelString', 'Bummer');
end


% --- Executes on button press in optRefNT.
function optRefNT_Callback(hObject, eventdata, handles)
% hObject    handle to optRefNT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.optRefNT.Value = 1;
handles.optRefAA.Value = 0;


% --- Executes on button press in optRefAA.
function optRefAA_Callback(hObject, eventdata, handles)
% hObject    handle to optRefAA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.optRefAA.Value = 1;
handles.optRefNT.Value = 0;
