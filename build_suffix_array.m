function SA = build_suffix_array(S, unique_SA)
%SA = BUILD_SUFFIX_ARRAY(S)
%  the resulting array [SA] will contain:
%  1st column: suffix start index in the string
%  2nd column: string ID (when given multiple strings in cell array [S])
%  3rd column: suffix frequency in given reference [S]
%  4th column: suffix start index relative to END of string
%  5th column: SA row indices sorted by 1st column value.
%  6th column: SA row indices sorted by 4th column value.
%  run setup.m to compile the required MEX file.
%
% Alon Diament, June 2015.
%
% changelog:
% 09/09/15: using MEX-file implementation to greatly reduce memory load.
%           calling merge_arrays() hierarchically (merge-sort style),
%           to speed things up.
% 01/2018:  4th column added.
% 04/2018:  position sorting in 5th-6th columns.

%% INPUT OPTIONS

if nargin < 2
    unique_SA = false;
end

if ~iscell(S)
    SA = build_single_array(S);
else
    % handle cells containing multiple strings
    % (concatanate multiple arrays)
    tic;
    nS = length(S);
    SA = cell(nS, 1);
    for i = 1:nS
        if isempty(S{i})
            continue  % 03/03/17 avoiding error in mex-file
        end
        SA{i} = build_single_array(S{i});
        SA{i}(:, 3) = i;  % string id
        SA{i} = SA{i}(:, [1, 3, 2]);  % position, string id, frequency
        if ~mod(i, 1000)
            fprintf('string %d:\t%.1f sec\n', i, toc);
        end
    end
    fprintf('SA merging:\t%.1f sec\n', toc);
    SA = SA(~cellfun(@isempty, SA));  % 03/03/17 same as above
    % merging arrays
    while length(SA) > 1
        nS = length(SA);
        merged_SA = cell(ceil(nS/2), 1);
        for i = 2:2:nS
            merged_SA{i/2} = merge_arrays(SA{i-1}, SA{i}, S, unique_SA);
        end
        if i/2 < length(merged_SA)
            merged_SA(end) = SA(end);
        end
        SA = merged_SA;
    end
    SA = SA{1};
    fprintf('SA coord sort:\t%.1f sec\n', toc);

    % 07/04/18 for position-specific chimera
    lens = cellfun(@length, S);
    SA(:, 4) = SA(:, 1) - lens(SA(:, 2)) - 1;  % position from STOP, -1 at end of string
    [~, SA(:, 5)] = sort(SA(:, 1));  % sorted start
    [~, SA(:, 6)] = sort(SA(:, 4));  % sorted end

    fprintf('SA complete:\t%.1f sec\n', toc);
end
end


function SA = build_single_array(S)
% [S] is a string, [SA] - its suffix array

nS = length(S);
SA = cell(nS, 3);
for i = 1:nS % place i in string
    SA{i,1} = S(i:end);
    SA{i,2} = i;
    SA{i,3} = 1;
end

[~, ind] = sort(SA(:, 1));
SA = cell2mat(SA(ind, 2:3));
end
