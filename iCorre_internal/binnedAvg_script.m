%%% Stitch registered files for QC post-hoc 
%
%   -needed, e.g., when specified stitched downsampled file is too big for TIF format (~4GB)
%
%---------------------------------------------------------------------------------------------------
clearvars;

root_dir = 'J:\Data & Analysis';
bin_width = 1;

[fname,fpath] = uigetfile(fullfile(root_dir),'Select One or More Files to Stitch and Downsample','MultiSelect', 'on');
for i = 1:numel(fname)
stack_path{i} = fullfile(fpath,fname{i});
end

%Load stack_info.mat
save_dir = fullfile(fpath,'..\.'); %Parent of fpath
stackInfo = load(fullfile(save_dir,'stack_info.mat'));

binnedAvg_batch(stack_path,save_dir,stackInfo,bin_width);