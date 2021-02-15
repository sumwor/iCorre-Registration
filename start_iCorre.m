%% start_iCorre
%
%Purpose: Script to process any number of 1- or 2-channel imaging stacks. Included
%   functions correct rigid and non-rigid movement artifacts using a flexible
%   recursive approach.
%
%Author: MJ Siniscalchi, 171212
%
%SETUP:
%       (1) To determine the needed files/MATLAB packages, run these three lines in the console:
%           [fList,pList] = matlab.codetools.requiredFilesAndProducts('start_batchProcessing.m');
%           {fList{:}}'
%           {pList.Name}'
%       (2) Download and install the necessary MATLAB components. Local toolboxes are included in 
%               the repository found at https://github.com/michaelsiniscalchi/iCorre-Registration
%       (3) Run this script to begin. Set hyperparameters using the dialog box.  
%
%   NOTES: 
%       *BATCH PROCESSING: iCorre_batch can register image stacks from multiple data directories back-to-back.
%       *Specify batch processing by entering full path to batch directory (parent dir. to data dirs.).
%       *Set SEARCH FILTER to specify a subset of data directories within batch directory.
%
%EDITS:
%   180709mjs Began rewrite for registration directly from raw TIF files 
%           (previous version required concatenation prior to movement correction)
%
%--------------------------------------------------------------------------
clearvars;

%% Set Hyperparameters

% Add Paths to Local Toolboxes
%addpath(genpath(pwd)); %Location of iCorre Registration directory

% Get User Input (See comments below for more explanation of parameter settings.)
path_settings = 'E:\labcode\iCorre-Registration-master\iCorre-Registration-master\user_settings.mat'; %Default user settings file; edit to specify eg a path within your data directory hierarchy
[root_dir, params] = getUserSettings(path_settings);

% If Batch Processing not Specified
if isempty(root_dir)
    data_dir = uigetdir(path_settings,'Please Specify Data Directory');
    [root_dir, data_dir, ~] = fileparts(data_dir); %Specify a "batch" of one data directory
    search_filter = ['*' data_dir '*']; 
end

%% Recursive Movement Correction (optional batch processing)

% Main iCorre Wrapper
[status,msg] = iCorre_batch(root_dir, params); %iCorre_batch(root_dir,search_filter,params)

