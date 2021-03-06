%%% getUserParams()
%
% PURPOSE: To get user inputs defining the hyperparameters for image registration in iCorre 
% AUTHOR: MJ Siniscalchi 200226
%
%---------------------------------------------------------------------------------------------------

function [ root_dir, params ] = getUserSettings(path_settings)

%% CHECK FOR SAVED SETTINGS

%If a path is specified, use those settings
if nargin < 1
    path_settings = string(fullfile(pwd,'user_settings'));
end

%% FIXED PARAMETERS
% maxRepSeed = 1; %Number of repetitions for ref. image (seed) fixed at 1


%% SET DEFAULTS AND EXTRACT SAVED SETTINGS

% These default settings worked well for most data acquired at 4 Hz with 256x256 pixels
settings = [...
    
    "Batch Directory (optional)",                                       "root_dir",         "";...
    "Batch Processing",                                                 "batch_pro",        "T";...
    "Save Settings As...",                                              "path_settings",    path_settings;...
    
    "ScanImage Version",                                                "scim_ver",         "5";...
    
    "Number of Frames to Average for Reference Image",                  "nFrames_seed",     "1000";...
    
    "Seed correction:",                                                 "maxRepSeed"        "0";...
    "Rigid Correction: Max. Number of Repeats",                         "maxRepRMC",        "1";...
    "Non-Rigid Correction: Max. Number of Repeats",                     "maxRepNRMC",       "0";...
        
    "Rigid Correction: Max. Shift (pixels)",                            "max_shift",        "10";...
    "Non-Rigid Correction: Max. Deviation from Rigid Shift (pixels)",   "max_dev",          "5";...
    "Non-Rigid Correction: Patch Width",                                "grid_size",        "16";...
        
    "Error Tolerance (pixels)",                                         "max_err",          "1";...
    "Reference Channel (0 for 1-color imaging)",                        "ref_channel",      "0";...
    "Follower Channel (0 for 1-color imaging)",                         "reg_channel",      "0";...
    "Save Concatenated Copy of Corrected Stacks? (T/F)",                "do_stitch",        "T";...
    "Downsample Factor for Concatenated Stack",                         "bin_width",        "20";...
    "Delete MAT files at Completion? (T/F)",                            "delete_mat",       "F";...
    
    ];

if exist(path_settings,'file')
    s = load(path_settings);
    settings = s.settings;
end

opts = struct('Resize','on','WindowStyle','modal','Interpreter','none');
settings(:,3) = string(inputdlg(settings(:,1),'*** iCORRE: Set Image Registration Parameters ***',1,settings(:,3),opts));

%% Save settings
save_path = settings(strcmp(settings(:,2),"path_settings"),3);
save(save_path,'settings');

%% Convert to correct MATLAB classes

% String or Numeric
for i = 1:size(settings,1)
    switch settings{i,2}
        case "root_dir" %Strings
            params.(settings{i,2}) = settings{i,3};
        case {"scim_ver","nFrames_seed","max_err","ref_channel","reg_channel","bin_width"} %Numeric
            params.(settings{i,2}) = str2double(settings{i,3});
        case {"grid_size","max_shift","max_dev"} %Numeric
            params.(settings{i,2}) = repmat(str2double(settings{i,3}),1,2);
        case "maxRepSeed"
            maxReps = str2double(settings(ismember(settings(:,2),["maxRepSeed","maxRepRMC","maxRepNRMC"]),3))';
            params.max_reps = maxReps; 
    end
    switch settings{i,3} %Logical
        case "T"
            params.(settings{i,2}) = true;
        case "F"
            params.(settings{i,2}) = false;
    end
end

root_dir        = params.root_dir;
