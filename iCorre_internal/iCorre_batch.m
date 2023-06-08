%%% iCorre_batch
%PURPOSE: Script for iterative implementation of NoRMCorre (non-rigid movement correction),
%   with option to apply correction based on a separate anatomical reference channel.
%AUTHOR: MJ Siniscalchi, from tools developed by EA Pnevmatikakis (Flatiron Institute, Simons Foundation)
%DATE: 180718
%
%INPUTS
%   string root_dir:        Main directory containing all data directories to process.
%   string search_filter:   Wildcard string used to define data directory names within root_dir, e.g., '*RuleSwitching*'.
%   struct params:          Iterative movement correction parameters.
%       fields {'max_reps','max_err','nFrames_seed'}(see <<Set Parameters>> below).
%
%OUTPUTS
%   logical status:     Logical mask for successfully processed data directories.
%   cell msg:           Any associated error messages, indexed according to data directory.
%
%-----------------------------------------------------------------------------------------------------------

function [ status, err_msg ] = iCorre_batch(root_dir,params)


%% Set Parameters (if not included as input args)
if nargin<2
    params.max_reps = [1,1,1]; %maximum number of repeats; [seed, rigid, non-rigid]
    params.max_err = 1; %threshold abs(dx)+abs(dy) per frame
    params.nFrames_seed = 1000; %nFrames to AVG for for initial ref image
end

%% Get list of data directories

% if batch_pro is true, process multiply sessions under one root directory
% at once
% if batch_pro is false, process one session only 
if params.batch_pro
    temp = dir(root_dir); %Edit to specify data directories
    data_dirs = {temp.name};
    temp = ~(strcmp(data_dirs,'.') | strcmp(data_dirs,'..'));
    data_dirs = data_dirs(temp); %remove '.' and '..' from directory list
    disp('Directories for movement correction:');
    disp(data_dirs');
else
    data_dirs{1} = '';
    disp('Directories for movement correction:');
    disp(root_dir);
end


%% Setup parallel pool for faster processing
if isempty(gcp('nocreate'))
    try
        parpool
    catch err
        warning(err.message);
    end
end

%% Get directories and filenames
for i=1:numel(data_dirs)
        status(i) = true;
        err_msg{i} = [];
        % Define & Create Subdirectories
        dirs.main = fullfile(root_dir,data_dirs{i});
        dirs.raw = fullfile(root_dir,data_dirs{i},'raw');
        dirs.QC = fullfile(root_dir,data_dirs{i},'QualityControl');
        dirs.mat = fullfile(root_dir,data_dirs{i},'mat');
        dirs.save = fullfile(root_dir,data_dirs{i},'registered'); %to save registered stacks as TIFF
        if params.ref_channel~=params.reg_channel && params.do_stitch %If two-color co-registration
            dirs.save_ref = fullfile(root_dir,data_dirs{i},'registered_ref_channel'); %Save dir for registered reference channel
        end
        
        % If No 'raw' Directory, Create It & Move the TIF Files There
        if ~exist(dirs.raw,'dir')
            mkdir(dirs.raw);
            movefile(fullfile(dirs.main,'*.tif'),dirs.raw);
        end
        % Create remaining directories
        field_names = fieldnames(dirs);
        for j=1:numel(field_names)
            if ~exist(dirs.(field_names{j}),'dir')
                mkdir(dirs.(field_names{j}));
            end
        end
        
        % Define All Filepaths Based on Data Filename
        temp = dir(fullfile(dirs.raw,'*.tif'));
        file_names = cell(numel(temp),1);
        frameStep = 1000; % stitch 1000 frames, save it in a single mat file
        matcount = 1;
        for j=1:numel(temp)
            file_names{j} = temp(j).name;
            paths.raw{j} = fullfile(dirs.raw,file_names{j});
            if mod(j, frameStep) == 1
                paths.mat{matcount} = fullfile(dirs.mat,[file_names{j}(1:end-10), num2str(matcount), '.mat']); %for .mat file
                matcount = matcount + 1;
            end
        end
        
        % Additional Paths for Metadat
        paths.regData = fullfile(root_dir,data_dirs{i},'reg_info.mat'); %Matfile containing registration data
        paths.stackInfo = fullfile(root_dir,data_dirs{i},'stack_info.mat'); %Matfile containing image header info and tag struct for writing to TIF
        
        %Display Paths
        disp(['Data directory: ' dirs.raw]);
        %disp({['Path for stacks as *.mat files: ' dirs.mat];...
        %    ['Path for info file: ' paths.regData]});
        
        clearvars temp;
        
        %% Load raw TIFs and convert to MAT for further processing, this is useful for aligning the imaging with behavior

             % check whether stack_info already exists
             if ~exist(paths.stackInfo) 
                 disp('Converting *.TIF files to *.MAT for movement correction...');
                 %stackInfo = get_imgInfo(paths.raw, params); %Extract header info from image stack (written by ScanImage)
                 stackInfo.tags = tiff2mat(paths.raw, paths.mat, params.ref_channel); %Batch convert all TIF stacks to MAT and get info.
                 save(paths.stackInfo,'-STRUCT','stackInfo','-v7.3');
             else
                 stackInfo = load(paths.stackInfo);
             end
        % no need to convert to mat file for suit2p data
        % but do need the image height and width info
        % read the first tiff file
        tiff0 = fullfile(dirs.raw, file_names(1));
        t = Tiff(tiff0{1}, 'r');
        warning('off','all')
        tiffdata = read(t);
        warning
        stiff = size(tiffdata);
        stackInfo.imageWidth = stiff(1);
        stackInfo.imageHeight = stiff(2);
        stackInfo.nFrames = length(file_names);

        close(t);
        
        %% Correct RIGID, then NON-RIGID movement artifacts iteratively
        % Set parameters
        RMC_shift = 0.5*params.max_shift; %Params value used for SEED (based on reference image); therefter, narrow the degrees freedom
        NRMC_shift = 0.5*params.max_shift; %[10,10,0]; max dev = [5,5,0] for 256x256
        overlap = params.grid_size/4; %[16,16] for 256x256
        
        options.seed = NoRMCorreSetParms('d1',stackInfo.imageHeight,'d2',stackInfo.imageWidth,...
            'max_shift',params.max_shift,...
            'boundary','zero','upd_template',false,'use_parallel',true,'output_type','mat'); %Initial rigid correct for drift
        options.RMC = NoRMCorreSetParms('d1',stackInfo.imageHeight,'d2',stackInfo.imageWidth,...
            'max_shift',RMC_shift,...
            'boundary','zero','upd_template',false,'use_parallel',true,'output_type','mat'); %Rigid Correct; avg whole stack for template
        options.NRMC = NoRMCorreSetParms('d1',stackInfo.imageHeight,'d2',stackInfo.imageWidth,......
            'grid_size',params.grid_size,'overlap_pre',overlap,'overlap_post',overlap,...
            'max_shift',NRMC_shift,'max_dev',params.max_dev,'correct_bidir',false,...
            'boundary','zero','upd_template',false,'use_parallel',true,'output_type','mat'); %Non-rigid Correct; {'correct_bidir',true} threw error in iCorre.m on some data sets.
        
        % try it later
        options_nonrigid = NoRMCorreSetParms('d1',stackInfo.imageHeight,'d2',stackInfo.imageHeight,'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);

        % Initialize .mat file
        save(paths.regData,'params','-v7.3'); %so that later saves in the loop can use -append
        

        %% if dirs.QC is empty, run motion correction
        QCcontent = dir(dirs.QC);
         % if there are no quality assessment already, run motion correction
        % Iterative movement correction
        % modify this function
        template = getRefImg(paths.raw,stackInfo,params.nFrames_seed); %Generate initial reference image to use as template
        
        options_label = fieldnames(options);
       
        for m = find(params.max_reps) %seed, rigid, non-rigid
             if numel(QCcontent) == 2 
            tic;
            disp(' ');
            disp([upper(options_label{m}) ' registration in progress...']);
            
            [template,nReps,err_mat.(options_label{m})] = ...
                iCorre(dirs.QC, paths.mat,options.(options_label{m}),options_label{m},template,...
                params.max_err,params.max_reps(m));
            
            template_out.(options_label{m}) = template;
            nRepeats.(options_label{m}) = nReps;
            run_times.(options_label{m}) = toc;
            
            
            
            
            disp([options_label{m} ' correction complete. Elapsed Time: ' num2str(run_times.(options_label{m})) 's']);
            save(paths.regData,'template_out','nRepeats','run_times','-append'); %save values
            
            %% plot quality control figures
            plot_qualCtrl(dirs.QC, options_label{m});
             end
        end
        
        % display movement correction result analysis
        
        %Update regInfo and save
        field_names = fieldnames(options);
        options = rmfield(options,field_names(~params.max_reps));
        save(paths.regData,'options','-append');
        
        %% Save registered stacks as TIFF
        tic;
        
        %Apply registration to second channel if needed
        if params.ref_channel~=params.reg_channel
            paths.save = applyShifts_batch(paths,dirs.save,stackInfo,params.reg_channel); %Apply shifts and save .TIF files
            if params.do_stitch
                disp('Getting global downsampled stack (binned avg.) and max projection from reference channel...');
                binnedAvg_batch(paths.mat,dirs.save_ref,stackInfo,params.bin_width); %Save binned avg and projection to ref-channel dir
                if params.delete_mat
                        rmdir(dirs.mat,'s'); %DELETE .MAT dir...
                end
                disp('Getting global downsampled stack (binned avg.) and max projection of co-registered frames...');
                binnedAvg_batch(paths.save,dirs.main,stackInfo,params.bin_width); %Save binned avg and projection to main data dir
            end
            
        else
            %Save registered stacks as .TIF
            for k = 1:numel(paths.mat)
                if exist(paths.mat{k})
                    S = load(paths.mat{k},'stack'); %load into struct to avoid eval(stack_names{k})
                    saveTiff(S.stack,stackInfo.tags,fullfile(dirs.save,[options_label{m} '_' file_names{k}])); %saveTiff(stack,img_info,save_path))
                end

            end
            if params.do_stitch %Generate global and summary stacks for quality control
                disp('Getting global downsampled stack (binned avg.) and max projection of registered frames...');
                %if exist(paths.mat{1})
                %    binnedAvg_batch(paths.mat,dirs.main,stackInfo,params.bin_width); %Save binned avg and projection to main data dir
                %else
                    regfiles = dir(fullfile(dirs.save,'*.tif'));
                    for rr = 1:length(regfiles)
                        paths.reg{rr} = fullfile(regfiles(rr).folder, regfiles(rr).name);
                    end
                    binnedAvg_batch(paths.reg,dirs.main,stackInfo,params.bin_width);
                %end
            end
            if params.delete_mat
                rmdir(dirs.mat,'s'); %DELETE .MAT dir...
            end
        end
        
        clearvars S;
        run_times.saveTif = toc;
        save(paths.regData,'run_times','-append'); %save parameters
        
        clearvars '-except' dirs root_dir data_dirs file_names dirs paths params stackInfo options_label run_times status msg i m;
       
    
    if ~exist('err','var')
        status(i) = true;
        err_msg{i} = [];
    else
        status(i) = false;
        err_msg{i} = err;
        if exist(paths.regData)
            save(paths.regData,'err_msg','-append'); %save error msg
        else
            save(paths.regData,'err_msg'); 
        end
    end
    
    clearvars '-except' dirs root_dir data_dirs params i status err_msg;
    
    %end
end %end <<for i=1:numel(data_dirs)>>

