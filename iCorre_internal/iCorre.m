%%% iCorre (Iterative movement correction)
%
%Author: MJ Siniscalchi, using NoRMCorre developed by EA Pnevmatikakis (Flatiron Institute, Simons Foundation)
%Purpose: Iterative implementation of NoRMCorre (non-rigid movement correction)
%
%INPUT ARGUMENTS
% path_names, options_in, options_label, template_in, max_err, max_reps
%
%
%OUTPUTS
%
%   template_out:   Final mean projection after registration
%   nReps:          Total number of repeats used.
%
%SAVED VARIABLES
%
%   sum_shifts: Structure used for applying all shifts to another channel (both rigid and non-rigid)
%   options
%
%Edits:
%       180710 Arg 1, 'path_names' is now a cell array containing complete paths to the .mat files for each stack.
%       180714 There seemed to be a problem with col_shift being fed in with the options struct for nReps>1.
%               (Instead, try saving once for apply_shifts().
%--------------------------------------------------------------------------

function [ template_out, nReps, err_mat ] = iCorre(path_QC, path_names, options_in, options_label, template_in, max_err, max_reps )

%Initialize variables
err = max_err+1; %initialize arbitrarily > than max_err for first iteration
nReps = 0; %initialize

local_avg = zeros(size(template_in,1),size(template_in,2),numel(path_names),'uint16'); %matrix of mean frames taken from each stack; used for new template (=grand avg).

%Disable parfor warning
w = warning; %get warning state
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

while (max(err) > max_err && nReps < max_reps) %continue if i<max_rep
    
    nReps = nReps+1; %Increment # Repeats
    
    if exist('h','var')
        close(h);
    end
    h = waitbar(0,[upper(options_label) ' registration, iteration #' num2str(nReps) ' (0%)'],'Name','Progress');
    
    err = []; %reset value before each iteration
    k = 0; %counter var for indexing global err values.
    
    QCcount = 1;
    saveLength = numel(path_names);
    tCPU = zeros(numel(path_names), 10); tGPU = zeros(numel(path_names),10);
    for i=1:numel(path_names)
        
        % variables for quality control, saved every 100 trials to solve
        % the memory outage
        if mod(i,saveLength) == 1
            nnYc = zeros(1, saveLength); mmYc =zeros(1, saveLength); cYc = cell(1,saveLength); mYc = zeros(size(template_in,1),size(template_in,2),saveLength); vYc = zeros(1,saveLength);
            cMRc = cell(1, saveLength); mMRc = zeros(size(template_in,1),size(template_in,2),saveLength); vMRc = zeros(1, saveLength);
            cMNRc = cell(1, saveLength); mMNRc = zeros(size(template_in,1),size(template_in,2),saveLength); vMNRc = zeros(1, saveLength);
            Tc = zeros(1, saveLength);
            shifts_rc = cell(1, saveLength);
            shifts_nrc = cell(1, saveLength); shifts_xc = cell(1, saveLength); shifts_yc = cell(1, saveLength);
            
        end
        %Display within iteration as percent
        temp = (i-1)/numel(path_names);
        msg = [upper(options_label) ' registration, iteration #' num2str(nReps) '  (' num2str(temp*100,2) '%)'];
        waitbar(temp,h,msg);
        
        %Movement correction
        S = load(path_names{i}); % contains variables: 'options','stack','sum_shifts'
        if isfield(S,'options')
            options = S.options; %Load struct so field can be appended for new registration type
        end
        
        %% check the speed of cpu/gpu computing
        
if gpuDevicecCount == 1
         [stack,shifts,~ ,options.(options_label),~] = ...
            normcorre_batch_gpu(S.stack,options_in,template_in); %use parallel processing toolbox (parfor loop)
            
else
        [stack,shifts,~ ,options.(options_label),~] = ...
            normcorre_batch(S.stack,options_in,template_in); %use parallel processing toolbox (parfor loop)

end
        local_avg(:,:,i) = mean(stack,3); %take mean of each stack for later grand avg frame to be used as new template.

        %% quality control
        % correlation coefficients - downsample ten times
        % compute metrics
        
        nnY = quantile(S.stack(:),0.005); nnYc(i) = nnY;
        mmY = quantile(S.stack(:),0.995); mmYc(i) = mmY;

        [cY,mY,vY] = motion_metrics(S.stack,10); 
        cYc{i} = cY; mYc(:,:,i) = mY; vYc(i) = vY;   % downsample
        if strcmp(options_label,'RMC')
            [cMR,mMR,vMR] = motion_metrics(stack,10);
            cMRc{i} = cMR; mMRc(:,:,i) = mMR; vMRc(i) = vMR; 
        elseif strcmp(options_label,'NRMC')
            [cMNR,mMNR,vMNR] = motion_metrics(stack,10);
            cMNRc{i} = cMNR; mMNRc(:,:,i) = mMNR; vMNRc(i) = vMNR; 
        end
        T = length(cY);
        Tc(i) = T;
        
        % shifts
        if strcmp(options_label,'RMC')
            shifts_r = squeeze(cat(3,shifts(:).shifts));
            shifts_rc{i} = shifts_r;
        elseif strcmp(options_label,'NRMC')
            shifts_nr = cat(ndims(shifts(1).shifts)+1,shifts(:).shifts);
            shifts_nr = reshape(shifts_nr,[],ndims(S.stack)-1,T);
            shifts_nrc{i} = shifts_nr;
            shifts_x = squeeze(shifts_nr(:,1,:))'; shifts_xc{i} = shifts_x;
            shifts_y = squeeze(shifts_nr(:,2,:))'; shifts_yc{i} = shifts_y;
        end


        % save the quality control varibles every 100 trials
        if mod(i,saveLength) == 0
            
            if strcmp(options_label,'RMC')
                QCfile = fullfile(path_QC,['QC_RM_',num2str(QCcount)','_',num2str(nReps),'.mat']);
                save(QCfile,'nnYc','mmYc','cYc','mYc','vYc','cMRc','mMRc','vMRc','Tc','shifts_rc','options_label');

            elseif strcmp(options_label,'NRMC')
                QCfile = fullfile(path_QC,['QC_NRM_',num2str(QCcount),'_',num2str(nReps),'.mat']);
                save(QCfile,'nnYc','mmYc','cYc','mYc','vYc','cMNRc','mMNRc','vMNRc','Tc','shifts_nrc','shifts_xc','shifts_yc','options_label');
            end
            QCcount = QCcount + 1;
        end
        
        if nReps==1 %Save col_shift for correction of a second chan using apply_shifts
            save(path_names{i},'options','-append'); %***Note: may need to sum col_shift from each iteration for apply_shifts()
        end
                
        %Calculate values for iteration criteria
        field_names = fieldnames(shifts);
        for j=1:numel(shifts)
            dx = shifts(j).shifts(:,:,:,1); %translations by frame, dim-1; err(j) = norm([dx dy])
            dy = shifts(j).shifts(:,:,:,2); %translations by frame, dim-2
            err(j+k) = max(sqrt(dx(:).^2 + dy(:).^2)); %translation distance; if NRMC, use MAX of grid
            if ~isfield(S,'sum_shifts')
                sum_shifts.(options_label) = shifts; %initialize 
            elseif ~isfield(S.sum_shifts,options_label)
                sum_shifts = S.sum_shifts;
                sum_shifts.(options_label) = shifts; %initialize for new registration type
            else
                sum_shifts = S.sum_shifts;
                sum_shifts.(options_label) = S.sum_shifts.(options_label);
                for n = 1:numel(field_names)
                    sum_shifts.(options_label)(j).(field_names{n}) = ...
                        sum_shifts.(options_label)(j).(field_names{n}) + shifts(j).(field_names{n});
                end
            end
            
        end
        k = j+k; % increment by nFrames
        
        save(path_names{i},'stack','sum_shifts','options','-append'); %save running sum of local shifts in MAT file.
        
        clearvars S stack shifts sum_shifts
    end
    
    % plot the quality control files & save the data mat
    
    if nReps==1
        err_mat = NaN(numel(err),max_reps); %initialize output var (size depends on cumulative length of all stacks)
    end
    
    err_mat(:,nReps) = err; %matrix of translation errors (nGrids x nReps)
    template_in = mean(local_avg,3); %take grand avg as template for next iteration
    clearvars local_templates
    
end %end While loop


if exist('err_mat','var')
    err_mat = err_mat(:,~isnan(err_mat(1,:))); %Matrix of errors (nFrames x nRepeats); trim off unpopulated columns
    template_out = template_in; %Final (global) mean projection after registration
else
    err_mat = [];
    template_out = [];
end

close(h); %Close waitbar
warning(w); %Restore warning state






