function binnedAvg_batch(stack_path,save_dir,stackInfo,bin_width)

tic; %Start timer

% Check data type
[pathname,filename,ext] = fileparts(stack_path{1});
stack_is_tif = strcmp(ext,'.tif');
stack_is_mat = strcmp(ext,'.mat');

% Initialize variables
curr = struct('stack',[]); %The current stack within 'stack_path'
next = struct('stack',[]); %The next stack within 'stack_path'

% Get indices for binned averaging
%Start frame for each segment to average
idx.DS_global = 1:bin_width:sum(stackInfo.nFrames); 
%Idx of first frame in each stack relative to first frame in session
idx.firstFrameGlobal = [0; cumsum(stackInfo.nFrames(1:end-1))]+1; 
%Get local bin indices for each stack
for i = 1:numel(idx.firstFrameGlobal)
    %Get idx of first frame in each averaging bin for the current stack
    firstFramelocal = idx.DS_global(find(idx.DS_global>=idx.firstFrameGlobal(i),1,'first'))...
        - idx.firstFrameGlobal(i) + 1; 
    idx.DS_local{i} = firstFramelocal:bin_width:stackInfo.nFrames(i); %Idx of each start frame relative to first frame in stack
end
idx.DS_local{end} = idx.DS_local{end}...
    (idx.DS_local{end}+bin_width-1 <= stackInfo.nFrames(end)); %remove last idx if < DS_factor from global last frame
idx.DS_global = idx.DS_global(1:numel([idx.DS_local{:}]));

%% Generate downsampled stack and max projection and save as TIF
stack_downsample = zeros(stackInfo.imageHeight,stackInfo.imageWidth,numel(idx.DS_global),'uint16'); %Initialize downsample stack

%Load first stack
if stack_is_tif
    fname = cell(numel(stack_path),1);
    for i = 1:numel(stack_path)
        [~,fname{i},~] = fileparts(stack_path{i});
    end
    curr.stack = loadtiffseq(pathname,strcat(fname{1},ext)); %Load first stack
elseif stack_is_mat
    curr = load(stack_path{1},'stack'); %Load first stack from *.MAT
else
    disp('Stack must be *.tif or saved as ''stack'' in a *.mat file.');
end

kk = 1; %Counter var for global idx 
for j = 1:numel(stack_path)
    disp(['Concatenating image frames and downsampling. Substack ' num2str(j) '/' num2str(numel(stack_path)) '...']);
    
    %Get local max for global max projection
    stack_max(:,:,j) = max(curr.stack,[],3);   
    
    if j<numel(stack_path) %Load next stack to pick up remainder of frames for averaging
        if stack_is_tif
            next.stack = loadtiffseq(pathname,strcat(fname{j+1},ext)); 
        else
            next = load(stack_path{j+1},'stack'); %Load next stack to pick up remainder of frames for averaging
        end
        curr.stack = cat(3,curr.stack,next.stack(:,:,1:bin_width)); %Concatenate with first n (=downsample_factor) frames of next_stack
    end
    
    for k = 1:numel(idx.DS_local{j})
        temp_idx = (idx.DS_local{j}(k) : idx.DS_local{j}(k) + bin_width - 1); %Segment to average
        stack_downsample(:,:,kk) = mean(curr.stack(:,:,temp_idx),3);
        kk = kk + 1;
    end
    
    curr = next;
end

saveTiff(stack_downsample,stackInfo.tags,fullfile(...
    save_dir,[filename(1:end-4) '_DS' num2str(bin_width) '.tif'])); %Save max projection for entire session
saveTiff(max(stack_max,[],3),stackInfo.tags,fullfile(save_dir,'reg_stackMax.tif')); %Save max projection for entire session

toc