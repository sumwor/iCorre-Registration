function save_path = applyShifts_batch(path_names, save_dir, stackInfo, chan_number)

w = warning; %get warning state
warning('off','all'); %TROUBLESHOOT invalid ImageDescription tag from ScanImage

%Console display
disp('Applying shifts to second channel...'); 
S = load(path_names.mat{1},'options');
field_names = fieldnames(S.options);
disp('Hyperparameter sets:');
disp(field_names);

%Progress bar
h = waitbar(0,'Applying shifts...','Name','Progress'); 

save_path = cell(numel(path_names.mat),1); %Initialize output var
for i = 1:numel(path_names.mat)
    
    %Update progress bar
    temp = (i-1)/numel(path_names.mat);
    msg = ['Applying shifts...  (' num2str(temp*100,2) '%)'];
    waitbar(temp,h,msg);    
    
    %Load tiff and get channel for co-registration 
    [path,filename,ext] = fileparts(path_names.raw{i});
    stack = loadtiffseq(path,[filename ext]); % load raw stack (.tif)
    if ~isempty(chan_number) %Check for correction based on structural channel (work-in-progress...)
        stack = stack(:,:,chan_number:2:end); %Get single channel out of interleaved frames
    end
    
    %Load .MAT file and apply shifts from master registration
    S = load(path_names.mat{i},'options','sum_shifts');
    for j=1:numel(field_names)
        stack = apply_shifts(stack,S.sum_shifts.(field_names{j}),S.options.(field_names{j})); %apply shifts: apply_shifts(stack,shifts,options)
    end
  
    save_path{i} = fullfile(save_dir,strcat(field_names{end},'_',filename,'.tif'));
    saveTiff(stack,stackInfo.tags,fullfile(save_path{i})); %saveTiff(stack,img_info,save_path))
    
end

close(h); %Close waitbar
warning(w); %Revert warning state

