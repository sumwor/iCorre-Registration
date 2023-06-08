% delete extra frames & restitch

root_dir = "Z:\HongliWang\Madeline\raw_imaging\JUV014";
temp = dir(root_dir); %Edit to specify data directories
data_dirs = {temp.name};
temp = ~(strcmp(data_dirs,'.') | strcmp(data_dirs,'..'));
data_dirs = data_dirs(temp); %remove '.' and '..' from directory list
disp('Directories for movement correction:');
disp(data_dirs');



%% Get directories and filenames
for i=1:numel(data_dirs)

    % Define & Create Subdirectories
    dirs.main = fullfile(root_dir,data_dirs{i});
    dirs.raw = fullfile(root_dir,data_dirs{i},'raw');
    dirs.save = fullfile(root_dir,data_dirs{i},'registered'); %to save registered stacks as TIFF
    dirs.stackInfo = fullfile(root_dir,data_dirs{i},'stack_info.mat');

    %% Load raw TIFs and convert to MAT for further processing, this is useful for aligning the imaging with behavior
    stackInfo = load(dirs.stackInfo);
    % check whether stack_info already exists

    nFrames = numel(dir(fullfile(dirs.raw,'*.tif')));
    stackInfo.nFrames=nFrames;
    stackInfo.imageWidth = stackInfo.tags.ImageWidth;
    stackInfo.imageHeight = stackInfo.tags.ImageLength;
    frames2save = mod(nFrames,1000);

    % load the last tiff in dirs.save
    if frames2save ~= 0
        saveFiles = dir(fullfile(dirs.save,'*.tif'));
        lastSave = fullfile(saveFiles(end).folder,saveFiles(end).name);

        stacktemp = loadtiffseq(saveFiles(end).folder,saveFiles(end).name);
        stack = stacktemp(:,:,1:frames2save);

        saveTiff(stack,stackInfo.tags,lastSave);

        %% stitch file
        disp('Getting global downsampled stack (binned avg.) and max projection of registered frames...');

        regfiles = saveFiles;
        for rr = 1:length(regfiles)
            paths.reg{rr} = fullfile(regfiles(rr).folder, regfiles(rr).name);
        end
        binnedAvg_batch(paths.reg,dirs.main,stackInfo,1);
    end
    clearvars '-except' dirs root_dir data_dirs file_names dirs stackInfo;

end


