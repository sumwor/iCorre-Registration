%%% getRefImg()
%PURPOSE: Generate average projection to be used as reference for batch 
%           registration of multi-stack calcium imaging data. Projection is 
%           generated from a segment of frames halfway through the batch of stacks.
%AUTHOR: MJ Siniscalchi (Yale University), 180718 
%
%INPUTS
%   mat_path: A cell array, each cell containing full path to a .MAT file 
%       containing an image stack stored as 'stack' (X x Y x nFrames uint16). 
%   stackInfo: A structure with fields,
%       'imageLength', in pixels, self-explanatory.
%       'imageWidth',  in pixels, self-explanatory.
%       'nFrames', Each element = number of frames in corresponding stack within the batch.
%   nFrames: Number of frames to average.
%
%OUTPUTS
%   ref_img: Mean projection from N (=nFrames) frames.
%
%--------------------------------------------------------------------------

function ref_img = getRefImg(mat_path,stackInfo,nFrames)

start_frame = round(0.5*sum(stackInfo.nFrames)-0.5*nFrames); %Get frames from ~halfway through session.
stackIdx(1) = find(cumsum(stackInfo.nFrames)>=start_frame,1,'first'); %idx of stack that contains first frame in segment for initial averaging
stackIdx(2) = find(cumsum(stackInfo.nFrames)>=start_frame+nFrames-1,1,'first'); %idx of stack that contains last frame
tempStack = zeros(stackInfo.imageHeight,stackInfo.imageWidth,nFrames); %to store all frames for averaging

f1_global = 1 + [0; cumsum(stackInfo.nFrames(1:end-1))]; %Global idx for first frame of each stack
jj = 1; %Frame counter for full stack used for averaging (var tempStack). 
for j = stackIdx(1):stackIdx(2)
    S = load(mat_path{j},'stack');
    frameIdx = [1 stackInfo.nFrames(j)];  %Local idx to frames for avg projection
    if j == stackIdx(1)
        frameIdx(1) = start_frame - f1_global(j) + 1;
    elseif j == stackIdx(2)
        frameIdx(2) = start_frame - f1_global(j) + nFrames;
    end
    tempStack(:,:,jj:jj+(frameIdx(2)-frameIdx(1))) = S.stack(:,:,frameIdx(1):frameIdx(2));
    jj = jj+(frameIdx(2)-frameIdx(1))+1;
end

ref_img = mean(tempStack,3); % initial template image


