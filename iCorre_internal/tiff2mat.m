%% tiff2mat()
%
% PURPOSE: To convert movement-corrected TIFF files into 3D arrays stored in MAT format.
%           Returns struct 'info' which contains selected info from ScanImage
%           header as well as the extracted tag struct for writing to TIF.
% AUTHOR: MJ Siniscalchi, 190826
%
% NOTES:
%       *Could try this Vidrio class for loading TIFs (might speed the process) 
%               % ScanImage Tiff reader
%               import ScanImageTiffReader.ScanImageTiffReader; %Requires download!
%               reader = ScanImageTiffReader(raw_path{1});
%               data = reader.data; %This is the image stack.
%--------------------------------------------------------------------------

function tags = tiff2mat(tif_paths, mat_paths, chan_number)

if nargin<3
   chan_number=[]; %For interleaved 2-color imaging; channel to convert.
end

tic;

w = warning; %get warning state
warning('off','all'); %TROUBLESHOOT: invalid ImageDescription tag from ScanImage

img_info = imfinfo(tif_paths{1}); %Copy info from first raw TIF
img_info = img_info(1);
fields_info = {'Height',      'Width',     'BitsPerSample','SamplesPerPixel'};
fields_tiff = {'ImageLength', 'ImageWidth','BitsPerSample','SamplesPerPixel'};
for i=1:numel(fields_info)
    tags.(fields_tiff{i}) = img_info.(fields_info{i}); %Assign selected fields to tag struct
end
tags.Photometric = Tiff.Photometric.MinIsBlack;
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tags.Software = 'MATLAB';

matcount = 0;

for i=1:numel(tif_paths)
    
    [pathname,filename,ext] = fileparts(tif_paths{i});
    source = [filename ext]; %Store filename of source file
    
    if mod(i, 1000) == 1
        % save mat files for every 1000 frame
        disp(['Converting ' source '...']);
        if i > floor(numel(tif_paths)/1000)*1000
            stack = zeros(img_info.Width, img_info.Height, numel(tif_paths)-floor(numel(tif_paths)/1000)*1000);
        else
            stack = zeros(img_info.Width, img_info.Height, 1000);
        end
        matcount = matcount + 1;
        framecount = 1;
    end    
   
    stacktemp = loadtiffseq(pathname,source); % load raw stack (.tif)
    if chan_number %Check for correction based on structural channel
        stacktemp = stacktemp(:,:,chan_number:2:end); %Just convert reference channel
    end
    stack(:,:,framecount) = stacktemp;
    framecount = framecount + 1;
    

    if mod(i, 1000) == 0 || i==numel(tif_paths)
        save(mat_paths{matcount},'stack','tags','source','-v7.3');
    end
end

%Console display
[pathname,~,~] = fileparts(mat_paths{1});

disp(['Stacks saved as .MAT in ' pathname]);
disp(['Time needed to convert files: ' num2str(toc) ' seconds.']);
warning(w); %revert warning state

