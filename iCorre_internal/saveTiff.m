%%saveTiff(img_stack,img_info,save_path)
%
%PURPOSE: Saves imaging stacks from MATLAB arrays to TIF.
%AUTHOR: MJ Siniscalchi, 180507
%LAST EDIT: 180718
%
%INPUT ARGS:
%   'img_stack': Y x X x Time numeric array representing a stack of image frames
%   'img_info': A structure containing the tags (e.g., ImageHeight & ImageWidth)
%       to be attached to all frames of the TIF file. If tiff2mat.m was used
%       to generate MATLAB array, then struct 'stackInfo' will contain all necessary tags.
%   'save_path': (char) full save path with file name.
%
%EDITS:
%   180718mjs   If imfinfo() was used to get img_info, returns 'Height' and 'Width'
%               instead of valid Tiff tags; added translation of image info into valid tags.
%
%--------------------------------------------------------------------------

function saveTiff(img_stack,img_info,save_path)

%Get tags from imfinfo() structure or from get_stackInfo.m
if ~isfield(img_info,'ImageLength') || ~isfield(img_info,'ImageWidth') 
    img_info.ImageLength = img_info.Height; %imfinfo returns 'Height' and 'Width' instead of valid Tiff tags
    img_info.ImageWidth  = img_info.Width;
end
field_names = fieldnames(img_info);
field_names = field_names(ismember(field_names,Tiff.getTagNames)); %Trim off fields that are not valid tags

for i=1:numel(field_names)
    tags.(field_names{i}) = img_info.(field_names{i}); %Assign selected fields to tif tags
end
tags.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tags.Photometric = Tiff.Photometric.MinIsBlack;
tags.Software = 'MATLAB';

%Save stack as TIF
disp(['Saving stack as ' save_path '...']);
t = Tiff(save_path,'w8'); %open/create tif for writing
t.setTag(tags); %tag and write first frame
t.write(uint16(img_stack(:,:,1)));
for j = 2:size(img_stack,3) %create write dir, tag, and write subsequent frames
    t.writeDirectory();
    t.setTag(tags);
    t.write(uint16(img_stack(:,:,j)));
end
t.close();

end


