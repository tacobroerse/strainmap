function [Image]=LoadImages(Image,Param,Op)
% load images that form the basis of the piv analysis
%
%   [Image]=LoadImages(Image,Param,Op)
%       Image, structure containing information on path of images:
%       Image.Dir
%       information on the file type:
%       Image.FileType
%       If a subselection of the available images is wanted, this range can
%       be given by:
%       Image.StartIndex and Image.EndIndex
%       For videofiles: 
%       Image.ProjectionFile can specify a geotif that contains a
%       projection, such that coordinates for the images can be constructed
%       For videofiles:
%       Image.EpochsFile can specify the file that contains dates for all
%       frames
%       and if Op.Synthetic = 1, a sythetic field is constructed.
%       Op, structure with  options. See SetDefaults
%       Op.FlipYImage specifies whether y-axis should be reversed
%       if the resolution of images should be reduced, specify this using:
%       Param.ImageScale (value < 1)
%
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H

% check whether input is video
if strcmp(Image.FileType,'mp4')
    VideoFile=1;
else
    VideoFile=0;
end
if ~isfolder(Image.Dir)
    Image.Dir
    error('image directory not found')
end

if ~isfield(Image,'FileType')
    disp('no file type set, set to default png')
    Image.FileType='png';
end

Image.Files=dir(fullfile(Image.Dir,strcat('*.',Image.FileType)));
if isempty(Image.Files)
    % try jpg
    error(strcat('no images found in folder:',Image.Dir))
end

if isfield( Image,'StartIndex')
    % start from this index on to load, to ensure that the same images are
    % loaded as in the PIV analysis
    disp(strcat('read images starting from index:',num2str(Image.StartIndex)))
else
    Image.StartIndex=1;
    disp('no start index provided, start to load images starting from the first image')
end

if isfield(Image,'nImages')
    % last image to load
    Image.EndIndex=Image.StartIndex+Image.nImages-1;
    disp(strcat('read image until index:',num2str(Image.EndIndex)))
else
    if (~VideoFile)
        disp('no number of images specified, assume all remaining images have to be used')
        nAllImages=length(Image.Files);
        Image.EndIndex=nAllImages;
        Image.nImages=nAllImages;
    end
end

if isfield(Image,'EpochIndex')
    if length(Image.EpochIndex)==Image.nImages
        
        disp('images will be assigned to the following epochs:')
        Image.EpochIndex
    else
        error('number of found images is not the same as length of Image.EpochIndex')
    end
else
    if isfield(Image,'nImages')
       Image.EpochIndex=[1:Image.nImages];
    end
end


% now load images
if Op.FlipYImage
    disp('flip images in y-direction')
end

if strcmp(Image.FileType,'mp4')
    % get images from video
    if length(Image.Files) > 1
        error('more than one video file found')
    end
    FileName=strcat(Image.Dir,'/',Image.Files.name);
    disp(strcat('read frames from video file:',FileName))
    video = VideoReader(FileName);
    
    
    % get number of frames
    nFrames=video.Duration;
    if ~isfield(Image,'EndIndex')
        Image.EndIndex=nFrames;
        Image.nImages=nFrames;
        Image.EpochIndex=[1:Image.nImages];
    end
    

    
    % read frame by frame
    for i=Image.StartIndex:Image.EndIndex
        j=Image.EpochIndex(i);
        Image1 = read(video,i);
        % rescale
        if Param.ImageScale>1
            disp('make sure Param.ImageScale is < 1')
            error('images should only be downscaled')
        elseif Param.ImageScale == 1
            if Op.FlipYImage
                % flip to make y-axis normal instead of reversed
                Image.Scaled{j}=flipdim(Image1,1);
            else
                Image.Scaled{j}=Image1;
            end
        else
            if Op.FlipYImage
                Image.Scaled{j}=flipdim(imresize(Image1,Param.ImageScale),1);
            else
                Image.Scaled{j}=imresize(Image1,Param.ImageScale);
            end
        end
        
        
        
    end
    
    % tiff image to retreive coordinates
    if isfield(Image,'ProjectionFile')
        % only store projection
        [sizeimy,sizeimx,~]=size(Image.Scaled{j});
        
        [~,R]=geotiffread(strcat(Image.Dir,'/',Image.ProjectionFile));
        
        % set coordinates using R
        stepx=R.RasterExtentInWorldX/sizeimx;
        stepy=R.RasterExtentInWorldY/sizeimy;
        % vector length is one smaller than raster size
        Image.xvec = [R.XWorldLimits(1)+stepx/2:stepx:R.XWorldLimits(2)-stepx/2];
        Image.yvec = [R.YWorldLimits(1)+stepy/2:stepy:R.YWorldLimits(2)-stepy/2];
    end
    
    if isfield(Image,'EpochsFile')
        % read epochs from file
        DatesStruct=load(strcat(Image.Dir,'/',Image.EpochsFile));
        Dates=struct2array(DatesStruct);
        for i=Image.StartIndex:Image.EndIndex
           j=Image.EpochIndex(i);
           Image.Dates(j,:)=Dates(i,:);
        end
        
        % check if table should be converted to array
        if istable(Image.Dates)
            Image.Dates=table2array(Image.Dates);
        end
        
    end
    
    
    
else
    % images in separate files
    j=0;
    for i=Image.StartIndex:Image.EndIndex
        j=Image.EpochIndex(i);
        %j=j+1; % separate counter
        FullFile=strcat(Image.Dir,Image.Files(i).name);
        ImCorrected{j}=imread(FullFile);
        % rescale
        if Param.ImageScale>1
            disp('make sure Param.ImageScale is < 1')
            error('images should only be downscaled')
        elseif Param.ImageScale == 1
            if Op.FlipYImage
                % flip to make y-axis normal instead of reversed
                Image.Scaled{j}=flipdim(ImCorrected{j},1);
            else
                Image.Scaled{j}=ImCorrected{j};
            end
        else
            if Op.FlipYImage
                Image.Scaled{j}=flipdim(imresize(ImCorrected{j},Param.ImageScale),1);
            else
                Image.Scaled{j}=imresize(ImCorrected{j},Param.ImageScale);
            end
        end
        % image size
        [ny,nx,~]=size(ImCorrected{1});
        Image.xRange=[1 nx];
        Image.yRange=[1 ny];
        
    end
end


