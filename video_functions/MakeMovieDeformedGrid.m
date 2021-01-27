function [Frames,VideoSettings] = MakeMovieDeformedGrid(Image,Points,Cells,Epochs,Op,Param,VideoSettings)
%MakeMovieDeformedGrid Make movie of the deforming grid.
% [Frames,VideoSettings] = MakeMovieDeformedGrid(Image,Points,Cells,Epochs,Op,Param,VideoSettings)
% Image (optional) contains images used for the PIV analysis
% Points contains the coordinates of vertices, as calculated by FollowPoint.
% Cells contains the connectivity and cell values as stretch, strain,
% strain type, etc.
% Epochs contains time information.
% Param contains parameters, with defaults set in SetDefaults.
% Param.GridColor specifies the grid color. Param.GridColor='none' hides
% the grid.
% Param.FigureWidth contains the figure width.
% Param.LabelFontSize specifies the font size for labels.
% Param.TitleFontSize specifies the font size for the title.
% Op, is an options structure. Defaults are set by SetDefaults
% Op.UseImages, to use or not use images in the background. If set to 1,
% images (when supplied) of the original object are shown below the
% deforming grid.
% VideoSettings, settings related to videos.
% VideoSettings.Quality sets the video quality, default = 100
% VideoSettings.FrameRate, frame rate, default = 4
%
%
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H


Frames=[];

%if ~isempty(Image.Dir)
% make video of deformed grid with images in the background


if strcmp(Param.GridColor,'none')
   disp('grid color is set to none, for this function it is changed to black')
   Param.GridColor='k';
end

if ~isfield(VideoSettings,'Quality')
    VideoSettings.Quality=100;
    
end


    % set to a default value
    if ~isfield(VideoSettings,'FrameRate')
        VideoSettings.FrameRate=4;
    end




nPoints = length(Points.u);


if Op.UseImages
    nImages = length(Image.Scaled);
    % do a check on the number of images and number of PIV solutions
    % npoints should be nimages minus 1
    if nPoints ~= nImages - 1
        nPoints
        nImages
        error('number of points not equal to number of images minus one')
    end
end
Debug=0;

% size of cell array

[ny,nx]=size(Cells.Midx);
aspectratio = ny/nx;
% make movies of separate images

clear Frame v

VideoName=strcat(Param.SaveDir,'/Deformed_grid_video','.mp4');

v = VideoWriter(VideoName,'MPEG-4');
% settings for video
v.FrameRate = VideoSettings.FrameRate;
v.Quality = VideoSettings.Quality;
open(v)

fig=figure;
set(fig,'Position',[0.1 0.1 Param.FigureWidth ceil(Param.FigureWidth*aspectratio)])


% firstplot=1;
pointsendx = Points.x{Epochs.Index(end)};
pointsendx(isinf(pointsendx)) = nan;
pointsendy = Points.y{Epochs.Index(end)};
pointsendy(isinf(pointsendy)) = nan;

Points.minx=min([pointsendx ; Points.x{Epochs.Index(1)}],[],'all');
Points.maxx=max([pointsendx ; Points.x{Epochs.Index(1)}],[],'all');
Points.miny=min([pointsendy ; Points.y{Epochs.Index(1)}],[],'all');
Points.maxy=max([pointsendy ; Points.y{Epochs.Index(1)}],[],'all');

% reference object
if Op.UseImages
    Image.RI = imref2d(size(Image.Scaled{1}));
    Image.RI.XWorldLimits = Image.xRange;
    Image.RI.YWorldLimits = Image.yRange;
end


% position for time stamp
xpostext = double((-Points.minx+Points.maxx)*0.025+Points.minx);
ypostext = double((-Points.miny+Points.maxy)*0.975+Points.miny);

for ii=1:nPoints
    figure(fig)
    clf
    % get time index
    itime=Epochs.Index(ii);
    
    if Op.UseImages
        % show image
        imshow(Image.Scaled{itime},Image.RI);
        hold on
    end
    %ishold
    
    if Op.UpwardYAxis==0
        set(gca,'YDir','reverse')
    elseif Op.UpwardYAxis==1
        set(gca,'YDir','normal')
    end
    
    
    % form grid vertices
    vertices=[Points.x{itime}(Cells.IndexingVertices) Points.y{itime}(Cells.IndexingVertices)];
    % plot deformed grid
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',NaN(length(Cells.Connectivity),1),'FaceColor','flat','EdgeColor',ones(3,1)*0.5)
    
    
    % and add time stamp
    if isfield(Epochs,'StartDate')
        timestr=strcat({'period: '},Epochs.StartDate{1},'-',Epochs.EndDate{itime});
    else
        timestr=strcat({'time: '},{num2str(Epochs.Time(itime),'%4.1f\n')},strcat({' '},Param.TimeUnit));
    end
    text((xpostext),(ypostext),timestr,'FontSize',Param.LabelFontSize)
    hold off
    
    % axes
    axis equal
    xlim([Points.minx Points.maxx]);
    ylim([Points.miny Points.maxy]);
    box on
    drawnow
    % make frame
    Frame = getframe(fig);
    % ii
    % (Frame)
    writeVideo(v,Frame);
    clf
    % proceed to next image
end
close(v)
close(fig)



end

