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

LineWidthGrid=0.2;
Frames=[];

%if ~isempty(Image.Dir)
% make video of deformed grid with images in the background


if strcmp(Param.GridColor,'none')
    warning('grid color is set to none')
%    Param.GridColor='k';
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
        warning('number of points not equal to number of images minus one')
        
        % checking which epochs from the grid correspond to which image
        % epochs
        disp('will now check which epochs of the images best fit with the points epochs')
        if isfield(Epochs,'EndDate')
           disp('using end dates, rather than start dates')
         %   UseEndDate=1;
            for itime=1:nPoints
            EpochsVec(itime)=datetime(str2num(Epochs.EndDate{itime}(1:4)),...
            str2num(Epochs.EndDate{itime}(5:6)),str2num(Epochs.EndDate{itime}(7:8)));

            end
        else
          %  UseEndDate=0;
            Epochs=Epochs.Time;
        end
        
        timethreshold = mean(diff(Epochs.Time))*0.5;
        for itime=1:nPoints-1
            difftime=1e10;
            Image.TimeCorrIndex(itime)=NaN;
            for jtime=1:nImages
                difftimeij=abs(Image.Dates(jtime)-EpochsVec(itime));
                if difftimeij < difftime
                    difftime=difftimeij;
                    if difftimeij < timethreshold
                        Image.TimeCorrIndex(itime)=jtime;
                    end
                end
            end
        end
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

if ~Op.UseFigureCoordLimits
    % use data coordinate limits
    Points.minx=min([pointsendx ; Points.x{Epochs.Index(1)}],[],'all');
    Points.maxx=max([pointsendx ; Points.x{Epochs.Index(1)}],[],'all');
    Points.miny=min([pointsendy ; Points.y{Epochs.Index(1)}],[],'all');
    Points.maxy=max([pointsendy ; Points.y{Epochs.Index(1)}],[],'all');
else
    % use specified coordinate limits
    Points.minx=Param.FigureXLim(1);
    Points.maxx=Param.FigureXLim(2);
    Points.miny=Param.FigureYLim(1);
    Points.maxy=Param.FigureYLim(2);
end


% reference object
if Op.UseImages
    Image.RI = imref2d(size(Image.Scaled{1}));
    if isfield(Image,'xRange')
        % pixel size
        Image.RI.XWorldLimits = Image.xRange;
        Image.RI.YWorldLimits = Image.yRange;
    elseif isfield(Image,'xvec')
        Image.RI.XWorldLimits = [Image.xvec(1) Image.xvec(end)];
        Image.RI.YWorldLimits = [Image.yvec(1) Image.yvec(end)];
    end
end


% position for time stamp
xpostext = double((-Points.minx+Points.maxx)*0.025+Points.minx);
ypostext = double((-Points.miny+Points.maxy)*0.975+Points.miny);

for ii=1:nPoints-1
    figure(fig)
    clf
    % get time index
    itime=Epochs.Index(ii);
    ii
    if Op.UseImages
        % show image
        if isfield(Image,'TimeCorrIndex')
            jtime=Image.TimeCorrIndex(itime);
        else
            jtime=itime;
        end
        
        imshow(Image.Scaled{jtime},Image.RI);
        hold on
    end
    %ishold
    
    if Op.UpwardYAxis==0
        set(gca,'YDir','reverse')
    elseif Op.UpwardYAxis==1
        set(gca,'YDir','normal')
    end
    
    
    % form grid vertices, take plus one for updated
    vertices=[Points.x{itime+1}(Cells.IndexingVertices) Points.y{itime+1}(Cells.IndexingVertices)];
    % plot deformed grid
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',NaN(length(Cells.Connectivity),1),'FaceColor','flat','EdgeColor',Param.GridColor,'LineWidth',LineWidthGrid)
    
    % plot every other x point with a red color
    steppoint=10;
    %plot(vertices([1:steppoint:end],1),vertices([1:steppoint:end],2),'.r')
    plot(Points.x{itime+1}(1:steppoint:end,1:steppoint:end),Points.y{itime+1}(1:steppoint:end,1:steppoint:end),'.r')
    % and add time stamp
    if isfield(Epochs,'StartDate')
        timestr=strcat({'period: '},Epochs.StartDate{1},'-',Epochs.EndDate{itime});
    else
        if isdatetime(Epochs.Time(itime))
            timegrid=datestr(Epochs.Time(itime));
        else
            timegrid=num2str(Epochs.Time(itime),'%4.1f\n');
        end
        timestr=['time grid: ' timegrid ' ' Param.TimeUnit];
        
        if isfield(Image,'TimeCorrIndex')
            if isdatetime(Image.Dates(jtime))
                timeimage=datestr(Image.Dates(jtime));
            else
                timeimage=num2str(Image.Dates(jtime),'%4.1f\n');
            end
            
            timestr=['time grid: ' timegrid ' time image: ' timeimage ' ' Param.TimeUnit];
            
        end
        
        
    end
    text((xpostext),(ypostext),timestr,'FontSize',Param.LabelFontSize)
    hold off
    
    % axes
    axis equal
    xlim([Points.minx Points.maxx]);
    ylim([Points.miny Points.maxy]);
    box on
    
    
    drawnow
    
    myFrame = getframe(fig);
    size(myFrame.cdata)
    
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

