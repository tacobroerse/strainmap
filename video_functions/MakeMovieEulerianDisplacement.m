function [Frames,VideoSettings] = MakeMovieEulerianDisplacement(PIVresults,Epochs,Param,Op,VideoSettings,PlotStrain)
%MakeMovieEulerianDisplacement Make movie of the eulerian displacements.
% [Frames,VideoSettings] = MakeMovieTensor(Points,Cells,Epochs,Param,Op,VideoSettings,PlotStrain)
% PIVresults contains the coordinates and eulerian displacements
% Epochs contains time information.
% Param contains parameters, with defaults set in SetDefaults.
% Param.ColorPercentile specifies the percentile of the deformation values
% that are used for the color map. Values in the range [>0,100]. A value of
% 100 corresponds to the maximum absolute values.
% Param.FigureWidth contains the figure width.
% Param.LabelFontSize specifies the font size for labels.
% Param.TitleFontSize specifies the font size for the title.
% VideoSettings, settings related to videos.
% VideoSettings.Quality sets the video quality, default = 100
% VideoSettings.FrameRate, frame rate, default = 4
% PlotStrain, type of tensor to plot
% - 'IncrementalDisplacement' incremental displacements
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

% make video of deformed grid with images in the background
Settings.FontSize=14;
% string interpreter
Interpreter = 'tex';

if ~isfield(VideoSettings,'Quality')
    VideoSettings.Quality=100;
end


if ~isfield(VideoSettings,'FrameRate')
    VideoSettings.FrameRate=4;
end




% color for missing data
NaNColor=0.25*[1 1 1];

nEpochs = Epochs.nfulltimes;


% make movies of separate images

% color map
[vik]=makescientificcolormap('vik');

% size of cell array
[ny,nx]=size(PIVresults.x{1});
aspectratio = ny/nx;

% axes maxima
Points.minx=min([PIVresults.x{Epochs.Index(1)}],[],'all');
Points.maxx=max([PIVresults.x{Epochs.Index(1)}],[],'all');
Points.miny=min([PIVresults.y{Epochs.Index(1)}],[],'all');
Points.maxy=max([PIVresults.y{Epochs.Index(1)}],[],'all');








clear Frame v
% string for saving figures
TensorTypeStr=PlotStrain;
VideoName=strcat(Param.SaveDir,TensorTypeStr,'_video','.mp4');
if ~isfolder(Param.SaveDir)
    mkdir Param.SaveDir
end
v = VideoWriter(VideoName,'MPEG-4');
% settings for video
v.FrameRate = VideoSettings.FrameRate;
v.Quality = VideoSettings.Quality;
open(v)

if aspectratio < 1
    m=2;n=1;
    FigureWidth=Param.FigureWidth;
    FigureHeight=ceil(Param.FigureWidth*aspectratio);
else
    m=1;n=2;
    FigureHeight=ceil(Param.FigureWidth*aspectratio/2);
    FigureWidth=Param.FigureWidth;
end

fig=figure;
set(fig,'Position',[0 0 FigureWidth FigureHeight])
%set(fig,'Visible','Off')

firstplot=1;
% whether to center color map at one

cmax=0;
for ii=1:nEpochs
    
    
    % get time index
    itime=Epochs.Index(ii);
    
    
    % selection of quantity to plot
    if strcmp(PlotStrain,'IncrementalDisplacement')
        
        % select displacements and apply mask
        if itime==1
            [fieldx,fieldy]=SelectDisplacements(PIVresults,Epochs,Op);
        end
        
        tensorstr={'\delta u','\delta v'};
        deformstr='Incremental displacement';
        
        cmaxi=prctile(abs([fieldx{itime}(:) ; fieldy{itime}(:)]),Param.ColorPercentile);
        if cmaxi > cmax
            % update cmax (color scale maximum)
            cmax=cmaxi;
        end
    else
        PlotStrain
        
        error('invalid strain option')
        
    end
end

% individual plot of principal strains
for ii=1:nEpochs
    clf % clear frame
    
    % get time index
    itime=Epochs.Index(ii);
    
    
    
    % coordinates
    xvec=unique(PIVresults.x{itime});
    yvec=unique(PIVresults.y{itime});
    
    
    for i=1:2
        
        colormap(vik)
        % strain xx
        subplot(m,n,i);
        if i==1
            field=flipud(fieldx{itime});
            cmaxi=cmax;
        elseif i==2
            field=flipud(fieldy{itime});
            cmaxi=cmax;
        end
        titlestr=tensorstr{i};
        % plot strain
        hold on
        box on
        % check empty cells
        imAlpha=ones(size(field));
        imAlpha(isnan(field))=0;
        
        imagesc(xvec,yvec,field,'AlphaData',imAlpha)
        if Op.UpwardYAxis==0
            set(gca,'YDir','reverse')
        else
            set(gca,'YDir','normal')
        end
        % box on
        axis equal
        xlim([Points.minx Points.maxx])
        ylim([Points.miny Points.maxy])
        if Op.UpwardYAxis==0
            set(gca,'YDir','reverse')
        end
        title(titlestr,'interpreter',Interpreter)
        colorbar
        caxis([-cmaxi  cmaxi])
        
        
        % set background color (for NANs)
        set(gca,'color',NaNColor);
        
        hold off
    end
    
    if firstplot
        % position for time stamp
        xpostext = 0.5;
        %  xpostext = (-Points.minx+Points.maxx)*0.025+Points.minx;
        %  ypostext = (-Points.miny+Points.maxy)*0.975+Points.miny;
        ypostext = 0.96;
        
        firstplot=0;
    end   % vertices=[Points.x{itime}(Cells.IndexingVertices) Points.y{itime}(Cells.IndexingVertices)];
    
    % and add time stamp
    if isfield(Epochs,'StartDate')
        timestr=strcat({'period: '},Epochs.StartDate{1},'-',Epochs.EndDate{itime});
    else
        timestr=strcat({'time: '},{num2str(Epochs.Time(itime),'%4.1f\n')},strcat({' '},Param.TimeUnitPlot));
    end
    
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    
    text(xpostext,ypostext,timestr,'FontSize',Settings.FontSize,'HorizontalAlignment','Center')
    
    drawnow
    % make frame
    Frame = getframe(fig);
    writeVideo(v,Frame);
    
    % proceed to next image
end
close(v)
close(fig)





end

