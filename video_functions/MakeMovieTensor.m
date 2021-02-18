function [Frames,VideoSettings] = MakeMovieTensor(Points,Cells,Epochs,Param,Op,VideoSettings,PlotStrain)
%MakeMovieTensor Make movie of a tensor field.
% [Frames,VideoSettings] = MakeMovieTensor(Points,Cells,Epochs,Param,Op,VideoSettings,PlotStrain)
% Points contains the coordinates of vertices, as calculated by FollowPoint.
% Cells contains the connectivity and cell values as stretch, strain,
% strain type, etc.
% Epochs contains time information.
% Param contains parameters, with defaults set in SetDefaults.
% Param.GridColor specifies the grid color. Param.GridColor='none' hides
% the grid.
% Param.ColorPercentile specifies the percentile of the deformation values
% that are used for the color map. Values in the range [>0,100]. A value of 
% 100 corresponds to the maximum absolute values.
% Param.FigureWidth contains the figure width.
% Param.LabelFontSize specifies the font size for labels.
% Param.TitleFontSize specifies the font size for the title.
% Op, is an options structure. Defaults are set by SetDefaults
% Op.Coordinates, whether to take the original coordinates
% (Op.Coordinates='Initial'), or the deformed coordinates
% (Op.Coordinates='Deformed').
% VideoSettings, settings related to videos.
% VideoSettings.Quality sets the video quality, default = 100
% VideoSettings.FrameRate, frame rate, default = 4
% PlotStrain, type of tensor to plot
% - 'IncrementalInfinitesimalStrain' incremental infinitesimal strain
% tensor, plus rotation angle from incremental rotation tensor
% - 'InfinitesimalStrain' infinitesimal strain
% tensor, plus rotation angle from infinitesimal rotation tensor
% - 'IncrementalGreenFiniteStrain' incremental finite Green strain tensor
% - 'GreenFiniteStrain' cumulative finite Green strain tensor
% PlotEpoch: epoch number for plot. Numeric or for the last epoch 'final'
% - 'LeftStretchV' left-stretch tensor and rotation angle from rotation
% tensor. Comes from decomposition F = VR (first rotation, then stretch)
% - 'RightStretchU' right-stretch tensor and rotation angle from rotation
% tensor. Comes from decomposition F = RU (first stretch, then rotation)
% - 'DisplacementGradient' incremental displacement gradient
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



% color for grid
EdgeColor=Param.GridColor;
% color for missing data
NaNColor=0.25*[1 1 1];

nEpochs = Epochs.nfulltimes;


% make movies of separate images

% color map
[bwg]=makecolormap('brownwhitegreen2');

% size of cell array
[ny,nx]=size(Cells.Midx);
aspectratio = ny/nx;
% axes maxima
pointsendx = Points.x{Epochs.Index(end)};
pointsendx(isinf(pointsendx)) = nan;
pointsendy = Points.y{Epochs.Index(end)};
pointsendy(isinf(pointsendy)) = nan;

Points.minx=min([pointsendx ; Points.x{Epochs.Index(1)}],[],'all');
Points.maxx=max([pointsendx ; Points.x{Epochs.Index(1)}],[],'all');
Points.miny=min([pointsendy ; Points.y{Epochs.Index(1)}],[],'all');
Points.maxy=max([pointsendy ; Points.y{Epochs.Index(1)}],[],'all');

CumulativeMeasure=0;
% check whether Cell quantities are Lagrangian or Eulerian
if ~isfield(Cells,'RefType')
    disp('no Cells.RefType found, assuming Lagrangian quantities')
    RefType='Lagrangian';
else
    if strcmp(Cells.RefType,'Lagrangian')
        RefType=Cells.RefType;
    elseif strcmp(Cells.RefType,'Eulerian')
        RefType=Cells.RefType;
    else
        error('unknown reference type')
    end
end

if strcmp(Op.Coordinates,'Initial')
    % nothing special, just take initial coordinates
    xvec=unique(Cells.Midx);
    yvec=unique(Cells.Midy);
elseif strcmp(Op.Coordinates,'Updated')
    % take the updated coordinates
    xvec=Points.x; % quadrilateral x coordinates
    yvec=Points.y; % quadrilateral y coordinates
end


clear Frame v
% string for saving figures
TensorTypeStr=PlotStrain;

    
if strcmp(RefType,'Eulerian')
    VideoName=strcat(Param.SaveDir,'Eulerian_',TensorTypeStr,'_video','.mp4');
else
    VideoName=strcat(Param.SaveDir,TensorTypeStr,'_video','.mp4');
end
if ~isfolder(Param.SaveDir)
    mkdir Param.SaveDir
end
v = VideoWriter(VideoName,'MPEG-4');
% settings for video
v.FrameRate = VideoSettings.FrameRate;
v.Quality = VideoSettings.Quality;
open(v)

fig=figure;
set(fig,'Position',[0.1 0.1 Param.FigureWidth ceil(Param.FigureWidth*aspectratio/1.3)])
%set(fig,'Visible','Off')
m=2;n=2;% size panels
firstplot=1;
% whether to center color map at one
CenterColorAtOne=0;

for ii=1:nEpochs
    clf
    
    % get time index
    itime=Epochs.Index(ii);
    
    
    % selection of quantity to plot
    if strcmp(PlotStrain,'IncrementalInfinitesimalStrain')
        fieldxx = squeeze(Cells.InfStrainIncrmt{itime}(1,1,:,:));
        fieldxy = squeeze(Cells.InfStrainIncrmt{itime}(1,2,:,:));
        fieldyy = squeeze(Cells.InfStrainIncrmt{itime}(2,2,:,:));
        fieldantixy = Cells.VorticityIncrmt{itime};
        tensorstr={'\epsilon_{xx}','\epsilon_{xy}','\omega','\epsilon_{yy}'};
        deformstr='incremental infinitesimal strain';
        if ii==1
            % color axes
            cmax=prctile(abs(Cells.InfStrainIncrmt{end}(:)),Param.ColorPercentile);
            cmaxantixy=prctile(abs(Cells.VorticityIncrmt{end}(:)),Param.ColorPercentile);
        end

    elseif strcmp(PlotStrain,'InfinitesimalStrain')
        fieldxx = squeeze(Cells.InfStrain{itime}(1,1,:,:));
        fieldxy = squeeze(Cells.InfStrain{itime}(1,2,:,:));
        fieldyy = squeeze(Cells.InfStrain{itime}(2,2,:,:));
        fieldantixy = Cells.Vorticity{itime};
        tensorstr={'\epsilon_{xx}','\epsilon_{xy}','\omega','\epsilon_{yy}'};
        deformstr='infinitesimal strain';
        if ii==1
            % color axes
            cmax=prctile(abs(Cells.InfStrain{end}(:)),Param.ColorPercentile);
            cmaxantixy=prctile(abs(Cells.Vorticity{end}(:)),Param.ColorPercentile);
        end
        CumulativeMeasure=1;
    elseif strcmp(PlotStrain,'IncrementalGreenFiniteStrain')
        fieldxx = squeeze(Cells.GreenStrainIncrmt{itime}(1,1,:,:));
        fieldxy = squeeze(Cells.GreenStrainIncrmt{itime}(1,2,:,:));
        fieldyy = squeeze(Cells.GreenStrainIncrmt{itime}(2,2,:,:));
        tensorstr={'E_{xx}','E_{xy}','','E_{yy}'};
        fieldantixy = NaN(size(fieldxx));
        deformstr='incremental Green finite strain';
        if ii==1
            % color axes
            cmax=prctile(abs(Cells.GreenStrainIncr{end}(:)),Param.ColorPercentile);
            cmaxantixy=0;
        end
    elseif strcmp(PlotStrain,'GreenFiniteStrain')
        fieldxx = squeeze(Cells.GreenStrain{itime}(1,1,:,:));
        fieldxy = squeeze(Cells.GreenStrain{itime}(1,2,:,:));
        fieldyy = squeeze(Cells.GreenStrain{itime}(2,2,:,:));
        fieldantixy = NaN(size(fieldxx));
        tensorstr={'E_{xx}','E_{xy}','','E_{yy}'};
        deformstr='Green finite strain';
        if ii==1
            % color axes
            cmax=prctile(abs(Cells.GreenStrain{end}(:)),Param.ColorPercentile);
            cmaxantixy=0;
        end
        CumulativeMeasure=1;
    elseif strcmp(PlotStrain,'LeftStretchV')
        fieldxx = squeeze(Cells.V{itime}(1,1,:,:));
        fieldxy = squeeze(Cells.V{itime}(1,2,:,:));
        fieldyy = squeeze(Cells.V{itime}(2,2,:,:));
        fieldantixy = Cells.FiniteRotAngle{itime};
        tensorstr={'V_{xx}','V_{xy}','\theta [deg]','V_{yy}'};
        deformstr='Left stretch tensor V';
        if ii==1
            % color axes
            for iii=1:nEpochs
            cmax(iii)=prctile(abs(Cells.V{iii}(:)),Param.ColorPercentile);
            cmaxantixy(iii)=prctile(abs(Cells.FiniteRotAngle{iii}(:)),Param.ColorPercentile);
            end
            cmax=max(cmax);
            cmaxantixy=max(cmaxantixy);
        end
        CumulativeMeasure=1;
        CenterColorAtOne=1;
    elseif strcmp(PlotStrain,'RightStretchU')
        fieldxx = squeeze(Cells.U{itime}(1,1,:,:));
        fieldxy = squeeze(Cells.U{itime}(1,2,:,:));
        fieldyy = squeeze(Cells.U{itime}(2,2,:,:));
        fieldantixy = Cells.FiniteRotAngle{itime};
        tensorstr={'U_{xx}','U_{xy}','\theta [deg]','U_{yy}'};
        deformstr='Right stretch tensor U';
        if ii==1
            % color axes
            cmax=prctile(abs(Cells.U{end}(:)),Param.ColorPercentile);
            cmaxantixy=prctile(abs(Cells.FiniteRotAngle{end}(:)),Param.ColorPercentile);
        end
        CenterColorAtOne=1;
        CumulativeMeasure=1;
    elseif strcmp(PlotStrain,'DisplacementGradient')
        fieldxx = squeeze(Cells.H{itime}(1,1,:,:));
        fieldxy = squeeze(Cells.H{itime}(1,2,:,:));
        fieldyy = squeeze(Cells.H{itime}(2,2,:,:));
        fieldantixy = squeeze(Cells.H{itime}(2,1,:,:));
        deformstr='Displacement gradient';
        tensorstr={'$\displaystyle\frac{du}{dx}$',...
            '$\displaystyle\frac{du}{dy}$',...
            '$\displaystyle\frac{dv}{dx}$',...
            '$\displaystyle\frac{dv}{dy}$'};
        Interpreter = 'latex';
        if ii==1
            % color axes
            cmax=prctile(abs(Cells.H{end}(:)),Param.ColorPercentile);
            cmaxantixy=cmax;
        end
    elseif strcmp(PlotStrain,'DisplacementGradientSum')
        fieldxx = squeeze(Cells.HSum{itime}(1,1,:,:));
        fieldxy = squeeze(Cells.HSum{itime}(1,2,:,:));
        fieldyy = squeeze(Cells.HSum{itime}(2,2,:,:));
        fieldantixy = squeeze(Cells.HSum{itime}(2,1,:,:));
        deformstr='Accumulated displacement gradient';
        tensorstr={'$\displaystyle\frac{du}{dx}$',...
            '$\displaystyle\frac{du}{dy}$',...
            '$\displaystyle\frac{dv}{dx}$',...
            '$\displaystyle\frac{dv}{dy}$'};
        Interpreter = 'latex';
        if ii==1
            % color axes
            cmax=prctile(abs(Cells.HSum{end}(:)),Param.ColorPercentile);
            cmaxantixy=cmax;
        end
        CumulativeMeasure=1;
    else
        PlotStrain
        error('invalid strain option')
    end
    
   
    
    
    
    for i=1:4
        box on
        colormap(bwg)
        % strain xx
        subplot(m,n,i);
        if i==1
            field=fieldxx;
            cmaxi=cmax;
        elseif i==2
            field=fieldxy;
            cmaxi=cmax;
        elseif i==4
            field=fieldyy;
            cmaxi=cmax;
        elseif i==3
            field=fieldantixy;
            cmaxi=cmaxantixy;
        end
        titlestr=tensorstr{i};
        % plot strain
        hold on
        % check empty cells
        imAlpha=ones(size(field));
        imAlpha(isnan(field))=0;
        if strcmp(Op.Coordinates,'Initial')
            imagesc(xvec,yvec,field,'AlphaData',imAlpha)
        elseif strcmp(Op.Coordinates,'Updated')
            % vertex array (2 columns)
            vertices=[xvec{itime}(Cells.IndexingVertices) yvec{itime}(Cells.IndexingVertices)];
            patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape(field',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor)
            
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

if CenterColorAtOne
    % only center diagonal terms
    if i==1 || i==4
        caxis([-cmaxi+1  cmaxi+1])
    else
    caxis([-cmaxi  cmaxi])
    end
else
    caxis([-cmaxi  cmaxi])
end
        
        
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
        if CumulativeMeasure
        timestr=strcat({'period: '},Epochs.StartDate{1},'-',Epochs.EndDate{itime});
        else
             timestr=strcat({'period: '},Epochs.StartDate{itime},'-',Epochs.EndDate{itime});
        end
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

