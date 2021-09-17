function [Frames,VideoSettings] = MakeMovieStrainSingleField(Points,Cells,Epochs,Param,Op,VideoSettings,PlotType)
%MakeMovieStrainSingleField Make movie of a scalar field.
% [Frames,VideoSettings] = MakeMovieStrainSingleField(Points,Cells,Epochs,Param,Op,VideoSettings,PlotType)
% Points contains the coordinates of vertices, as calculated by FollowPoint.
% Cells contains the connectivity and cell values as stretch, strain,
% strain type, etc.
% Epochs contains time information.
% Param contains parameters, with defaults set in SetDefaults.
% Param.nColorsStrainType specifies the number of discrete colors used for
% the strain type plot.
% Param.PercentileTresholdTypeStrainPlot specifies the strain magnitude
% below which strain types are reduced in transparency.
% Param.GridColor specifies the grid color. Param.GridColor='none' hides
% the grid.
% Param.FigureWidth contains the figure width.
% Param.LabelFontSize specifies the font size for labels.
% Param.TitleFontSize specifies the font size for the title.
% Op, is an options structure. Defaults are set by SetDefaults
% Op.IncludeTypeExpansionAndContraction: include biaxial extension or
% shortening, Op.IncludeTypeExpansionAndContraction == 0: color map from
% uniaxial extension to uniaxial shortening over strike-slip.
% Op.IncludeTypeExpansionAndContraction == 1; color map from biaxial to
% uniaxial extension to uniaxial to biaxial shortening over strike-slip.
% Op.Coordinates, whether to take the original coordinates
% (Op.Coordinates='Initial'), or the deformed coordinates
% (Op.Coordinates='Deformed').
% VideoSettings, settings related to videos.
% VideoSettings.Quality sets the video quality, default = 100
% VideoSettings.FrameRate, frame rate, default = 4
% PlotType, string: 'StrainType' (cumulative strain type) or 'IncrmtStrainType'
% (incremental strain type).
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


% frame rate, set to a default value
if ~isfield(VideoSettings,'FrameRate')
    VideoSettings.FrameRate=4;
end


% color for grid
EdgeColor=Param.GridColor;
% color for missing data
NaNColor=1*[1 1 1];

nEpochs = Epochs.nfulltimes;

Debug=0;
% make movies of separate images


% colors
% blue green orange red
if Op.IncludeTypeExpansionAndContraction
    nColors=Param.nColorsStrainType*2-1;
    
else
    nColors=Param.nColorsStrainType;
end

% colormaps
[bwg]=makecolormap('brownwhitegreen');
roma=makescientificcolormap('roma',nColors);
gwb=flipud(bwg);
romaO=makescientificcolormap('romaO');


% size of cell array

[ny,nx]=size(Cells.Midx);
aspectratio = ny/nx;
% % axes maxima
% Points.minx=min([Points.x{Epochs.Index(end)} ; Points.x{Epochs.Index(1)}],[],'all');
% Points.maxx=max([Points.x{Epochs.Index(end)} ; Points.x{Epochs.Index(1)}],[],'all');
% Points.miny=min([Points.y{Epochs.Index(end)} ; Points.y{Epochs.Index(1)}],[],'all');
% Points.maxy=max([Points.y{Epochs.Index(end)} ; Points.y{Epochs.Index(1)}],[],'all');

% axes maxima
pointsendx = Points.x{Epochs.Index(end)};
pointsendx(isinf(pointsendx)) = nan;
pointsendy = Points.y{Epochs.Index(end)};
pointsendy(isinf(pointsendy)) = nan;

Points.minx=min([pointsendx ; Points.x{Epochs.Index(1)}],[],'all');
Points.maxx=max([pointsendx ; Points.x{Epochs.Index(1)}],[],'all');
Points.miny=min([pointsendy ; Points.y{Epochs.Index(1)}],[],'all');
Points.maxy=max([pointsendy ; Points.y{Epochs.Index(1)}],[],'all');

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
TensorTypeStr=PlotType;
VideoName=strcat(Param.SaveDir,TensorTypeStr,'_video','.mp4');

v = VideoWriter(VideoName,'MPEG-4');
% settings for video
v.FrameRate = VideoSettings.FrameRate;
v.Quality = VideoSettings.Quality;
open(v)

fig=figure;
set(fig,'Position',[0.1 0.1 Param.FigureWidth ceil(Param.FigureWidth*aspectratio)])
%set(fig,'Visible','Off')



firstplot=1;
invertcolor = 0;
CenterColorAtOne=0;
CyclicColor=0;

for ii=1:nEpochs
    % clear frame
    clf
    
    box on
    daspect([1 1 1])
    
    
    % get time index
    itime=Epochs.Index(ii);
    
    
    % selection of quantity to plot
    if strcmp(PlotType,'StrainType')
        field = (Cells.StrainType{itime});
        fieldmagnitude = Cells.MagnitudeStrain{itime};
        fieldmagnitudefinal = Cells.MagnitudeStrain{Epochs.Index(end)};
        if ii==1
            deformstr='type of strain';
            titlestr='type of strain';
            % color axes
            if Op.IncludeTypeExpansionAndContraction
                nColors=Param.nColorsStrainType*2-1;
                cmax=pi/2;
            else
                nColors=Param.nColorsStrainType;
                cmax=pi/4;
            end
        end
    elseif strcmp(PlotType,'IncrmtStrainType')
        field = (Cells.StrainTypeIncrmt{itime});
        fieldmagnitude = Cells.MagnitudeStrainIncrmt{itime};
        fieldmagnitudefinal = Cells.MagnitudeStrainIncrmt{Epochs.Index(end)};
        if ii==1
            deformstr='type of strain';
            titlestr='type of strain';
            % color axes
            if Op.IncludeTypeExpansionAndContraction
                nColors=Param.nColorsStrainType*2-1;
                cmax=pi/2;
            else
                nColors=Param.nColorsStrainType;
                cmax=pi/4;
            end
        end
    elseif strcmp(PlotType,'MeanRotation')
        
        field = squeeze(Cells.FiniteRotAngle{itime});
        if ii==1
            titlestr='average material rotation [deg]';
            deformstr='rotation';
            cmax=max(abs(squeeze(Cells.FiniteRotAngle{end}(:))));
        end
        
    elseif strcmp(PlotType,'Dilatation')
        CenterColorAtOne=1;
        field = real(Cells.Dilatation{itime});
     
        if ii==1
            titlestr='dilatation';
            deformstr='dilatation';
            cmax=prctile(abs(squeeze(Cells.Dilatation{end}(:))-1),Param.ColorPercentile);
        end
        
        invertcolor = 1;
    elseif strcmp(PlotType,'StrainDirection')
        field = (Cells.DominantPrincplAngle{itime});
        fieldmagnitude = (Cells.MagnitudeStrain{itime});
        fieldmagnitudefinal = Cells.MagnitudeStrain{Epochs.Index(end)};
        
        if ii==1
            deformstr='dominant strain direction';
            titlestr='dominant strain direction';
            % color axes
            cmax=90;
            CyclicColor=1;
        end
    elseif strcmp(PlotType,'MagnitudeStrain')
       
        field = (Cells.MagnitudeStrain{itime});
       
        if ii==1
            deformstr='strain magnitude';
            titlestr='largest principal Hencky strain';
            % color axes
            cmax=prctile(abs(squeeze(Cells.MagnitudeStrain{end}(:))),Param.ColorPercentile);
           
        end
    else
        PlotType
        error('invalid strain option')
    end
    
    
    
    
    % plot
    % check empty cells
    if strcmp(PlotType,'StrainType') || strcmp(PlotType,'IncrmtStrainType') || strcmp(PlotType,'StrainDirection')
        % normalisation factor, using percentile. Everything above a certain
        % treshold will be opaque, below will be gradually transparant
        if ~isfield(Param,'PercentileTresholdTypeStrainPlot')
            Param.PercentileTresholdTypeStrainPlot=95;
        end
        
        PercentileTreshold=Param.PercentileTresholdTypeStrainPlot;
        % normalisation strain (based on final epoch)
        normStrain=prctile(abs(fieldmagnitudefinal(~isnan(fieldmagnitudefinal))),PercentileTreshold,'all');
        
        if firstplot
            disp(strcat('reducing transparancy for largest strain eigenvalues smaller than:',num2str(normStrain,2)))
        end
        % make transparancy mask
        
        AlphaPower=1;
        imAlpha=(abs(fieldmagnitude)/normStrain).^AlphaPower;
        % set everything above the treshold to 1
        imAlpha(imAlpha>1)=1;
        imAlpha(1,1)=1; % otherwise it does not work
    else
        % only use nans
        imAlpha=ones(size(field));
        imAlpha(isnan(field))=0;
        
    end
    %hold on
    if strcmp(Op.Coordinates,'Initial')
        imagesc(xvec,yvec,field,'AlphaData',imAlpha)
    elseif strcmp(Op.Coordinates,'Updated')
        % vertex array (2 columns)
        vertices=[xvec{itime}(Cells.IndexingVertices) yvec{itime}(Cells.IndexingVertices)];
        patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape((field)',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor,'FaceVertexAlphaData',reshape(imAlpha',nx*ny,1),'FaceAlpha','flat')
        
    end
    
    
    xlim([Points.minx Points.maxx])
    ylim([Points.miny Points.maxy])
    
    
    if Op.UpwardYAxis==0
        set(gca,'YDir','reverse')
    end
    title(titlestr,'interpreter',Interpreter,'FontSize',Param.TitleFontSize)
    
    % color map
    if strcmp(PlotType,'StrainType') || strcmp(PlotType,'IncrmtStrainType')
        colormap(roma)
    elseif CyclicColor
        colormap(romaO)
    else
        if invertcolor
            colormap(gwb)
        else
            colormap(bwg)
        end
    end
    % strain type has different ticks
    if strcmp(PlotType,'StrainType') || strcmp(PlotType,'IncrmtStrainType')
        colorbar('Ticks',[ -pi/4 , 0, pi/4],...
            'TickLabels',{'Shortening','Strike-Slip','Extension'},...
            'FontSize',Param.LabelFontSize)
    else
        colorbar
    end
    
    if CenterColorAtOne
        % center colorscale at 1 (for i.e. dilatation)
        caxis([-cmax+1 cmax+1])
    else
        caxis([-cmax cmax])
    end
    
    
    % set background color (for NANs)
    set(gca,'color',NaNColor);
    hold off
    
    if firstplot
        % position for time stamp
        xpostext = 0.5;
        ypostext = 0.96;
        
        firstplot=0;
    end
    
    % and add time stamp
    
    if isempty(Param.TimeUnit)
        timestr=strcat({' time: '},{num2str(Epochs.Time(itime))});
    else
        if isdatetime(Epochs.Time(itime))
            timestr= strcat({datestr(Epochs.Time(itime))});
        else
            timestr=strcat({' time: '},{num2str(Epochs.Time(itime),'%4.1f\n')},{' '},{Param.TimeUnit});
        end
    end
    
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    
    text(xpostext,ypostext,timestr,'FontSize',Settings.FontSize,'HorizontalAlignment','Center','FontSize',Param.TitleFontSize)
    
    
    
    drawnow
    % make frame
    Frame = getframe(fig);
    
    writeVideo(v,Frame);
    
    % proceed to next image
end
close(v)
close(fig)





end

