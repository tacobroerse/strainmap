function PlotStrainType(Cells,Param,Op,Points,PlotType,PlotEpoch,Epochs,Fault,dxfstruct)
% PlotStrainType(Cells,Param,Op,Points,PlotType,PlotEpoch,Epochs,Fault,dxfstruct)
%PlotStrainType plots the type of strain (shortening, strike-slip,
% extension) from the logarithmic principal stretches, i.e. the Hencky
% strains, the logarithm of the principal stretches V.
% PlotStrainType(Cells,Param,Op,Points,PlotType,PlotEpoch,Epochs)
% Cells contains the strain type and connectivity.
% Param contains parameters, with defaults set in SetDefaults.
% Param.nColorsStrainType specifies the number of discrete colors used for
% the plot.
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
% Points contains the coordinates of vertices, as calculated by FollowPoint.
% PlotType, string: 'StrainType' (cumulative strain type) or 'IncrmtStrainType'
% (incremental strain type).
% PlotEpoch: epoch number for plot. Numeric or for the last epoch 'final'
% Epochs contains time information.
%PlotStrainType(Cells,Param,Op,Points,PlotType,PlotEpoch,Epochs,Fault,dxfstruct)
% allows to overlay fault lines
% Fault.xlim and Fault.ylim specify the coordinate limits that correspond
% to the dxfstruct containing all lines.
%
%
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H

if isempty(PlotEpoch)
    PlotEpoch = 'final';
end
if nargin < 8
    OverlayPlot = 0;
else
    OverlayPlot= 1;
    if nargin == 9
        DXF=1;
    else
        DXF=0;
    end
end

% blue green orange red
if Op.IncludeTypeExpansionAndContraction
    nColors=Param.nColorsStrainType*2-1;
    colorlims=[-pi/2 pi/2];
else
    nColors=Param.nColorsStrainType;
    colorlims=[-pi/4 pi/4];
end

% use scientific color map roma
roma=makescientificcolormap('roma',nColors);
vikO=makescientificcolormap('vikO');
romaO=makescientificcolormap('romaO');
cyclic_kovesi_c1=makescientificcolormap('cyclic_kovesi_c1');
cyclicmap=romaO;

% color for grid
OpColorbar=1;
EdgeColor=Param.GridColor;

% color for missing data
NaNColor=1*[1 1 1];

% size of cell array
indexval=find(~isnan(Points.x{end}) & ~isinf(Points.x{end}));
xlims=[min([Points.x{end}(indexval); Points.x{1}(indexval)]) max([Points.x{end}(indexval); Points.x{1}(indexval)])];
ylims=[min([Points.y{end}(indexval); Points.y{1}(indexval)]) max([Points.y{end}(indexval); Points.y{1}(indexval)])];

xrange=xlims(2)-xlims(1);
yrange=ylims(2)-ylims(1);
aspectratio = xrange/yrange;
% number of grid cells
[nx,ny]=size(Cells.Midx);

% selection of epoch of field
if ~isnumeric(PlotEpoch)
    if strcmp(PlotEpoch,'final')
        fileepochstr='final';
        itime=length(Cells.F);
        epochstr=' final epoch';
    else
        disp(strcat('Plot Epoch:',PlotEpoch))
        disp('choose from: final (default) or numeric epoch')
        
        error('unrecognized option')
    end
else
    itime=PlotEpoch;
    fileepochstr=strcat('epoch_',num2str(itime));
    epochstr=cell2mat(strcat({' epoch '},{num2str(itime)}));
end

% coordinate vectors
if strcmp(Op.Coordinates,'Initial')
    % nothing special, just take initial coordinates
    xvec=unique(Cells.Midx);
    yvec=unique(Cells.Midy);
elseif strcmp(Op.Coordinates,'Updated')
    % take the updated coordinates
    xvec=Points.x{itime+1}; % quadrilateral x coordinates, itime plus 1 to have fully updated coordinates
    yvec=Points.y{itime+1}; % quadrilateral y coordinates
end

% string interpreter
Interpreter = 'tex';
% string for saving figures
PlotTypeStr=PlotType;
% selection of quantity to plot
if strcmp(PlotType,'StrainType')
    fieldtype = (Cells.StrainType{itime});
    fieldmagnitude = (Cells.MagnitudeStrain{itime});
    fieldmagnitudefinal = Cells.MagnitudeStrain{Epochs.Index(end)};
    deformstr='strain type';
elseif strcmp(PlotType,'IncrmtStrainType')
    fieldtype = (Cells.StrainType{itime});
    fieldmagnitude = (Cells.MagnitudeStrain{itime});
    fieldmagnitudefinal = Cells.MagnitudeStrain{Epochs.Index(end)};
    deformstr='incremental strain type';
elseif strcmp(PlotType,'StrainDirection')
    fieldtype = (Cells.DominantPrincplAngle{itime});
    fieldmagnitude = (Cells.MagnitudeStrain{itime});
    fieldmagnitudefinal = Cells.MagnitudeStrain{Epochs.Index(end)};
    deformstr='dominant strain direction';
else
    PlotType
    error('invalid strain option')
end

titlestr=strcat({deformstr},{epochstr});




% field figure
fig=figure;hold on
set(fig,'Position',[0 0 Param.FigureWidth ceil(Param.FigureWidth/aspectratio)])
daspect([1 1 1])

box on




% normalisation factor, using percentile. Everything above a certain
% treshold will be opaque, below will be gradually transparant
if ~isfield(Param,'PercentileTresholdTypeStrainPlot')
    Param.PercentileTresholdTypeStrainPlot=90;
end

% set transparancy

% using a threshold
PercentileTreshold=Param.PercentileTresholdTypeStrainPlot;

normStrain=prctile(abs(fieldmagnitudefinal(~isnan(fieldmagnitudefinal))),PercentileTreshold,'all');

disp(strcat('reducing transparancy for largest strain eigenvalues smaller than:',num2str(normStrain)))
% make transparancy mask

% reduction based on magnitude
AlphaPower=1;
imAlpha=(abs(fieldmagnitude)/normStrain).^AlphaPower;
imAlpha(1,1)=1;
% set everything above the treshold to 1
imAlpha(imAlpha>1)=1;


if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,fieldtype,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape((fieldtype)',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor,'FaceVertexAlphaData',reshape(imAlpha',nx*ny,1),'FaceAlpha','flat')
end

box on


xlim(xlims)
ylim(ylims)


if ~Op.ShowCoordinates
    % remove coordinates
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end

if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end

title(titlestr,'interpreter',Interpreter,'FontSize',Param.TitleFontSize)

if OpColorbar
    if strcmp(PlotType,'StrainType') || strcmp(PlotType,'IncrmtStrainType')
        % colors
        colormap(roma)
        caxis(colorlims)
        cb=colorbar('Ticks',[ -pi/4 , 0, pi/4],...
            'TickLabels',{'Shortening','Strike-Slip','Extension'});
    elseif strcmp(PlotType,'StrainDirection')
        % colors
        colormap(cyclicmap)

        caxis([-90 90]);
        cb=colorbar;
        set(cb,'Ticks',[-90:45:90])
        cb.Label.String = '[deg]';
    end
    set(cb,'FontSize',Param.LabelFontSize)
    %
    
end
% set background color (for NANs)
set(gca,'color',NaNColor);

if OverlayPlot
    
    % change axes limits, to comply with overlay image
    xlim(Fault.xlims)
    ylim(Fault.ylims)
    % make second axes for overlay
    h_ax=gca;
    h_ax_line = axes('position', get(h_ax, 'position')); % Create a new axes in the same position as the first one, overlaid on top
    axis(h_ax_line,'off')
    drawnow
    if DXF
        plotdxf(dxfstruct)
    else
        h=imshow(Fault.Lines);
        set(h, 'AlphaData', Fault.Mask);
    end
    
    box off
    h_ax_line.YTick=[];
    h_ax_line.XTick=[];
    
    
end

% save image
if Op.SaveFigures
    if ~isfolder(Param.SaveDir)
        mkdir(Param.SaveDir)
        disp(strcat('making folder:',Param.SaveDir))
    end
    if OverlayPlot
        SaveFigName=strcat(Param.SaveDir,'/',PlotTypeStr,'_with_faults_',fileepochstr);
    else
        SaveFigName=strcat(Param.SaveDir,'/',PlotTypeStr,'_',fileepochstr);
    end
    savefig(fig,SaveFigName)
    if Op.SavePng
        % save as png
        print(fig,SaveFigName,'-painters','-dpng',Param.FigureResolution)
        
    end
end



end

