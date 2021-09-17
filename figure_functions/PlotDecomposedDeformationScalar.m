function PlotDecomposedDeformationScalar(Cells,CellsLaplace,Param,Op,Points,PlotType,PlotEpoch,Fault,dxfstruct)
%PlotDeformationScalar plots scalar measures of deformation.
% PlotDeformationScalar(Cells,Param,Op,Points,PlotType,PlotEpoch,Epochs)
% Cells contains the deformation measures and connectivity.
% Param contains parameters, with defaults set in SetDefaults.
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
% PlotType, string: 'Dilatation' (relative area change) or 'MeanRotation'
% (average rotation angle from polar decomposition rotation tensor R).
% PlotEpoch: epoch number for plot. Numeric or for the last epoch 'final'
% Epochs contains time information.
% PlotDeformationScalar(Cells,Param,Op,Points,PlotType,PlotEpoch,Epochs,Fault,dxfstruct)
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

% if epoch is not specified, take last
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

% color for grid
OpColorbar=1;
EdgeColor=Param.GridColor;
CenterColorAtOne = 0;
LineWidth=0.5;

% color for missing data
NaNColor=1*[1 1 1];

% size of cell array
[nx,ny]=size(Cells.Midx);


% aspect ratio images
aspectratio=nx/ny;

% colors
[bwg]=makecolormap('brownwhitegreen');
vikO=makescientificcolormap('vikO');
romaO=makescientificcolormap('romaO');
gwb=flipud(bwg);

invertcolor = 0;
circularcolor = 0;
colorstartat0 = 0;

% selection of epoch of field
if ~isnumeric(PlotEpoch)
    if strcmp(PlotEpoch,'final')
        fileepochstr='final';
         itime=length(Cells.F);
        epochstr=' final epoch';
    else
        disp(strcat('Plot Epoch:',PlotEpoch))
        disp('choose from: final (default), last recorded or epoch')
        disp('epoch needs an additional (4th) argument with the epoch (index)')
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
    
    if strcmp(PlotType,'DecompositionDirection') || strcmp(PlotType,'DecompositionDirectionandMagnitude')
     % take updated coordinates, and determine mid points
    xmidvec=0.25*(xvec(1:end-1,1:end-1)+xvec(2:end,1:end -1)+xvec(1:end-1,2:end)+xvec(2:end,2:end));
    ymidvec=0.25*(yvec(1:end-1,1:end-1)+yvec(2:end,1:end -1)+yvec(1:end-1,2:end)+yvec(2:end,2:end));
    end
end

% string interpreter
Interpreter = 'tex';
% string for saving figures
PlotTypeStr=PlotType;
% selection of quantity to plot
if strcmp(PlotType,'DecompositionDirection')
   
    % rotation angle from rotation tensor R (polar decomposition)
    if strcmp(Op.Coordinates,'Initial')
    fieldvec = squeeze(CellsLaplace.RefAngle{itime})*180/pi;
    elseif strcmp(Op.Coordinates,'Updated')
        % rotate reference orientation by rotation during deformation
        fieldvec = squeeze(CellsLaplace.RefAngle{itime}).*squeeze(CellsLaplace.theta{itime})*180/pi;
    end
    % scalar field is the same as decomposition direction
    field=fieldvec;
circularcolor=1;
colorstartat0=1;
    deformstr='decomposition direction [deg]';  
elseif strcmp(PlotType,'DecompositionDirectionandMagnitude')
    % plot strain magnitude and arrows of decomposition direction
    if strcmp(Op.Coordinates,'Initial')
        fieldvec = squeeze(CellsLaplace.RefAngle{itime})*180/pi;
    elseif strcmp(Op.Coordinates,'Updated')
        % rotate reference orientation by rotation during deformation
        fieldvec = squeeze(CellsLaplace.RefAngle{itime}).*squeeze(CellsLaplace.theta{itime})*180/pi;
    end
    % scalar field is strain magnitude
    field=Cells.MagnitudeStrain{itime};
    circularcolor=0;
    
    deformstr='finite strain magnitude [deg]'; 
%    
%     % rotation angle from rotation tensor R (polar decomposition)
%     field = squeeze(Cells.FiniteRotAngle{itime});
% 
%     deformstr='average rotation [deg]';  
% elseif strcmp(PlotType,'Dilatation')
%     % dilatation
%     CenterColorAtOne=1;
%     field = real(Cells.Dilatation{itime});
% 
%     deformstr='dilatation';
%     invertcolor = 1;
% elseif  strcmp(PlotType,'PrincStretchAngle')
%     field = squeeze(Cells.PrincplAngleDegV{itime});
% 
%     deformstr='left stretch largest principal strain direction [deg]';  
%     circularcolor = 1;
else
    PlotType
    error('invalid strain option')
end

titlestr=strcat({deformstr},{epochstr});


if CenterColorAtOne
    cmaxfield=prctile(abs(field(:)-1),Param.ColorPercentile);
    % dilatation should be centered at one
else
    cmaxfield=prctile(abs(field(:)),Param.ColorPercentile);
end


% field figure

fig=figure;hold on
set(fig,'Position',[0 0 Param.FigureWidth ceil(Param.FigureWidth*aspectratio*1.3) ])


box on
if invertcolor
    colormap(gwb);
elseif circularcolor
    colormap(romaO);
else
    colormap(bwg)
end
% field 

imAlpha=ones(size(field));
imAlpha(isnan(field))=0;
if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,field,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    % vertex array (2 columns)
    % itime + 1 has the updated coordinates
    vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
    patch('vertices',vertices,'faces',Cells.Connectivity,'FaceVertexCData',reshape(field',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor,'LineWidth',LineWidth)
    if strcmp(PlotType,'DecompositionDirection') || strcmp(PlotType,'DecompositionDirectionandMagnitude')
       % plot arrows as well 
       u=cosd(fieldvec);
       v=sind(fieldvec);
       quiver(xmidvec,ymidvec,u,v,'k')
       quiver(xmidvec,ymidvec,-u,-v,'k')
    end
end

box on
     daspect([1 1 1])
     xlim([min([Points.x{end}(:); Points.x{1}(:)]) max([Points.x{end}(:); Points.x{1}(:)])])
     ylim([min([Points.y{end}(:); Points.y{1}(:)]) max([Points.y{end}(:); Points.y{1}(:)])])
    
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end

if ~Op.ShowCoordinates
    % remove coordinates
   set(gca,'XTick',[]) 
   set(gca,'YTick',[])
end

title(titlestr,'interpreter',Interpreter,'FontSize',Param.TitleFontSize)
if OpColorbar
cb=colorbar;
set(cb,'FontSize',Param.LabelFontSize)

if  strcmp(PlotType,'PrincStretchAngle')
    % set range to  [-90 90]
    cmaxfield=90;
end

if CenterColorAtOne 
    caxis([-cmaxfield+1 cmaxfield+1])
elseif colorstartat0
    caxis([0 cmaxfield])
else
    caxis([-cmaxfield cmaxfield])
end
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
        print(fig,SaveFigName,'-dpng',Param.FigureResolution)
    end
end



end

