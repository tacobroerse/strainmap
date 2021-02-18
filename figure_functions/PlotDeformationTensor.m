function PlotDeformationTensor(Cells,Param,Op,Points,PlotStrain,PlotEpoch)
%PlotDeformationTensor plot deformation tensors, cumulative or incremental
% PlotDeformationTensor(Cells,Param,Op,Points,PlotStrain,PlotEpoch)
% Cells contains the displacements and connectivity.
% Param contains parameters, with defaults set in SetDefaults.
% Param.ColorPercentile specifies the percentile of the deformation values
% that are used for the color map. Values in the range [>0,100]. A value of 
% 100 corresponds to the maximum absolute values.
% Param.GridColor specifies the grid color. Param.GridColor='none' hides
% the grid.
% Param.FigureWidth contains the figure width.
% Param.LabelFontSize specifies the font size for labels.
% Param.TitleFontSize specifies the font size for the title.
% Op, is an options structure. Defaults are set by SetDefaults
% Op.Coordinates, whether to take the original coordinates
% (Op.Coordinates='Initial'), or the deformed coordinates
% (Op.Coordinates='Deformed').
% Points contains the coordinates of vertices, as calculated by FollowPoint.
% PlotStrain, string that specifies the type of tensor to show:
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
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H 

% if epoch is not specified, take last
if nargin < 6
    PlotEpoch = 'final';
end

% color for grid
EdgeColor=Param.GridColor;
LineWidth=0.5;
% color for missing data
NaNColor=0.25*[1 1 1];

% size of cell array
[nx,ny]=size(Cells.Midx);

% aspect ratio images
aspectratio=nx/ny;

% make color map (ColorBrewer)
[bwg]=makecolormap('brownwhitegreen2');
%[rwg]=makecolormap('redwhitegreen');
%cork=makescientificcolormap('cork');

CenterColorAtOne=0;

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
end

% string interpreter
Interpreter = 'tex';
% string for saving figures
TensorTypeStr=PlotStrain;

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

% selection of quantity to plot
if strcmp(PlotStrain,'IncrementalInfinitesimalStrain')
    fieldxx = squeeze(Cells.InfStrainIncrmt{itime}(1,1,:,:));
    fieldxy = squeeze(Cells.InfStrainIncrmt{itime}(1,2,:,:));
    fieldyy = squeeze(Cells.InfStrainIncrmt{itime}(2,2,:,:));
    fieldantixy = Cells.VorticityIncrmt{itime};
    tensorstr={'\epsilon_{xx}','\epsilon_{xy}','\omega','\epsilon_{yy}'};
    deformstr='incremental infinitesimal strain';  
elseif strcmp(PlotStrain,'IncrementalInfinitesimalStrain')
    fieldxx = squeeze(Cells.InfStrainIncrmt{itime}(1,1,:,:));
    fieldxy = squeeze(Cells.InfStrainIncrmt{itime}(1,2,:,:));
    fieldyy = squeeze(Cells.InfStrainIncrmt{itime}(2,2,:,:));
    fieldantixy = Cells.VorticityIncrmt{itime};
    tensorstr={'\epsilon_{xx}','\epsilon_{xy}','\omega','\epsilon_{yy}'};
    deformstr='incremental infinitesimal strain';  
elseif strcmp(PlotStrain,'InfinitesimalStrain')
    fieldxx = squeeze(Cells.InfStrain{itime}(1,1,:,:));
    fieldxy = squeeze(Cells.InfStrain{itime}(1,2,:,:));
    fieldyy = squeeze(Cells.InfStrain{itime}(2,2,:,:));
    fieldantixy = Cells.Vorticity{itime};
    tensorstr={'\epsilon_{xx}','\epsilon_{xy}','\omega','\epsilon_{yy}'};
    deformstr='infinitesimal strain';
elseif strcmp(PlotStrain,'IncrementalGreenFiniteStrain')
    fieldxx = squeeze(Cells.GreenStrainIncrmt{itime}(1,1,:,:));
    fieldxy = squeeze(Cells.GreenStrainIncrmt{itime}(1,2,:,:));
    fieldyy = squeeze(Cells.GreenStrainIncrmt{itime}(2,2,:,:));
    tensorstr={'E_{xx}','E_{xy}','','E_{yy}'};
    fieldantixy = NaN(size(fieldxx));
    deformstr='incremental Green finite strain';
elseif strcmp(PlotStrain,'GreenFiniteStrain')
    fieldxx = squeeze(Cells.GreenStrain{itime}(1,1,:,:));
    fieldxy = squeeze(Cells.GreenStrain{itime}(1,2,:,:));
    fieldyy = squeeze(Cells.GreenStrain{itime}(2,2,:,:));
    fieldantixy = NaN(size(fieldxx));
    tensorstr={'E_{xx}','E_{xy}','','E_{yy}'};
    deformstr='Green finite strain';
elseif strcmp(PlotStrain,'LeftStretchV')
    fieldxx = squeeze(Cells.V{itime}(1,1,:,:));
    fieldxy = squeeze(Cells.V{itime}(1,2,:,:));
    fieldyy = squeeze(Cells.V{itime}(2,2,:,:));
    fieldantixy = Cells.FiniteRotAngle{itime};
    tensorstr={'V_{xx}','V_{xy}','\theta [deg]','V_{yy}'};
    deformstr='Left stretch tensor V';
    CenterColorAtOne=1;
elseif strcmp(PlotStrain,'RightStretchU')
    fieldxx = squeeze(Cells.U{itime}(1,1,:,:));
    fieldxy = squeeze(Cells.U{itime}(1,2,:,:));
    fieldyy = squeeze(Cells.U{itime}(2,2,:,:));
    fieldantixy = Cells.FiniteRotAngle{itime};
    tensorstr={'U_{xx}','U_{xy}','\theta [deg]','U_{yy}'};
    deformstr='Right stretch tensor U';
    CenterColorAtOne=1;
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

else
    PlotStrain
    error('invalid strain option')
end

titlestr=strcat({deformstr},{epochstr});


% color scale
if CenterColorAtOne
    cmaxfield=prctile(abs([fieldxx(:);fieldxy(:);fieldyy(:)]-1),Param.ColorPercentile);
else
    cmaxfield=prctile(abs([fieldxx(:);fieldxy(:);fieldyy(:)]),Param.ColorPercentile);
end

% field figure
m=2;n=2;
fig=figure;
set(fig,'Position',[0 0 Param.FigureWidth ceil(Param.FigureWidth*aspectratio/1.3)])


box on
colormap(bwg)
% field xx
subplot(m,n,1)
imAlpha=ones(size(fieldxx));
imAlpha(isnan(fieldxx))=0;
if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,fieldxx,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    % vertex array (2 columns)
    % itime + 1 has the updated coordinates
    vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape(fieldxx',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor,'LineWidth',LineWidth)
    
end
box on
axis equal
axis tight
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end

title(tensorstr{1},'interpreter',Interpreter)

colorbar

if CenterColorAtOne
    % center color map at 1
    caxis([-cmaxfield cmaxfield]+1)
else
    caxis([-cmaxfield cmaxfield])
end

% set background color (for NANs)
set(gca,'color',NaNColor);


% field xy
subplot(m,n,2)
imAlpha=ones(size(fieldxy));
imAlpha(isnan(fieldxy))=0;
if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,fieldxy,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape(fieldxy',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor,'LineWidth',LineWidth)
    
end
box on
axis equal
axis tight
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end
title(tensorstr{2},'interpreter',Interpreter)

caxis([-cmaxfield cmaxfield])
colorbar

% set background color (for NANs)
set(gca,'color',NaNColor);

if ~isempty(find(~isnan(fieldantixy)))
    
    % vorticity
    subplot(m,n,3)
    imAlpha=ones(size(fieldantixy));
    imAlpha(isnan(fieldantixy))=0;
    if strcmp(Op.Coordinates,'Initial')
        imagesc(xvec,yvec,fieldantixy,'AlphaData',imAlpha)
    elseif strcmp(Op.Coordinates,'Updated')
        patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape(fieldantixy',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor,'LineWidth',LineWidth)
        
    end
    box on
    axis equal
    axis tight
    if Op.UpwardYAxis==0
        set(gca,'YDir','reverse')
    end
    title(tensorstr{3},'interpreter',Interpreter,'FontSize',Param.TitleFontSize)
    colorbar

        cmaxvort=prctile(abs([fieldantixy(:)]),Param.ColorPercentile);

    caxis([-cmaxvort cmaxvort])
    % set background color (for NANs)
    set(gca,'color',NaNColor);
end

% field yy
subplot(m,n,4)
imAlpha=ones(size(fieldyy));
imAlpha(isnan(fieldyy))=0;
if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,fieldyy,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape(fieldyy',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor,'LineWidth',LineWidth)
end
box on
axis equal
axis tight
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end
title(tensorstr{4},'interpreter',Interpreter,'FontSize',Param.TitleFontSize)
if CenterColorAtOne
    % center color map at 1
    caxis([-cmaxfield cmaxfield]+1)
else
    caxis([-cmaxfield cmaxfield])
end
colorbar

% set background color (for NANs)
set(gca,'color',NaNColor);

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,titlestr,'HorizontalAlignment','Center','FontSize',Param.TitleFontSize)

% save image
if Op.SaveFigures
    if ~isfolder(Param.SaveDir)
        mkdir(Param.SaveDir)
        disp(strcat('making folder:',Param.SaveDir))
    end
    if strcmp(RefType,'Eulerian')
        % add eulerian to file name
        SaveFigName=strcat(Param.SaveDir,'/','Eulerian_',TensorTypeStr,'_',fileepochstr);
    else
    SaveFigName=strcat(Param.SaveDir,'/',TensorTypeStr,'_',fileepochstr);
    end
    savefig(fig,SaveFigName)
    if Op.SavePng
        % save as png
        print(fig,SaveFigName,'-dpng',Param.FigureResolution)
    end
end


end

