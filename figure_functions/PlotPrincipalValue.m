function PlotPrincipalValue(Cells,Points,Param,Op,PlotStrain,PlotEpoch,Epochs)
%PlotPrincipalValue plots principal values of strain or stretch tensors.
% PlotPrincipalValue(Cells,Points,Param,Op,PlotStrain,PlotEpoch,Epochs)
% Cells contains the deformation measures and connectivity.
% Points contains the coordinates of vertices, as calculated by FollowPoint.
% Param contains parameters, with defaults set in SetDefaults.
% Param.GridColor specifies the grid color. Param.GridColor='none' hides
% the grid.
% Param.FigureWidth contains the figure width.
% Param.LabelFontSize specifies the font size for labels.
% Param.TitleFontSize specifies the font size for the title.
% Op, is an options structure. Defaults are set by SetDefaults
% Op.Coordinates, whether to take the original coordinates
% (Op.Coordinates='Initial'), or the deformed coordinates
% (Op.Coordinates='Deformed').
% Op.PlotLogarithmicColorMapPrincplStretch use logarithmic color map
% Op.PlotPrincplStretchAsStrain subtract 1 from stretch to get principal
% strain
% Op.ShowCoordinates, whether or not to show coordinates
% PlotStrain, string:
% 'IncrementalInfinitesimalStrain', incremental infinitesimal strain
% 'InfinitesimalStrain', infinitesimal strain
% 'IncrementalGreenFiniteStrain', incremental finite Green-Lagrangian
% strain
% 'GreenFiniteStrain', finite Green-Lagrangian strain
% 'LeftStretchV', finite left stretch
% 'RightStretchU', finite right stretch
% PlotEpoch: epoch number for plot. Numeric or for the last epoch 'final'
% Epochs contains time information.
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
if nargin == 7
    HasEpochs = 1;
else
    HasEpochs = 0;
end

% color for grid
EdgeColor=Param.GridColor;
%   plot strain tensor
NaNColor=0.25*[1 1 1];

[ny,nx]=size(Cells.Midx);
aspectratio = ny/nx;

% color maps
% brown white green
%cork=makescientificcolormap('cork');
bwg=makecolormap('brownwhitegreen2');

PlotLogarithmicColorMapPrincplStretch =0;

% check on options
if Op.PlotLogarithmicColorMapPrincplStretch && Op.PlotPrincplStretchAsStrain
    warning('it is not possible to plot stretch as strain, and plot stretches logarithmically simultaneously')
    warning('switching off plotting stretch as strain')
    Op.PlotPrincplStretchAsStrain=0;
end


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
    if HasEpochs  && isfield(Epochs,'StartDate')
        fileepochstr=strcat(Epochs.StartDate{1},'_',Epochs.EndDate{itime});
        epochstr=strcat({' period: '},Epochs.StartDate{1},{' '},Epochs.EndDate{itime});
    else
        fileepochstr=strcat('epoch_',num2str(itime));
        epochstr=strcat({' epoch '},{num2str(itime)});
    end
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

% string for saving figures
TensorTypeStr=PlotStrain;
% selection of quantity to plot
if strcmp(PlotStrain,'IncrementalInfinitesimalStrain')
    
    fieldmax = Cells.PrincplStrainMaxInfStrainIncr{itime};
    fieldmin = Cells.PrincplStrainMinInfStrainIncr{itime};
    
    tensorstr={'\epsilon_{max}','\epsilon_{min}'};
    deformstr='Principal strains incremental infinitesimal strain';
    
elseif strcmp(PlotStrain,'InfinitesimalStrain')
    fieldmax = Cells.PrincplStrainMaxInfStrain{itime};
    fieldmin = Cells.PrincplStrainMinInfStrain{itime};
    tensorstr={'\epsilon_{max}','\epsilon_{min}'};
    deformstr='Principal strains infinitesimal strain';
elseif strcmp(PlotStrain,'IncrementalGreenFiniteStrain')
    fieldmax = Cells.PrincplStrainMaxGreenStrainIncr{itime};
    fieldmin = Cells.PrincplStrainMinGreenStrainIncr{itime};
    tensorstr={'E_{max}','E_{min}'};
    deformstr='Principal strains incremental Green finite strain';
elseif strcmp(PlotStrain,'GreenFiniteStrain')
    fieldmax = Cells.PrincplStrainMaxGreenStrain{itime};
    fieldmin = Cells.PrincplStrainMinGreenStrain{itime};
    tensorstr={'E_{max}','E_{min}'};
    deformstr='Principal strains Green finite strain';
elseif strcmp(PlotStrain,'LeftStretchV')
    
    if Op.PlotPrincplStretchAsStrain
        fieldstrainmax = Cells.PrincplStretchMaxV{itime}-1;
        fieldstrainmin = Cells.PrincplStretchMinV{itime}-1;
        fieldmax=fieldstrainmax;
        fieldmin=fieldstrainmin;
        tensorstr={'\lambda_{max}','\lambda_{min}'};
        deformstr='Principal strains from left stretch tensor V';
    else
        fieldmax = Cells.PrincplStretchMaxV{itime};
        fieldmin = Cells.PrincplStretchMinV{itime};
        tensorstr={'\lambda_{max}','\lambda_{min}'};
        deformstr='Principal stretch \lambda (V)';
        if Op.PlotLogarithmicColorMapPrincplStretch
            PlotLogarithmicColorMapPrincplStretch  = 1;
        end
    end
    
elseif strcmp(PlotStrain,'RightStretchU')
    if Op.PlotPrincplStretchAsStrain
        fieldstrainmax = Cells.PrincplStretchMaxU{itime}-1;
        fieldstrainmin = Cells.PrincplStretchMinU{itime}-1;
        fieldmax=fieldstrainmax;
        fieldmin=fieldstrainmin;
        tensorstr={'\lambda_{max}','\lambda_{min}'};
        deformstr='Principal strains from right stretch tensor U';
    else
        fieldmax = Cells.PrincplStretchMaxU{itime};
        fieldmin = Cells.PrincplStretchMinU{itime};
        tensorstr={'U_{max}','U_{min}'};
        deformstr='Principal stretches right stretch tensor U';
        if Op.PlotLogarithmicColorMapPrincplStretch
            PlotLogarithmicColorMapPrincplStretch  = 1;
        end
    end
    
else
    PlotStrain
    
    error('invalid strain option')
    
end
titlestr=char(strcat({deformstr},epochstr));




% individual plot of principal strains

cmax=prctile(abs(fieldmax(:)),Param.ColorPercentile);
if Op.PlotPrincplStretchAsStrain
    clims=[-cmax cmax];
else
    if PlotLogarithmicColorMapPrincplStretch
        clims=[0 cmax];
    else
        clims=[-cmax cmax];
    end
end
if aspectratio < 1
    m=2;n=1;
    FigureWidth=Param.FigureWidth;
    FigureHeight=ceil(Param.FigureWidth*aspectratio);
else
    m=1;n=2;
    FigureHeight=ceil(Param.FigureWidth*aspectratio);
    FigureWidth=Param.FigureWidth;
end

fig2=figure;
set(fig2,'Position',[0 0 FigureWidth FigureHeight])
colormap(bwg)
% strain 11
subplot(m,n,1)
imAlpha=ones(size(fieldmax));
imAlpha(isnan(fieldmax))=0;
if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,fieldmax,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape(fieldmax',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor)
end
box on
axis equal tight
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end
title(tensorstr{1},'FontSize',Param.TitleFontSize)
cb1=colorbar;
set(cb1,'FontSize',Param.LabelFontSize)
if PlotLogarithmicColorMapPrincplStretch
    % plot logarithmic colormapping
    set(gca,'ColorScale','log')
    logcmax=log10(cmax);
    caxis([10^(-logcmax) 10^(logcmax)])
    roundcmax=ceil(logcmax);
    nticks=5;
    tickintervals = round(logspace(-roundcmax,roundcmax,2*nticks*roundcmax+1),1,'significant');
    cb1.Ticks = tickintervals;
else
    caxis(clims)
end
if ~Op.ShowCoordinates
    % remove coordinates
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end

% set background color (for NANs)
set(gca,'color',NaNColor);

% strain 22
subplot(m,n,2)
imAlpha=ones(size(fieldmin));
imAlpha(isnan(fieldmin))=0;
if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,fieldmin,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape(fieldmin',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor)
end
box on
axis tight
axis equal
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end
title(tensorstr{2},'FontSize',Param.TitleFontSize)
cb2=colorbar;
set(cb2,'FontSize',Param.LabelFontSize)
if PlotLogarithmicColorMapPrincplStretch
    % plot logarithmic colormapping
    set(gca,'ColorScale','log')
    logcmax=log10(cmax);
    caxis([10^(-logcmax) 10^(logcmax)])
    cb2.Ticks = tickintervals;
else
    caxis(clims)
end

if ~Op.ShowCoordinates
    % remove coordinates
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
% set background color (for NANs)
set(gca,'color',NaNColor);

% set general title
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,titlestr,'HorizontalAlignment','Center')

% save image
if Op.SaveFigures
    if ~isfolder(Param.SaveDir)
        mkdir(Param.SaveDir)
        disp(strcat('making folder:',Param.SaveDir))
    end
    if Op.PlotLogarithmicColorMapPrincplStretch
        ColorStr='_log_colormap';
    else
        ColorStr=[];
    end
    SaveFigName=strcat(Param.SaveDir,'/',TensorTypeStr,'_Principal_values_epoch',fileepochstr,ColorStr);
    savefig(fig2,SaveFigName)
    if Op.SavePng
        % save as png
        print(fig2,SaveFigName,'-dpng',Param.FigureResolution)
    end
end




end


