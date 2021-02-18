function PlotEulerian(PIVresults,Param,Op,PlotStrain,PlotEpoch,Epochs)
%PlotEulerian plots eulerian displacements
% PlotPrincipalValue(Cells,Points,Param,Op,PlotStrain,PlotEpoch,Epochs)
% PIVresults contains the Eulerian displacements and coordinates.
% Param contains parameters, with defaults set in SetDefaults.
% Param.GridColor specifies the grid color. Param.GridColor='none' hides
% the grid.
% Param.FigureWidth contains the figure width.
% Param.LabelFontSize specifies the font size for labels.
% Param.TitleFontSize specifies the font size for the title.
% PlotStrain, string:
% 'IncrementalDisplacement', incremental displacement
% PlotEpoch: epoch number for plot. Numeric or for the last epoch 'final'
% Epochs contains time information.
%
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2021
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

[ny,nx]=size(PIVresults.x{1});
aspectratio = ny/nx;

% color map
[vik]=makescientificcolormap('vik');






% selection of epoch of field
if ~isnumeric(PlotEpoch)
    if strcmp(PlotEpoch,'final')
        fileepochstr='final';
        itime=length(PIVresults.x);
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

xvec=unique(PIVresults.x{itime});
yvec=unique(PIVresults.y{itime});



% string for saving figures
TensorTypeStr=PlotStrain;
% selection of quantity to plot
if strcmp(PlotStrain,'IncrementalDisplacement')
    
    % select displacements and apply mask
    [fieldx,fieldy]=SelectDisplacements(PIVresults,Epochs,Op);
    fieldx=flipud(fieldx{itime});
    fieldy=flipud(fieldy{itime});
    tensorstr={'\delta u','\delta v'};
    deformstr='Incremental displacement';
    
else
    PlotStrain
    
    error('invalid strain option')
    
end
titlestr=char(strcat({deformstr},epochstr));




% individual plot of principal strains
cmax=prctile(abs([fieldx(:) ; fieldy(:)]),Param.ColorPercentile);

clims=[-cmax cmax];

if aspectratio < 1
    m=2;n=1;
    FigureWidth=Param.FigureWidth;
    FigureHeight=ceil(Param.FigureWidth*aspectratio);
else
    m=1;n=2;
    FigureHeight=ceil(Param.FigureWidth*aspectratio/2);
    FigureWidth=Param.FigureWidth;
end

fig2=figure;
set(fig2,'Position',[0 0 FigureWidth FigureHeight])
colormap(vik)
% field x
subplot(m,n,1)
imAlpha=ones(size(fieldx));
imAlpha(isnan(fieldx))=0;
imagesc(xvec,yvec,fieldx,'AlphaData',imAlpha)

box on
axis equal tight
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
else
    set(gca,'YDir','normal')
end

title(tensorstr{1},'FontSize',Param.TitleFontSize)
cb1=colorbar;
set(cb1,'FontSize',Param.LabelFontSize)

caxis(clims)

if ~Op.ShowCoordinates
    % remove coordinates
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end

% set background color (for NANs)
set(gca,'color',NaNColor);

% field y
subplot(m,n,2)
imAlpha=ones(size(fieldy));
imAlpha(isnan(fieldy))=0;
imagesc(xvec,yvec,fieldy,'AlphaData',imAlpha)

box on
axis tight
axis equal
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
else
    set(gca,'YDir','normal')
end
title(tensorstr{2},'FontSize',Param.TitleFontSize)
cb2=colorbar;
set(cb2,'FontSize',Param.LabelFontSize)

caxis(clims)

if ~Op.ShowCoordinates
    % remove coordinates
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end
% set background color (for NANs)
set(gca,'color',NaNColor);

% set general title
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,titlestr,'HorizontalAlignment','Center','FontSize',Param.TitleFontSize)

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
    SaveFigName=strcat(Param.SaveDir,'/',TensorTypeStr,'_Eulerian_epoch',fileepochstr,ColorStr);
    savefig(fig2,SaveFigName)
    if Op.SavePng
        % save as png
        print(fig2,SaveFigName,'-dpng',Param.FigureResolution)
    end
end




end


