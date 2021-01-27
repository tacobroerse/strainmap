function PlotDisplacement(Cells,Param,Op,Points,PlotDeformation,PlotEpoch)
%PlotDisplacement plot displacement, cumulative or incremental
% PlotDisplacement(Cells,Param,Op,Points,PlotDeformation,PlotEpoch)
% Cells contains the displacements and connectivity.
% Param contains parameters, with defaults set in SetDefaults.
% Param.ColorPercentile specifies the percentile of the displacement values
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
% PlotDeformation, string: 'CumulativeDisplacement' or
% 'IncrementalDisplacement'.
% PlotEpoch: epoch number for plot. Numeric or for the last epoch 'final'
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
if nargin < 6
    PlotEpoch = 'final';
end

% color for grid
EdgeColor=Param.GridColor;
LineWidth=0.5;
% color for missing data
NaNColor=0.25*[1 1 1];

% size of cell array
[nx,ny]=size(Points.x{1});

% aspect ratio images
aspectratio=nx/ny;

% color map
[vik]=makescientificcolormap('vik');


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
    xvec=unique(Points.x{1});
    yvec=unique(Points.y{1});
elseif strcmp(Op.Coordinates,'Updated')
    % take the updated coordinates
    xvec=Points.x{itime+1}; % quadrilateral x coordinates, itime plus 1 to have fully updated coordinates
    yvec=Points.y{itime+1}; % quadrilateral y coordinates
end

% string interpreter
Interpreter = 'tex';
% string for saving figures
TensorTypeStr=PlotDeformation;
% selection of quantity to plot
if strcmp(PlotDeformation,'IncrementalDisplacement')
    fieldx = squeeze(Points.u{itime}(:,:));
    fieldy = squeeze(Points.v{itime}(:,:));
   
    
    tensorstr={'\delta u','\delta v'};
    deformstr='incremental displacement';  
elseif strcmp(PlotDeformation,'CumulativeDisplacement')
    fieldx = squeeze(Points.utot{itime}(:,:));
    fieldy = squeeze(Points.vtot{itime}(:,:));
   
    
    tensorstr={'u','v'};
    deformstr='cumulative displacement';  

else
    PlotDeformation
    error('invalid strain option')
end

titlestr=strcat({deformstr},{epochstr});

cmaxfield=prctile(abs([fieldx(:);fieldy(:)]),Param.ColorPercentile);
% end

% field figure
m=1;n=2;
fig=figure;
set(fig,'Position',[0 0 Param.FigureWidth ceil(Param.FigureWidth*aspectratio)])


box on
colormap(vik)
% field xx
subplot(m,n,1)
imAlpha=ones(size(fieldx));
imAlpha(isnan(fieldx))=0;
if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,fieldx,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    % vertex array (2 columns)
    % itime + 1 has the updated coordinates
    vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape(fieldx',nx*ny,1),'FaceColor','interp','EdgeColor',EdgeColor,'LineWidth',LineWidth)
    
end
box on
axis equal
axis tight
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end

title(tensorstr{1},'interpreter',Interpreter,'FontSize',Param.TitleFontSize)

colorbar
caxis([-cmaxfield cmaxfield])

% set background color (for NANs)
set(gca,'color',NaNColor);


% field y
subplot(m,n,2)
imAlpha=ones(size(fieldy));
imAlpha(isnan(fieldy))=0;
if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,fieldxy,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape(fieldy',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor,'LineWidth',LineWidth)
    
end
box on
axis equal
axis tight
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end
title(tensorstr{2},'interpreter',Interpreter,'FontSize',Param.TitleFontSize)
caxis([-cmaxfield cmaxfield])
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
    SaveFigName=strcat(Param.SaveDir,'/',TensorTypeStr,'_',fileepochstr);
    savefig(fig,SaveFigName)
    if Op.SavePng
        % save as png
        print(fig,SaveFigName,'-dpng',Param.FigureResolution)
    end
end



end

