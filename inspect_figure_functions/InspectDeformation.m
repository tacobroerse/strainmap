function InspectDeformation(Cells,Param,Op,Points,PlotType,PlotEpoch,Epochs)
%InspectDeformation Show interactive figure of model and allow for
% inspecting time evolution of stretch/strain of single cells.
% click 1 to zoom out, click 2 to zoom in. Click on cell location to see
% the temporal evolution of a single cell.
% InspectDeformation(Cells,Param,Op,Points,PlotType,PlotEpoch)
% Cells contains the connectivity and cell values as stretch, strain,
% strain type, etc.
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
% Points contains the coordinates of vertices, as calculated by FollowPoint.
% PlotType, string: 
% 'StrainType' cumulative strain type
% 'IncrmtStrainType' incremental strain type
% 'MeanRotation' mean rotation (from polar decomposition)
% 'IncrmtRotation' incremental rotation (small strain definition)
% 'Dilatation' dilatation, relative area change
%
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H

%  inspect deformation
% USAGE: InspectDeformation(Cells,Param,Op,Points,PlotType,PlotEpoch)

% if epoch is not specified, take last
if nargin < 5
    PlotEpoch = 'final';
end


if exist('Epochs','var')
   
else
    Epochs=[];
end
% color for grid
Op.Colorbar=1;
EdgeColor=Param.GridColor;
CenterColorAtOne = 0;
LineWidth=0.5;
% color for missing data
NaNColor=0.25*[1 1 1];

% size of cell array
[nx,ny]=size(Cells.Midx);


% aspect ratio images
aspectratio=nx/ny;

% blue green orange red
if Op.IncludeTypeExpansionAndContraction
    nColors=Param.nColorsStrainType*2-1;
    colorlims=[-pi/2 pi/2];
else
    nColors=Param.nColorsStrainType;
    colorlims=[-pi/4 pi/4];
end

% colormaps
[bwg]=makecolormap('brownwhitegreen');
roma=makescientificcolormap('roma',nColors);

% selection of epoch of field
if ~isnumeric(PlotEpoch)
    if strcmp(PlotEpoch,'final')
        %fileepochstr='final';
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
  %  fileepochstr=strcat('epoch_',num2str(itime));
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

% selection of quantity to plot
if strcmp(PlotType,'MeanRotation')
    
    field = (Cells.FiniteRotAngle{itime});
    
    deformstr='rotation [deg]';
elseif strcmp(PlotType,'StrainType')
    
    field = (Cells.StrainType{itime});
    fieldmagnitude = (Cells.MagnitudeStrain{itime});
    deformstr='strain type';
elseif strcmp(PlotType,'IncrmtStrainType')
    
    field = (Cells.StrainTypeIncrmt{itime});
    fieldmagnitude = (Cells.MagnitudeStrainIncrmt{itime});
    deformstr='incremental strain type';
elseif strcmp(PlotType,'IncrmntRotation')
    
    field = (Cells.VorticityIncrmt{itime});
    
    deformstr='incremental rotation';
elseif strcmp(PlotType,'Dilatation')
    CenterColorAtOne=1;
    field = real(Cells.Dilatation{itime});
    
    deformstr='dilatation';
elseif strcmp(PlotType,'MagnitudeStrain')
    field=Cells.MagnitudeStrain{itime};
    deformstr='strain magnitude';
else
    PlotType
    error('invalid strain option')
end

titlestr=strcat({deformstr},{epochstr});

% maximum colormap values
if strcmp(PlotType,'StrainType') || strcmp(PlotType,'IncrmtStrainType')
    if Op.IncludeTypeExpansionAndContraction
        cmaxfield=pi/2;
    else
        cmaxfield=pi/4;
    end
else
    if CenterColorAtOne
        cmaxfield=prctile(abs(field(:)-1),Param.ColorPercentile);
        % dilatation should be centered at one
    else
        cmaxfield=prctile(abs(field(:)),Param.ColorPercentile);
    end
end
% end

% field figure

fig=figure;hold on
set(fig,'Position',[0 0 Param.FigureWidth ceil(Param.FigureWidth*aspectratio*1.3) ])


box on
if strcmp(PlotType,'StrainType')
    Cmap=roma;
   
else
    Cmap=bwg;
   
end
colormap(Cmap)
% field
if strcmp(PlotType,'StrainType') ||  strcmp(PlotType,'IncrmtStrainType') 
    % normalisation factor, using percentile. Everything above a certain
    % treshold will be opaque, below will be gradually transparant
    if ~isfield(Param,'PercentileTresholdTypeStrainPlot')
        Param.PercentileTresholdTypeStrainPlot=95;
    end
    
    PercentileTreshold=Param.PercentileTresholdTypeStrainPlot;
    normStrain=prctile(abs(fieldmagnitude(~isnan(fieldmagnitude))),PercentileTreshold,'all');
    disp(strcat('reducing transparancy for largest strain eigenvalues smaller than:',num2str(normStrain,2)))
    % make transparancy mask
    
    AlphaPower=1;
    imAlpha=(abs(fieldmagnitude)/normStrain).^AlphaPower;
    % set everything above the treshold to 1
    imAlpha(imAlpha>1)=1;
else
    % only use nans
    imAlpha=ones(size(field));
    imAlpha(isnan(field))=0;
    
end


if strcmp(Op.Coordinates,'Initial')
    imagesc(xvec,yvec,field,'AlphaData',imAlpha)
elseif strcmp(Op.Coordinates,'Updated')
    % vertex array (2 columns)
    % itime + 1 has the updated coordinates
    vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape((field)',nx*ny,1),'FaceColor','flat','EdgeColor',EdgeColor,'FaceVertexAlphaData',reshape(imAlpha',nx*ny,1),'FaceAlpha','flat')
end



box on
axis equal
axis tight
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end

title(titlestr,'interpreter',Interpreter)

% color bar
if Op.Colorbar
    if strcmp(PlotType,'StrainType') ||  strcmp(PlotType,'IncrmtStrainType') 
        % strain type has different color bar
        colorbar('Ticks',[ -pi/4 , 0, pi/4],...
            'TickLabels',{'Shortening','Strike-Slip','Extension'})
        
        caxis([-cmaxfield cmaxfield])
    else
        
        colorbar
        
        if CenterColorAtOne
            caxis([-cmaxfield+1 cmaxfield+1])
        else
            caxis([-cmaxfield cmaxfield])
        end
    end
end
% set background color (for NANs)
set(gca,'color',NaNColor);

% now allow user to select a grid cell



% mid points of cells
xvectransp=xvec';
yvectransp=yvec';
xmid = mean(xvectransp(Cells.Connectivity),2);
ymid = mean(yvectransp(Cells.Connectivity),2);
while 0<1
    figure(fig)
    % allow for clicking multiple times
    disp('click on cell to inspect')
    disp('type 1 to zoom out')
    disp('type 2 to zoom in')
    
    while 0<1
        [x,y,b] = ginput(1);
        if isempty(b)
            break; % get out of while loop
        elseif b==49 % (key 1)
            % zoom out
            ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
            axis([x-width/2 x+width/2 y-height/2 y+height/2]);
            zoom(1/2);
        elseif b==50 % key 2
            % zoom in
            ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
            axis([x-width/2 x+width/2 y-height/2 y+height/2]);
            zoom(2);
        else
            % area clicked, proceed
            % get corresponding cell from coordinates
            [ix,iy]=getCell(x,y,Cells.Connectivity,xmid,ymid,xvectransp,yvectransp);
            if isempty(ix)
                disp('did not find a corresponding cell, try again')
            else
                
                % make figure
                MakeFigureRotationInTime(ix,iy,Cells,PlotType,xvec,yvec,field,imAlpha,Op,Cmap,Param,Epochs)
            end
            break
        end
    end
    
    prompt = 'choose another cell? (y/n)';
    
    str = input(prompt,'s');
    if strcmp(str,'y')
        % stay in while loop, choose another point
    elseif strcmp(str,'n')
        % leave while loop
        break
    else
        disp('unknown option')
    end
    
    
end
end

function [ix,iy] = getCell(xq,yq,connectivity,midx,midy,verticesx,verticesy)
% find cell that belongs to given coordinates

ix=[];iy=[];
[ny,nx]=size(verticesx);
% reduce by one to go from vertices to cells
ny=ny-1;
nx=nx-1;

% find 4 closest cells next to coordinate
distance=sqrt((midx(:)-xq).^2+(midy(:)-yq).^2);
% sort
[~,indicessorted]=sort(distance);
icell =0;
CellFound=0;
while icell <= 8  && ~CellFound
    icell=icell+1;
    % check whether point falls within cell
    indexcell=indicessorted(icell);
    
    % indices are corresponding to cells, which is different from the
    % indices of its vertices
    
    % get points of current cell
    [VerticesNr]=connectivity(indexcell,:);
    % point locations
    
    xVertices = verticesx(VerticesNr);
    yVertices = verticesy(VerticesNr);
    
    in = inpolygon(xq,yq,xVertices,yVertices);
    if in
        % found cell
        CellFound=1;
 
        % find corresponding indices
        [ix,iy]=ind2sub([ny nx],indexcell);
        % test
        plot([xVertices xVertices(1)],[yVertices yVertices(1)],'-k','LineWidth',2)
        
    end
end
end

function MakeFigureRotationInTime(ix,iy,Cells,PlotType,xvec,yvec,field,imAlpha,Op,Cmap,Param,Epochs)
% make detail figure of rotation/strain type/dilatation in time

alltimes=length(Cells.F);

if strcmp(PlotType,'MeanRotation')
    vec=zeros(alltimes,1);
    for itime=1:alltimes
        vec(itime) = (Cells.FiniteRotAngle{itime}(iy,ix));
    end
    deformstr='rotation [deg]';
    usesubplot=0;
elseif strcmp(PlotType,'StrainType') || strcmp(PlotType,'Dilatation') || strcmp(PlotType,'IncrmtStrainType') || strcmp(PlotType,'MagnitudeStrain')
    vec=zeros(alltimes,1);
    vecincrmt=zeros(alltimes,1);
    vecmagn=zeros(alltimes,1);
    vecdilatation=zeros(alltimes,1);
    vecmaxprincplstrainV=zeros(alltimes,1);
    vecminprincplstrainV=zeros(alltimes,1);
    for itime=1:alltimes
        vecincrmt(itime) = (Cells.StrainTypeIncrmt{itime}(iy,ix));
        vec(itime) = (Cells.StrainType{itime}(iy,ix));
        % magnitude
    
        vecmagn(itime) = real(Cells.MagnitudeStrain{itime}(iy,ix));
        vecmaxprincplstrainV(itime) = real(Cells.PrincplStretchMaxV{itime}(iy,ix)) ;
        vecminprincplstrainV(itime) = real(Cells.PrincplStretchMinV{itime}(iy,ix)) ;
        % dilatation
        vecdilatation(itime) = Cells.Dilatation{itime}(iy,ix);
    end
    usesubplot=1;
    deformstr='strain type';
elseif strcmp(PlotType,'IncrmtRotation')
    vec=zeros(alltimes,1);
    for itime=1:alltimes
        vec(itime) = (Cells.VorticityIncrmt{itime}(iy,ix));
    end
    deformstr='incremental rotation';
    usesubplot=0;
else
    PlotType
    error('invalid strain option')
end

% plot rotations
fig=figure;%('Position',[0.2 0.2 0.6 0.8]);
fig.Position=[500 500 800 600];
hold on


ax1=gca;
cmap = ax1.ColorOrder;
if usesubplot
    pos1=[0.15 0.55 0.35 0.35];
    ax1=subplot('Position',pos1);
end
hold on
box on
grid on

if isempty(Epochs)
    timevec=[1:alltimes];
else
    timevec=Epochs.Time(1:alltimes);
    
end
h1(1)=plot(timevec,vec,'.-','Color',cmap(1,:));
if strcmp(PlotType,'StrainType') || strcmp(PlotType,'Dilatation') || strcmp(PlotType,'IncrmtStrainType')
    legendstr1{1}='mean strain type';
    legendstr1{2}='incremental strain type';
    legendstr1{3}='dilatation';
    
    % incremental value
    h1(2)=plot(timevec,vecincrmt,'.-','Color',cmap(2,:));
    axl=gca;
    % plot final value
    plot(timevec,ones(size(timevec))*Cells.StrainType{end}(iy,ix),'k:')
    set(axl,'YTick',[-pi/2, -pi/4 , 0, pi/4 , pi/2],...
        'YTickLabel',{'Biaxial \newline Shortening','Uniaxial \newline Shortening','Strike-Slip','Uniaxial \newline Extension','Biaxial \newline Extension'})
    ylim([-pi/2 pi/2])
    
    % dilatation
    yyaxis right
    axr = gca;
    axr.YColor = 'k';
    h1(3) = plot(timevec,vecdilatation,'-','Color',cmap(3,:),'LineWidth',1);
    ylabel('dilatation []')
    maxy=max(abs(vecdilatation-1));
    if maxy == 0
        maxy=0.1;
    end
    ylim([-maxy maxy]+1)
    
    yyaxis left
    xlim([timevec(1) timevec(end)])
    % magnitude of mean strain
    pos2=[0.15 0.1 0.35 0.35];
    ax2=subplot('Position',pos2);
    hold on; box on
    grid on
    h2(1)=plot(timevec,vecmagn,'-','Color',cmap(1,:),'LineWidth',3);
    % magnitude of incremental strain
    h2(2)=plot(timevec,(log(vecmaxprincplstrainV)),'-','Color',cmap(2,:),'LineWidth',2);
    h2(4)=plot(timevec,(log(vecminprincplstrainV)),'-','Color',cmap(3,:),'LineWidth',2);
    % absolute values
    h2(3)=plot(timevec,abs(log(vecmaxprincplstrainV)),'--','Color',cmap(2,:),'LineWidth',1);
    h2(5)=plot(timevec,abs(log(vecminprincplstrainV)),'--','Color',cmap(3,:),'LineWidth',1);
    
    legendstr2{1}='mean strain magnitude';
    legendstr2{2}='log principal stretch 1';
    legendstr2{4}='log principal stretch 2';
    legendstr2{3}='log absolute principal stretch 1';
    legendstr2{5}='log absolute principal stretch 2';
    
    ylabel('strain magnitude []')
    
    xlim([timevec(1) timevec(end)])
    xlabel('epoch')
    
    % plot field
    pos3=[0.65 0.1 0.30 0.8];
    ax3=subplot('Position',pos3);
    hold on; box on
    % size of cell array
    [nx,ny]=size(Cells.Midx);
    
    vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
    % plot
    patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape((field)',nx*ny,1),'FaceColor','flat','EdgeColor','k','FaceVertexAlphaData',reshape(imAlpha',nx*ny,1),'FaceAlpha','flat')
    box on
    axis equal
    % now highlight selected cell, but first get
    indexcell=sub2ind([ny nx],ix,iy);
    % get points of current cell
    [VerticesNr]=Cells.Connectivity(indexcell,:);
    % vertex locations
    xvectransp=xvec';
    yvectransp=yvec';
    xVertices = xvectransp(VerticesNr);
    yVertices = yvectransp(VerticesNr);
    % plot cell
    plot([xVertices xVertices(1)],[yVertices yVertices(1)],'-k','LineWidth',2)
    % axis tight
    % set background color (for NANs)
    NaNColor=0.25*[1 1 1];
    set(gca,'color',NaNColor);
    % coordinate limits
    indexnonan=~isnan(xvec)&~isinf(xvec);
    xsize=max(xvec(indexnonan),[],'all')-min(xvec(indexnonan),[],'all');
    windowsize=xsize/10;
    xlim(mean(xVertices)+[-1 1]*windowsize/2);
    ylim(mean(yVertices)+[-1 1]*windowsize/2*0.8/0.3);
   daspect([1 1 1])
    if Op.UpwardYAxis==0
        set(gca,'YDir','reverse')
    end
    %     end
    colormap(Cmap)
    if strcmp(PlotType,'StrainType') || strcmp(PlotType,'IncrmtStrainType')
    % strain type has different color bar
    cb=colorbar('Ticks',[ -pi/4 , 0, pi/4],...
        'TickLabels',{'Shortening','Strike-Slip','Extension'});
    
   %  cb.Location='southoutside';
    if Op.IncludeTypeExpansionAndContraction
        cmaxfield=pi/2;
    else
        cmaxfield=pi/4;
    end
        caxis([-cmaxfield cmaxfield])
    else 
        % dilatation, center at one
        cmaxfield=max(abs(field(:)-1));
        caxis([-cmaxfield+1 cmaxfield+1])
        cb=colorbar;
    end
    
    
    %  plot(timevec,abs(log(vecmaxprincplstrainV))./abs(log(vecminprincplstrainV)))
    
    subplot(ax1)
elseif strcmp(PlotType,'MeanRotation')
    ylabel('mean rotation [deg]')
elseif strcmp(PlotType,'IncrmtRotation')
    ylabel('incremental rotation [deg]')

end

title(strcat('point [ix:iy] [',num2str(ix),':',num2str(iy),']'))
xlabel('epoch')

%xlim([1 alltimes])

if usesubplot
    subplot(ax1)
    LegendLocation='Best';
    legend(h1,legendstr1,'Location',LegendLocation,'Box','off')
    
    subplot(ax2)
    LegendLocation='Best';
    legend(h2,legendstr2,'Location',LegendLocation,'Box','off')
    
end

prompt = 'save figure? (y/n)';

str = input(prompt,'s');
if strcmp(str,'y')
    % save figure
    SaveDir=strcat(Param.SaveDir,'/Timeseries/');
    if ~isfolder(SaveDir)
        mkdir(SaveDir)
        disp(strcat('making folder:',SaveDir))
    end
    SaveFigName=strcat(SaveDir,'/','timeseries_ix_',num2str(ix),'_iy_',num2str(iy));
    savefig(fig,SaveFigName)
    if Op.SavePng
        % save as png
        print(fig,SaveFigName,'-dpng',Param.FigureResolution)
    end
end








end



