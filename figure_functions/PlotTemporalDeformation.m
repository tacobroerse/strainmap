function PlotTemporalDeformation(Cells,Param,Op,Points,Epochs,PlotType,GridCells,PlotEpoch)
%  plot temporal deformation
% USAGE: InspectDeformation(Cells,Param,Op,Points,PlotType,PlotIDS,PlotEpoch)

% if epoch is not specified, take last
if nargin < 7
    PlotEpoch = 'final';
end

% relative size of the zoom windows
relwindowsize=0.05;
% color for grid
Op.Colorbar=1;
EdgeColor=Param.GridColor;
CenterColorAtOne = 0;
AxisLineWidth=1;

nColors=Param.nColorsStrainType;
% color for missing data
NaNColor=[1 1 1];

% size of cell array
[ny,nx]=size(Cells.Midx);
alltimes=length(Cells.F);
%timevec=[1:alltimes];
if strcmp(Param.TimeUnit,'sec') && strcmp(Param.TimeUnitPlot,'min')
    timescale=1/60;
elseif strcmp(Param.TimeUnit,Param.TimeUnitPlot)
    % same units
    timescale=1;
else
    error('time units not implemented')
end
timevec=Epochs.Time*timescale;
timestr=strcat('time [',Param.TimeUnitPlot,']');

% aspect ratio images
%aspectratio=ny/nx;
% wrong, since there is rotation


% colormaps
bwg=makecolormap('brownwhitegreen2');
roma=makescientificcolormap('roma',nColors);
cmap=makecolormap('qualitative1');

cmap2=cmap(5:end,:);
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

if ~Op.UseFigureCoordLimits
    % use data coordinate limits
    minx=min(xvec(:));
    maxx=max(xvec(:));
    miny=min(yvec(:));
    maxy=max(yvec(:));
else
    % use specified coordinate limits
    minx=Param.FigureXLim(1);
    maxx=Param.FigureXLim(2);
    miny=Param.FigureYLim(1);
    maxy=Param.FigureYLim(2);
end
xrange=maxx-minx;
yrange=maxy-miny;
aspectratio=yrange/xrange;

% circles and deformed circles for legend
% set circle
radius=pi/10;
angle=linspace(0,2,40)*pi;
% circle coordinates
xc=radius*cos(angle);
yc=radius*sin(angle);

% deformation gradients
Fbext = [1.5 0 ; 0 1.5];
Fext = [1 0 ; 0 1.5];
Fss = [1 1 ; 0 1];
Fsh = [1 0 ; 0 .5];
Fbsh = [0.5 0 ; 0 0.5];

% biaxial extension
temp = Fbext*[xc ; yc];
x_bext = temp(1,:);
y_bext = temp(2,:);

% uniaxial extension
temp = Fext*[xc ; yc];
x_ext = temp(1,:);
y_ext = temp(2,:);

% strike slip 
temp = Fss*[xc ; yc];
x_ss = temp(1,:);
y_ss = temp(2,:);

% uniaxial shortening
temp = Fsh*[xc ; yc];
x_sh = temp(1,:);
y_sh = temp(2,:);

% biaxial shortening
temp = Fbsh*[xc ; yc];
x_bsh = temp(1,:);
y_bsh = temp(2,:);


% string interpreter
Interpreter = 'tex';


% selection of quantity to plot
if strcmp(PlotType,'AllRotation')
    
    field = (Cells.FiniteRotAngle{itime});
    
    deformstr='rotation [deg]';
elseif strcmp(PlotType,'StrainType')
    
    field = (Cells.StrainType{itime});
    fieldmagnitude = (Cells.MagnitudeStrain{itime});
    deformstr='strain type';
elseif strcmp(PlotType,'Dilatation')
    CenterColorAtOne=1;
    field = real(Cells.Dilatation{itime});
    
    deformstr='dilatation';
    
else
    PlotType
    error('invalid strain option')
end

titlestr=strcat({deformstr},{epochstr});

% maximum colormap values
if strcmp(PlotType,'StrainType')
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

% panel sizes
nmaxrowspage=4;% max number of rows per page

dside = 0.05;
dsideleft = dside*3/4;
dxmain = 0.225;
dxsub=.18;

% main figure (only one if number of grid cells to plot is less than
% nmaxrowspage
fig(1)=figure;hold on
set(fig(1),'Position',[0 0 Param.FigureWidth ceil(Param.FigureWidth*sqrt(2))])


% field figure
% make subplot

% plot in top left
% x y dx dy
x = dsideleft;
dx = dxmain;
dy = dx*aspectratio;
y = 1-dside -dy;

pos1=[x y dx dy];
mainsub=subplot('Position',pos1);
mainsub.LineWidth=AxisLineWidth;
box on
if strcmp(PlotType,'StrainType')
    colormap(roma)
else
    colormap(bwg)
end
% field
if strcmp(PlotType,'StrainType')
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

if ~Op.ShowCoordinates
    % remove coordinates
    set(gca,'XTick',[])
    set(gca,'YTick',[])
end

box on
% set aspect ratio (equal)
daspect([1 1 1])
% set coordinate limits
xlim([minx maxx])
ylim([miny maxy])
%axis tight
if Op.UpwardYAxis==0
    set(gca,'YDir','reverse')
end


title(titlestr,'interpreter',Interpreter,'FontSize',Param.TitleFontSize)

% color bar
if Op.Colorbar
    if strcmp(PlotType,'StrainType')
        % strain type has different color bar
        cb=colorbar('Ticks',[ -pi/4 , 0, pi/4],...
            'TickLabels',{'Shortening','Strike-Slip','Extension'});
        set(cb,'Location','SouthOutside')
        set(cb,'TickDirection','out')
        set(cb,'TickLength',0.03)
        set(cb,'FontSize',Param.LabelFontSize)
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

% x y dx dy
x = dside;
dx = dxsub/3;
dy = dxsub;
y = dxsub;
% legend for 3rd panel
posleg=[x y dx dy];
        
axleg=subplot('Position',posleg);
axis off
hold on

set(axleg,'Clipping','off')

daspect([1 1 1])

% plot 5 circles
plot(xc+0.5,yc+pi,':k')
plot(xc+0.5,yc+pi/2,':k')
plot(xc+0.5,yc,':k')
plot(xc+0.5,yc-pi/2,':k')
plot(xc+0.5,yc-pi,':k')

% plot distorted circles
plot(x_bext+0.5,y_bext+pi,'k')
plot(x_ext+0.5,y_ext+pi/2,'k')
plot(x_ss+0.5,y_ss,'k')
plot(x_sh+0.5,y_sh-pi/2,'k')
plot(x_bsh+0.5,y_bsh-pi,'k')

% text (reverse order)
xtext=1.2;

text(0,pi*1.35,'strain type legend','FontSize',Param.LabelFontSize)

text(xtext,-pi,'biaxial shortening','FontSize',Param.LabelFontSize)
text(xtext,-pi/2,'uniaxial shortening','FontSize',Param.LabelFontSize)
text(xtext,0,'strike-slip','FontSize',Param.LabelFontSize)
text(xtext,pi/2,'uniaxial extension','FontSize',Param.LabelFontSize)
text(xtext,pi,'biaxial extension','FontSize',Param.LabelFontSize)

ylim([-pi pi])

%%

nrstr=['a','b','c','d','e','f','g','h','i','j','k'];

nGridCells=length(GridCells.name);

% now panels of individual elements
if Op.UseSameYLimitsMultiPanel
maxstrain=[];
maxdilatation=[];

for isub=1:nGridCells
    ix=GridCells.ix(isub);
    iy=GridCells.iy(isub);
    % determine maximum values
    for itime=1:alltimes
    maxstrain = max([maxstrain abs(log(Cells.PrincplStretchMaxV{itime}(iy,ix))) abs((log(Cells.PrincplStretchMinV{itime}(iy,ix))))]);
    
    % determine maximum dilatation values
  %  maxdilatation = max([maxdilatation max(Cells.Dilatation{itime}(iy,ix))]);
    end
end
maxstrain=ceil(maxstrain);
end
%maxdilatation=maxdilatation*1.1;
maxdilatation=2;
k=0;
for isub=1:nGridCells
    ix=GridCells.ix(isub);
    iy=GridCells.iy(isub);
    titlestr=[nrstr(isub) ': ' GridCells.name{isub}];
    % three columns
    for ipanel = 1:3
        k=k+1;
        % figure number
        fignum=floor((isub-1)/nmaxrowspage)+1;
        % set figure nr
        % check if new figure has to be made
        if (isub-1)/nmaxrowspage == round((isub-1)/nmaxrowspage) && ipanel == 1 && isub ~=1
            % make new figure
            fig(fignum)=figure;
            set(fig(fignum),'Position',[0 0 Param.FigureWidth ceil(Param.FigureWidth*sqrt(2))])   
        end
        
        % position and size of panel
        x = dside*(0.25+ipanel)+dxmain+(ipanel-1)*dxsub;
        dx = dxsub;
        dy = dxsub;
        iisub=mod(isub,nmaxrowspage);  % start anew at new figure
        if iisub==0 
            iisub=nmaxrowspage;
        end
        y = 1-(dside+dxsub)*(iisub);
        
        % position of panel
        pos{k}=[x y dx dy];
        % new subplot
        axsub(k)=subplot('Position',pos{k});
        axsub(k).LineWidth=AxisLineWidth;
        if ipanel == 1
            % plot zoom of strain type
            hold on; box on
            if ~Op.ShowCoordinates
                % remove coordinates
                set(gca,'XTick',[])
                set(gca,'YTick',[])
            end
            % size of cell array
            
            
            vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
            % plot
            patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',reshape((field)',nx*ny,1),'FaceColor','flat','EdgeColor','k','FaceVertexAlphaData',reshape(imAlpha',nx*ny,1),'FaceAlpha','flat')
            box on
            daspect([1 1 1])
            % axis equal
            % [nx,ny]=size(Cells.Midx);
            % now highlight selected cell, but first get
            indexcell=sub2ind([nx ny],ix,iy);
            % get points of current cell
            [VerticesNr]=Cells.Connectivity(indexcell,:);
            % vertex locations
            xvectransp=xvec';
            yvectransp=yvec';
            xVertices = xvectransp(VerticesNr);
            yVertices = yvectransp(VerticesNr);
            % plot cell outlines
            plot([xVertices xVertices(1)],[yVertices yVertices(1)],'-r','LineWidth',1)
            
            
            
            set(gca,'color',NaNColor);
            % coordinate limits, zoom in on element
            xsize=max(xvec(:))-min(xvec(:));
            windowsize=xsize*relwindowsize;
            xlims=[mean(xVertices)+[-1 1]*windowsize/2];
            ylims=[mean(yVertices)+[-1 1]*windowsize/2];
            xlim(xlims);
            ylim(ylims);
            
            
            %   ylim(mean(yVertices)+[-1 1]*windowsize/2*0.8/0.3);
            if Op.UpwardYAxis==0
                set(gca,'YDir','reverse')
            end
            
            
            %     end
            if strcmp(PlotType,'StrainType')
                colormap(axsub(k),roma)
            else
                colormap(axsub(k),bwg)
            end
            
            % max values
            caxis([-cmaxfield cmaxfield])
            
            % title
            ylabel(titlestr,'FontSize',Param.LabelFontSize)
            
            % title on top
            if iisub == 1
                title('strain type zoom','FontSize',Param.TitleFontSize)
            end
            
            % plot area in main plot
            figure(fig(1))
            subplot(mainsub)
            hold on
            % plot zoom outlines
            plot([xlims(1) xlims(2) xlims(2) xlims(1) xlims(1)],[ylims(1) ylims(1) ylims(2) ylims(2) ylims(1)],'-k','LineWidth',1)
            % add number
            text(xlims(2)+0.3*abs(xlims(2)-xlims(1)),ylims(2)+0.3*abs(xlims(2)-xlims(1)),nrstr(isub),'FontSize',Param.LabelFontSize)
            figure(fig(fignum))
            
        elseif ipanel==2
            % plot stretch in time
            hold on; box on
            clear h2
            
            % make vectors of strain with time
            
            vecmagn=zeros(alltimes,1);
            vecmaxprincplstrainV=zeros(alltimes,1);
            vecminprincplstrainV=zeros(alltimes,1);
            
            vec1princplstrainV=zeros(alltimes,1);
            vec2princplstrainV=zeros(alltimes,1);
            for itime=1:alltimes
                vecmagn(itime) = real(Cells.MagnitudeStrain{itime}(iy,ix));
                vecmaxprincplstrainV(itime) = real(Cells.PrincplStretchMaxV{itime}(iy,ix)) ;
                vecminprincplstrainV(itime) = real(Cells.PrincplStretchMinV{itime}(iy,ix)) ;
       
                % instead of min max, take largest and smallest amplitude
                [~,imin] = min([abs(log(Cells.PrincplStretchMaxV{itime}(iy,ix))) abs(log(Cells.PrincplStretchMinV{itime}(iy,ix)))]);
                [~,imax] = max([abs(log(Cells.PrincplStretchMaxV{itime}(iy,ix))) abs(log(Cells.PrincplStretchMinV{itime}(iy,ix)))]);
  
                
                if imin == 1
                   vec1princplstrainV(itime) = Cells.PrincplStretchMaxV{itime}(iy,ix);
                else
                   vec1princplstrainV(itime) = Cells.PrincplStretchMinV{itime}(iy,ix); 
                end
                if imax == 1
                   vec2princplstrainV(itime) = Cells.PrincplStretchMaxV{itime}(iy,ix);
                else
                   vec2princplstrainV(itime) = Cells.PrincplStretchMinV{itime}(iy,ix); 
                end
                
            end
            
            grid on
           % h2(1)=plot(timevec,vecmagn,'-','Color',cmap(1,:),'LineWidth',3);
            % magnitude of incremental strain
            h2(1)=plot(timevec,(log(vecmaxprincplstrainV)),'-','Color',cmap(2,:),'LineWidth',2);
            h2(2)=plot(timevec,(log(vecminprincplstrainV)),'-','Color',cmap(3,:),'LineWidth',2);
            % absolute values
            %   h2(3)=plot(timevec,abs(log(vecmaxprincplstrainV)),'--','Color',cmap(2,:),'LineWidth',2);
            %   h2(5)=plot(timevec,abs(log(vecminprincplstrainV)),'--','Color',cmap(3,:),'LineWidth',2);
            % values as used in strain type definition
      %      h2(3) = plot(timevec,(log(vec1princplstrainV)),'.','Color',cmap(1,:),'LineWidth',3);
       %     h2(4) = plot(timevec,(log(vec2princplstrainV)),'.','Color',cmap(4,:),'LineWidth',3);
     
            
            legendstr{1}='strain magnitude';
            legendstr{1}='ln(\lambda_{max})';
            legendstr{2}='ln(\lambda_{min})';
            %   legendstr{3}='|log(\lambda_{max})|';
            %   legendstr{5}='|log(\lambda_{min})|';
         %   legendstr{3}='ln(\lambda_{a})';
         %   legendstr{4}='ln(\lambda_{b})';
            
            if iisub==1
                legend(h2,legendstr,'Location','NorthWest','FontSize',Param.LabelFontSize)
                legend('boxoff')
                
                title('deformation magnitude','FontSize',Param.TitleFontSize)
            end
            
           
            if ~Op.UseSameYLimitsMultiPanel
                % otherwise, use same maximum strain for all plots
                maxstrain=max([abs(log(vecmaxprincplstrainV)) ; abs(log(vecminprincplstrainV))]);
            end
            if maxstrain > 1
                maxstrain=ceil(maxstrain);
            else
                % round of to .5s
                maxstrain=ceil(maxstrain*2)/2;
            end
          
         axsub(k).YMinorGrid='on';
         %   if maxstrain > 1
                axsub(k).YTick=[-maxstrain 0 maxstrain];
                
         %   else
          %      axsub(k).YTick=[-maxstrain:.25:maxstrain];
         %   end
          ylim([-maxstrain maxstrain])
          
          xlim([0 max(timevec)])
          % labels
          xlabel(timestr,'FontSize',Param.LabelFontSize)
          ylbl=ylabel('Hencky strain []','FontSize',Param.LabelFontSize);
          %
          ylbl.VerticalAlignment='baseline';
         % ylbl.Position=[-33 0 -1];
         set(ylbl, 'Units', 'Normalized', 'Position', [-0.12, 0.5, 0]);
           
        elseif ipanel == 3
            hold on
            grid on
            box on
            % make vectors of strain with time
            vec=zeros(alltimes,1);
            vecincrmt=zeros(alltimes,1);
            vecdilatation=zeros(alltimes,1);
            for itime=1:alltimes
                vecincrmt(itime) = (Cells.StrainTypeIncrmt{itime}(iy,ix));
                vec(itime) = (Cells.StrainType{itime}(iy,ix));
                
                %   % dilatation
                vecdilatation(itime) = Cells.Dilatation{itime}(iy,ix);
            end
            clear legendstr
            legendstr{1}='cumulative strain type';
            legendstr{2}='incremental strain type';
            legendstr{3}='dilatation';
            yyaxis right
            % strain type
            h1(1)=plot(timevec,vec,'-','Color',cmap2(3,:),'LineWidth',2);
            % incremental value
            h1(2)=plot(timevec,vecincrmt,'.-','Color',cmap2(2,:));
            axr=gca;
            axr.YColor = 'k';
            % plot final value
           
            set(axr,'YTick',[-pi/2, -pi/4 , 0, pi/4 , pi/2]);%,...
            set(axr,'YTickLabel',[]);
           
            ylim([-pi/2 pi/2])
            
            % dilatation
            yyaxis left
            axl = gca;
            axl.YColor = 'k';
            %axl.YTick=[0:1:maxdilatation];
            h1(3) = plot(timevec,vecdilatation,'-','Color',cmap2(1,:),'LineWidth',2);
            
            ylbl=ylabel('dilatation []','FontSize',Param.LabelFontSize);
            set(ylbl,'VerticalAlignment','middle')
           
            ylim([0 maxdilatation]);
            xlim([0 max(timevec)])
            xlabel(timestr,'FontSize',Param.LabelFontSize)
            
            yyaxis right
            
            
            if iisub==1
                legend(h1,legendstr,'Location','NorthWest')
                legend('boxoff')
                
                title('dilatation / strain type','FontSize',Param.TitleFontSize)
            end
            
            % add additional panel for explanation of right y-axis
            posextra=[x+dx*1.05 y dx/3 dy];
            
            axextra(isub)=subplot('Position',posextra);
            hold on
            axis off
            set(axextra(isub),'Clipping','off')
            
            daspect([1 1 1])
            
            
            % plot 5 circles
            plot(xc+0.5,yc+pi,':k')
            plot(xc+0.5,yc+pi/2,':k')
            plot(xc+0.5,yc,':k')
            plot(xc+0.5,yc-pi/2,':k')
            plot(xc+0.5,yc-pi,':k')
            
            % plot distorted circles
            plot(x_bext+0.5,y_bext+pi,'k')
            plot(x_ext+0.5,y_ext+pi/2,'k')
            plot(x_ss+0.5,y_ss,'k')
            plot(x_sh+0.5,y_sh-pi/2,'k')
            plot(x_bsh+0.5,y_bsh-pi,'k')
            
            
            
            ylim([-pi pi])
            
        end
        
        
        % reset panel sizes to same size
        OuterPosition=axsub(k).OuterPosition;
        dxsub2=OuterPosition(4);
        axsub(k).OuterPosition=[OuterPosition(1) OuterPosition(2) dxsub2 dxsub2];
        
    end
    
end


for ifig=1:length(fig)
% save image
if Op.SaveFigures
    if ~isfolder(Param.SaveDir)
        mkdir(Param.SaveDir)
        disp(strcat('making folder:',Param.SaveDir))
    end
    
        SaveFigName=strcat(Param.SaveDir,'/',strcat('analogue_model_nr',num2str(ifig),'_temporal_strain_and_type'));
    
    savefig(fig(ifig),SaveFigName)
    if Op.SavePng
        % save as png
        print(fig(ifig),SaveFigName,'-painters','-dpng',Param.FigureResolution)

    end
end
end
end





