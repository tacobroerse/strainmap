function [Stats,StatArea]=AreaStatistics(Cells,Param,Op,Points,PlotType,PlotEpoch,StatArea)
%AreaStatistics Show interactive figure of model and allow for
% doing statistics on displacement and stretch in a specified area
% click 1 to zoom out, click 2 to zoom in. Click on cell location to see
% the temporal evolution of a single cell.
% [Stats,StatArea]=AreaStatistics(Cells,Param,Op,Points,PlotType)
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
% 'StrainMagnitude' largests absolute Hencky strain (log of finite stretch)
% [Stats,StatArea]=AreaStatistics(Cells,Param,Op,Points,PlotType,PlotEpoch)
% PlotEpoch specifies the epoch for plotting. Final is defaults.
% [Stats,StatArea]=AreaStatistics(Cells,Param,Op,Points,PlotType,PlotEpoch,StatArea)
% StatArea, optional, provides the area of the area of interest
% Stats return statistics for the chosen area
%
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H

% do statistics on area
% USAGE: [Stats,StatArea]=AreaStatistics(Cells,Param,Op,Points,PlotType,PlotEpoch,StatArea)

% if epoch is not specified, take last
if nargin < 5
    PlotEpoch = 'final';
end

% initialise
Stats=[];
if ~exist('StatArea','var')
    StatArea=[];
    SelectArea=1;
else
    SelectArea=0;
end

% color for grid
Op.Colorbar=1;
EdgeColor=Param.GridColor;
CenterColorAtOne = 0;
% color for missing data
NaNColor=0.25*[1 1 1];

% size of cell array
[nx,ny]=size(Cells.Midx);


% aspect ratio images
aspectratio=nx/ny;

% blue green orange red
if Op.IncludeTypeExpansionAndContraction
    nColors=Param.nColorsStrainType*2-1;
else
    nColors=Param.nColorsStrainType;
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
elseif strcmp(PlotType,'StrainMagnitude')
    field = (Cells.MagnitudeStrain{itime});
    
    deformstr='strain magnitude (largest absolute Hencky strain)';
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

figmain=figure;hold on
set(figmain,'Position',[0 0 Param.FigureWidth ceil(Param.FigureWidth*aspectratio*1.3) ])


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


% now allow user to select an area

% mid points of cells
xvectransp=xvec';
yvectransp=yvec';
xmid = mean(xvectransp(Cells.Connectivity),2);
ymid = mean(yvectransp(Cells.Connectivity),2);
while 0<1
    figure(figmain)
    % allow for clicking multiple times
    disp('click in figure to select region of interest')
    disp('type 1 to zoom out')
    disp('type 2 to zoom in')
    
    while 0<1
        if ~SelectArea
            % no need to manually select area, as area has been provided as
            % argument
            x=[];y=[];b=[];
        else
            [x,y,b] = ginput(1);
        end
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
            if SelectArea
                % point clicked, proceed
                % draw area
                hROI = drawrectangle;
                % get corresponding cell from coordinates
                % first lower-left corner
                xll=hROI.Position(1);
                yll=hROI.Position(2);
                % upper right corner
                xur=hROI.Position(1)+hROI.Position(3);
                yur=hROI.Position(2)+hROI.Position(4);
                
                % get indices
                [ix1,iy1]=getCell(xll,yll,Cells.Connectivity,xmid,ymid,xvectransp,yvectransp);
                [ix2,iy2]=getCell(xur,yur,Cells.Connectivity,xmid,ymid,xvectransp,yvectransp);
                
                if isempty(ix1) || isempty(ix2)
                    warning('did not find a corresponding cell, try again')
                else
                    % save for output
                    StatArea.ix=[ix1 ix2];
                    StatArea.iy=[iy1 iy2];
                    StatArea.xlim=[xll xur];
                    StatArea.ylim=[yll yur];
                    
                    % MakeFigureRotationInTime(ix,iy,Cells,PlotType,xvec,yvec,field,imAlpha,Op,Cmap,Param)
                end
                % do statistics on this area
                [Stats]=StatisticsROI(StatArea,Cells,Points);
            end
            break
        end
        
    end
    
    prompt = 'choose another area? (y/n)';
    
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

if Op.SaveFigures
    while 0<1
        prompt = 'save figure? (y/n)';
        
        str = input(prompt,'s');
        if strcmp(str,'y')
            % save figure
            SaveDir=strcat(Param.SaveDir,'/statistics/');
            if ~isfolder(SaveDir)
                mkdir(SaveDir)
                disp(strcat('making folder:',SaveDir))
            end
            
            SaveFigName=strcat(SaveDir,'/area_statistics');
          
            savefig(figmain,SaveFigName)
            if Op.SavePng
                % save as png
                print(figmain,SaveFigName,'-painters','-dpng',Param.FigureResolution)
                
            end
            break
        elseif strcmp(str,'n')
            % break out while loop
            break
        else
            % stay in loop and ask again
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

    function [Stats]=StatisticsROI(StatArea,Cells,Points)
        % do statistics on displacements and strain
        Stats=[];
        ntimes=length(Points.u);
        iX1=StatArea.ix(1);
        iX2=StatArea.ix(2);
        iY1=StatArea.iy(1);
        iY2=StatArea.iy(2);
        
        % check for order of indices
        if iY1 > iY2
            iY2old=iY2;
            iY2=iY1;
            iY1=iY2old;
        end
        if iX1 > iX2
            iX2old=iX2;
            iX2=iX1;
            iX1=iX2old;
        end
        
        % displacement in time
        for iTime=1:ntimes
            CumUt(iTime,:,:)=Points.utot{iTime}(iY1:iY2+1,iX1:iX2+1);
            CumVt(iTime,:,:)=Points.vtot{iTime}(iY1:iY2+1,iX1:iX2+1);
            IncrUt(iTime,:,:)=Points.u{iTime}(iY1:iY2+1,iX1:iX2+1);
            IncrVt(iTime,:,:)=Points.v{iTime}(iY1:iY2+1,iX1:iX2+1);
        end
        
        
        fig=figure;
        subplot(1,2,1)
        plot(CumUt(:,1,1),CumVt(:,1,1)); hold on
        plot(IncrUt(:,1,1),IncrVt(:,1,1),'.');
        title('particle path for left lower corner point')
        legend('cumulative displacement','incremental displacement')
        daspect([1 1 1])
        subplot(1,2,2)
        plot(CumUt(:,1,3),CumVt(:,1,3)); hold on
        plot(IncrUt(:,1,3),IncrVt(:,1,3),'.');
        title('particle path for neighbour of left lower corner point')
        legend('cumulative displacement','incremental displacement')
        daspect([1 1 1])
         % save figure
         if Op.SaveFigures
             SaveDir=strcat(Param.SaveDir,'/statistics/');
             if ~isfolder(SaveDir)
                 mkdir(SaveDir)
                 disp(strcat('making folder:',SaveDir))
             end
             
             SaveFigName=strcat(SaveDir,'/particle_path_example');
             
             savefig(fig,SaveFigName)
             if Op.SavePng
                 % save as png
                 print(fig,SaveFigName,'-painters','-dpng',Param.FigureResolution)
                 
             end
         end
        
        % make histogram of displacements
        
        fig=figure;
        subplot(1,2,1);hold on
        title('PDF u displacement (x) \newline in area of interest')
        histogram(IncrUt(:),'Normalization','pdf')
        
        pdu=fitdist(IncrUt(:),'Normal');
        xu=linspace(-3*pdu.sigma,3*pdu.sigma,100);
        yu = pdf(pdu,xu);
        plot(xu,yu,'LineWidth',2)
        axis square
        
        subplot(1,2,2); hold on
        histogram(IncrVt(:),'Normalization','pdf')
        pdv=fitdist(IncrVt(:),'Normal');
        xv=linspace(-3*pdv.sigma,3*pdv.sigma,100);
        yv = pdf(pdv,xv);
        plot(xv,yv,'LineWidth',2)
        title('PDF v displacement (y) \newline in area of interest')
        axis square
         % save figure
         if Op.SaveFigures
             SaveFigName=strcat(SaveDir,'/PDF_displacement');    
             savefig(fig,SaveFigName)
             if Op.SavePng
                 % save as png
                 print(fig,SaveFigName,'-painters','-dpng',Param.FigureResolution)    
             end
         end
        
        
        % save for output
        Stats.usigma=pdu.sigma;
        Stats.umu=pdu.mu;
        Stats.vsigma=pdv.sigma;
        Stats.vmu=pdv.mu;
        
        
        % now strain in time
        for iTime=1:ntimes
            % finite
            V11(iTime,:,:)=Cells.V{iTime}(1,1,iY1:iY2,iX1:iX2);
            V12(iTime,:,:)=Cells.V{iTime}(1,2,iY1:iY2,iX1:iX2);
            V22(iTime,:,:)=Cells.V{iTime}(2,2,iY1:iY2,iX1:iX2);
            % incremental
            eps11(iTime,:,:)=Cells.InfStrainIncrmt{iTime}(1,1,iY1:iY2,iX1:iX2);
            eps12(iTime,:,:)=Cells.InfStrainIncrmt{iTime}(1,2,iY1:iY2,iX1:iX2);
            eps22(iTime,:,:)=Cells.InfStrainIncrmt{iTime}(2,2,iY1:iY2,iX1:iX2);
        end
        
        % make figure of cumulative strain
        fig=figure;
        subplot(2,1,1);
        iplotx=1;iploty=1;
        plot([1:ntimes],V11(:,iplotx,iploty)-1); hold on
        plot([1:ntimes],V12(:,iplotx,iploty));
        plot([1:ntimes],V22(:,iplotx,iploty)-1);
        
        
        plot([1:ntimes],eps11(:,iplotx,iploty),'.'); hold on
        plot([1:ntimes],eps12(:,iplotx,iploty),'.');
        plot([1:ntimes],eps22(:,iplotx,iploty),'.');
        ylabel('strain []')
        xlabel('epochs')
        legend('V11-1','V12','V22-1','\delta\epsilon_{11}','\delta\epsilon_{12}','\delta\epsilon_{22}')
        title('cumulative finite strain and incremental strain \newline for a single point')
        
               subplot(2,1,2);
        iplotx=1;iploty=1;
        plot([1:ntimes],V11(:,iplotx,iploty+2)-1); hold on
        plot([1:ntimes],V12(:,iplotx,iploty+2));
        plot([1:ntimes],V22(:,iplotx,iploty+2)-1);
        
        
        plot([1:ntimes],eps11(:,iplotx,iploty+2),'.'); hold on
        plot([1:ntimes],eps12(:,iplotx,iploty+2),'.');
        plot([1:ntimes],eps22(:,iplotx,iploty+2),'.');
        ylabel('strain []')
        xlabel('epochs')
     %   legend('V11-1','V12','V22-1','\delta\epsilon_{11}','\delta\epsilon_{12}','\delta\epsilon_{22}')
     %   title('cumulative finite strain and incremental strain \newline for a single point')
                 % save figure
         if Op.SaveFigures
             SaveFigName=strcat(SaveDir,'/cumulative_strain_example');    
             savefig(fig,SaveFigName)
             if Op.SavePng
                 % save as png
                 print(fig,SaveFigName,'-painters','-dpng',Param.FigureResolution)    
             end
         end
        
        % make histograms for incremental strain
        fig=figure;
        subplot(1,3,1);hold on
        histogram(eps11(:),'Normalization','pdf')
        
        title('PDF \epsilon_{xx} \newline in area of interest')
        pd{1}=fitdist(eps11(:),'Normal');
        x1=linspace(-3*pd{1}.sigma,3*pd{1}.sigma,100);
        y1 = pdf(pd{1},x1);
        plot(x1,y1,'LineWidth',2)
        axis square
        
        
        subplot(1,3,2);hold on
        histogram(eps12(:),'Normalization','pdf')
        title('PDF \epsilon_{xy} \newline in area of interest')
        pd{2}=fitdist(eps12(:),'Normal');
        x1=linspace(-3*pd{2}.sigma,3*pd{2}.sigma,100);
        y1 = pdf(pd{2},x1);
        plot(x1,y1,'LineWidth',2)
        axis square
        subplot(1,3,3);hold on
        histogram(eps22(:),'Normalization','pdf')
        title('PDF \epsilon_{yy} \newline in area of interest')
        pd{3}=fitdist(eps22(:),'Normal');
        x1=linspace(-3*pd{3}.sigma,3*pd{3}.sigma,100);
        y1 = pdf(pd{3},x1);
        plot(x1,y1,'LineWidth',2)
        axis square
        
                 % save figure
         if Op.SaveFigures
             SaveFigName=strcat(SaveDir,'/PDF_incremental_strain');    
             savefig(fig,SaveFigName)
             if Op.SavePng
                 % save as png
                 print(fig,SaveFigName,'-painters','-dpng',Param.FigureResolution)    
             end
         end
        
        % save for output
        Stats.epsxxsigma=pd{1}.sigma;
        Stats.epsxxmu=pd{1}.mu;
        Stats.epsxysigma=pd{2}.sigma;
        Stats.epsxymu=pd{2}.mu;
        Stats.epsyysigma=pd{3}.sigma;
        Stats.epsyymu=pd{3}.mu;
        
        
        % histograms finite stretch
        fig=figure;
        subplot(1,3,1);hold on
        histogram(V11(:)-1,'Normalization','pdf')
        title('PDF V_{xx} -1 \newline in area of interest')
        pd{1}=fitdist(V11(:)-1,'Normal');
        x1=linspace(-3*pd{1}.sigma,3*pd{1}.sigma,100);
        y1 = pdf(pd{1},x1);
        plot(x1,y1,'LineWidth',2)
        axis square
        subplot(1,3,2);hold on
        histogram(V12(:),'Normalization','pdf')
        title('PDF V_{xy} \newline in area of interest')
        pd{2}=fitdist(V12(:),'Normal');
        x1=linspace(-3*pd{2}.sigma,3*pd{2}.sigma,100);
        y1 = pdf(pd{2},x1);
        plot(x1,y1,'LineWidth',2)
        axis square
        subplot(1,3,3);hold on
        histogram(V22(:)-1,'Normalization','pdf')
        title('PDF V_{yy} -1 \newline in area of interest')
        pd{3}=fitdist(V22(:)-1,'Normal');
        x1=linspace(-3*pd{3}.sigma,3*pd{3}.sigma,100);
        y1 = pdf(pd{3},x1);
        plot(x1,y1,'LineWidth',2)
        axis square
        
                 % save figure
         if Op.SaveFigures
             SaveFigName=strcat(SaveDir,'/PDF_finite_strain');    
             savefig(fig,SaveFigName)
             if Op.SavePng
                 % save as png
                 print(fig,SaveFigName,'-painters','-dpng',Param.FigureResolution)    
             end
         end
        
        % save for output
        Stats.Vxxsigma=pd{1}.sigma;
        Stats.Vxxmu=pd{1}.mu+1;
        Stats.Vxysigma=pd{2}.sigma;
        Stats.Vxymu=pd{2}.mu;
        Stats.Vyysigma=pd{3}.sigma;
        Stats.Vyymu=pd{3}.mu+1;
    end

end