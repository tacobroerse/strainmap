function PlotPrincipalVector(Cells,Points,Param,Op,PlotStrain,PlotEpoch,Epochs)
%PlotPrincipalVector plots principal vectors of strain or stretch tensors.
% PlotPrincipalVector(Cells,Points,Param,Op,PlotStrain,PlotEpoch,Epochs)
% Cells contains the deformation measures and connectivity.
% Points contains the coordinates of vertices, as calculated by FollowPoint.
% Param contains parameters, with defaults set in SetDefaults.
% Param.GridColor specifies the grid color. Param.GridColor='none' hides
% the grid.
% Param.FigureWidth contains the figure width.
% Param.LabelFontSize specifies the font size for labels.
% Param.TitleFontSize specifies the font size for the title.
% Param.VectorPercentileTresholdTypeStrainPlot, only show vectors above a
% certain percentile. Values in range [>0 100]. Default = .1, i.e. .1%. 
% Param.NumVectorsPrincipalStrains, if
% Op.ShowPrincipalVectors='gridnearest', number of vectors in x direction.
% Op, is an options structure. Defaults are set by SetDefaults
% Op.Coordinates, whether to take the original coordinates
% (Op.Coordinates='Initial'), or the deformed coordinates
% (Op.Coordinates='Deformed').
% Op.ColorPrincipalVectors, colouring of principal vectors. 'binary' uses
% two colors, depending on the largest or smallest eigen value. 'colormap'
% colors the principal vectors according to the strain type.
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
% plot principal strains



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
LineWidth=0.5;


[ny,nx]=size(Cells.Midx);
aspectratio = ny/nx;

% color maps
% blue green orange red
nColors=Param.nColorsStrainType;
roma=makescientificcolormap('roma',nColors);


% % check on options
if Op.PlotHenckyStrainInsteadOfStretch && Op.PlotPrincplStretchAsStrain
     warning('it is not possible to plot stretch as strain, and plot stretches logarithmically simultaneously')
     warning('switching off plotting stretch as strain')
     Op.PlotPrincplStretchAsStrain=0;
end

% arrow defaults
LineWidthVec=1;
MaxHeadSize=30;

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
    if HasEpochs && isfield(Epochs,'StartDate')
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
    
    
    % take updated coordinates, and determine mid points
    xmidvec=0.25*(xvec(1:end-1,1:end-1)+xvec(2:end,1:end -1)+xvec(1:end-1,2:end)+xvec(2:end,2:end));
    ymidvec=0.25*(yvec(1:end-1,1:end-1)+yvec(2:end,1:end -1)+yvec(1:end-1,2:end)+yvec(2:end,2:end));
    
end


file_str=[];
% selection of quantity to plot
if strcmp(PlotStrain,'IncrementalInfinitesimalStrain')
    TypeDef='strain';
    fieldstrainmax = Cells.PrincplStrainMaxInfStrainIncr{itime};
    fieldstrainmin = Cells.PrincplStrainMinInfStrainIncr{itime};
    vecmax = Cells.PrincplVecMaxInfStrainIncr{itime};
    vecmin = Cells.PrincplVecMinInfStrainIncr{itime};
   % tensorstr={'\epsilon_{max}','\epsilon_{min}'};
    deformstr='Principal strains incremental infinitesimal strain';

elseif strcmp(PlotStrain,'InfinitesimalStrain')
    TypeDef='strain';
    fieldstrainmax = Cells.PrincplStrainMaxInfStrain{itime};
    fieldstrainmin = Cells.PrincplStrainMinInfStrain{itime};
    vecmax = Cells.PrincplVecMaxInfStrain{itime};
    vecmin = Cells.PrincplVecMinInfStrain{itime};
   % tensorstr={'\epsilon_{max}','\epsilon_{min}'};
    deformstr='Principal strains infinitesimal strain';
elseif strcmp(PlotStrain,'IncrementalGreenFiniteStrain')
    TypeDef='strain';
    fieldstrainmax = Cells.PrincplStrainMaxGreenStrainIncr{itime};
    fieldstrainmin = Cells.PrincplStrainMinGreenStrainIncr{itime};
    vecmax = Cells.PrincplVecMaxGreenStrainIncr{itime};
    vecmin = Cells.PrincplVecMinGreenStrainIncr{itime};
   % tensorstr={'E_{max}','E_{min}'};
    deformstr='Principal strains incremental Green finite strain';
elseif strcmp(PlotStrain,'GreenFiniteStrain')
    TypeDef='strain';
    fieldstrainmax = Cells.PrincplStrainMaxGreenStrain{itime};
    fieldstrainmin = Cells.PrincplStrainMinGreenStrain{itime};
    vecmax = Cells.PrincplVecMaxGreenStrain{itime};
    vecmin = Cells.PrincplVecMinGreenStrain{itime};
   % tensorstr={'E_{max}','E_{min}'};
    deformstr='Principal strains Green finite strain';
elseif strcmp(PlotStrain,'LeftStretchV')
    
    if Op.PlotPrincplStretchAsStrain
        TypeDef='strain';
        fieldstrainmax = Cells.PrincplStretchMaxV{itime}-1;
        fieldstrainmin = Cells.PrincplStretchMinV{itime}-1;
        fieldmax=fieldstrainmax;
        fieldmin=fieldstrainmin;
        vecmax = Cells.PrincplVecMaxV{itime};
        vecmin = Cells.PrincplVecMinV{itime};
     %   tensorstr={'V_{max}','V_{min}'};
        file_str='_strains_';
        deformstr='Principal strains from left stretch tensor V';
    elseif Op.PlotHenckyStrainInsteadOfStretch
        TypeDef='strain';
        fieldmax = (Cells.PrincplStretchMaxV{itime});
        fieldmin = (Cells.PrincplStretchMinV{itime});
        fieldstrainmax=log(fieldmax);
        fieldstrainmin=log(fieldmin);
        vecmax = (Cells.PrincplVecMaxV{itime});
        vecmin = (Cells.PrincplVecMinV{itime});
       % tensorstr={'V_{max}','V_{min}'};
       file_str='_Hencky_strains_';
        deformstr='ln\lambda (V)';
    else
        TypeDef='stretch';
        fieldmax = Cells.PrincplStretchMaxV{itime};
        fieldmin = Cells.PrincplStretchMinV{itime};
        fieldstrainmax=fieldmax-1;
        fieldstrainmin=fieldmin-1;
        vecmax = Cells.PrincplVecMaxV{itime};
        vecmin = Cells.PrincplVecMinV{itime};
       % tensorstr={'V_{max}','V_{min}'};
        deformstr='\lambda (V)';
      %  deformstr='Principal stretches left stretch tensor V';
    end
    
elseif strcmp(PlotStrain,'RightStretchU')
    if Op.PlotPrincplStretchAsStrain
        TypeDef='strain';
        fieldstrainmax = Cells.PrincplStretchMaxU{itime}-1;
        fieldstrainmin = Cells.PrincplStretchMinU{itime}-1;
        fieldmax=fieldstrainmax;
        fieldmin=fieldstrainmin;
        vecmax = Cells.PrincplVecMaxU{itime};
        vecmin = Cells.PrincplVecMinU{itime};
     %   tensorstr={'U_{max}','U_{min}'};
        deformstr='Principal strains from right stretch tensor U';
         file_str='_strains_';
    elseif Op.PlotHenckyStrainInsteadOfStretch
                TypeDef='strain';
        fieldmax = (Cells.PrincplStretchMaxU{itime});
        fieldmin = (Cells.PrincplStretchMinU{itime});
        fieldstrainmax=log(fieldmax);
        fieldstrainmin=log(fieldmin);
        vecmax = (Cells.PrincplVecMaxU{itime});
        vecmin = (Cells.PrincplVecMinU{itime});
     %   tensorstr={'U_{max}','U_{min}'};
        deformstr='ln\lambda (U)';
         file_str='_Hencky_strains_';
    else

        TypeDef='stretch';
        fieldmax = Cells.PrincplStretchMaxU{itime};
        fieldmin = Cells.PrincplStretchMinU{itime};
        fieldstrainmax=fieldmax-1;
        fieldstrainmin=fieldmin-1;
        vecmax = Cells.PrincplVecMaxU{itime};
        vecmin = Cells.PrincplVecMinU{itime};
      %  tensorstr={'U_{max}','U_{min}'};
        deformstr='Principal stretches right stretch tensor U';
    end
    
else
    PlotStrain
    
    error('invalid strain option')
 
end
% string for saving figures
TensorTypeStr=strcat(PlotStrain,file_str);

titlestr=char(strcat({deformstr},epochstr));

% split in x and y components
EigVecMaxx=squeeze(vecmax(1,:,:));
EigVecMaxy=squeeze(vecmax(2,:,:));
EigVecMinx=squeeze(vecmin(1,:,:));
EigVecMiny=squeeze(vecmin(2,:,:));

%titlestr=strcat({deformstr},{epochstr});
Debug=0;
if ~Debug
    FigureWidth=Param.FigureWidth;
    FigureHeight=ceil(Param.FigureWidth*aspectratio);
    % first figure, show principal strains as orthogonal vectors
    fig1=figure; hold on
    set(fig1,'Position',[0 0 FigureWidth FigureHeight])
   % axes=gca;
    if ~strcmp(Param.GridColor,'none')
        % plot grid
        vertices=[xvec(Cells.IndexingVertices) yvec(Cells.IndexingVertices)];
        % plot empty grid
        patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',zeros(nx*ny,1),'FaceColor','none','EdgeColor',EdgeColor,'LineWidth',LineWidth)
    end
    
    if strcmp(Op.ShowPrincipalVectors,'all')
        if numel(Cells.Midx) > 5000
            warning('plotting all principal vectors make take a while due to the large number of cells')
            warning('consider option Op.ShowPrincipalVectors=gridnearest')
        end
        % show vectors for all cells
        xRange=max(Cells.Midx(:))-min(Cells.Midx(:));
        MaxVec=max(fieldstrainmax(:));
        % scaling of vectors, such that vectors do not overlap
        if strcmp(TypeDef,'strain')
            VectorScaling=xRange/nx/MaxVec;
        elseif strcmp(TypeDef,'stretch')
           
            dd=Param.ElementSize;
            VectorScaling=dd/2; % divided by 2 since it starts at the origin
        end
        if VectorScaling == 0
            error('vectors will be scaled with 0, please review definitions')
        end
        % make indices of which maximum eigen values are negative and which are
        % positive
        if strcmp(TypeDef,'strain')
            indexposmax=fieldstrainmax(:)>0;
            indexnegmax=fieldstrainmax(:)<0;
            
            % smallest eigenvalue
            
            indexposmin=fieldstrainmin(:)>0;
            indexnegmin=fieldstrainmin(:)<0;
        elseif strcmp(TypeDef,'stretch')
            indexposmax=fieldmax(:)>0;
            indexnegmax=[];
            
            % smallest eigenvalue
            
            indexposmin=fieldmin(:)>0;
            indexnegmin=[];
        end
    elseif strcmp(Op.ShowPrincipalVectors,'gridnearest')
        % only show a subset of the principal vectors, to prevent too crowded
        % figures
        disp('find cells nearest to an evenly distributed grid')
        nCells=numel(Cells.Midx);
        
        xRange=max(Cells.Midx(:))-min(Cells.Midx(:));
        % space between vectors
        dx = xRange/Param.NumVectorsPrincipalStrains;
        % make homogeneous grid where the vectors should be plotted
        xPosVec = [min(Cells.Midx(:)):dx:max(Cells.Midx(:))];
        yPosVec = [min(Cells.Midy(:)):dx:max(Cells.Midy(:))];
        [xPosGrid,yPosGrid]=meshgrid(xPosVec,yPosVec);
        xPosGrid=xPosGrid(:);
        yPosGrid=yPosGrid(:);
        % find nearest points in cells to this grid
        
        iiPoint=0;
        for iPoint=1:length(xPosGrid)
            distances=sqrt((xmidvec-xPosGrid(iPoint)).^2+(ymidvec-yPosGrid(iPoint)).^2);
            % find nearest
            [MinDist,IndexGridi]=min(reshape(distances,nCells,1));
            if MinDist <= dx/sqrt(2)
                iiPoint = iiPoint + 1;
                IndexGrid(iiPoint)=IndexGridi;
            end
        end
        % largest eigen value
        % make first a vector
        PrincipStrainMaxVecAll=fieldstrainmax(:);
        PrincipStrainMaxVec=NaN(nCells,1);
        PrincipStrainMaxVec(IndexGrid)=PrincipStrainMaxVecAll(IndexGrid);
        % find indices with positive/negative values of largest eigenvalue
        indexposmax=PrincipStrainMaxVec>0;
        indexnegmax=PrincipStrainMaxVec<0;
        
        % smallest eigenvalue
        % make first a vector
        PrincipStrainMinVecAll=fieldstrainmin(:);
        PrincipStrainMinVec=NaN(nCells,1);
        PrincipStrainMinVec(IndexGrid)=PrincipStrainMinVecAll(IndexGrid);
        % find indices with positive/negative values of smallest eigenvalue
        indexposmin=PrincipStrainMinVec>0;
        indexnegmin=PrincipStrainMinVec<0;
        
        
        % scaling of vectors, such that vectors do not overlap
        MaxVec=max(abs([PrincipStrainMaxVec ; PrincipStrainMinVec]));
       
        VectorScaling=1*dx/MaxVec;
    else
        error('unknown option for Op.ShowPrincipalVectors')
    end
    
    if strcmp(TypeDef,'stretch')
        % vector components
        umax = EigVecMaxx.*(fieldmax)*VectorScaling;
        vmax = EigVecMaxy.*(fieldmax)*VectorScaling;
        
        umin = EigVecMinx.*(fieldmin)*VectorScaling;
        vmin = EigVecMiny.*(fieldmin)*VectorScaling;
    elseif strcmp(TypeDef,'strain')
        umax = EigVecMaxx.*(fieldstrainmax)*VectorScaling;
        vmax = EigVecMaxy.*(fieldstrainmax)*VectorScaling;
        
        umin = EigVecMinx.*(fieldstrainmin)*VectorScaling;
        vmin = EigVecMiny.*(fieldstrainmin)*VectorScaling;
    end
    
    % filter out components with too low values, in range [>0 100]
    if ~isfield(Param,'VectorPercentileTresholdTypeStrainPlot')
        Param.VectorPercentileTresholdTypeStrainPlot=.1;
    end
    
    PercentileTreshold=Param.VectorPercentileTresholdTypeStrainPlot;
    % first check that not all values are the same
    minStrain=min(abs(fieldstrainmax(:)));
    maxStrain=max(abs(fieldstrainmax(:)));
    if maxStrain > 5*minStrain
        normStrain=prctile(abs(fieldstrainmax(~isnan(fieldstrainmax))),PercentileTreshold,'all');
    else
        % in case all strains are equal, do not exclude strains
        normStrain=0;
    end
    disp(strcat('omitting vectors for largest strain eigenvalues smaller than:',num2str(normStrain)))
    
    uvmax=sqrt(umax.^2+vmax.^2);
    disp(strcat('percentage of cells excluded:',num2str(numel(find(uvmax<normStrain))/numel(uvmax)*100,3)))
    % set values below treshold to zero
    umax(uvmax<normStrain)=0;
    vmax(uvmax<normStrain)=0;
    umin(uvmax<normStrain)=0;
    vmin(uvmax<normStrain)=0;
    
    if strcmp(Op.ColorPrincipalVectors,'binary')
        colorMax='r';
        colorMin='k';
        
        TypeArrow='arrow';
        %TypeArrow='quiver';
        
        if strcmp(TypeArrow,'quiver')
            % positives (outward vectors)
            quiver(xmidvec(indexposmax),ymidvec(indexposmax),umax(indexposmax),vmax(indexposmax),0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            quiver(xmidvec(indexposmax),ymidvec(indexposmax),-umax(indexposmax),-vmax(indexposmax),0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            % negatives (inward vectors)
            quiver(xmidvec(indexnegmax)+umax(indexnegmax),ymidvec(indexnegmax)+vmax(indexnegmax),-umax(indexnegmax),-vmax(indexnegmax),0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            quiver(xmidvec(indexnegmax)-umax(indexnegmax),ymidvec(indexnegmax)-vmax(indexnegmax),umax(indexnegmax),vmax(indexnegmax),0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            
            % smallest strain component
            
            % positives (outward vectors)
            quiver(xmidvec(indexposmin),ymidvec(indexposmin),umin(indexposmin),vmin(indexposmin),0,'Color',colorMin,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            quiver(xmidvec(indexposmin),ymidvec(indexposmin),-umin(indexposmin),-vmin(indexposmin),0,'Color',colorMin,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            % negatives (inward vectors)
            quiver(xmidvec(indexnegmin)+umin(indexnegmin),ymidvec(indexnegmin)+vmin(indexnegmin),-umin(indexnegmin),-vmin(indexnegmin),0,'Color',colorMin,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            quiver(xmidvec(indexnegmin)-umin(indexnegmin),ymidvec(indexnegmin)-vmin(indexnegmin),umin(indexnegmin),vmin(indexnegmin),0,'Color',colorMin,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            
        else
            daspect([1 1 1])
            headwidth=max(uvmax(:))/2;
            headlength=headwidth;
            % positives (outward vectors)
            if ~isempty(indexposmax)
                arrow3t([xmidvec(indexposmax) ymidvec(indexposmax)],[xmidvec(indexposmax) ymidvec(indexposmax)]+[umax(indexposmax)  vmax(indexposmax)],'k-1',headwidth,headlength)%,0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
                arrow3t([xmidvec(indexposmax) ymidvec(indexposmax)],[xmidvec(indexposmax) ymidvec(indexposmax)]-[umax(indexposmax) vmax(indexposmax)],'k-1',headwidth,headlength)%,0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            end
            % negatives (inward vectors)
            if ~isempty(indexnegmax)
                arrow3t([xmidvec(indexnegmax) ymidvec(indexnegmax)]+[umax(indexnegmax) vmax(indexnegmax)],[xmidvec(indexnegmax) ymidvec(indexnegmax)],'k-1',headwidth,headlength)%,0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
                arrow3t([xmidvec(indexnegmax) ymidvec(indexnegmax)]-[umax(indexnegmax) vmax(indexnegmax)],[xmidvec(indexnegmax) ymidvec(indexnegmax)],'k-1',headwidth,headlength)%,0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            end
            
            % smallest strain component
            % positives (outward vectors)
            if ~isempty(indexposmin)
                arrow3t([xmidvec(indexposmin) ymidvec(indexposmin)],[xmidvec(indexposmin) ymidvec(indexposmin)]+[umin(indexposmin) vmin(indexposmin)],'r-1',headwidth,headlength)%,0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
                arrow3t([xmidvec(indexposmin) ymidvec(indexposmin)],[xmidvec(indexposmin) ymidvec(indexposmin)]-[umin(indexposmin) vmin(indexposmin)],'r-1',headwidth,headlength)%,0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            end
            % negatives (inward vectors)
            if ~isempty(indexnegmin)
                arrow3t([xmidvec(indexnegmin) ymidvec(indexnegmin)]+[umin(indexnegmin) vmin(indexnegmin)],[xmidvec(indexnegmin) ymidvec(indexnegmin)],'r-1',headwidth,headlength)%,0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
                arrow3t([xmidvec(indexnegmin) ymidvec(indexnegmin)]-[umin(indexnegmin) vmin(indexnegmin)],[xmidvec(indexnegmin) ymidvec(indexnegmin)],'r-1',headwidth,headlength)%,0,'Color',colorMax,'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            end
          
        end
        axis normal
        axis equal
        
    elseif strcmp(Op.ColorPrincipalVectors,'colormap')
        % use colormap to map compression - transform - extension on vectors
        
        % can only be done in quiver using loops, one color at a time
        
        ColorArray=roma;
        % determine where in the [-1 1] range the straintype is, use discretize
        % to bin
        % first cut off at 1
        StrainTypeCutOff=Cells.StrainType{itime};
        StrainTypeCutOff((StrainTypeCutOff)>1)=1;
        StrainTypeCutOff((StrainTypeCutOff)<-1)=-1;
        
        % make edges first, this is needed when teh StrainTypeCutOff does not
        % contain the full range [-1 1]
        [~,EdgesBins]=discretize([-1 1],nColors);
        ColorRow=discretize(StrainTypeCutOff,EdgesBins);
        

        
        % now loop
        % positives (outward vectors)
        for iCell=find(indexposmax')
            %   if umax(iCell)~=0 && vmax(iCell) ~=0
            quiver(xmidvec(iCell),ymidvec(iCell),umax(iCell),vmax(iCell),0,'Color',ColorArray(ColorRow(iCell),:),'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            quiver(xmidvec(iCell),ymidvec(iCell),-umax(iCell),-vmax(iCell),0,'Color',ColorArray(ColorRow(iCell),:),'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            %  end
        end
        for iCell=find(indexnegmax')
            
            %  if umax(iCell)~=0 && vmax(iCell) ~=0
            % negatives (inward vectors)
            quiver(xmidvec(iCell)+umax(iCell),ymidvec(iCell)+vmax(iCell),-umax(iCell),-vmax(iCell),0,'Color',ColorArray(ColorRow(iCell),:),'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            quiver(xmidvec(iCell)-umax(iCell),ymidvec(iCell)-vmax(iCell),umax(iCell),vmax(iCell),0,'Color',ColorArray(ColorRow(iCell),:),'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            % end
        end
        % positives (outward vectors)
        for iCell=find(indexposmin')
            % if umax(iCell)~=0 && vmax(iCell) ~=0
            quiver(xmidvec(iCell),ymidvec(iCell),umin(iCell),vmin(iCell),0,'Color',ColorArray(ColorRow(iCell),:),'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            quiver(xmidvec(iCell),ymidvec(iCell),-umin(iCell),-vmin(iCell),0,'Color',ColorArray(ColorRow(iCell),:),'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            % end
        end
        for iCell=find(indexnegmin')
            % if umax(iCell)~=0 && vmax(iCell) ~=0
            % negatives (inward vectors)
            quiver(xmidvec(iCell)+umin(iCell),ymidvec(iCell)+vmin(iCell),-umin(iCell),-vmin(iCell),0,'Color',ColorArray(ColorRow(iCell),:),'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            quiver(xmidvec(iCell)-umin(iCell),ymidvec(iCell)-vmin(iCell),umin(iCell),vmin(iCell),0,'Color',ColorArray(ColorRow(iCell),:),'LineWidth',LineWidthVec,'MaxHeadSize',MaxHeadSize)
            %  end
        end
        
        colormap(ColorArray)
        colorbar
        colorbar('Ticks',[-1+(1/nColors),0,1-(1/nColors)],...
            'TickLabels',{'Shortening','Strike-slip','Extension'})
        caxis([-1 1])
    end
 
     title(titlestr,'FontSize',Param.TitleFontSize)
     
     % axes
     
     
     % set axis limits based on coordinate range from first and last epoch
     box on
     
     %axis tight
     %axis equal
     daspect([1 1 1])
     xlim([min([Points.x{end}(:); Points.x{1}(:)]) max([Points.x{end}(:); Points.x{1}(:)])])
     ylim([min([Points.y{end}(:); Points.y{1}(:)]) max([Points.y{end}(:); Points.y{1}(:)])])
    
     if ~Op.ShowCoordinates
         % remove coordinates
         set(gca,'XTick',[])
         set(gca,'YTick',[])
     end
     if Op.UpwardYAxis==0
         set(gca,'YDir','reverse')
     end
     
     if Op.SaveFigures
         if ~isfolder(Param.SaveDir)
             mkdir(Param.SaveDir)
             disp(strcat('making folder:',Param.SaveDir))
        end
        SaveFigName=strcat(Param.SaveDir,'/',TensorTypeStr,'_Principal_directions_',Op.ColorPrincipalVectors,'_epoch',fileepochstr);
        savefig(fig1,SaveFigName)
        if Op.SavePng
            % save as png
            print(fig1,SaveFigName,'-dpng',Param.FigureResolution)
        end
    end
    
    
end


end


