function [RegGrid,xGrid,yGrid,IndexField]=InterpolateToRegularGrid(IrregGrid,Param,Op,Points,Epochs,xGrid,yGrid,IndexField)
% SaveAsRegularGrid
% interpolate Cell output on irregular grid to regularly spaced grids and
% save output
% SaveAsRegularGrid(Cells,Param,Op,Points,Epochs,'Quantity')
% Cells contains the connectivity and cell values as stretch, strain,
% strain type, etc.
% Param, is a parameter structure. Defaults are set by SetDefaults
% Op, is an options structure. Defaults are set by SetDefaults
% Points contains the coordinates of vertices, as calculated by FollowPoint.
% Epochs contains timing information
% Quantity is selected quantity from Cells structure
% RegGrids, output on regular grid
% xGrid and yGrid, x and y coordinate arrays
% IndexField provides the indices of the grid where there is a solution in
% the irregular grid, to prevent extrapolation on locations where there is
% no solution.
% xGrid, yGrid and IndexField can be reused, and optionally used as input
% for a next run
%
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2021
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H

%  inspect deformation
% USAGE: [RegGrid,xGrid,yGrid,IndexField]=InterpolateToRegularGrid(IrregGrid,Param,Op,Points,Epochs,xGrid,yGrid,IndexField)


% check if selected quantity exists
nnearbypoints=100;

% % if isfield(Param,'MaxDistGridInterpolation')
% %     scalefactorstrain=Param.MaxDistGridInterpolation;
% % else
% % scalefactorstrain=5;
% % end
%disp(['assume points lie within ',num2str(scalefactorstrain),' times the original grid size from quadrilateral centers'])
% check number of epochs
alltimes=length(IrregGrid);

% set output grid
xmin=[];
xmax=[];
ymin=[];
ymax=[];
for itime = 1:alltimes
    % check first nan values
    indexnonan=~isnan(Points.x{itime+1})&~isinf(Points.x{itime+1});
    % find minimum values
    xmin=min([xmin ; min(Points.x{itime+1}(indexnonan),[],'all')]);
    xmax=max([xmax ; max(Points.x{itime+1}(indexnonan),[],'all')]);
    ymin=min([ymin ; min(Points.y{itime+1}(indexnonan),[],'all')]);
    ymax=max([ymax ; max(Points.y{itime+1}(indexnonan),[],'all')]);
end

if  ~exist('xGrid','var') || ~exist('yGrid','var')
    % make vector using the maximum coordinates
    xq=[xmin:Points.dx:xmax];
    yq=[ymin:Points.dy:ymax];
    
    [xGrid,yGrid]=ndgrid(xq,yq);
end

if ~isa(xGrid,'double')
    xGrid=double(xGrid);
    yGrid=double(yGrid);
end

% check whether indexfield should be made first
if ~exist('IndexField','var')
    NoIndexFieldSupplied=1;
else
    NoIndexFieldSupplied=0;
end
% loop on time
for itime = 1:alltimes
    % get field
    irregfield=IrregGrid{itime};
    
    % get coordinates of cell quantity
    % take the updated coordinates
    xarray=Points.x{itime+1}; % quadrilateral x coordinates, itime plus 1 to have fully updated coordinates
    yarray=Points.y{itime+1}; % quadrilateral y coordinates
    
    % take updated coordinates, and determine mid points
    xmidarray=0.25*(xarray(1:end-1,1:end-1)+xarray(2:end,1:end -1)+xarray(1:end-1,2:end)+xarray(2:end,2:end));
    ymidarray=0.25*(yarray(1:end-1,1:end-1)+yarray(2:end,1:end -1)+yarray(1:end-1,2:end)+yarray(2:end,2:end));
    
    if ~isa(xmidarray,'double')
        xmidarray=double(xmidarray);
        ymidarray=double(ymidarray);
    end
    
    
    % find invalid values
    indexnonan=find(~isnan(irregfield));
    
    % make interpolation function
    F = scatteredInterpolant(xmidarray(indexnonan),ymidarray(indexnonan),irregfield(indexnonan),'linear','none');
    
    % find locations where cells have solution
    if NoIndexFieldSupplied
        IndexField{itime}=zeros(numel(xGrid),1);
        % find which points of the regular grid fall inside the
        % quadrilaterals with a valid solution
        for iPoint=1:numel(xGrid)
            % for each point, check whether it falls within one of the
            % quadrilaterals
            
            if xGrid(iPoint) <= xmax && xGrid(iPoint) >= xmin && yGrid(iPoint) <= ymax && yGrid(iPoint) >= ymin
                
                
                % first check nearby quadrilaterals
                distance=sqrt((xmidarray-xGrid(iPoint)).^2+(ymidarray-yGrid(iPoint)).^2);
                
                % check nearby cells
                
                [mindist,imin]=mink(distance(:),nnearbypoints);
                k=1;
             %   if mindist < Points.dx * scalefactorstrain
                    while IndexField{itime}(iPoint) == 0 && k<=nnearbypoints
                        iQuad=imin(k);
                        % check if grid point falls within quadrilateral
                        [i,j]=ind2sub(size(irregfield),iQuad);
                        x1=Points.x{itime+1}(i,j);
                        x2=Points.x{itime+1}(i+1,j);
                        x4=Points.x{itime+1}(i,j+1);
                        x3=Points.x{itime+1}(i+1,j+1);
                        y1=Points.y{itime+1}(i,j);
                        y2=Points.y{itime+1}(i+1,j);
                        y4=Points.y{itime+1}(i,j+1);
                        y3=Points.y{itime+1}(i+1,j+1);
                        
                        if inpolygon(xGrid(iPoint),yGrid(iPoint),[x1 x2 x3 x4],[y1 y2 y3 y4])
                            % check for nonan
                            if find(iQuad==indexnonan)
                                IndexField{itime}(iPoint)=1;
                                % break out of while loop
                            end
                            
                            
                        end
                        k=k+1;
                    end % end loop on nearby points
                    
             %   end % end if point is close by enough
            end % end for loop on all points
            
        end
    end
    
    % interpolate using the regular coordinates
    
    RegGrid{itime}=NaN(size(xGrid));
    % only use grid points that fall within cells with solutions
    indexesInterp=find(IndexField{itime});
    if ~isempty(indexesInterp)
        RegGrid{itime}(indexesInterp)=F(xGrid(indexesInterp),yGrid(indexesInterp));
    end
    
    
    
    
end
end