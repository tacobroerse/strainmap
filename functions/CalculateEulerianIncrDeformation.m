function [CellsEulerian]=CalculateEulerianIncrDeformation(PIVresults,Epochs,Op)
%CalculateEulerianIncrDeformation calculates incremental strain from pixel displacements
% pixel displacements are in an Eulerian frame and thus denote space
% referenced displacement. For Lagrangian strains, use CalculateDeformation
% instead.
%
% [CellsEulerian]=CalculateEulerianIncrDeformation(PIVresults,Epochs,Op) 
% PIVresults contains the
% coordinates and incremental spatial displacements. 
% Epochs contains time information.
% Op contains options. For this function it should contain:
% Op.DisplacementGradient: 'Simple','MidPoint', or 'ShapeFunctions', where
% only ShapeFunctions is recommended. ShapeFunctions uses bilinear shape
% functions for quadrilaterals. From the displacement gradient the
% deformation gradient F is computed, from which all other deformation
% measures are derived. 
%
% In every step the infinitesimal strain or Green-Lagrangian strain is
% computed for a quadrilateral in between 4 points, which is a square for the regular grid. Mind that
% infinitesimal strain is based on the assumption of small rotations, and will not be valid
% for large rotations or shears > 0.1 .
% The incremental deformation gradient tensor can be computed using the displacement at 3
% points (the standard definition) or using 4 points (better) using linear
% shape functions (similar as FEM computations). These linear shape functions provide a bi-linear
% interpolation of the displacements (u and v, in x resp y direction). When shape functions are used the
% strain is evaluated in the middle of the cell.
%
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2021
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H
%

% time information
nEpochs = length(PIVresults.x);
alltimes=Epochs.Index;
starttime = alltimes(1);

% size of grid
[ny,nx]=size(PIVresults.x{starttime});

%% make grid coordinates and connectivity

% x y is meshgrid
% x is ordered
%
% x1 x2 .. xn
% x1 x2 .. xn
% x1 x2 .. xn
% x1 x2 .. xn

% y is ordered
%
% y1 y1 y1 y1
% y2 y2 y2 y2
% ....
% yn yn yn yn

% determine grid spacing
PIVresults.dx=abs(diff(PIVresults.x{1}(1,1:2)));
PIVresults.dy=abs(diff(PIVresults.y{1}(1:2,1)));

% select the right type of displacement (original, filtered, smoothed)
 [u,v]=SelectDisplacements(PIVresults,Epochs,Op);

switch Op.DisplacementGradient
    
    case 'Simple'
        disp('Op.DisplacementGradient == Simple is not recommended')
        % middle points of deforming cells (at the first epoch)
        CellsEulerian.Midx = PIVresults.x{starttime}(1:ny-1,1:nx-1)+0.5*PIVresults.dx;
        CellsEulerian.Midy = PIVresults.y{starttime}(1:ny-1,1:nx-1)+0.5*PIVresults.dy;
    case 'MidPoint'
        disp('Op.DisplacementGradient == MidPoint is not recommended')
        % middle points of deforming cells are identical to point locations
        % (at the first epoch), without the points at the boundaries
        CellsEulerian.Midx = PIVresults.x{starttime}(2:ny-1,2:nx-1);
        CellsEulerian.Midy = PIVresults.y{starttime}(2:ny-1,2:nx-1);
    case 'ShapeFunctions'
        
        % middle points of deforming cells (at the first epoch)
        CellsEulerian.Midx = PIVresults.x{starttime}(1:ny-1,1:nx-1)+0.5*PIVresults.dx;
        CellsEulerian.Midy = PIVresults.y{starttime}(1:ny-1,1:nx-1)+0.5*PIVresults.dy;
        
    otherwise
        error(strcat('invalid option for Op.DisplacementGradient:',Op.DisplacementGradient))
end



% make polygons of coordinates, for plotting

% compute connectivity between xvec and yvec (which are the vertices of
% the quadrilaterals)

switch Op.DisplacementGradient
    case {'Simple','ShapeFunctions'}
        
        Connectivity=zeros(nx-1,ny-1,4);
        
        for i=1:nx-1
            for j=1:ny-1
                %   k=k+1;
                
                % ordering of vertices
                % 4--3
                % |  |
                % 1--2
                %
                Connectivity(i,j,1:4) = [ nx*(j)+i nx*(j)+i+1 nx*(j-1)+i+1 nx*(j-1)+i];
            end
        end
        
        % reshape from ny x nx x 4 to (nx*ny) x 4 array
        CellsEulerian.Connectivity=reshape(Connectivity,(nx-1)*(ny-1),4);
        
        % coordinates differ for each epoch, but indices do not
        % this vector can be used to convert an array into a vector
        % Vec = Array(Index)
        CellsEulerian.IndexingVertices=reshape(reshape([1:nx*ny],ny,nx)',nx*ny,1);
        
    case 'MidPoint'
        
        Connectivity=zeros(nx-2,ny-2,4);
        
        for i=1:nx-2 % loop over x
            for j=1:ny-2 % loop over y
                %   k=k+1;
                
                % ordering of vertices
                %    1
                %   / \
                %  4-o-2
                %   \ /
                %    3
                % midpoint is not part of element, but only the center
                %
                Connectivity(i,j,1:4) = [nx*(j-1)+i+1 nx*(j)+i+2 nx*(j+1)+i+1 nx*(j)+i];
            end
        end
        
        % reshape
        CellsEulerian.Connectivity=reshape(Connectivity,(ny-2)*(nx-2),4);
        
        % coordinates differ for each epoch, but indices do not
        % this vector can be used to convert an array into a vector
        % Vec = Array(Index)
        CellsEulerian.IndexingVertices=reshape(reshape([1:nx*ny],ny,nx)',nx*ny,1);
end

%% determine displacement gradient tensor F

switch Op.DisplacementGradient
    
    case 'Simple'
        % simple: uses not all displacements of neighbouring nodes
        %    4-----3
        %    |     |
        %    |     |
        %    1-----2   but uses only the difference in displacements between
        %    1-2 and 1-4. i.e., no information of the displacement in 3 is
        %    incorporated
        %
        %
        %
        % in case of large deformation this approach may not be optimal
        
        % loop over time
        for itime=alltimes
            
            % width of cells, using current widths of cells
            % current points are the points before the displacement of the
            % current time step, which is what is needed for the definition
            % of the computation of the displacement gradient
            DxMat = abs(diff(PIVresults.x{itime}(1:ny-1,:),1,2));
            
            DyMat = abs(diff(PIVresults.y{itime}(:,1:nx-1),1,1));
            
       
            % compute delta displacement
            % x direction is dimension 2
            % y direction is dimension 1
            dudx = diff(u{itime}(1:ny-1,:),1,2)./ DxMat; % du in x direction
            dudy = diff(u{itime}(:,1:nx-1),1,1)./ DyMat; % du in y direction
            dvdx = diff(v{itime}(1:ny-1,:),1,2)./ DxMat; % dv in x direction
            dvdy = diff(v{itime}(:,1:nx-1),1,1)./ DyMat; % dv in y direction
            
            % displacement gradient tensor H
            CellsEulerian.H{itime}(1,1,:,:) = dudx;
            CellsEulerian.H{itime}(1,2,:,:) = dudy;
            CellsEulerian.H{itime}(2,1,:,:) = dvdx;
            CellsEulerian.H{itime}(2,2,:,:) = dvdy;
        end
        
    case 'MidPoint'
        % midpoint: computes strain not in middle of cell, but at a
        % vertex
        %          i,j+1
        %            |
        %            |
        % i-1,j-----i,j-----i+1,j
        %            |
        %            |
        %          i,j-1
        %
        % in case of large deformation this approach may not be optimal
        
        
        
        % loop over time
        for itime=alltimes
           
            % two times width of cells, using current widths of cells
            % current points are the points before the displacement of the
            % current time step, which is what is needed for the definition
            % of the computation of the displacement gradient
            
             Dx2Mat = abs(diff(PIVresults.x{itime}(1:ny-2,:),1,2))+abs(diff(PIVresults.x{itime}(2:ny-1,:),1,2));
             Dy2Mat = (abs(diff(PIVresults.y{itime}(:,1:nx-2),1,1)) + abs(diff(PIVresults.y{itime}(:,2:nx-1),1,1)));
  
            
            dudx = (diff(u{itime}(1:ny-2,:),1,2) + diff(u{itime}(2:ny-1,:),1,2) )...
                ./ Dx2Mat; % du in x direction
            dudy = (diff(u{itime}(:,1:nx-2),1,1) + diff(u{itime}(:,2:nx-1),1,1) )...
                ./ Dy2Mat; % du in y direction
            dvdx = (diff(v{itime}(1:ny-2,:),1,2) + diff(v{itime}(2:ny-1,:),1,2)) ...
                ./ Dx2Mat; % dv in x direction
            dvdy = (diff(v{itime}(:,1:nx-2),1,1) + diff(v{itime}(:,2:nx-1),1,1)) ...
                ./ Dy2Mat; % dv in y direction
            
            % displacement gradient tensor H
            Cells.H{itime}(1,1,:,:) = dudx;
            Cells.H{itime}(1,2,:,:) = dudy;
            Cells.H{itime}(2,1,:,:) = dvdx;
            Cells.H{itime}(2,2,:,:) = dvdy;
        end
        
    case 'ShapeFunctions'
        
        % use linear shape functions for quadrilaterals (deformed rectangles)
        % t  4-----3
        % ^  |     |
        % |  |     |
        %    1-----2  -> s
        %
        % t and s are local coordinates of the four points of a quadrilateral
        % grid (initially square, but deforms as time progresses)
        % (both are in the range [-1 1] (FEM standards)
        %
        % global coordinates x and y are related to local coordinates t and s by shape functions:
        % N = [N1, N2, N3, N4] = 1/4 [(1 - s)(1 - t), (1 + s)(1 - t), (1 + s)(1 + t), (1 - s)(1 + t)]
        %
        % [x,y] = N [X_l Y_l]'
        %
        % similarly for the displacements
        % u = N d, where d are nodal displacements, u = [u1 u2 u3 u4], with
        % 1..4 denoting node number
        %
        %
        %
        % spatial derivatives of u and v in x and y directions
        %
        % | du/dx  du/dy | = | du/ds  du/dt | | ds/dx  ds/dy | = | du/ds  du/dt | | dx/ds  dx/dt |-1
        % | dv/dx  dv/dy |   | dv/ds  dv/dt | | dt/dx  dt/dy |   | dv/ds  dv/dt | | dy/ds  dy/dt |
        %
        % can thus be written as derivates of t and s
        %
        % d/dt = 1/4 [-(1-s) ,-(1+s) , (1+s) , (1-s)  ]
        % d/ds = 1/4 [-(1-t) ,(1-t)  ,(1+t) , -(1+t)  ]
        %
        % evaluate at s=0 and t=0, which is the centroid of the cell
        %
        % then simply multiply with nodal coordinates x, y or nodal
        % displacements u,v to compute the full derivatives
        % finally take the inverse of the second matrix
        
        ddt = 0.25 * [-1 -1 1 1];
        dds = 0.25 * [-1 1 1 -1];
        
        % loop over time
        %         iitime=0; % additional counter
        for itime=alltimes
            
            %             iitime=iitime+1;
            %             dt=Epochs.TimeSteps(iitime);
            CellsEulerian.u{itime}=NaN(ny-1,nx-1);
            CellsEulerian.v{itime}=NaN(ny-1,nx-1);
            
            % loop over cells
            
            for ix=1:nx-1
                for iy=1:ny-1
                    clear Node indexx indexy
                    Node.x=zeros(4,1);Node.y=zeros(4,1);Node.u=zeros(4,1);Node.v=zeros(4,1);
                    
                    % make indices
                    if Op.UpwardYAxis
                        indexy(1)=iy+1;
                        indexy(2)=iy+1;
                        indexy(3)=iy;
                        indexy(4)=iy;
                        
                        indexx(1)=ix;
                        indexx(2)=ix+1;
                        indexx(3)=ix+1;
                        indexx(4)=ix;
                    else
                        indexy(1)=iy;
                        indexy(2)=iy;
                        indexy(3)=iy+1;
                        indexy(4)=iy+1;
                        
                        indexx(1)=ix;
                        indexx(2)=ix+1;
                        indexx(3)=ix+1;
                        indexx(4)=ix;
                    end
                    % make nodal coordinate vectors and displacement vectors
                    for ivertex=1:4
                        Node.x(ivertex) = PIVresults.x{itime}(indexy(ivertex),indexx(ivertex));
                        Node.y(ivertex) = PIVresults.y{itime}(indexy(ivertex),indexx(ivertex));
                        Node.u(ivertex) =  u{itime}(indexy(ivertex),indexx(ivertex));
                        Node.v(ivertex) =  v{itime}(indexy(ivertex),indexx(ivertex));
                    end
                    
                    
                    
                    % check for NaN (due to mask or when point moved outside
                    % area of interest)
                    
                    if ~isempty(find(isnan(Node.x),1)) || ~isempty(find(isnan(Node.u),1))
                        Du=NaN(2);
                    else
                        
                        % spatial derivative of displacement
                        Matduvdts = [ dds * Node.u  ddt*Node.u ; ...
                            dds * Node.v  ddt*Node.v ];
                        
                        Matdxydtds = [ dds * Node.x  ddt*Node.x ; ...
                            dds * Node.y  ddt*Node.y ];
                        %  Du = Matduvdts * inv(Matdxydtds)
                        Du = Matduvdts/(Matdxydtds);
                        
                        % keep track when the cell has a solution
                        CellsEulerian.LastEpochData(iy,ix)=itime;
                        
                        % calculate velocity in middle of element
                        CellsEulerian.u{itime}(iy,ix) = 0.25 * ( Node.u(1) + Node.u(2) + Node.u(3) + Node.u(4));
                        CellsEulerian.v{itime}(iy,ix) = 0.25 * ( Node.v(1) + Node.v(2) + Node.v(3) + Node.v(4));
                        
                        
                    end
                    
                    % displacement gradient tensor
                    CellsEulerian.H{itime}(:,:,iy,ix) = Du;
               
                end
            end % end loop over grids
            
        end % end time loop
end % end switch


%% now setup deformation gradient tensor F
% F_ij = u_i,j + delta_ij
% in the remainder of this function, all tensors will be derived from F
% per time step is F 4 - dimensional
% dimension 1 and 2 represent the tensor dimensions, dimension 3 is the
% location in the cell in x direction, dimension 4 is the location of the
% cell in y direction
for itime=alltimes
    % incremental F
    CellsEulerian.FIncrmt{itime}(1,1,:,:) = CellsEulerian.H{itime}(1,1,:,:) + 1;
    CellsEulerian.FIncrmt{itime}(1,2,:,:) = CellsEulerian.H{itime}(1,2,:,:) ;
    CellsEulerian.FIncrmt{itime}(2,1,:,:) = CellsEulerian.H{itime}(2,1,:,:) ;
    CellsEulerian.FIncrmt{itime}(2,2,:,:) = CellsEulerian.H{itime}(2,2,:,:) + 1;
    
    
    
end % end loop time

%% compute stretch and strain tensors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% infinitesimal strain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CellsEulerian.RefType='Eulerian';
% compute incremental infinitesimal strain
if Op.InfinitesimalStrainIncrmt
    
    for itime=alltimes
        for ix=1:nx-1
            for iy=1:ny-1
                CellsEulerian.InfStrainIncrmt{itime}(:,:,iy,ix) = 0.5 * (CellsEulerian.FIncrmt{itime}(:,:,iy,ix) +  CellsEulerian.FIncrmt{itime}(:,:,iy,ix)' ) - eye(2,2);
                
            end % end loop iy
        end % end loop ix
        % vorticity (not a tensor)
        CellsEulerian.VorticityIncrmt{itime} = squeeze(0.5*(CellsEulerian.H{itime}(2,1,:,:) - CellsEulerian.H{itime}(1,2,:,:)));
   
    end % end loop time
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% finite strain %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Op.GreenFiniteStrainIncrmt
    % incremental finite strain
    
    for itime=alltimes
        for ix=1:nx-1
            for iy=1:ny-1
                CellsEulerian.GreenStrainIncrmt{itime}(:,:,iy,ix) = 0.5 * (CellsEulerian.FIncrmt{itime}(:,:,iy,ix)'*CellsEulerian.FIncrmt{itime}(:,:,iy,ix) - eye(2,2));
            end % end loop iy
        end % end loop ix
    end % end loop time
end



for itime=alltimes
    % initialize
    CellsEulerian.DilatationIncrmt{itime}=NaN(ny-1,nx-1);
    
    for ix=1:nx-1
        for iy=1:ny-1
            F=CellsEulerian.FIncrmt{itime}(:,:,iy,ix);
            
            % %
            % jacobian
            J = det(F);
            % %                         % % dilatation
            CellsEulerian.DilatationIncrmt{itime}(iy,ix) = J^(1/2);
        end % end iy
    end % end ix
end

end








