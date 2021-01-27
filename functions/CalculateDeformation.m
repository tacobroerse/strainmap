function [Cells]=CalculateDeformation(Points,Epochs,Op)
%CalculateDeformation calculate deformation from pixel displacements
% pixel displacements need to be converted from an Eulerian frame to a
% Lagrangian frame using FollowPoint.
%
% [Cells]=CalculateDeformation(Points,Epochs,Op) Points contains the
% coordinates and incremental material displacements, as calculated by
% FollowPoint. The displacements at epoch t lead to the coordinates at
% epoch t+1.
% Epochs contains time information.
% Op contains options. For this function it should contain:
% Op.DisplacementGradient: 'Simple','MidPoint', or 'ShapeFunctions', where
% only ShapeFunctions is recommended. ShapeFunctions uses bilinear shape
% functions for quadrilaterals. From the displacement gradient the
% deformation gradient F is computed, from which all other deformation
% measures are derived. These measures include:
% Op.InfinitesimalStrainIncrmt, logical for incremental infinitesimal
% strain
% Op.InfinitesimalStrain, logical for computing infinitesimal strain, only
% valid for small rotations and shear (<0.1).
% Op.FiniteStrainIncrmt, logical for computing incremental Green-Lagrangian finite strain
% Op.FiniteStrain, logical for computing Green-Lagrangian finite strain
% Op.FiniteStretchU and Op.FiniteStretchV are logicals for computing the
% left-stretch and right-stretch using polar decomposition of F.
%
% In every step the infinitesimal strain or Green-Lagrangian strain is
% computed for a quadrilateral in between 4 points, which originally starts out as
% a square, but evolves as a quadrilateral due to deformation. Mind that
% infinitesimal strain is based on the assumption of small rotations, and will not be valid
% for large rotations or shears > 0.1 .
% Green-Lagrangian strain is not sensitive to rotations, but is not
% linearly related to stretch.
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
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H
%

% time information
nEpochs = length(Points.u);
alltimes=Epochs.Index;
starttime = alltimes(1);

% size of grid
[ny,nx]=size(Points.x{starttime});

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

switch Op.DisplacementGradient
    
    case 'Simple'
        disp('Op.DisplacementGradient == Simple is not recommended')
        % middle points of deforming cells (at the first epoch)
        Cells.Midx = Points.x{starttime}(1:ny-1,1:nx-1)+0.5*Points.dx;
        Cells.Midy = Points.y{starttime}(1:ny-1,1:nx-1)+0.5*Points.dy;
    case 'MidPoint'
        disp('Op.DisplacementGradient == MidPoint is not recommended')
        % middle points of deforming cells are identical to point locations
        % (at the first epoch), without the points at the boundaries
        Cells.Midx = Points.x{starttime}(2:ny-1,2:nx-1);
        Cells.Midy = Points.y{starttime}(2:ny-1,2:nx-1);
    case 'ShapeFunctions'
        
        % middle points of deforming cells (at the first epoch)
        Cells.Midx = Points.x{starttime}(1:ny-1,1:nx-1)+0.5*Points.dx;
        Cells.Midy = Points.y{starttime}(1:ny-1,1:nx-1)+0.5*Points.dy;
        
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
        Cells.Connectivity=reshape(Connectivity,(nx-1)*(ny-1),4);
        
        % coordinates differ for each epoch, but indices do not
        % this vector can be used to convert an array into a vector
        % Vec = Array(Index)
        Cells.IndexingVertices=reshape(reshape([1:nx*ny],ny,nx)',nx*ny,1);
        
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
        Cells.Connectivity=reshape(Connectivity,(ny-2)*(nx-2),4);
        
        % coordinates differ for each epoch, but indices do not
        % this vector can be used to convert an array into a vector
        % Vec = Array(Index)
        Cells.IndexingVertices=reshape(reshape([1:nx*ny],ny,nx)',nx*ny,1);
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
            DxMat = abs(diff(Points.x{itime}(1:ny-1,:),1,2));
            
            DyMat = abs(diff(Points.y{itime}(:,1:nx-1),1,1));
            
            %    dt=Epochs.TimeSteps(itime);
            % compute delta displacement
            % x direction is dimension 2
            % y direction is dimension 1
            dudx = diff(Points.u{itime}(1:ny-1,:),1,2)./ DxMat; % du in x direction
            dudy = diff(Points.u{itime}(:,1:nx-1),1,1)./ DyMat; % du in y direction
            dvdx = diff(Points.v{itime}(1:ny-1,:),1,2)./ DxMat; % dv in x direction
            dvdy = diff(Points.v{itime}(:,1:nx-1),1,1)./ DyMat; % dv in y direction
            
            % displacement gradient tensor H
            Cells.H{itime}(1,1,:,:) = dudx;
            Cells.H{itime}(1,2,:,:) = dudy;
            Cells.H{itime}(2,1,:,:) = dvdx;
            Cells.H{itime}(2,2,:,:) = dvdy;
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
            Dx2Mat=zeros(ny,nx);
            Dy2Mat=zeros(ny,nx);
            % two times width of cells, using current widths of cells
            % current points are the points before the displacement of the
            % current time step, which is what is needed for the definition
            % of the computation of the displacement gradient
            
             Dx2Mat = abs(diff(Points.x{itime}(1:ny-2,:),1,2))+abs(diff(Points.x{itime}(2:ny-1,:),1,2));
            Dy2Mat = (abs(diff(Points.y{itime}(:,1:nx-2),1,1)) + abs(diff(Points.y{itime}(:,2:nx-1),1,1)));
%             
%             Dx2Mat(2:ny-1,:) = abs(diff(Points.x{itime}(:,:),1,2))+abs(diff(Points.x{itime}(:,:),1,2));
%             Dx2Mat(
%             Dy2Mat(2:ny-1,:) = (abs(diff(Points.y{itime}(:,1:nx-2),1,1)) + abs(diff(Points.y{itime}(:,2:nx-1),1,1)));
%             
            
            dudx = (diff(Points.u{itime}(1:ny-2,:),1,2) + diff(Points.u{itime}(2:ny-1,:),1,2) )...
                ./ Dx2Mat; % du in x direction
            dudy = (diff(Points.u{itime}(:,1:nx-2),1,1) + diff(Points.u{itime}(:,2:nx-1),1,1) )...
                ./ Dy2Mat; % du in y direction
            dvdx = (diff(Points.v{itime}(1:ny-2,:),1,2) + diff(Points.v{itime}(2:ny-1,:),1,2)) ...
                ./ Dx2Mat; % dv in x direction
            dvdy = (diff(Points.v{itime}(:,1:nx-2),1,1) + diff(Points.v{itime}(:,2:nx-1),1,1)) ...
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
            Cells.u{itime}=NaN(ny-1,nx-1);
            Cells.v{itime}=NaN(ny-1,nx-1);
            
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
                        Node.x(ivertex) = Points.x{itime}(indexy(ivertex),indexx(ivertex));
                        Node.y(ivertex) = Points.y{itime}(indexy(ivertex),indexx(ivertex));
                        Node.u(ivertex) = Points.u{itime}(indexy(ivertex),indexx(ivertex));
                        Node.v(ivertex) = Points.v{itime}(indexy(ivertex),indexx(ivertex));
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
                        Cells.LastEpochData(iy,ix)=itime;
                        
                        % calculate velocity in middle of element
                        Cells.u{itime}(iy,ix) = 0.25 * ( Node.u(1) + Node.u(2) + Node.u(3) + Node.u(4));
                        Cells.v{itime}(iy,ix) = 0.25 * ( Node.v(1) + Node.v(2) + Node.v(3) + Node.v(4));
                        
                        
                    end
                    
                    % displacement gradient tensor
                    Cells.H{itime}(:,:,iy,ix) = Du;
                    %                     % different elements
                    %                     Cells.dudx{itime}(iy,ix) = Du(1,1); % du in x direction
                    %                     Cells.dudy{itime}(iy,ix) = Du(1,2); % du in y direction
                    %                     Cells.dvdx{itime}(iy,ix) = Du(2,1); % dv in x direction
                    %                     Cells.dvdy{itime}(iy,ix) = Du(2,2); % dv in y direction
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
    Cells.FIncrmt{itime}(1,1,:,:) = Cells.H{itime}(1,1,:,:) + 1;
    Cells.FIncrmt{itime}(1,2,:,:) = Cells.H{itime}(1,2,:,:) ;
    Cells.FIncrmt{itime}(2,1,:,:) = Cells.H{itime}(2,1,:,:) ;
    Cells.FIncrmt{itime}(2,2,:,:) = Cells.H{itime}(2,2,:,:) + 1;
    
    % full cumulative deformation, achieved by product of all
    % previous FIncrmt, i.e F{n} = FIncrm{n}*FIncrm{n-1}*..*FIncr{1})
    % or F{n} = FIncrm{n}*F{n-1}
    if itime==alltimes(1)
        Cells.F{itime}=Cells.FIncrmt{itime};
    else
        for ix=1:nx-1
            for iy=1:ny-1
                Cells.F{itime}(:,:,iy,ix)=Cells.FIncrmt{itime}(:,:,iy,ix)*Cells.F{itime-1}(:,:,iy,ix);
            end
        end
        
    end
    
end % end loop time

%% compute stretch and strain tensors

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% infinitesimal strain %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute incremental infinitesimal strain
if Op.InfinitesimalStrainIncrmt
    
    for itime=alltimes
        for ix=1:nx-1
            for iy=1:ny-1
                Cells.InfStrainIncrmt{itime}(:,:,iy,ix) = 0.5 * (Cells.FIncrmt{itime}(:,:,iy,ix) +  Cells.FIncrmt{itime}(:,:,iy,ix)' ) - eye(2,2);
                
            end % end loop iy
        end % end loop ix
        % vorticity (not a tensor)
        Cells.VorticityIncrmt{itime} = squeeze(0.5*(Cells.H{itime}(2,1,:,:) - Cells.H{itime}(1,2,:,:)));
   
    end % end loop time
end

% compute current infinitesimal strain, only has a physical meaning when rotations are small

if Op.InfinitesimalStrain
    for itime=alltimes
        for ix=1:nx-1
            for iy=1:ny-1
                Cells.InfStrain{itime}(:,:,iy,ix) = 0.5 * (Cells.F{itime}(:,:,iy,ix) +  Cells.F{itime}(:,:,iy,ix)' ) - eye(2,2);
                
                % vorticity (not a tensor)
                RotTensor = 0.5 * (Cells.F{itime}(:,:,iy,ix) - Cells.F{itime}(:,:,iy,ix)' );
                Cells.Vorticity{itime}(iy,ix) = RotTensor(2,1);
                
            end % end loop iy
        end % end loop ix
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% finite strain %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Op.GreenFiniteStrainIncrmt
    % incremental finite strain
    
    for itime=alltimes
        for ix=1:nx-1
            for iy=1:ny-1
                Cells.GreenStrainIncrmt{itime}(:,:,iy,ix) = 0.5 * (Cells.FIncrmt{itime}(:,:,iy,ix)'*Cells.FIncrmt{itime}(:,:,iy,ix) - eye(2,2));
            end % end loop iy
        end % end loop ix
    end % end loop time
end

if Op.GreenFiniteStrain
    % compute current finite strain
    for itime=alltimes
        for ix=1:nx-1
            for iy=1:ny-1
                Cells.GreenStrain{itime}(:,:,iy,ix) = 0.5 * (Cells.F{itime}(:,:,iy,ix)'*Cells.F{itime}(:,:,iy,ix) - eye(2,2));
            end % end loop iy
        end % end loop ix
    end % end loop time
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% finite stretch and rotation %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Op.FiniteStretchU || Op.FiniteStretchV
    for itime=alltimes
        % initiate
        Cells.FiniteRotAngle{itime}=NaN(ny-1,nx-1);
        Cells.R{itime}=NaN(2,2,ny-1,nx-1);
        if Op.FiniteStretchU
            Cells.U{itime}=NaN(2,2,ny-1,nx-1);
        end
        if Op.FiniteStretchV
            Cells.V{itime}=NaN(2,2,ny-1,nx-1);
        end
        for ix=1:nx-1
            for iy=1:ny-1
                % polar decomposition
                if ~isnan(Cells.F{itime}(:,:,iy,ix)) % only if F has a value
                    %[R,U,V,~,RotAngleDeg]=PolarDecompositionNumeric(Cells.F{itime}(:,:,iy,ix));
                    [R,U,V,~,RotAngleDeg]=PolarDecompositionAnalytic(Cells.F{itime}(:,:,iy,ix));
                    Cells.R{itime}(1:2,1:2,iy,ix) = R;
                    if Op.FiniteStretchU
                        Cells.U{itime}(1:2,1:2,iy,ix) = U;
                    end
                    if Op.FiniteStretchV
                        Cells.V{itime}(1:2,1:2,iy,ix) = V;
                    end
                    Cells.FiniteRotAngle{itime}(iy,ix)=RotAngleDeg;
                end
            end % end loop iy
        end % end loop ix
    end % end loop time
end

% % finite strain rate
% if Op.FiniteStrainRate
%     for itime=alltimes
%         % initiate
%         Cells.FiniteStrainRate{itime}=NaN(2,2,ny-1,nx-1);
%
%         for ix=1:nx-1
%             for iy=1:ny-1
%                 % finite strain rate
%                 if itime==alltimes(1)
%                     F2=Cells.F{itime+1}(:,:,iy,ix);
%                     F1=Cells.F{itime}(:,:,iy,ix);
%                     dt=Epochs.TimeSteps(itime);
%                 elseif itime==nEpochs
%                     F2=Cells.F{itime}(:,:,iy,ix);
%                     F1=Cells.F{itime-1}(:,:,iy,ix);
%                     dt=Epochs.TimeSteps(itime-1);
%                 else
%                     F2=Cells.F{itime+1}(:,:,iy,ix);
%                     F1=Cells.F{itime-1}(:,:,iy,ix);
%                     dt=Epochs.TimeSteps(itime)+Epochs.TimeSteps(itime-1);
%                 end
%                 % numerical derivative, using midpoint rule
%                 Fdot=(F2-F1)/dt;
%                 if ~isnan(Fdot) % only if F has a value
%                     % strain rate E = 0.5*(F'*Fdot + Fdot'*F)
%                     F=Cells.F{itime}(:,:,iy,ix);
%                     Cells.FiniteStrainRate{itime}(1:2,1:2,iy,ix) = 0.5*(F'*Fdot+Fdot'*F);
%
%
%                 end
%             end % end loop iy
%         end % end loop ix
%     end % end loop time
% end

for itime=alltimes
    % initialize
    Cells.Dilatation{itime}=NaN(ny-1,nx-1);
    
    for ix=1:nx-1
        for iy=1:ny-1
            F=Cells.F{itime}(:,:,iy,ix);
            
            % %
            % jacobian
            J = det(F);
            % %                         % % dilatation
            Cells.Dilatation{itime}(iy,ix) = J^(1/2);
        end % end iy
    end % end ix
end

end

function [R,U,V,RotAngleRad,RotAngleDeg]=PolarDecompositionNumeric(F)
% polar decomposition using numerical algorithm
% based on svd (singular value decomposition
% returns rotation matrix R
% and right stretch tensor U (stretch in initial state)
% and left stretch tensor V (stretch in final state)
% that relate to F as
% F = R*U = V*R

[u,s,v] = svd(F);
% rotation matrix
R = u*v';
% right stretch matrix in original frame
U = v*s*v';
% left stretch matrix in final frame
V = R*U*R';

% F = V*R
% R = u*v'
% V = u*s*u' => V = R*U*R' = u*v' * v*s*v' * v*u' = u*I*s*I*u'

% as a bonus, the rotation angle
RotAngleRad = atan2(R(2,1),R(1,1)); % returns angle between pi and -pi
RotAngleDeg = RotAngleRad * 180/pi; % angle between -180 and 180 degrees
end


function [R,U,V,RotAngleRad,RotAngleDeg]=PolarDecompositionAnalytic(F)
% % polar decomposition using the Hoger & Carlson Method
%
% % step 1 compute symmetric right Cauchy-Green deformation tensor C
%
% % Cauchy-Green deformation tensor C
C = F'*F;

% % step 2 find eigenvalues of C

mu1 = 0.5 * (C(1,1) + C(2,2) + sqrt(4*C(1,2)*C(2,1) + (C(1,1) - C(2,2))^2) );
mu2 = 0.5 * (C(1,1) + C(2,2) - sqrt(4*C(1,2)*C(2,1) + (C(1,1) - C(2,2))^2) );

% step 3 compute invariants of Cauchy Green tensor C. eq 5.1 from Hoger &
% Carlson
IC = mu1 + mu2;
IIC = mu1*mu2;

% step 4 compute invariants of right stretch tensor U eq 5.2
IU = sqrt(IC + 2*sqrt(IIC));
IIU = sqrt(IIC);

% step 5 compute U, eq 3.3
U = (C+IIU*eye(2))/IU;

% step 6 compute U inverse, eq 4.1
Uinv = -IU/(IIU*(IIU+IC)+IIC)*(C-(IIU+IC)*eye(2));

% step 7 compute R (rotation matrix)
R = F*Uinv;

% step 8 compute V (left stretch tensor)
%V = R*U*R';
V = F*R';
% as a bonus, the rotation angle
RotAngleRad = atan2(R(2,1),R(1,1));
RotAngleDeg = RotAngleRad * 180/pi;

end
%






