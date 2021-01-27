function [CellsLaplace]=DecomposeLaplaceStretch(Cells,Points,Epochs,Op)
%DecomposeLaplaceStretch decomposes the deformation gradient into Laplace
%Stretch, an upper triangular stretch tensor with an easy to interpret
% stretch and shear components, and a rotation matrix.
% [CellsLaplace]=DecomposeLaplaceStretch(Cells,Epochs,Op)
% Cells contains principal stretches, computed by
% CalculateDeformation (stretch) and PrincipalStrain (principal stretch).
% Points contains information on the last epoch with data
% Epochs contains time information.
% Op, is an options structure.
% Op.DirectionType contains the option for choosing a reference
% orientation. As the Laplace stretch depends on the orientation of the
% reference system, we need additional information how the reference
% orientation should be. Ideally the x-axis should be aligned with the
% direction of simple shear, i.e. the direction of a shear zone or fault.
% Op.DirectionType='SteadyStateIncrementalStretch' takes the (spatially
% varying) reference orientation that leads to the smallest changes in
% incremental Laplace Stretch. This should lead to a direction that
% coincides with shear zone or fault directions.
% CellsLaplace is a structure that contains the decomposed deformation.
%
%
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2021

% size of cell array
[ny,nx]=size(Cells.Midx);

alltimes=Epochs.Index;
ntimes=length(alltimes);

% initialise
for itime = alltimes
    CellsLaplace.Stretch{itime}=NaN(2,2,ny,nx);
    CellsLaplace.R{itime}=NaN(2,2,ny,nx);
    CellsLaplace.RefAngle{itime}=NaN(ny,nx);
    CellsLaplace.a{itime}=NaN(ny,nx);
    CellsLaplace.b{itime}=NaN(ny,nx);
    CellsLaplace.gamma{itime}=NaN(ny,nx);
    CellsLaplace.theta{itime}=NaN(ny,nx);
end

% function for rotation matrix
RotMat = @(theta) [cos(theta) -sin(theta) ; sin(theta) cos(theta)];

% loop on cells

ixx=110;
iyy=275;

ixx=107;
ixx=117;
iyy=265;
iyy=170;
%ixx=214;
%iyy=142;

nxx=20;
nyy=140;
%for ix=ixx-nxx:ixx+nxx
  %  for iy=iyy-nyy:iyy+nyy
for ix=1:nx
    for iy=1:ny
        
        
        % deformation gradient
        ntimesi=Points.LastEpochData(iy,ix);
        if ntimesi >= 3
            % at least 3 time steps are needed
            F=cell(ntimesi,1);
            
            for itime=1:ntimesi
                
                F{itime} = Cells.F{itime}(1:2,1:2,iy,ix);
            end
            if strcmp(Op.DirectionType,'SteadyStateIncrementalStretch')
                % search for reference angle that optimizes steady state
                RefAngleMin=0;
                RefAngleMax=pi;
                % find angle that minimizes the change in incremental Laplace
                % Stretch
              %  RefAngle=fminbnd(@(RefAngle) ObjectiveFunction(RefAngle,F,ntimesi),RefAngleMin,RefAngleMax);
               % FOptions =  optimset('PlotFcns',@optimplotfval);
                [RefAngle,FVal,ExitFlag]=fminbnd(@(RefAngle) ObjectiveFunction(RefAngle,F,ntimesi),RefAngleMin,RefAngleMax);%,FOptions);
       
                % save for output
                CellsLaplace.RefAngle{itime}(iy,ix)=RefAngle;
                
                
            else
                error('not implemented option for Op.DirectionType')
            end
            
            % find Laplace Stretch using this angle for all epochs
            for itime=1:ntimesi
                % rotate to chosen reference frame
                F2 = RotMat(RefAngle) * F{itime} * RotMat(RefAngle)';
                
                % QR decomposition
                [U,R,a,b,gamma,theta]=QRdecompositionUpTriang(F2);
                
                % determine corotating reference direction
                
                
                % save for output
                CellsLaplace.Stretch{itime}(:,:,iy,ix)=U;
                CellsLaplace.R{itime}(:,:,iy,ix)=R;
                
                CellsLaplace.a{itime}(iy,ix)=a;
                CellsLaplace.b{itime}(iy,ix)=b;
                CellsLaplace.gamma{itime}(iy,ix)=gamma;
                CellsLaplace.theta{itime}(iy,ix)=theta;
            end
            
            
        end
    end% end iy
    
    
    
end % end ix
end

function f=ObjectiveFunction(RefAngle,F,ntimes)

% objective function to minimize. Here we minimize the change in
% incremental Laplace Stretch


% function for rotation matrix
RotMat = @(theta) [cos(theta) -sin(theta) ; sin(theta) cos(theta)];

U=cell(ntimes,1);
Uinv=cell(ntimes,1);
IncrU=cell(ntimes-1,1);
aIncr=zeros(ntimes-1,1);
bIncr=zeros(ntimes-1,1);
gammaIncr=zeros(ntimes-1,1);
 thetaincr=zeros(ntimes-1,1);
for itime=1:ntimes
    
    % rotate to a different reference frame orientation
    
    Frot=RotMat(RefAngle)*F{itime}*RotMat(RefAngle)';
    
    % now do a QR decomposition in this reference orientation
    
   % [U{itime},Uinv{itime}]=QRdecompositionUpTriangCompact(Frot);
    [~,theta]=QRdecompositionUpTriangCompact2(Frot);
    % make incremental Laplace Stretch
   
%     % check whether the first epoch can also be used 
     if itime>1
%         %         % U(i) = Uincr(i)
%         %         IncrU{itime}=eye(2)*Uinv{iref,itime-1};
%         %     else
%         % U(i)=Uincr(i)*U(i-1)
%         IncrU{itime-1}=U{itime}*Uinv{itime-1};
%         
%         % incremental stretch parameters
%         aIncr(itime-1)=IncrU{itime-1}(1,1);
%         bIncr(itime-1)=IncrU{itime-1}(2,2);
%         % incremental shear
%         gammaIncr(itime-1)=IncrU{itime-1}(1,2)/aIncr(itime-1);

         thetaincr(itime-1)=theta-thetaprev;
     end
    thetaprev=theta;
    
end
% change in incremental Laplace Stretch parameters
% aIncrDiff=diff(aIncr);
% bIncrDiff=diff(bIncr);
% gammaIncrDiff=diff(gammaIncr);
% 
% % objective function value, this may be improved later on
%StretchChange=abs(aIncrDiff)+abs(bIncrDiff)+abs(gammaIncrDiff);
% take the mean
%f=mean(StretchChange);

 f=sqrt(mean((thetaincr).^2));

 
end





function  [U,R,a,b,gamma,theta,Uinv]=QRdecompositionUpTriang(F)
% do upper triangular QR decomposition of deformation gradient into :
% a distortion V and rotation matrix R, via F=R*U
% see Freed  & Srinivasa 2015 (Logarithmic strain and its material
% derivative for a QR decompisition of the deformation gradient)
% notation follows Freed et al. 2020 Laplace Stretch Eulerian and
% Lagrangian formulations

U=zeros(2,2);

% right cauchy green deformation tensor
C = F'*F;

%  QT decomposition
% distortion U
U(1,1) = sqrt(C(1,1));
U(1,2) = C(1,2)/U(1,1);
U(2,2) = sqrt(C(2,2)-U(1,2)^2);

% get extensions
% this assumes an order of deformation
% first shear, then extensions
% i.e. F = R * Lambda * gamma
% Lambda = [a 0 ; 0 b]
% gamma = [1 gamma ; 0 1]
a=U(1,1);
b=U(2,2);
gamma=U(1,2)/U(1,1);

% rotation matrix R = F * inv(U)
Uinv=[1/U(1,1) -U(1,2)/(U(1,1)*U(2,2)) ;...
    0          1/U(2,2) ];
R=F*Uinv;


% rotation angle
theta = atan2(R(2,1),R(1,1));

end

function  [U,Uinv]=QRdecompositionUpTriangCompact(F)
% do upper triangular QR decomposition of deformation gradient into :
% a distortion V and rotation matrix R, via F=R*U
% see Freed  & Srinivasa 2015 (Logarithmic strain and its material
% derivative for a QR decompisition of the deformation gradient)
% notation follows Freed et al. 2020 Laplace Stretch Eulerian and
% Lagrangian formulations
% compact version for iterations

U=zeros(2,2);

% right cauchy green deformation tensor
C = F'*F;

%  QT decomposition
% distortion U
U(1,1) = sqrt(C(1,1));
U(1,2) = C(1,2)/U(1,1);
U(2,2) = sqrt(C(2,2)-U(1,2)^2);


% inverse of U
Uinv=[1/U(1,1) -U(1,2)/(U(1,1)*U(2,2)) ;...
    0          1/U(2,2) ];


end

function  [U,theta]=QRdecompositionUpTriangCompact2(F)
% do upper triangular QR decomposition of deformation gradient into :
% a distortion V and rotation matrix R, via F=R*U
% see Freed  & Srinivasa 2015 (Logarithmic strain and its material
% derivative for a QR decompisition of the deformation gradient)
% notation follows Freed et al. 2020 Laplace Stretch Eulerian and
% Lagrangian formulations
% compact version for iterations

U=zeros(2,2);

% right cauchy green deformation tensor
C = F'*F;

%  QT decomposition
% distortion U
U(1,1) = sqrt(C(1,1));
U(1,2) = C(1,2)/U(1,1);
U(2,2) = sqrt(C(2,2)-U(1,2)^2);


% inverse of U
Uinv=[1/U(1,1) -U(1,2)/(U(1,1)*U(2,2)) ;...
    0          1/U(2,2) ];
R=F*Uinv;


% rotation angle
theta = atan2(R(2,1),R(1,1));

end
