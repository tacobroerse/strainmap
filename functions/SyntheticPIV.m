function PIVresults=SyntheticPIV(Op,SynOp)
%SyntheticPIV creates displacement data for synthetic cases
%
%   [PIVresults] = SyntheticPIV(Op,SynOp)
%       Op, structure with  options. See SetDefaults
%       Op.TypeSynthetic provides options for a number of synthetic cases
%       Op.TypeSynthetic='elongationx' elongation in x direction
%       Op.TypeSynthetic='elongationy' elongation in y direction
%       Op.TypeSynthetic='shearzoney' shear zone in y direction
%       Op.TypeSynthetic='shearzonex' shear zone in x direction
%       Op.TypeSynthetic='rotation' rotation
%       Op.TypeSynthetic='simpleshear' simple shear
%       Op.TypeSynthetic='smoothrotationshearzonex' rotating shear zone in x direction
%       Op.TypeSynthetic='smoothrotationshearzoney' rotating shear zone in y direction
%       Op.TypeSynthetic='selection_slip' selection of strain types, single epoch
%       SynOp, structure with options for syntetic fields
%       SynOp.nx is number of x stepts, SynOp.ny is number of y steps,
%       SynOp.totaldispl is total displacement, SynOp.totalangle is total
%       rotation angle, SynOp.ntimes is total nr of epochs
%       PIVresults returns a structure with coordinates and displacement
%       data 
%
%
% STRAINMAP
% 
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory: 
% https://doi.org/10.31223/X5FS3H

if nargin == 1
    ntimes = 50;
    nx = 50;
    ny = 50;
    totaldispl=5;
    totalangle=0;
    
else
    ntimes=SynOp.ntimes;
    nx=SynOp.nx;
    ny=SynOp.ny;
    totaldispl=SynOp.totaldispl;
    totalangle=SynOp.totalangle;
   
end

% make initial grid

PIVresults.dxy=1;
for itime=1:ntimes
    
    [PIVresults.x{itime},PIVresults.y{itime}]=ndgrid(linspace(0,nx,nx),linspace(0,ny,ny));
    
end

% we use ndgrids 
%
% x is ordered
% 
% x1 x1 x1 x1
% x2 x2 x2 x2
% ....
% xn xn xn xn
% 
% y is ordered
% 
% y1 y2 .. yn
% y1 y2 .. yn
% y1 y2 .. yn
% y1 y2 .. yn

PIVresults.units = '[px] respectively [px/frame]';

for itime=1:ntimes
    % displacements
    exponent=1;
    
    if strcmp(Op.TypeSynthetic,'elongationx')
        % u function of x
      %  normfactor=2/sum([1:ntimes])/(nx^exponent);
        normfactor=totaldispl/ntimes;
        PIVresults.u_filtered{itime} = normfactor*(itime/ntimes*PIVresults.x{1}).^exponent;
        PIVresults.v_filtered{itime} = zeros(size(PIVresults.y{1}));
    elseif strcmp(Op.TypeSynthetic,'elongationy')
        % v function of y
       % normfactor=2/sum([1:ntimes])/(ny^exponent);
       normfactor=totaldispl/ntimes;
        PIVresults.v_filtered{itime} = normfactor*(PIVresults.y{1});
        PIVresults.u_filtered{itime} = zeros(size(PIVresults.x{1}));
    elseif strcmp(Op.TypeSynthetic,'shearzoney')
        % shear zone
        width=0.8;
        normfactor=totaldispl;
        % normfactor=totaldispl/sum(1:ntimes);
        PIVresults.v_filtered{itime} = normfactor/ntimes*(erf((PIVresults.x{1}-nx/2)/(nx/2)*pi/width));
        PIVresults.u_filtered{itime} = zeros(size(PIVresults.x{1}));
        
    elseif strcmp(Op.TypeSynthetic,'shearzonex')
        % shear zone
        width=0.8;
        normfactor=totaldispl;
        %      normfactor=totaldispl/sum(1:ntimes);
        PIVresults.u_filtered{itime} = normfactor/ntimes*(erf((PIVresults.y{1}-ny/2)/(ny/2)*pi/width));
        PIVresults.v_filtered{itime} = zeros(size(PIVresults.x{1}));
        
    elseif strcmp(Op.TypeSynthetic,'rotation')
        % rotation
        % angle per time step
        theta=totalangle/ntimes;
        
        if itime==1
            xcentered=PIVresults.x{1}-nx/2;
            ycentered=PIVresults.y{1}-ny/2;
        end
        
        PIVresults.u_filtered{itime}=(xcentered * cosd(theta) + ycentered * -sind(theta))-xcentered;
        PIVresults.v_filtered{itime}=(xcentered * sind(theta) + ycentered * cosd(theta))-ycentered;
        
        
    elseif strcmp(Op.TypeSynthetic,'simpleshear')
        
        
        
        % centered coordinates
        if itime==1
            xcentered=PIVresults.x{1}-nx/2;
            ycentered=PIVresults.y{1}-ny/2;
        end
        
        % top goes to right with totaldispl, bottom goes to left with totaldispl, all linear
        normfactor = totaldispl/ny/2/ntimes;
        
        PIVresults.u_filtered{itime}=normfactor*ycentered;
        PIVresults.v_filtered{itime}=ycentered*0;
        
    elseif strcmp(Op.TypeSynthetic,'smoothrotationshearzonex')
        %
        %       % rotating shear zone
        % angle per time step
        theta=totalangle/ntimes;
        % current angle
        thetatot=theta*itime;
        % shear zone width
        width=0.8;
        % center coordinates
        if itime==1
            xcentered=PIVresults.x{1}-nx/2;
            ycentered=PIVresults.y{1}-ny/2;
        end
        
        
        
        % calculate rotation needed: new location (rotated) minus old
        rotation_u=(xcentered * cosd(theta) + ycentered * -sind(theta))-xcentered;
        rotation_v=(xcentered * sind(theta) + ycentered * cosd(theta))-ycentered;
        
        
        % shear zone in eulerian reference frame
        
        normfactor=totaldispl;
        
        % coordinates are rotated, so rotate to find local y
        %xrotated=cosd(thetatot)*xcentered + sind(thetatot)*ycentered;
        yrotated=-sind(thetatot)*xcentered + cosd(thetatot)*ycentered;
        
        % displacement is function of y rotated
        displ=normfactor/ntimes*(erf((yrotated)/(ny/2)*pi/width));
        
        % rotate shear distplacements back and add
        PIVresults.u_filtered{itime}=displ*cosd(thetatot)+rotation_u;
        PIVresults.v_filtered{itime}=displ*sind(thetatot)+rotation_v;
        
    elseif strcmp(Op.TypeSynthetic,'smoothrotationshearzoney')
        %
        %       % rotating shear zone
        % angle per time step
        theta=totalangle/ntimes;
        % current angle
        thetatot=theta*itime;
        % shear zone width
        width=1.2;
        % center coordinates
        if itime==1
            xcentered=PIVresults.x{1}-nx/2;
            ycentered=PIVresults.y{1}-ny/2;
        end
        
        
        % calculate rotation needed: new location (rotated) minus old
        rotation_u=(xcentered * cosd(theta) + ycentered * -sind(theta))-xcentered;
        rotation_v=(xcentered * sind(theta) + ycentered * cosd(theta))-ycentered;
        
        % shear zone in eulerian reference frame
        
        normfactor=totaldispl;
        
        % coordinates are rotated, so rotate to find local x
        xrotated=(cosd(thetatot)*xcentered + sind(thetatot)*ycentered);
        yrotated=-(sind(thetatot)*xcentered + cosd(thetatot)*ycentered);
        
        % displacement is function of y rotated
        displ=normfactor/ntimes*(erf((xrotated)/(nx/2)*pi/width));
        
        %          % rotate shear distplacements : u = Q*[0 ; displ] + 
        PIVresults.u_filtered{itime}=-displ*sind(thetatot)+rotation_u;
        PIVresults.v_filtered{itime}=displ*cosd(thetatot)+rotation_v;
    elseif strcmp(Op.TypeSynthetic,'rotationshearzoney')
        %
        %       % rotating shear zone
        % angle per time step
        theta=totalangle/ntimes;
        % current angle
        thetatot=theta*itime;
        % shear zone width
        width=0.8;
        % center coordinates
        if itime==1
            xcentered=PIVresults.x{1}-nx/2;
            ycentered=PIVresults.y{1}-ny/2;
        end
        
        
        % calculate rotation needed: new location (rotated) minus old
        rotation_u=(xcentered * cosd(theta) + ycentered * -sind(theta))-xcentered;
        rotation_v=(xcentered * sind(theta) + ycentered * cosd(theta))-ycentered;
        
        % shear zone in eulerian reference frame
        
        normfactor=totaldispl;
        
        % coordinates are rotated, so rotate to find local x
        xrotated=(cosd(thetatot)*xcentered + sind(thetatot)*ycentered);
        yrotated=-(sind(thetatot)*xcentered + cosd(thetatot)*ycentered);
        
        % displacement is function of y rotated
        displ = normfactor/ntimes*(((xrotated)/(nx/2)/width));
        
        %          % rotate shear distplacements
        PIVresults.u_filtered{itime}=-displ*sind(thetatot)+rotation_u;
        PIVresults.v_filtered{itime}=displ*cosd(thetatot)+rotation_v;
        
    elseif strcmp(Op.TypeSynthetic,'selection_slip')
        ntimes = 1;
        % make selection of slip
        
        
        
        % make new grid
        PIVresults=rmfield(PIVresults,'x');
        PIVresults=rmfield(PIVresults,'y');
        
        nSlip = linspace(-1,1,SynOp.nSlip);
        RotAngles = linspace(-SynOp.totalangle,SynOp.totalangle,SynOp.nRotAngles)*pi/180;
        RotMat = @(theta) [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
        
        
        ny=2*(2*SynOp.nRotAngles+1);
        nx=2*(2*(SynOp.nSlip+1)+1);
        PIVresults.x{itime}=NaN(ny,nx);
        PIVresults.y{itime}=NaN(ny,nx);
        PIVresults.u_filtered{itime}=NaN(ny,nx);
        PIVresults.v_filtered{itime}=NaN(ny,nx);
        
        for irotation = 1:length(RotAngles)
            for islip = 1:length(nSlip)+2
                % make coordinates with separate cells
                scalegrid=4;
                [tempx,tempy]=ndgrid(scalegrid*linspace(islip-1/scalegrid,islip,2),scalegrid*linspace(irotation-1/scalegrid,irotation,2));
                
                index1=scalegrid*[irotation-1/scalegrid irotation]-2;
                index2=scalegrid*[islip-1/scalegrid islip]-2;
                PIVresults.x{itime}(index1,index2)=tempx;
                PIVresults.y{itime}(index1,index2)=tempy;
                
                % centered coordinates
                xcentered=tempx-mean(tempx(:));
                ycentered=tempy-mean(tempy(:));
                
                % current rotation angle
                theta=RotAngles(irotation);
                
                if islip <= length(nSlip) + 1
                    
                    
                    %                 % make displacements
                    orientation = 0; % orientation of fault
                    % slip is composed from a strike slip and a stretch component (standing for
                    % normal or thrust slip under an arbitrary angle)
                    if islip <= length(nSlip)
                        strikeslipmagnitude=1;
                        dipslipmagnitude=1;
                        % strike slip contribution (linear)
                        strikeslip = (1-abs(nSlip(islip)))*strikeslipmagnitude;
                        % dip slip contribution (cannot be made linear)
                        dipslip = (1+dipslipmagnitude)^nSlip(islip);
                        % in local frame
                        deformF = [1 strikeslip ; 0 dipslip] ;
                    else
                        
                        % 1st dip slip contribution (cannot be made linear)
                        dipslip1 = 0.5;
                        % dip slip contribution (cannot be made linear)
                        dipslip2 = 2;
                        % in local frame
                        deformF = [dipslip1 0 ; 0 dipslip2] ;
                    end
                    % rotate to the current rotated frame
                    % deformFrot = RotMat(theta)*deformF*RotMat(theta)';
                    if orientation > pi/4 || orientation < -pi/4
                        error('make sure orientation of fault lies bewteen -pi/4 and pi/4')
                    else
                        % loop on vertices of cell
                        for iy=1:2
                            for ix=1:2
                                xy=[xcentered(iy,ix) ; ycentered(iy,ix)];
                                [uv] = deformF*xy-xy;
                                % add deformation displacement to the rotation
                                % displacements
                                PIVresults.u_filtered{itime}(index1(iy),index2(ix))=uv(1,:);
                                PIVresults.v_filtered{itime}(index1(iy),index2(ix))=uv(2,:);
                            end
                        end
                        
                    end
                
                    
                else
                    % finally, cells with no deformation, just rotation
                    for iy=1:2
                        for ix=1:2
                            PIVresults.u_filtered{itime}(index1(iy),index2(ix))=0;
                            PIVresults.v_filtered{itime}(index1(iy),index2(ix))=0;
                        end
                    end
                end
                % calculate rotation contribution: new location (rotated) minus old
                for iy=1:2
                    for ix=1:2
                        xy=[xcentered(iy,ix)+PIVresults.u_filtered{itime}(index1(iy),index2(ix)) ; ycentered(iy,ix)+PIVresults.v_filtered{itime}(index1(iy),index2(ix))];
                        [uv] = RotMat(theta) * xy - xy;
                        PIVresults.u_filtered{itime}(index1(iy),index2(ix))=PIVresults.u_filtered{itime}(index1(iy),index2(ix))+uv(1,:);
                        PIVresults.v_filtered{itime}(index1(iy),index2(ix))=PIVresults.v_filtered{itime}(index1(iy),index2(ix))+uv(2,:);
                    end
                end
                
                
            end % end loop on slip types
            
            
            
        end
        
        
        % break out of for loop
        
        break
    else
        error('unrecognized strain option')
    end
end


% show result
itime=ntimes;
%itime=1;
figure;
quiver(PIVresults.x{1},PIVresults.y{1},PIVresults.u_filtered{itime},PIVresults.v_filtered{itime},0)
title(strcat(Op.TypeSynthetic,' at time ',':',num2str(itime)))
% upward definition of y-axis
PIVresults.DirectionYaxis='Upward';
axis equal
axis tight

PIVresults.u_original=PIVresults.u_filtered;
PIVresults.v_original=PIVresults.v_filtered;

for itime=1:ntimes
    PIVresults.u_smoothed{itime}=[];
    PIVresults.v_smoothed{itime}=[];
end
% mask
for itime=1:ntimes
    PIVresults.typevector_original{itime}=ones(size(PIVresults.x{1}));
end
end

