function [Cells]=GetTypeStrain(Cells,Epochs,Op)
%GetTypeStrain determines the type of strain (shortening, strike-slip,
% extension) from the logarithmic principal stretches, i.e. the Hencky
% strains, the logarithm of the principal stretches V.
% [Cells]=GetTypeStrain(Cells,Epochs,Op)
% Cells contains principal stretches, computed by
% CalculateDeformation (stretch) and PrincipalStrain (principal stretch).
% Epochs contains time information.
% Op, is an options structure. To compute the Hencky strain, the
% Finite Stretch V or U should be available, thus Op.FiniteStretchV == 1 or
% Op.FiniteStretchU == 1. To compute the incremental strain type, the Green
% Finite strain should be availbe, thus Op.GreenFiniteStrainIncrmt == 1.
% Cells output is the Cells.StrainType, in the range [-pi/2 pi/2]. For interpretations
% see PlotStrainType. Cells.MagnitudeStrain is the largest absolute
% logarithmic principal stretch. For the incremental strain type we have
% outputs Cells.StrainTypeIncrmt, and Cells.MagnitudeStrainIncrmt.
% Cells.HenckyStrainPrincplAngle provides the dominant deformation
% direction.
% 
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H

% size of cell array
[ny,nx]=size(Cells.Midx);


alltimes=Epochs.Index;


% initialise
for itime = alltimes
    Cells.StrainType{itime}=NaN(ny,nx);
    Cells.MagnitudeStrain{itime}=NaN(ny,nx);
    Cells.StrainTypeIncrmt{itime}=NaN(ny,nx);
    Cells.MagnitudeStrainIncrmt{itime}=NaN(ny,nx);
    Cells.DominantPrincplAngle{itime}=NaN(ny,nx);
end




% check availability of data
% get time vectors of maximum and minimum finite
% stretch
if Op.FiniteStretchV
    % use left-stretch tensor V
    CalcFiniteStrainType=1;
    PrincplStretchMax = Cells.PrincplStretchMaxV(~cellfun(@isempty, Cells.PrincplStretchMaxV));
    PrincplStretchMin = Cells.PrincplStretchMinV(~cellfun(@isempty, Cells.PrincplStretchMinV));
elseif Op.FiniteStretchU
    % use right-stretch tensor U
    CalcFiniteStrainType=1;
    PrincplStretchMax = Cells.PrincplStretchMaxU(~cellfun(@isempty, Cells.PrincplStretchMaxU));
    PrincplStretchMin = Cells.PrincplStretchMinU(~cellfun(@isempty, Cells.PrincplStretchMinU));
else
    CalcFiniteStrainType=0;
    warning(strcat('no Principal stretches available, set Op.FiniteStretchV and/or Op.FiniteStretchU', ...
        ' to 1 before running CalculateDeformation and PrincipalStrain'))
end
% incremental

if Op.GreenFiniteStrainIncrmt
    CalcIncrStrainType=1;
    PrincplStrainMaxGreenStrainIncr = Cells.PrincplStrainMaxGreenStrainIncr(~cellfun(@isempty, Cells.PrincplStrainMaxGreenStrainIncr));
    PrincplStrainMinGreenStrainIncr = Cells.PrincplStrainMinGreenStrainIncr(~cellfun(@isempty, Cells.PrincplStrainMinGreenStrainIncr));
else
    CalcIncrStrainType=0;
    warning(strcat('no incremental Green strain available, set Op.GreenFiniteStrainIncrmt', ...
        ' to 1 before running CalculateDeformation and PrincipalStrain'))
end

% loop on cells
for ix=1:nx
    for iy=1:ny
        % determine type of strain
        theta=NaN;
       
        if CalcFiniteStrainType
            % get time vectors of logarithm of maximum and minimum finite
            % stretch
            MaxStrain(alltimes)=cellfun(@(c) c(iy,ix),PrincplStretchMax);
            MinStrain(alltimes)=cellfun(@(c) c(iy,ix),PrincplStretchMin);
            
            
            
            % log
            MaxStrain=log(MaxStrain);
            MinStrain=log(MinStrain);
            
            
        end
        
        if CalcIncrStrainType
            % incremental
            MaxStrainIncr(alltimes)=cellfun(@(c) c(iy,ix),PrincplStrainMaxGreenStrainIncr);
            MinStrainIncr(alltimes)=cellfun(@(c) c(iy,ix),PrincplStrainMinGreenStrainIncr);
        end
        
        % loop over all epochs
        for itime = alltimes
            if CalcFiniteStrainType
                
                if abs(MinStrain(itime)) > abs(MaxStrain(itime))
                    % negative Hencky strain is larger than positive Hencky strain
                    phi=atan2(MaxStrain(itime),MinStrain(itime));
                else
                    % determine finite angle in 2D plane
                    phi=atan2(MinStrain(itime),MaxStrain(itime));
                end
                theta=nan;
                
                % convert to number in [-pi/2 pi/2] range
                % first rotate to get the strike angle (-pi/4) at 0 and pi
                if phi <= pi/4 && phi >= -pi/4
                    theta = phi + pi/4;
                elseif phi < -pi/4
                    theta = phi + pi/4;
                elseif phi >= pi/4
                    theta = phi -7/4*pi;
                end
                
                
                % mirror in y-axis
                if ~isnan(theta)
                    if abs(theta) > pi/2  % actually: theta < -pi/2
                        theta = -theta - pi;
                    end
                end
                Cells.StrainType{itime}(iy,ix)=theta;
                % the strain magnitude is here defined as the largest
                % absolute logarithmic strain
                
                [Cells.MagnitudeStrain{itime}(iy,ix)] = max([abs(log(Cells.PrincplStretchMaxV{itime}(iy,ix))) abs(log(Cells.PrincplStretchMinV{itime}(iy,ix)))]);
                
                % PrincipalStrain already calculates the direction of
                % largest principal stretch. However these will generally
                % be positive stretches (extensions), while as strains
                % (stretch - 1) these are better visible
                [~,imax]=max([abs(Cells.PrincplStretchMaxV{itime}(iy,ix)-1) abs(Cells.PrincplStretchMinV{itime}(iy,ix)-1)]);
                if imax == 1
                    % maximum principal stretch is also largest absolute
                    % strain
                    Cells.DominantPrincplAngle{itime}(iy,ix) = Cells.PrincplAngleDegV{itime}(iy,ix);
                else
                    % minimum principal strain is largest absolute strain
                    % 
                    % angle defined wrt x axis, counterclockwise
                    Cells.DominantPrincplAngle{itime}(iy,ix) = Cells.PrincplAngleDegV{itime}(iy,ix)-90;
                    if Cells.DominantPrincplAngle{itime}(iy,ix) < -90
                        % angle should be in range [-90 90]
                        Cells.DominantPrincplAngle{itime}(iy,ix)=Cells.DominantPrincplAngle{itime}(iy,ix)+180;
                    end
                end
                
            end % end CalcFiniteStrainType
            
            if CalcIncrStrainType
                
                % also incremental type
                
                % determine angle in complex plane
                phi=atan2(MinStrainIncr(itime),MaxStrainIncr(itime));
                theta=nan;
               
                % convert to number in [-pi/2 pi/2] range
                % first rotate to get the strike angle (-pi/4) at 0 and pi
                if phi <= pi/4 && phi >= -pi/4
                    theta = phi + pi/4;
                elseif phi < -pi/4
                    theta = phi + pi/4;
                elseif phi >= pi/4
                    theta = phi -7/4*pi;
                end
                
                % mirror in y-axis
                if ~isnan(theta)
                    if abs(theta) > pi/2  % actually: theta < -pi/2
                        theta = -theta - pi;
                    end
                end
                
                Cells.StrainTypeIncrmt{itime}(iy,ix) = theta;
                % magnitude: the largest absolute strain
                Cells.MagnitudeStrainIncrmt{itime}(iy,ix) = max([abs(MinStrainIncr(itime)) abs(MaxStrainIncr(itime))]);

            end % end CalcIncrStrainType
        end % end time
        %   end % end switch
        
    end % end iy
end % end ix




