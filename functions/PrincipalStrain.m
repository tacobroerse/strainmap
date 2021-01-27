function [Cells]=PrincipalStrain(Cells,Epochs,Op)
%PrincipalStrain calculates principal strains and stretches 
% 
% [Cells]=PrincipalStrain(Cells,Epochs,Op) with Cells the structure that
% contains Stretch and Strain tensors, which are calculated by
% CalculateDefomration. 
% Epochs contains time information, such as set by SortPIVresults
% PrincipalStrain calculates eigen values and vectors for tensors depending on what
% has been set in the options structure Op:
% Op.InfinitesimalStrainIncrmt, logical for incremental infinitesimal
% strain
% Op.InfinitesimalStrain, logical for computing infinitesimal strain, only
% valid for small rotations and shear (<0.1).
% Op.FiniteStrainIncrmt, logical for computing incremental Green-Lagrangian finite strain
% Op.FiniteStrain, logical for computing Green-Lagrangian finite strain
% Op.FiniteStretchU and Op.FiniteStretchV are logicals for computing the
% left-stretch and right-stretch using polar decomposition of F.
%
% Principal (eigen) values and vectors are computed using Matlab's function
% eig. The PrincipalStrain function returns principal values, vectors and the angle of the
% largest eigen vector. This function sorts principal values according to
% the largest and smallest absolute values.
%
% Note that infinitesimal strains should only be used for small
% deformations and rotations. 
%
% Principal values of the Green-Lagrangian finite strain are related to
% finite principal stretches, but are not proportional. This function
% outputs the finite principal stretches, derived from the Green-Lagrangian
% principal values, by Principal stretch = sqrt(1+2*EigValue)-1.
%
% For large deformation, the most useful values are derived from the FiniteStretches, where the
% principal stretches are equal to the relative length change. Principal
% stretches are equal to finite principal strains + 1.
%
% Questions/bugs -> d.b.t.broerse@uu.nl
%
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H
%
[ny,nx]=size(Cells.Midx);

%
alltimes=Epochs.Index;
for itime=alltimes
    
    for ix=1:nx
        for iy=1:ny
            
     
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% infinitesimal strain %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % compute incremental infinitesimal strain eigen values
            if Op.InfinitesimalStrainIncrmt
                if isempty(find(isnan(Cells.InfStrainIncrmt{itime}(:,:,iy,ix))))
                    % determine principal strains
                    [EigValueMax,EigValueMin,EigVecMax,EigVecMin,PrincplAngle]=EigenVectors(Cells.InfStrainIncrmt{itime}(:,:,iy,ix));
                else
                    % no strains available for this cell
                    EigValueMax=NaN;EigValueMin=NaN;EigVecMax=NaN(1:2);EigVecMin=NaN(1:2);PrincplAngle=NaN;
                end
              
                % save to Cells structure
                Cells.PrincplStrainMaxInfStrainIncr{itime}(iy,ix) = EigValueMax;
                Cells.PrincplStrainMinInfStrainIncr{itime}(iy,ix) = EigValueMin;
                Cells.PrincplVecMaxInfStrainIncr{itime}(1:2,iy,ix) = EigVecMax;
                Cells.PrincplVecMinInfStrainIncr{itime}(1:2,iy,ix) = EigVecMin;
                Cells.PrincplAngleDegInfStrainIncr{itime}(iy,ix)=PrincplAngle;
                
                
            end
            
           
            % infinitesimal strain 
            if Op.InfinitesimalStrain
                if isempty(find(isnan(Cells.InfStrain{itime}(:,:,iy,ix))))
                    % determine principal strains
                    [EigValueMax,EigValueMin,EigVecMax,EigVecMin,PrincplAngle]=EigenVectors(Cells.InfStrain{itime}(:,:,iy,ix));
                else
                    % no strains available for this cell
                    EigValueMax=NaN;EigValueMin=NaN;EigVecMax=NaN(1:2);EigVecMin=NaN(1:2);PrincplAngle=NaN;
                end
              
                % save to Cells structure
                Cells.PrincplStrainMaxInfStrain{itime}(iy,ix) = EigValueMax;
                Cells.PrincplStrainMinInfStrain{itime}(iy,ix) = EigValueMin;
                Cells.PrincplVecMaxInfStrain{itime}(1:2,iy,ix) = EigVecMax;
                Cells.PrincplVecMinInfStrain{itime}(1:2,iy,ix) = EigVecMin;
                Cells.PrincplAngleDegInfStrain{itime}(iy,ix)=PrincplAngle;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% finite strain %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if Op.GreenFiniteStrainIncrmt
                % incremental finite strain

                if isempty(find(isnan(Cells.GreenStrainIncrmt{itime}(:,:,iy,ix))))
                    % determine principal strains
                    [EigValueMax,EigValueMin,EigVecMax,EigVecMin,PrincplAngle]=EigenVectors(Cells.GreenStrainIncrmt{itime}(:,:,iy,ix));
                else
                    % no strains available for this cell
                    EigValueMax=NaN;EigValueMin=NaN;EigVecMax=NaN(1:2);EigVecMin=NaN(1:2);PrincplAngle=NaN;
                end
                
                % save to Cells structure
                % conversion to strain see Malvern eq 4.5.18
                Cells.PrincplStrainMaxGreenStrainIncr{itime}(iy,ix) = sqrt(1+2*EigValueMax)-1;
                Cells.PrincplStrainMinGreenStrainIncr{itime}(iy,ix) = sqrt(1+2*EigValueMin)-1;
                Cells.PrincplVecMaxGreenStrainIncr{itime}(1:2,iy,ix) = EigVecMax;
                Cells.PrincplVecMinGreenStrainIncr{itime}(1:2,iy,ix) = EigVecMin;
                Cells.PrincplAngleDegGreenStrainIncr{itime}(iy,ix)=PrincplAngle;
            end
            
            if Op.GreenFiniteStrain
                %     % compute current finite strain principal strains
                if isempty(find(isnan(Cells.GreenStrain{itime}(:,:,iy,ix))))
                    % determine principal strains
                    [EigValueMax,EigValueMin,EigVecMax,EigVecMin,PrincplAngle]=EigenVectors(Cells.GreenStrain{itime}(:,:,iy,ix));
                else
                    % no strains available for this cell
                    EigValueMax=NaN;EigValueMin=NaN;EigVecMax=NaN(1:2);EigVecMin=NaN(1:2);PrincplAngle=NaN;
                end

                % save to Cells structure
                % conversion to strain see Malvern eq 4.5.18
                Cells.PrincplStrainMaxGreenStrain{itime}(iy,ix) = sqrt(1+2*EigValueMax)-1;
                Cells.PrincplStrainMinGreenStrain{itime}(iy,ix) = sqrt(1+2*EigValueMin)-1;
                Cells.PrincplVecMaxGreenStrain{itime}(1:2,iy,ix) = EigVecMax;
                Cells.PrincplVecMinGreenStrain{itime}(1:2,iy,ix) = EigVecMin;
                Cells.PrincplAngleDegGreenStrain{itime}(iy,ix)=PrincplAngle;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% finite stretch and rotation %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if Op.FiniteStretchU || Op.FiniteStretchV
                %     % compute current finite stretch principal
                %     stretches: right stretch tensor U
                if Op.FiniteStretchU
                    if isempty(find(isnan(Cells.U{itime}(:,:,iy,ix))))
                        % determine principal strains
                        [EigValueMax,EigValueMin,EigVecMax,EigVecMin,PrincplAngle]=EigenVectors(Cells.U{itime}(:,:,iy,ix));
                    else
                        % no strains available for this cell
                        EigValueMax=NaN;EigValueMin=NaN;EigVecMax=NaN(1:2);EigVecMin=NaN(1:2);PrincplAngle=NaN;
                    end
                    
                    % save to Cells structure
                    Cells.PrincplStretchMaxU{itime}(iy,ix) = EigValueMax;
                    Cells.PrincplStretchMinU{itime}(iy,ix) = EigValueMin;
                    Cells.PrincplVecMaxU{itime}(1:2,iy,ix) = EigVecMax;
                    Cells.PrincplVecMinU{itime}(1:2,iy,ix) = EigVecMin;
                    Cells.PrincplAngleDegU{itime}(iy,ix)=PrincplAngle;
                end
                if Op.FiniteStretchV
                    %     % compute current finite stretch principal
                    %     stretches: left stretch tensor V
                    if isempty(find(isnan(Cells.V{itime}(:,:,iy,ix))))
                        % determine principal strains
                        [EigValueMax,EigValueMin,EigVecMax,EigVecMin,PrincplAngle]=EigenVectors(Cells.V{itime}(:,:,iy,ix));
                    else
                        % no strains available for this cell
                        EigValueMax=NaN;EigValueMin=NaN;EigVecMax=NaN(1:2);EigVecMin=NaN(1:2);PrincplAngle=NaN;
                    end
                    
                    % save to Cells structure
                    Cells.PrincplStretchMaxV{itime}(iy,ix) = EigValueMax;
                    Cells.PrincplStretchMinV{itime}(iy,ix) = EigValueMin;
                    Cells.PrincplVecMaxV{itime}(1:2,iy,ix) = EigVecMax;
                    Cells.PrincplVecMinV{itime}(1:2,iy,ix) = EigVecMin;
                    Cells.PrincplAngleDegV{itime}(iy,ix)=PrincplAngle;
                end
            end
            
%             % finite strain rate
%             if Op.FiniteStrainRate
%                 if isempty(find(isnan(Cells.FiniteStrainRate{itime}(:,:,iy,ix))))
%                     % determine principal strains
%                     [EigValueMax,EigValueMin,EigVecMax,EigVecMin,PrincplAngle]=EigenVectors(Cells.FiniteStrainRate{itime}(:,:,iy,ix));
%                 else
%                     % no strains available for this cell
%                     EigValueMax=NaN;EigValueMin=NaN;EigVecMax=NaN(1:2);EigVecMin=NaN(1:2);PrincplAngle=NaN;
%                 end
%                 
%                 % save to Cells structure
%                 Cells.PrincplStrainRateMax{itime}(iy,ix) = EigValueMax;
%                 Cells.PrincplStrainRateMin{itime}(iy,ix) = EigValueMin;
%                 Cells.PrincplStrainRateVecMax{itime}(1:2,iy,ix) = EigVecMax;
%                 Cells.PrincplStrainRateVecMin{itime}(1:2,iy,ix) = EigVecMin;
%                 Cells.PrincplStrainRateAngleDeg{itime}(iy,ix)=PrincplAngle;
%             end
            
        end
    end
end


end

function [EigValueMax,EigValueMin,EigVecMax,EigVecMin,PrincplAngle]=EigenVectors(Tensor)
% This function determines eigen values, vectors and angles. It sorts
% according to the absolute eigen values.
EigValueMax=NaN;
EigValueMin=NaN;
EigVecMax=NaN;
EigVecMin=NaN;
PrincplAngle=NaN;
% compute eigen values and axes, and sort
if ~any(isnan(Tensor(:))) && ~any(isinf(Tensor(:)))
    % determine eigen vector using matlab's eig function
    [EigVector,EigValue]=eig(Tensor);
    
    % keep only diagonal
    EigValue=diag(EigValue);
    % check largest (absolute) eigen value
    [~,indexmax]=max(abs(EigValue));
    [~,indexmin]=min(abs(EigValue));
    
    % exception, for simple shear both values are identical
    if indexmax == indexmin
        % skip the absolute value and take the maximum positive
        % value
        [~,indexmax]=max((EigValue));
        [~,indexmin]=min((EigValue));
    end
    
    % store eigen value
    EigValueMax=EigValue(indexmax);
    EigValueMin=EigValue(indexmin);
    %
    % store eigen vectors
    EigVecMax(1:2)=squeeze(EigVector(1:2,indexmax));
    EigVecMin(1:2)=squeeze(EigVector(1:2,indexmin));
    
    % principal axis orientation, defined wrt x axis, counterclockwise
    PrincplAngle = atand(EigVecMax(2)/EigVecMax(1));
end
end

