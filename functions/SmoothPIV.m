function [PIVresults]=SmoothPIV(PIVresults,Epochs,Param,Op)
%SmoothPIV smooths the eulerian displacement field. Uses 3rd party function
%by D. Garcia, smoothn2 for robust spline smoothing. 
%
% [PIVresults]=SmoothPIV(PIVresults,Epochs,Param,Op)
%       PIVresults, structure with coordinates and displacement
%       data. 
%       Op, structure with  options. See SetDefaults
%       Op.Smooth==1 enables smoothing
%       Param, structure with default parameters
%       Param.Smooth contains smoothing parameter. See help smoothn2 for
%       explanation
%       Epochs, time structure. Epochs.Index contains indices of the data
%       that will be used in StrainMap.
%
% more info on smoothing see help smoothn2. Garcia D, Robust smoothing of gridded data in one and higher
%   dimensions with missing values. Computational Statistics & Data
%   Analysis, 2010;54:1167-1178. 
%
% STRAINMAP
% 
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory: 
% https://doi.org/10.31223/X5FS3H
%

if Op.Smooth
    % select type of displacement
    if strcmp(Op.StepType,'filtered')
        u=PIVresults.u_filtered;
        v=PIVresults.v_filtered;
    elseif strcmp(Op.StepType,'original')
        u=PIVresults.u_original;
        v=PIVresults.v_original;
    elseif strcmp(Op.StepType,'smoothed')
        u=PIVresults.u_smoothed;
        v=PIVresults.v_smoothed;
    else
        error('invalid option for Op.StepType')
    end
    disp(strcat('smooth displacements using smoothn, with parameter:',num2str(Param.Smooth)))
    for itime=Epochs.Index
        
        if Op.Mask
            
            % make mask
            Mask{itime}=ones(size(u{1}));
            Mask{itime}(PIVresults.typevector_original{itime}==0)=NaN;
            
            % mask out mask
            u{itime}=u{itime}.*Mask{itime};
            v{itime}=v{itime}.*Mask{itime};
            
        end
        
        % use smoothn for moving average filter
        [Z]=smoothn2({u{itime},v{itime}},Param.Smooth);
        PIVresults.u_smoothedlocally{itime} = Z{1};
        PIVresults.v_smoothedlocally{itime} = Z{2};
    end
else
    
    disp('smoothing disabled, continue without smoothing')
    
end
end

