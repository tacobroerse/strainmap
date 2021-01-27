function [u,v]=SelectDisplacements(PIVresults,Epochs,Op)
%SelectDisplacements selects the correct displacements depending on the settings in Op
%
% [u,v]=SelectDisplacements(PIVresults,Epochs,Op)
%       PIVresults, structure with coordinates and displacement
%       data. 
%       Epochs, time structure. 
%       Op, structure with  options. See SetDefaults
%       u and v are displacements, selected from PIVresults
%
%
% STRAINMAP
% 
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory: 
% https://doi.org/10.31223/X5FS3H
% 
Verbose=Op.Verbose;
if Op.Smooth
    if isfield(PIVresults,'u_smoothedlocally')
        % smoothed field produces by SmoothPIV (instead of smoothing in
        % PIVlab)
        u=PIVresults.u_smoothedlocally;
        v=PIVresults.v_smoothedlocally;
        if Verbose
            disp('select locally smoothed displacements')
        end
    else
        disp('WARNING: there should be a locally smoothed field, but that was not found')
        
        disp('please first run SmoothPIV first, or set Op.Smooth to false')
        disp('continue with unsmoothed data')
        Op.Smooth = 0;
    end
end
if ~Op.Smooth
    if strcmp(Op.DisplacementType,'filtered')
        u=PIVresults.u_filtered;
        v=PIVresults.v_filtered;
        if Verbose
            disp('select filtered displacements')
        end
    elseif strcmp(Op.DisplacementType,'original')
        u=PIVresults.u_original;
        v=PIVresults.v_original;
        if Verbose
            disp('select original displacements')
        end
    elseif strcmp(Op.DisplacementType,'smoothed')
        u=PIVresults.u_smoothed;
        v=PIVresults.v_smoothed;
        if Verbose
            disp('select PIVlab smoothed displacements')
        end
    else
        error('invalid option for Op.DisplacementType')
    end
end

%% use the same mask as used in PIVlab

alltimes = Epochs.Index;

% mask can be time variable
if Op.Mask
    if Verbose
        disp('apply mask')
    end
    
    for itime = alltimes
        % make mask
        Mask{itime}=ones(size(u{1}));
        Mask{itime}(PIVresults.typevector_original{itime}==0)=NaN;
      
        % mask out mask
        u{itime}=u{itime}.*Mask{itime};
        v{itime}=v{itime}.*Mask{itime};
    end
end

end

