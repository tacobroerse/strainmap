function [PIVresults,Epochs,PIVresultsorig]=SortPIVresults(PIVresults,Param,Op,Epochs,Area)
%SortPIVresults sorts displacement data and converts to displacements from
%velocities if necessary. Input data can be limited to a smaller area, or
%time range. Furthermore, the reference frame of the data can be rotated.
%
% [PIVresults,Epochs]=SortPIVresults(PIVresults,Param,Op)
%       PIVresults, structure with coordinates and displacement
%       data (or velocities). Output PIVresults will have displacements
%       only.
%       Op, structure with  options. See SetDefaults
%       Op.UpwardYAxis [0,1] defines whether y-axis should be pointing upwards, checked against data
%       Param, structure with default parameters
%       Param.RotAngle contains a rotation angle under which results can be
%       rotated if Op.Rotate==1, in the range [90,180,270,0]
%       Epochs, time structure. Epochs.Index will set to all epochs.
% [PIVresults,Epochs]=SortPIVresults(PIVresults,Param,Op,Epochs)
%       Epochs.Index contains all epochs that
%       should be included in the analysis. If supplied, PIVresults will be
%       cropped accordingly.
% [PIVresults,Epochs]=SortPIVresults(PIVresults,Param,Op,Epochs,Area)
%       if Op.ConstrainArea==1, Area.xlim and Area.ylim (2x1 vector) contain
%       the Area limits to which the data will be cropped.
%
%
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H
%


if ~isfield(PIVresults,'ConvertedVelocityToDisplacements')
    PIVresults.ConvertedVelocityToDisplacements=0;
end
if ~isfield(PIVresults,'DirectionYaxis')
    PIVresults.DirectionYaxis='Downward'; % default
end
if nargin < 4
    Epochs.Index=[1:length(PIVresults.x)];
end
if ~isfield(Epochs,'Index')
    Epochs.Index=[1:length(PIVresults.x)];
end
if nargin < 5
    Area=[];
end

if strcmp(Op.Quantity,'displacements')
    PIVresults.ConvertedVelocityToDisplacements=1;
end


% do a check on the index
dIndex=abs(diff(Epochs.Index));
if ~isempty(find(dIndex>1))
    error('index has to be continuous')
end

% make time vector
ntimes=length(PIVresults.x);
nobs=ntimes;


%% copy only the epochs of interest

PIVresultsorig=PIVresults;


% check whether velocities have already been converted to displacements
if PIVresults.ConvertedVelocityToDisplacements==0
    disp('convert velocities to displacements')
    disp(strcat('current units:',PIVresults.units))
    
    if strcmp(PIVresults.units,'[m] respectively [m/s]')
        if strcmp(Param.TimeUnit,'min')
            TimeScale=Param.TimeStep*60;
        elseif strcmp(Param.TimeUnit,'sec')
            TimeScale=Param.TimeStep;
        else
            error('not implemented time step unit')
        end
        % now convert all velocities to displacements
        for i=1:ntimes
            
            PIVresults.u_original{i}=PIVresultsorig.u_original{i}*TimeScale;
            PIVresults.v_original{i}=PIVresultsorig.v_original{i}*TimeScale;
            if isfield(PIVresultsorig,'u_filtered')
                PIVresults.u_filtered{i}=PIVresultsorig.u_filtered{i}*TimeScale;
                PIVresults.v_filtered{i}=PIVresultsorig.v_filtered{i}*TimeScale;
            end
            if isfield(PIVresultsorig,'u_smoothed')
                PIVresults.u_smoothed{i}=PIVresultsorig.u_smoothed{i}*TimeScale;
                PIVresults.v_smoothed{i}=PIVresultsorig.v_smoothed{i}*TimeScale;
            end
            
        end
        
        PIVresults.ConvertedVelocityToDisplacements=1;
        % change units
        
        PIVresults.units='[m] respectively [m]';
        PIVresults.unitsCoords='[m]';
        PIVresults.unitsDisplacement='[m]';
        disp(strcat('new units:',PIVresults.units))
    elseif strcmp(PIVresults.units,'[m] respectively [km/yr]')
        TimeScale=PIVresults.DiffTime;
        VelocityScale=1e3/365.24;
        % now convert all velocities to displacements
        for i=1:ntimes
            
            PIVresults.u_original{i}=PIVresultsorig.u_original{i}*VelocityScale*TimeScale;
            PIVresults.v_original{i}=PIVresultsorig.v_original{i}*VelocityScale*TimeScale;
            if isfield(PIVresultsorig,'u_filtered')
                PIVresults.u_filtered{i}=PIVresultsorig.u_filtered{i}*VelocityScale*TimeScale;
                PIVresults.v_filtered{i}=PIVresultsorig.v_filtered{i}*VelocityScale*TimeScale;
            end
            if isfield(PIVresultsorig,'u_smoothed')
                PIVresults.u_smoothed{i}=PIVresultsorig.u_smoothed{i}*VelocityScale*TimeScale;
                PIVresults.v_smoothed{i}=PIVresultsorig.v_smoothed{i}*VelocityScale*TimeScale;
            end
            
            
        end
        
        PIVresults.ConvertedVelocityToDisplacements=1;
        % change units
        
        PIVresults.units='[m] respectively [m]';
        PIVresults.unitsCoords='[m]';
        PIVresults.unitsDisplacement='[m]';
        disp(strcat('new units:',PIVresults.units))
    elseif strcmp(PIVresults.units,'[m] respectively [m/d]')
        
        VelocityScale=1;
        % now convert all velocities to displacements
        for i=1:ntimes
            TimeScale=PIVresults.DiffTime(i);
            PIVresults.u_original{i}=PIVresultsorig.u_original{i}*VelocityScale*TimeScale;
            PIVresults.v_original{i}=PIVresultsorig.v_original{i}*VelocityScale*TimeScale;
            if isfield(PIVresultsorig,'u_filtered')
                PIVresults.u_filtered{i}=PIVresultsorig.u_filtered{i}*VelocityScale*TimeScale;
                PIVresults.v_filtered{i}=PIVresultsorig.v_filtered{i}*VelocityScale*TimeScale;
            end
            if isfield(PIVresultsorig,'u_smoothed')
                PIVresults.u_smoothed{i}=PIVresultsorig.u_smoothed{i}*VelocityScale*TimeScale;
                PIVresults.v_smoothed{i}=PIVresultsorig.v_smoothed{i}*VelocityScale*TimeScale;
            end
            
            
        end
        
        PIVresults.ConvertedVelocityToDisplacements=1;
        % change units
        
        PIVresults.units='[m] respectively [m]';
        PIVresults.unitsCoords='[m]';
        PIVresults.unitsDisplacement='[m]';
        disp(strcat('new units:',PIVresults.units))
    elseif strcmp(PIVresults.units,'[px] respectively [px/frame]')
        disp('velocities are equal to displacements, no need for conversion')
        disp(strcat('current units:',PIVresults.units))
    else
        
        error('unknown units, exiting')
    end
    
    
else
    disp('velocities are already converted to incremental displacements')
    disp(strcat('current units:',PIVresults.units))
end

%% reverse y axis (because Images have a flipped y axis in matlab)
if Op.UpwardYAxis && strcmp(PIVresults.DirectionYaxis,'Downward')
    disp('change y-axis to pointing upward')
    
    miny=min(PIVresults.y{1},[],'all');
    maxy=max(PIVresults.y{1},[],'all');
    %
    %     % function for flipping
    reversey = @(y) -y+maxy+miny;
    %     % flip y and v
    for itime = 1:ntimes
        PIVresults.y{itime}= reversey(PIVresults.y{itime});
        % flip displacements in y direction: v
        PIVresults.v_original{itime}=-PIVresults.v_original{itime};
        if isfield(PIVresultsorig,'u_filtered')
            PIVresults.v_filtered{itime}=-PIVresults.v_filtered{itime};
        end
        if isfield(PIVresultsorig,'u_smoothed')
            PIVresults.v_smoothed{itime}=-PIVresults.v_smoothed{itime};
        end
    end
    % change status
    PIVresults.DirectionYaxis='Upward';
elseif Op.UpwardYAxis && strcmp(PIVresults.DirectionYaxis,'Upward')
    disp('y-axis already pointing upward')
end

%% rotate results when asked

if Op.Rotate
    disp(strcat('rotate results by:',num2str(Param.RotAngle),' degrees clockwise'))
    
    PIVresultsRot=PIVresults;
    if ~isfield(PIVresults,'Rotation')
        switch Param.RotAngle
            case -90
                PIVresultsRot.Rotation='rotated_90_counter_clockwise';
                %   maxx=max(PIVresults.x{1},[],'all'); % use this to prevent negative y coordinates
                for itime = 1:ntimes
                    % rotations are needed to keep ndgriddedness
                    % -y becomes x
                    PIVresultsRot.x{itime} = rot90(-PIVresults.y{itime},1);
                    % x becomes y
                    PIVresultsRot.y{itime} = rot90(PIVresults.x{itime},1);
                    % -v becomes u
                    PIVresultsRot.u_smoothed{itime} = rot90(-PIVresults.v_smoothed{itime},1);
                    PIVresultsRot.u_original{itime} = rot90(-PIVresults.v_original{itime},1);
                    PIVresultsRot.u_filtered{itime} = rot90(-PIVresults.v_filtered{itime},1);
                    % u becomes v
                    PIVresultsRot.v_smoothed{itime} = rot90(PIVresults.u_smoothed{itime},1);
                    PIVresultsRot.v_original{itime} = rot90(PIVresults.u_original{itime},1);
                    PIVresultsRot.v_filtered{itime} = rot90(PIVresults.u_filtered{itime},1);
                    % rotate typevector (mask)
                    PIVresultsRot.typevector_original{itime} = rot90(PIVresults.typevector_original{itime},-1);
                    
                end
            case 90
                error('not implemented yet')
            case 180
                error('not implemented yet')
            otherwise
                error('only angles of a multiple of 90 are possible')
        end
        % write back
        PIVresults=PIVresultsRot;
        clear PIVresultsRot
    else
        disp('rotation already performed, skipping')
    end
end

if Op.FlipuDirection
    disp('flip sign of u (x-displacements')
    % flip sign of u
    for i=1:ntimes
        PIVresults.u_original{i}=-PIVresults.u_original{i};
        if isfield(PIVresultsorig,'u_filtered')
            PIVresults.u_filtered{i}=-PIVresults.u_filtered{i};
        end
        if isfield(PIVresultsorig,'u_smoothed')
            PIVresults.u_smoothed{i}=-PIVresults.u_smoothed{i};
        end
    end
end

if Op.FlipvDirection
    disp('flip sign of v (y-displacements')
    % flip sign of v
    for i=1:ntimes
        PIVresults.v_original{i}=-PIVresults.v_original{i};
        if isfield(PIVresultsorig,'u_filtered')
            PIVresults.v_filtered{i}=-PIVresults.v_filtered{i};
        end
        if isfield(PIVresultsorig,'u_smoothed')
            PIVresults.v_smoothed{i}=-PIVresults.v_smoothed{i};
        end
    end
end

% select subset of data using based on location
if Op.ConstrainArea
    if isempty(Area)
        error('to trim area an Area has to be supplied')
    end
    
    Trimy=find(PIVresults.x{1}(1,:)<=Area.xlim(2) & PIVresults.x{1}(1,:)>=Area.xlim(1));
    Trimx=find(PIVresults.y{1}(:,1)<=Area.ylim(2) & PIVresults.y{1}(:,1)>=Area.ylim(1));
    
    disp('trim original area to area of interest')
    % trim all structures
    for i=1:ntimes
        PIVresults.x{i}=PIVresults.x{i}(Trimx(1):Trimx(end),Trimy(1):Trimy(end));
        PIVresults.y{i}=PIVresults.y{i}(Trimx(1):Trimx(end),Trimy(1):Trimy(end));
        
        if isfield(PIVresultsorig,'typevector_original')
            PIVresults.typevector_original{i}=PIVresults.typevector_original{i}(Trimx(1):Trimx(end),Trimy(1):Trimy(end));
        end
        
        PIVresults.u_original{i}=PIVresults.u_original{i}(Trimx(1):Trimx(end),Trimy(1):Trimy(end));
        PIVresults.v_original{i}=PIVresults.v_original{i}(Trimx(1):Trimx(end),Trimy(1):Trimy(end));
        if isfield(PIVresultsorig,'u_filtered')
            PIVresults.u_filtered{i}=PIVresults.u_filtered{i}(Trimx(1):Trimx(end),Trimy(1):Trimy(end));
            PIVresults.v_filtered{i}=PIVresults.v_filtered{i}(Trimx(1):Trimx(end),Trimy(1):Trimy(end));
        end
        if isfield(PIVresultsorig,'u_smoothed')
            PIVresults.u_smoothed{i}=PIVresults.u_smoothed{i}(Trimx(1):Trimx(end),Trimy(1):Trimy(end));
            PIVresults.v_smoothed{i}=PIVresults.v_smoothed{i}(Trimx(1):Trimx(end),Trimy(1):Trimy(end));
        end
    end
    
end

% make all times (using original number of time steps)
SingleTimeStep=1;
if isempty(Param.TimeStep)
    if numel(Param.TimeStep>1)
        SingleTimeStep=0;
    else
        SingleTimeStep=1;
    end
end
if numel(Param.TimeStep)>1
    SingleTimeStep=0;
    
end

if SingleTimeStep
    AllTimes = linspace(1,nobs,nobs)*Param.TimeStep;
else
    for itime=1:length(PIVresults.StartDate)
        AllTimes(itime)=datetime(str2num(PIVresults.StartDate{itime}(1:4)),...
            str2num(PIVresults.StartDate{itime}(5:6)),str2num(PIVresults.StartDate{itime}(7:8)));
        
    end
end

% select epochs
Epochs.Time=AllTimes;
Epochs.SelectedTimes=AllTimes(Epochs.Index);

% times
if isfield(PIVresults,'StartDate')
    Epochs.StartDate=PIVresults.StartDate;
    Epochs.EndDate=PIVresults.EndDate;
end

% make timesteps (first is always default)
if isfield(PIVresults,'DiffTime')
    Epochs.TimeSteps=PIVresults.DiffTime;
else
    Epochs.TimeSteps = [Param.TimeStep diff(Epochs.Time)];
end


% all existing times
Epochs.nfulltimes = length(PIVresults.x);

end

