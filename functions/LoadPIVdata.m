function [PIVresults] = LoadPIVdata(PIV,Op)
%LoadPIVdata loads displacement (or velocity) data from PIVlab or geotiff
%
%   [PIVresults] = LoadPIVdata(PIV,Op)
%       PIV, structure with containing PIV.File containing the file name
%       including path (PIVlab data), or PIV.Dir (geotiff) containing the
%       path to the data.
%       and if Op.Synthetic = 1, a sythetic field is constructed.
%       Op, structure with  options. See SetDefaults
%       PIVresults returns a structure with coordinates and displacement
%       data (or velocities)
%
%
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H

PIVresults = [];
if Op.Synthetic
    % load synthetic test
    warning('Op.Synthetic set to true, stopping')
    return
else
    
end


if strcmp(Op.FileType,'Pivlab')
    % do not load other components from PIV lab because these are not used
    PIVresults=load(PIV.File,'u_original','v_original','x','y','u_smoothed','v_smoothed','u_filtered','v_filtered','units','typevector_original');
    
elseif strcmp(Op.FileType,'geotiff')
    if ~Op.MultipleInputFiles
        error('not yet implemented')
        %[fields,R] = geotiffread(PIV.File);
    else
        % set units
        PIVresults.units = strcat('[',Op.PositionUnits,'] respectively [',Op.VelocityUnits,']');
        disp(strcat('current displacement/velocity units:',PIVresults.units))
        % check all files that are available in the directory
        files=dir(strcat(PIV.File,'*'));
        nfiles=length(files);
        % open all files
        ktime=0;
        for ifile=1:nfiles
            % get dates from string
            date1 = files(ifile).name(PIV.date1index);
            date2 = files(ifile).name(PIV.date2index);
            % check whether the date is a number
            if ~isnumeric(str2num(date1))
                error('date is not numeric, please change PIV.date1index')
            end
            if ~isnumeric(str2num(date2))
                error('date is not numerica, please change PIV.date2index')
            end
            % check whether dates are
            datenum1 = datenum(date1,'yyyymmdd');
            datenum2 = datenum(date2,'yyyymmdd');
            
            % the total time represented by the data
            timespan = datenum2 - datenum1;
            
            % notification
            disp(' ')
            disp(['loading:',files(ifile).name])
            disp(['dates:',num2str(date1),':',num2str(date2),' period [days] ',num2str(timespan)])
            % only for images this is downward
            PIVresults.DirectionYaxis='Upward';
            if ifile == 1
                [temp,R]=geotiffread(strcat(PIV.Dir,files(ifile).name));
                includedata = 1;
                % use initial time difference
                % timespandefault = timespan;
                
            else
                % check whether the current date follows on the previous and
                % has the right time span
                % check whether dates are sequential
                if datenum1==datenum2prev+1
                    
                    
                    
                    % load all files
                    [temp,R]=geotiffread(strcat(PIV.Dir,files(ifile).name));
                    includedata = 1;
                    
                    
                    %                 elseif strcmp(PIV.Settings,'12_day_intervals')
                    %                     if strcmp(date1,date2prev) && timespandefault == timespan
                    %
                    %                         [temp,R]=geotiffread(strcat(PIV.Dir,files(ifile).name));
                    %                         includedata = 1;
                    %                         disp(strcat('loading:',files(ifile).name))
                    %                     else
                    %                         % disp(strcat('skip loading:',files(ifile).name))
                    %                         includedata = 0;
                    %                     end
                    %                 end
                    
                else
                    warning('dates are not sequential')
                    warning(['start date this file:',num2str(date1)])
                    warning(['previous date last file:',num2str(date2prev)])
                    error('quitting')
                end
            end
            
            if includedata
                ktime=ktime+1;
                PIVresults.DiffTime(ktime) = timespan;
                % store times
                PIVresults.StartDate{ktime} = files(ktime).name(PIV.date1index);
                PIVresults.EndDate{ktime} = files(ktime).name(PIV.date2index);
                
                
                % set velocities (conversion to displacement occurs later,
                % in SortPIVresults.m)
                PIVresults.u_original{ktime} = temp(:,:,1);
                PIVresults.v_original{ktime} = temp(:,:,2);
                if ktime == 1
                    [ny,nx]=size(PIVresults.u_original{1});
                end
                % set coordinates using R
                sizex=(R.CellExtentInWorldX);
                sizey=(R.CellExtentInWorldY);
                xvec = [R.XWorldLimits(1)+sizex/2:sizex:R.XWorldLimits(2)-sizex/2];
                yvec = [R.YWorldLimits(1)+sizey/2:sizey:R.YWorldLimits(2)-sizey/2];
                if numel(xvec) ~= nx
                    error('length x vector not equal to number of x in velocity field')
                end
                if numel(yvec) ~= ny
                    error('length y vector not equal to number of y in velocity field')
                end
                % make x and y arrays
                
                [PIVresults.y{ktime},PIVresults.x{ktime}]=ndgrid(yvec,xvec);
                
                
                % store previous date
                date2prev=date2;
                datenum2prev=datenum2;
                
            end
        end
    end
else
    error('unknown file type')
    
end
end


