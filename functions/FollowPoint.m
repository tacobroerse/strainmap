function [Points,GridCorrected]=FollowPoint(PIVresults,Op,Param,Epochs)
%FollowPoint creates the material paths using, interpolated, Eulerian grids of
%displacement.
%
% [Points,GridCorrected]=FollowPoint(PIVresults,Op,Param,Epochs)
%       PIVresults, structure with coordinates and displacement
%       data.
%       Op, structure with  options. See  SetDefaults
%       Param, structure with parameters. See  SetDefaults
%       Epochs, time structure.
%
% Follow points by tracking points throughout model evolution
% input are PIV eulerian displacements. Output is displacements in a
% Lagrangian frame (i.e., points follow the material, so that strain
% can be calculated). PIVresults are ordered as in PIVlab (Thielicke
% software). Because PIV results are provided in a non-deforming grid,
% interpolation is used to asses the displacements at material points.
% bi-linear, cubic and spline interpolation is possible (see
% GriddedInterpolant for more options). Bi-linear provides the largest
% dataset in case gaps are present in the displacement field.
%
% Outlier detection is provided by checking the time series in the original
% eulerian frame. I.e., we assume a relatively smooth eulerian displacement
% in time.
%
%
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H
%


Points=[];
GridCorrected=[];

Debug=0;
% ntimes
fulltimes = Epochs.nfulltimes; % all saved times
alltimes = Epochs.Index;
starttime = Epochs.Index(1);

%% grid properties
xmin=min(PIVresults.x{1}(:));
xmax=max(PIVresults.x{1}(:));
ymin=min(PIVresults.y{1}(:));
ymax=max(PIVresults.y{1}(:));


% determine type of input grid: MESHGRID vs. NDGRID
if PIVresults.x{1}(1,1) == PIVresults.x{1}(2,1)
    GridType='Meshgrid';
    % number of y and x grid points
    [ny,nx]=size(PIVresults.x{1});
else
    GridType='NDgrid';
    % number of y and x grid points
    [nx,ny]=size(PIVresults.x{1});
end

disp(strcat('type of grid:',GridType))

% select correct displacements depending on options
[u,v]=SelectDisplacements(PIVresults,Epochs,Op);


%% check outliers in original displacements (eulerian)

% matrix, for later use, which makes it easier to draw time series

% leave open whether grid is a meshgrid or NDgrid
[na,nb]=size(u{1});
umat=NaN(fulltimes,na,nb);
vmat=NaN(fulltimes,na,nb);


for itime = alltimes
    
    umat(itime,:,:)=u{itime};
    vmat(itime,:,:)=v{itime};
end


%% outlier check
% compare differences through time
% loop over points
GridCorrected.u=u;
GridCorrected.v=v;
GridCorrected.Outlier=zeros(fulltimes,na,nb);
GridCorrected.OutlierThisEpoch=zeros(fulltimes,1);

if Op.OutlierDetection
    % time
    epochs = alltimes;
    fullepochs = [1:fulltimes];
    for ia=1:na
        for ib=1:nb
            % outlier detection using pairs of double differences within a
            % sliding window
            
            % make time series for current (lagrangian) point
            utimeseries=umat(:,ia,ib);
            vtimeseries=vmat(:,ia,ib);
            
            
            
            % do check on NaNs
            indexavail = find(~isnan(utimeseries));
            
            % check if this point has data
            if ~isempty(indexavail)
                
                
                % minimum four points
                if length(indexavail) > 3
                    % check for outliers in both u and v
                    [OutlierVecu,utimeseriesFiltered,~] = FindReplaceOutliers(fullepochs,utimeseries,Param);
                    [OutlierVecv,vtimeseriesFiltered,~] = FindReplaceOutliers(fullepochs,vtimeseries,Param);
                    
                    % both u and v outlier vectors
                    % OutlierVec=OutlierVecu & OutlierVecv;
                    % any of u or v outlier vectors
                    OutlierVec = OutlierVecu & OutlierVecv;
                    % OutlierVec=isoutlier(uv(:,iy,ix),'movmedian',Param.WidthMovMedian,'ThresholdFactor',Param.ThresholdFactor);
                    EpochsOutlier=find(OutlierVec);
                    if ~isempty(EpochsOutlier)
                        % there is an outlier
                        GridCorrected.Outlier(EpochsOutlier,ia,ib)=1;
                        GridCorrected.OutlierThisEpoch(EpochsOutlier)=1;
                        
                        
                        % replace original values
                        
                        for itime = find(OutlierVec')
                            
                            GridCorrected.u{itime}(ia,ib) = utimeseriesFiltered(itime);
                            GridCorrected.v{itime}(ia,ib) = vtimeseriesFiltered(itime);
                        end
                    end
                end
            end
        end
    end
    
    
    
    %% plot corrections to initial eulerian grid
    
    for itime=alltimes
        if GridCorrected.OutlierThisEpoch(itime)
            outliers=squeeze(GridCorrected.Outlier(itime,:,:));
            disp(strcat('nr of outliers for epoch:',num2str(itime),':',num2str(length(find(outliers)))))
            
            if Debug
                figure; hold on
                scale=Param.VectorScale;
                % plot grid corrected for outliers
                quiver(PIVresults.x{itime},PIVresults.y{itime},GridCorrected.u{itime}*scale,GridCorrected.v{itime}*scale,0)
                % plot outliers (originals)
                quiver(PIVresults.x{itime}(outliers==1),PIVresults.y{itime}(outliers==1),u{itime}(outliers==1)*scale,v{itime}(outliers==1)*scale,0,'k','LineWidth',3)
                axis equal
                set(gca,'YDir','reverse') % because PIV results are in reverse y direction
                title(strcat('outliers (thick lines) in eulerian displacements, epoch:',num2str(itime)))
            end
        else
            disp(strcat('no outliers for epoch:',num2str(itime)))
        end
    end
    
    clear umat vmat u v
end

%% interpolate all fields
if alltimes(end) > 1
    for itime = alltimes
        
        
        % make gridded interpolant
        if strcmp(GridType,'Meshgrid')
            % transpose needed because an ND grid is needed for
            % griddedInterpolant
            if Op.UpwardYAxis == 0
                InterPolant.u{itime}=griddedInterpolant(PIVresults.x{itime}',PIVresults.y{itime}',GridCorrected.u{itime}',Op.InterpolantType,Op.ExtrapolationMethod);
                InterPolant.v{itime}=griddedInterpolant(PIVresults.x{itime}',PIVresults.y{itime}',GridCorrected.v{itime}',Op.InterpolantType,Op.ExtrapolationMethod);
            elseif Op.UpwardYAxis == 1
                
                % flip y-axis, because it has to be monitonically increasing
                InterPolant.u{itime}=griddedInterpolant(fliplr(PIVresults.x{itime}'),fliplr(PIVresults.y{itime}'),fliplr(GridCorrected.u{itime}'),Op.InterpolantType,Op.ExtrapolationMethod);
                InterPolant.v{itime}=griddedInterpolant(fliplr(PIVresults.x{itime}'),fliplr(PIVresults.y{itime}'),fliplr(GridCorrected.v{itime}'),Op.InterpolantType,Op.ExtrapolationMethod);
            else
                Op.UpwardYAxis
                error('invalid option for Op.UpwardYAxis')
            end
        elseif strcmp(GridType,'NDgrid')
            InterPolant.u{itime}=griddedInterpolant(PIVresults.x{itime},PIVresults.y{itime},GridCorrected.u{itime},Op.InterpolantType,Op.ExtrapolationMethod);
            InterPolant.v{itime}=griddedInterpolant(PIVresults.x{itime},PIVresults.y{itime},GridCorrected.v{itime},Op.InterpolantType,Op.ExtrapolationMethod);
        end
    end
end

%%%%%%%%%%%%%%%%
%% follow points
%%%%%%%%%%%%%%%%

%% initialise new grid (Points) (which may be denser or less dense than the original grid

if Param.DxScale==1
    % take the same grid as from PIV
    if strcmp(GridType,'Meshgrid')
        % take transpose, otherwise griddedInterpolant gets mad
        Points.x{starttime}=PIVresults.x{starttime}';
        Points.y{starttime}=PIVresults.y{starttime}';
        Points.u{starttime}=GridCorrected.u{starttime}';
        Points.v{starttime}=GridCorrected.v{starttime}';
    elseif strcmp(GridType,'NDgrid')
        % already in correct grid format
        Points.x{starttime}=PIVresults.x{starttime};
        Points.y{starttime}=PIVresults.y{starttime};
        Points.u{starttime}=GridCorrected.u{starttime};
        Points.v{starttime}=GridCorrected.v{starttime};
        
    end
    
    
    Points.dx=(xmax-xmin)/(nx-1);
    Points.dy=(ymax-ymin)/(ny-1);
    
    % final size of grid
    nx2=nx;
    ny2=ny;
else
    %make new grid
    
    % new mesh
    % number of points
    nx2=ceil((nx-1)*abs(Param.DxScale)+1);
    ny2=ceil((ny-1)*abs(Param.DxScale)+1);
    % x and y coordinates
    x=linspace(xmin,xmax,nx2);
    y=linspace(ymin,ymax,ny2);
    
    Points.dx=(xmax-xmin)/(nx2-1);
    Points.dy=(ymax-ymin)/(ny2-1);
    
    % mesh, using ndgrid
    [Points.x{starttime},Points.y{starttime}]=ndgrid(x,y);
    
    
end



%% update of coordinates and interpolation of displacement field
% add incremental displacement to coordinates
Points.LastEpochData=zeros(nx2,ny2);
for itime=alltimes
    if itime > starttime
        % update coordinates
        Points.x{itime}=Points.x{itime-1}+Points.u{itime-1};
        Points.y{itime}=Points.y{itime-1}+Points.v{itime-1};
    end
    
    % if itime > starttime % we have em already for the starttime
    % use these coordinates for interpolation
    Points.u{itime}=InterPolant.u{itime}(Points.x{itime},Points.y{itime});
    Points.v{itime}=InterPolant.v{itime}(Points.x{itime},Points.y{itime});
    
    % keep track of last epoch with data for each point
    isdata=~isnan(Points.u{itime});
    
    Points.LastEpochData(isdata)=itime;
    %  end
    
    
    
    
    
    % add additional coordinates to Points.x and y for final location
    if itime==alltimes(end)
        Points.x{itime+1}=Points.x{itime}+Points.u{itime};
        Points.y{itime+1}=Points.y{itime}+Points.v{itime};
    end
end


%% convert from ND grid to meshgrid, since that works better for plotting

% ND grid
% x is ordered
%
% x1 x1 x1 x1
% x2 x2 x2 x2
% ....
% xn xn xn xn

% y is ordered
%
% y1 y2 .. yn
% y1 y2 .. yn
% y1 y2 .. yn
% y1 y2 .. yn

% meshgrid


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


for itime=alltimes
    Points.x{itime}=Points.x{itime}';
    Points.y{itime}=Points.y{itime}';
    
    Points.u{itime}=Points.u{itime}';
    Points.v{itime}=Points.v{itime}';
    if itime==alltimes(end)
        Points.x{itime+1}=Points.x{itime+1}';
        Points.y{itime+1}=Points.y{itime+1}';
    end
end



Points.LastEpochData=Points.LastEpochData';
%% total displacement
% total displacement
%     Points.utot=zeros(size(Points.x{1}));
%     Points.vtot=zeros(size(Points.x{1}));
itime=starttime;
Points.utot{itime}=Points.u{itime};
Points.vtot{itime}=Points.v{itime};
for itime=alltimes(2:end)
    % update total displacement
    Points.utot{itime}=Points.utot{itime-1}+Points.u{itime};
    Points.vtot{itime}=Points.vtot{itime-1}+Points.v{itime};
end

fig=figure; hold on
for itime=alltimes
    % in scale
    quiver(Points.x{itime},Points.y{itime},Points.u{itime},Points.v{itime},0)
    % total displacement vector
    quiver(Points.x{starttime},Points.y{starttime},Points.utot{alltimes(end)},Points.vtot{alltimes(end)},0,'Color',ones(3,1)*0.5)
end

if strcmp(Op.UpwardYAxis,'Downward')
    set(gca,'YDir','reverse')
end

axis equal
title('incremental displacements (true scale)')
drawnow % this should cause matlab to wait until everything has been added to the plot


if Op.SaveFigures
    if ~isfolder(Param.SaveDir)
        mkdir(Param.SaveDir)
        disp(strcat('making folder:',Param.SaveDir))
    end
    SaveFigName=strcat(Param.SaveDir,'/incremental_lagrangian_displacements');
    savefig(fig,SaveFigName)
    % figures are rather large, so currently no saving
end

end

