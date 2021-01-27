function [Outlierfullvec,yinterp,maddifffullvec,Interpolated,MovMedianThreshold] = FindReplaceOutliers(xorig,yorig,Param,Op)
%FindReplaceOutliers checks for outliers and replaces these with
%interpolated values. Also interpolates NaN. No extrapolation.
%
% [Outlierfullvec,yinterp,maddifffullvec,Interpolated] = FindReplaceOutliers(xorig,yorig,Param,Op)
%       xorig is the epoch vector
%       yorig is the vector of displacement at a single location
%       in space
%       Param, a structure that should contain :
%       Param.WidthWindow provides the sliding windows under which the 
%       differences are calculated
%       Param.Threshold: the threshold for outlier detection is the MAD of the
%       differences times this threshold
%
%       Outlierfullvec is a logical vector with the indexes of the outliers
%       yinterp is the interpolated vector with outliers removed
%       maddifffullvec is a vector with the median absolute difference
%       Interpolated is a logical vector indicating interpolated values
%       when NaNs are in input vector
%       MovMedianThreshold is a vector with the local threshold value to
%       determine outliers, based on moving median difference window
%       (Param.WidthWindow) and scale value (Param.Threshold)
%       
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


Outlierfullvec=zeros(size(yorig));
maddifffullvec=zeros(size(yorig));
yinterp=yorig;
Interpolated=zeros(size(yorig));

if exist('Op','var')
else
    Op.Debug=0;
end

% clean data by removing continuous NaN at the beginning and end of y
% do check on NaNs
indexnonan = find(~isnan(yorig));
if isempty(indexnonan)
    return
end
if length(indexnonan) < 4
    % at least four points needed
    return
end

% remove nan at start and end of y
y=yorig(indexnonan);
x=xorig(indexnonan);
%y = yorig([indexnonan(1):indexnonan(end)]);
xinterp = xorig([indexnonan(1):indexnonan(end)]);

% index for interpolated values
Interpolated([indexnonan(1):indexnonan(end)])=1; % slightly difficult way to do it
Interpolated(indexnonan)=0;

% check if there are still nan somewhere in between real values
isnaninbetween = ~isempty(find(isnan(x)));


% initialise matrix
diffmat=NaN(length(y));

% compute all double differences between all pairs 
for i=1:length(y)
    for j=1:length(y)
        if i~=j && abs(i-j) <= Param.WidthWindow
          diffmat(i,j)=(y(i)-y(j))/(i-j);
        end
    end
end

% calculate the median and the MAD (median absolute difference)
diffrow=cell(1,length(y));
medianabsdiff=NaN(1,length(y));
maddiff=NaN(1,length(y));


for i=1:length(y)
    diffrow{i}=diffmat(i,~isnan(diffmat(i,:)));
    medianabsdiff(i)=median(abs(diffrow{i}));
    maddiff(i)= median(abs(diffrow{i} - medianabsdiff(i)));
end
    
% determine outliers on the basis of the median of MAD times a factor
MovMedianThreshold=movmedian(maddiff,Param.WidthWindow)*Param.ThresholdFactor;
Outlier = maddiff > MovMedianThreshold ;



% if there are outliers, or gaps, interpolate
if ~isempty(find(Outlier, 1)) || isnaninbetween

    % replace outliers using shape-preserving piecewise cubic interpolation
    xx=x(~Outlier);
    yy=y(~Outlier);
    if length(xx) > 2
        y2 = interp1(xx,yy,xinterp,'pchip');
        % return to original size, in case of nans
        yinterp([indexnonan(1):indexnonan(end)])=y2;
       
    end
end

% put nan back
Outlierfullvec(indexnonan)=Outlier;
maddifffullvec(indexnonan)=maddiff;

if Op.Debug
maxc=max(abs(diffmat(:)));
figure;imagesc(diffmat);colorbar;caxis([-maxc maxc])
title('all differences')
end


end

