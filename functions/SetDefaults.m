function [Op,Param,VideoSettings,Epochs,SynOp] = SetDefaults
%SetDefaults set defaults for the computation of lagrangian strain from eulerian
% displacements
%
%   [Op,Param,VideoSettings,Epochs] = SetDefaults returns
%       Op, structure with default options
%       Param, structure with default parameters
%       VideoSettings, structure with default settings for making videos
%       Epochs, empty structure, that will later be used to define time
%
%
% STRAINMAP
% 
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory: 
% https://doi.org/10.31223/X5FS3H

%%%%%%%%%%%%%%%%%%%
%%% general options
%%%%%%%%%%%%%%%%%%%
Op.Verbose=0; % providing extra information on steps taken in subroutines


%%%%%%%%%%%%%%%%%%%
% options related to input data
%%%%%%%%%%%%%%%%%%%

Op.Synthetic = 0; % whether or not to use synthetic fields
Op.UseImages=0; % load the original images that where used to construct displacement field
Param.ImageScale=1/2; % reduce image resolution
Op.PlotDisplacement=1; % plot du/dx du/dy dv/dx dv/dy
Op.UpwardYAxis = 1; % use upward y-axis instead of PIVdefault with downward axis
Op.FlipuDirection = 0; % flip direction of u
Op.FlipvDirection = 0; % flip direction of v
Op.Rotate = 0;% rotate field
Param.RotAngle = 0; % rotation angle, defined clockwise
Op.ConstrainArea = 0; % constrain area to smaller set
Op.FileType='Pivlab'; % type of input file, choose ['Pivlab','geotiff']
Op.MultipleInputFiles=0; % in case of 'geotiff' separate files may contain u and v displacements
Op.Mask=0; % use a mask 

% type of displacements
Op.DisplacementType='filtered'; % use filtered data (PIVLAB) [filtered,original,smoothed]
Op.Quantity='displacement';% input [displacements,velocities]


%%%%%%%%%%%%%%%%%%%
% time step information
Param.TimeStep = 1; % time step length
Param.TimeUnit='min';% time units of input [sec,min,days,years]
Param.TimeUnitPlot='min';% time units for output [sec,min,days,years]

% epoch information
Epochs = [];



%%%%%%%%%%%%%%%%%%%%%
% optional smoothing of input data
%%%%%%%%%%%%%%%%%%%%%
% displacement smoothing
Op.Smooth=0; % apply smoothing
Param.Smooth=1e1; % smoothing parameter

%%%%%%%%%%%%%%%%%%%
% outlier detection
%%%%%%%%%%%%%%%%%%%
Op.OutlierDetection=0; % experimental outlier detection
Param.WidthWindow=4; % sliding window (half)width
Param.ThresholdFactor=4; % threshold for outlier detection

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interpolation of input data fields to be able to construct material paths
%%%%%%%%%%%%%%%%%%%%%%%%%
Op.InterpolantType='linear'; % options see griddedInterpolant i.e. [linear,spline,cubic,makima]
Op.ExtrapolationMethod = 'none';% extrapolation for griddedInterpolant 
Param.DxScale=1; % increase in number of gridpoints wrt original


%%%%%%%%%%%%%%%%%%%
% type of strain definition
%%%%%%%%%%%%%%%%%%%
Op.DisplacementGradient = 'ShapeFunctions'; % [ShapeFunctions] FEM binlinear shape functions [ShapeFunctions, Simple, Midpoint]. Only ShapeFunctions is advised
Op.InfinitesimalStrainIncrmt = 1; % incremental infinitesimal
Op.InfinitesimalStrain=0; % infinitesimal (only valid for small shear/rotations)
Op.GreenFiniteStrainIncrmt=1; % Green finite strain incremental
Op.GreenFiniteStrain=0;% Green finite strain 
Op.FiniteStretchV=1; % finite left-stretch (polar decomposition)
Op.FiniteStretchU=0;% finite right-stretch (polar decomposition)




%

%%%%%%%%%%%%%%%%%%%
% figures settings
%%%%%%%%%%%%%%%%%%%

% general figure settings
Op.Coordinates='Updated'; % update coordinates in plots [Updated,original]
Op.SaveFigures=1; % saving figures 
Op.SavePng=1; % save figures as png
Param.FigureResolution='-r300';
Param.FigureWidth=900;
Param.LabelFontSize=15;
Param.TitleFontSize=15;
Param.ColorPercentile=99.5; % percentile for color maximum
Param.VectorScale=5;% scale vectors in plot
Param.GridColor='none'; % grid color
Op.ShowOriginalGrid=0;% show original grid in some plots (ViewDeformedGrid)
Op.MaxColorFigure='all_epochs'; % maximum color based on all epochs [all_epochs] or current epoch [all]
Op.ShowCoordinates=1;% show coordinates in plots
Op.UseFigureCoordLimits=0;% use different limits for coordinates in figures
Param.FigureXLim=[];% coordinate limits for figures
Param.FigureYLim=[];


% figures principal vectors / values
Op.ShowPrincipalVectors='gridnearest'; % where to show principal vectors, everywhere, or on a grid [gridnearest,all]
Param.NumVectorsPrincipalStrains=40; % number of points along x axis when Op.ShowPrincipalVectors='gridnearest'
Op.ColorPrincipalVectors='colormap';% type of colors for principal vectors [colormap,binary]
Op.PlotPrincplStretchAsStrain=0; % show principal stretches as principal strain
Op.PlotLogarithmicColorMapPrincplStretch=1;% show principal stretch using logarithmic color map
Op.PlotHenckyStrainInsteadOfStretch=0;% plot logarithmic principal stretches (Hencky strain) instead of finite principal stretch

% figures strain type
Op.IncludeTypeExpansionAndContraction=0; % include biaxial extension or shortening, 0=color map from uniaxial extension to uniaxial shortening over strike-slip
% 1=color map from biaxial to uniaxial extension to uniaxial to biaxial shortening over strike-slip
Param.PercentileTresholdTypeStrainPlot=95; % specifies the strain magnitude below which strain types are reduced in
% transparency.
Param.nColorsStrainType=13; % number of discrete colors for strain type plot

% image enhancement
Op.EnhanceContrastImage=1; % enhancement of image contrast
Param.ColorRangePerc=95; % value in range [0,100], percentile of color intensities that are stretched to the full RGB range

% multipanel figure temporal deformation
Op.UseSameYLimitsMultiPanel=1;

% video settings
VideoSettings=[];
VideoSettings.Quality=100;% video quality

% synthetic displacement field options
SynOp=[];
end

