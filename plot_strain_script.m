%%%%%%%%%%%%%%%%%%%%%%
% STRAIN MAP %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
% 
%
% The program StrainMap uses Eulerian, space-fixed, 2D displacements in time to
% construct Lagrangian, material-fixed, displacements, by tracing the
% material paths in time
% From the diplacement of points in time, StrainMap computes the 2D finite
% deformation of elements defined by the displaced points.
% Based on the principal stretch components, StrainMap classifies strain in
% a tectonically meaningfull way and qualitatively determines material
% undergoing: extension, strike-slip (shear), shortening, and oblique
% strike-slip, intermediates of strike-slip with perpendicular extension
% (transtension) or the intermediate of strike-slip with perpendicular
% shortening (transpression).
%
% This script is an example how to run StrainMap as script.
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020, d.b.t.broerse@uu.nl
% for more information on the theory: 
% https://doi.org/10.31223/X5FS3H

clear all
clc
close all
%% set paths
% main path to where the program is installed
cd '/Users/tacobroerse/MATLAB/work/strainmap/StrainMap_to_publish'

restoredefaultpath
userpath('reset')
% add the following to the path
addtopath={'./functions/','./figure_functions','./video_functions','./3rd_party_functions/','./additional_functions/','./inspect_figure_functions/'};
for i=1:length(addtopath)
    path(path,addtopath{i})
end




%% set defaults
[Op,Param,VideoSettings,Epochs,SynOp]=SetDefaults;
%% settings on input data and time information
%%%%%%%%%%%%%%%%%%%
% name of the experiment
Model.Name='Experiment1';

% input data file
PIV.File=strcat('PIVlab_',Model.Name);
% input data directory
PIV.Dir='yourdatalocationpathhere';
% time step information
Param.TimeUnit='min';
Param.TimeStep = 45/60; % 45 seconds 

%Epoch.Index=[1:10]; % smaller selection of epochs, than what is available


% if synthetic model is chosen, set a certain type of synthetic model
if Op.Synthetic
    
    Op.TypeSynthetic='elongationx';
    Op.TypeSynthetic='elongationy';
    Op.TypeSynthetic='shearzoney';
    Op.TypeSynthetic='shearzonex';
    Op.TypeSynthetic='rotation';
    Op.TypeSynthetic='simpleshear';
    Op.TypeSynthetic='smoothrotationshearzonex';
    Op.TypeSynthetic='smoothrotationshearzoney';
    
    Op.TypeSynthetic='selection_slip';
    Model.Name=strcat('synthetic_',Op.TypeSynthetic);
    SynOp=[];
end


Param.DxScale=2; % increase in number of gridpoints wrt original

% % special limits for figures
% Op.UseFigureCoordLimits=1;
Op.Rotate = 1; % rotate field
Param.RotAngle = -90;
% Param.FigureXLim=[-0.5 -0.2277];
% Param.FigureYLim=[0.35 0.875];

% settings for images 
Image.Dir=[];
Image.StartIndex=1;
Image.nImages=0; 


if ~Op.Synthetic
% file name
    PIV.File=strcat(PIV.Dir,PIV.File);
else
    PIV.Dir=strcat('synthetic_',Op.TypeSynthetic);
end
%% construct save directory
SmoothStr=[];GridRefinementStr=[];EpochStr=[];
if Op.Smooth
    SmoothStr=strcat('_smoothing_param_',num2str(Param.Smooth));
end
if Param.DxScale~=1
    GridRefinementStr=strcat('_grid_refinement_',num2str(Param.DxScale));
end

if ~isempty(Epochs) && ~isempty(Epochs.Index)
    if length(Epochs.Index) < Epochs.nfulltimes
       EpochStr=strcat('_epochs_',num2str(Epochs.Index(1)),'-',num2str(Epochs.Index(end))); 
    end
end
Param.Dir = strcat(PIV.Dir,Model.Name,'_Strain_results_',SmoothStr,GridRefinementStr,EpochStr);
Param.SaveDir = strcat(Param.Dir,'_Grid_',Op.DisplacementGradient,'/');
Param.SaveFile=strcat(Param.SaveDir,'saved_results.mat');


%% load PIV results
if Op.Synthetic
    PIVresults=SyntheticPIV(Op,SynOp);
else
    [PIVresults] = LoadPIVdata(PIV,Op);
end
%% load images
if Op.UseImages
    [Image]=LoadImages(Image,Param) ;
else 
    Image=[];
end

%% sort data change velocities into displacements
[PIVresults,Epochs]=SortPIVresults(PIVresults,Param,Op,Epochs);

%% smooth piv data using smoothn (same algorithm as in pivlab)
[PIVresults]=SmoothPIV(PIVresults,Epochs,Param,Op);

%% show difference with original
if Op.Smooth
    itime=20;
    Mask{itime}=ones(size(PIVresults.x{1}));
    Mask{itime}(PIVresults.typevector_original{itime}==0)=NaN;
    figure;
    set(gca,'YDir','reverse')
    hold on
    quiver(PIVresults.x{itime},PIVresults.y{itime},PIVresults.u_smoothedlocally{itime}.*Mask{itime}*Param.VectorScale,PIVresults.v_smoothedlocally{itime}.*Mask{itime}*Param.VectorScale,0,'Color','k','LineWidth',2)
    quiver(PIVresults.x{itime},PIVresults.y{itime},PIVresults.u_filtered{itime}.*Mask{itime}*Param.VectorScale,PIVresults.v_filtered{itime}.*Mask{itime}*Param.VectorScale,0,'Color','r','LineWidth',1)
end
%% show PIV displacement results
ViewPIVresults(PIVresults,Epochs,Param,Op)

%% follow points and detect outliers
[Points,PIVcorrected]=FollowPoint(PIVresults,Op,Param,Epochs);

%% show PIV displacement results with outliers removed
ViewPIVresults(PIVresults,Epochs,Param,Op,PIVcorrected)

%% compute deformations tensors
[Cells]=CalculateDeformation(Points,Epochs,Op);

%% compute principal strain
[Cells]=PrincipalStrain(Cells,Epochs,Op);

%% get type of strain using finite strain
[Cells]=GetTypeStrain(Cells,Epochs,Op);

%% save results
save(Param.SaveFile,'Epochs','Points','PIVresults','PIVcorrected','Param','Op','Cells','-v7.3')
%% INSPECTION FIGURES

%% inspect deformation
Param.GridColor='k';
InspectDeformation(Cells,Param,Op,Points,'IncrmtStrainType','final')
Param.GridColor='none';

%% FIGURES
%% plot displacement
PlotDisplacement(Cells,Param,Op,Points,'CumulativeDisplacement','final')
%% plot final stretch
PlotDeformationTensor(Cells,Param,Op,Points,'LeftStretchV','final')
%% plot incremental strain
PlotDeformationTensor(Cells,Param,Op,Points,'IncrementalInfinitesimalStrain','final')
%% plot strain type
PlotStrainType(Cells,Param,Op,Points,'StrainType','final',Epochs)
%% plot strain direction
PlotStrainType(Cells,Param,Op,Points,'StrainDirection','final',Epochs)
%% plot dilatation (area change)
PlotDeformationScalar(Cells,Param,Op,Points,'Dilatation','final')
%% plot average rotation
PlotDeformationScalar(Cells,Param,Op,Points,'MeanRotation','final')
%% plot stretch direction
PlotDeformationScalar(Cells,Param,Op,Points,'PrincStretchAngle','final')

%% plot principal values
PlotPrincipalValue(Cells,Points,Param,Op,'LeftStretchV','final',Epochs)
%% plot principal values
PlotPrincipalValue(Cells,Points,Param,Op,'IncrementalInfinitesimalStrain','final',Epochs)
%% plot principal vectors
PlotPrincipalVector(Cells,Points,Param,Op,'LeftStretchV','final',Epochs)
%% plot principal vectors
PlotPrincipalVector(Cells,Points,Param,Op,'IncrementalInfinitesimalStrain','final',Epochs)

%% MOVIES
%% make video of strain type
MakeMovieStrainSingleField(Points,Cells,Epochs,Param,Op,VideoSettings,'StrainType')
%% make video of incremental strain type
MakeMovieStrainSingleField(Points,Cells,Epochs,Param,Op,VideoSettings,'IncrmtStrainType')
%% make video of dominant strain direction
Param.PercentileTresholdTypeStrainPlot=80
MakeMovieStrainSingleField(Points,Cells,Epochs,Param,Op,VideoSettings,'StrainDirection')
Param.PercentileTresholdTypeStrainPlot=95;
%% make video of rotation
MakeMovieStrainSingleField(Points,Cells,Epochs,Param,Op,VideoSettings,'MeanRotation')
%% make video of dilatation
MakeMovieStrainSingleField(Points,Cells,Epochs,Param,Op,VideoSettings,'Dilatation')
%% make video of deformed grid
MakeMovieDeformedGrid(Image,Points,Cells,Epochs,Op,Param,VideoSettings)
%% make video of strain tensor: stretch
Param.GridColor='none';
MakeMovieTensor(Points,Cells,Epochs,Param,Op,VideoSettings,'LeftStretchV')
%% make video of strain tensor: incremental strain
MakeMovieTensor(Points,Cells,Epochs,Param,Op,VideoSettings,'IncrementalInfinitesimalStrain')

%% show deformed grid on top of images
% this has to be tested first on a data set with images
% if Op.UseImages
% ViewDeformedGrid(Points,Image,Cells,Epochs,Param,Op,PIVresults)
% end

%% strain type in time, for selected locations
close all
clear GridCells
Param.GridColor='none';
Op.ShowCoordinates=0;

   % define points, by selecting x and y indexes
    ipanel=1;
    GridCells.name{ipanel}='basal thrust';
    GridCells.ix(ipanel)=189;GridCells.iy(ipanel)=115;
    
    ipanel=2;
    GridCells.name{ipanel}='northern slope thrust wedge';
    GridCells.ix(ipanel)=210;GridCells.iy(ipanel)=100;
   
    
    ipanel=3;
    GridCells.name{ipanel}='first back-thrust';
    GridCells.ix(ipanel)=194;GridCells.iy(ipanel)=71;
    
    ipanel=4;
    GridCells.name{ipanel}='later back-thrust';
    GridCells.ix(ipanel)=193;GridCells.iy(ipanel)=39;


PlotTemporalDeformation(Cells,Param,Op,Points,Epochs,'StrainType',GridCells,'final')


