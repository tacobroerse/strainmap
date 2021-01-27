function [mapout]=makescientificcolormap(namemap,nsteps)
%makescientificcolormap
% easy script for loading scientific colormaps from F. Crameri.
% [mapout]=makescientificcolormap(namemap,nsteps)
% namemap: 'roma' (more can be added by adding .mat files in
% directory additional_functions
% [mapout]=makescientificcolormap(namemap,nsteps)
% nsteps: number of discrete steps to interpolate colormap.
%
% written by: Taco Broerse, 2020

if nargin == 1
    InterPolate = 0;
elseif nargin == 2
    InterPolate = 1;
end

Colormapfile=strcat('additional_functions/colormaps/',namemap,'.mat');
mapout=load(Colormapfile,'-mat');
mapout=struct2array(mapout);

% resample to desired number of steps
if InterPolate
    
    % interpolate colormap 
    mapout=interp1(linspace(0,1,size(mapout,1)),mapout,linspace(0,1,nsteps));
end