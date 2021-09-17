function ExportGridAsGeoTiff(xGrid,yGrid,Field,Epochs,Param,filenamestr)
%ExportGridAsGeoTiff exports a regular grid to geotiff. First run
%InterpolateToRegularGrid to interpolate results on irregular grid to a
%regular grid.
% ExportGridAsGeoTiff(xGrid,yGrid,Field,Epochs,filenamestr)
% xGrid contains the x coordinate grid
% yGrid contains the y coordinate grid
% Field contains the field to be saved
% Epochs contains time information
% filenamestr contains a string for file names
% Param, structure with parameters. Param.SaveDir specifies directory for
% saving outputs
%
%
% Questions/bugs -> d.b.t.broerse@uu.nl
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H


alltimes=length(Field);

for itime=1:alltimes
    % create file name
    SaveDir=strcat(Param.SaveDir,'/output_geotiff/');
    if ~isdir(SaveDir)
        mkdir(SaveDir)
    end
    if isfield(Epochs,'Time')
        FileName=strcat(SaveDir,filenamestr,'_',datestr(Epochs.Time(itime)));
    else
        FileName=strcat(SaveDir,filenamestr,'_',num2str(itime));
    end
    
    if itime==1
        % create geographic raster reference object
        xlims=[min(xGrid,[],'all') max(xGrid,[],'all')];
        ylims=[min(yGrid,[],'all') max(yGrid,[],'all')];
        
        R = maprefpostings;
        R.XWorldLimits = xlims;
        R.YWorldLimits = ylims;
        R.ColumnsStartFrom = 'north';
        R.RasterSize = size(Field{itime});
        if ~isfield(Param,'CoordRefSysCode')
            disp('no CoordRefSysCode found in structure Params')
            disp('for codes, see matlabroot/toolbox/map/mapproj/projdata/epsg_csv')
            error(' ')
        end
    end
    % save as geotiff
    geotiffwrite(FileName,Field{itime},R,'CoordRefSysCode',Param.CoordRefSysCode)
end

end