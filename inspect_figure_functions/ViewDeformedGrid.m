function ViewDeformedGrid(Points,Image,Cells,Epochs,Op,PIVresults)
% show deformed grid on top of original images
if Op.UseImages
    nPoints = length(Points.x);
    nImages = length(Image.Scaled);
    nepochs=nPoints-1;
    % npoints should be nimages minus 1
    if nPoints ~= nImages - 1
        nPoints
        nImages
        disp('number of points not equal to number of images minus one')
        
        % checking which epochs from the grid correspond to which image
        % epochs
        disp('will now check which epochs of the images best fit with the points epochs')
        
        timethreshold = mean(diff(Epochs.Time))*0.5;
        for itime=1:nPoints-1
            difftime=1e10;
            Image.TimeCorrIndex(itime)=NaN;
            for jtime=1:nImages
                difftimeij=abs(Image.Dates(jtime)-Epochs.Time(itime));
                if difftimeij < difftime
                    difftime=difftimeij;
                    if difftimeij < timethreshold
                        Image.TimeCorrIndex(itime)=jtime;
                    end
                end   
            end
        end
        
    end
    % take selection of data
    PointsR.x=Points.x;
    PointsR.y=flipud(Points.y);
    
    PointsR.minx=min(Points.x{Epochs.Index(1)},[],'all');
    PointsR.maxx=max(Points.x{Epochs.Index(1)},[],'all');
    PointsR.miny=min(Points.y{Epochs.Index(1)},[],'all');
    PointsR.maxy=max(Points.y{Epochs.Index(1)},[],'all');
    
    rangex=PointsR.maxx-PointsR.minx;
    rangey=PointsR.maxy-PointsR.miny;
    
    % add extra buffer
    PointsR.minx=PointsR.minx-0.05*rangex;
    PointsR.maxx=PointsR.maxx+0.05*rangex;
    PointsR.miny=PointsR.miny-0.05*rangey;
    PointsR.maxy=PointsR.maxy+0.05*rangey;
    
    CellsR.IndexingVertices=Cells.IndexingVertices;
    CellsR.Connectivity=Cells.Connectivity;
    
    % copy information
    ImageR.Im=Image.Scaled;
    ImageR.TimeCorrIndex=Image.TimeCorrIndex;
    
    ImageR.Dates=Image.Dates;
    % reference object
    ImageR.RI = imref2d(size(Image.Scaled{1}));
    if isfield(Image,'xRange')
        % pixel size
      ImageR.RI.XWorldLimits = Image.xRange;
      ImageR.RI.YWorldLimits = Image.yRange;
    elseif isfield(Image,'xvec')
      ImageR.RI.XWorldLimits = [Image.xvec(1) Image.xvec(end)];
      ImageR.RI.YWorldLimits = [Image.yvec(1) Image.yvec(end)];
    end
 
    
    if Op.ShowOriginalGrid
        PIVresultsR.x = PIVresults.x;
        PIVresultsR.y= PIVresults.y;
        PIVresultsR.mask = PIVresults.typevector_original;
    else
        PIVresultsR=[];
    end
    
    % y-axis direction
    Opts.UpwardYAxis= Op.UpwardYAxis;
    
    % make figure
    fig=figure;
    set(fig,'Position',[0.1 0.1 900 900])
    set(fig,'Name','Image','Toolbar','figure',...
        'NumberTitle','off')
    
    % Create an axes to plot in
    axes('Position',[.05 .05 .9 .95]);
    
    % sliders for epsilon and lambda
    slider1_handle=uicontrol(fig,'Style','slider','Max',nepochs,'Min',1,...
        'Value',1,'SliderStep',[1/(nepochs-1) 10/(nepochs-1)],...
        'Units','normalized','Position',[.02 .02 .14 .05]);
    uicontrol(fig,'Style','text','Units','normalized','Position',[.02 .07 .14 .04],...
        'String','Choose frame');
    
    % Set up callbacks
    vars=struct('slider1_handle',slider1_handle,'Points',PointsR,'Image',ImageR,'Cells',CellsR,'Epochs',Epochs,'PIVresults',PIVresultsR,'Op',Opts);
    set(slider1_handle,'Callback',{@slider1_callback,vars});
    % call plot function
    plotterfcn(vars)
else
    disp('skip, as no images have been supplied')
end

% End of main file

% Callback subfunctions to support UI actions
function slider1_callback(~,~,vars)
% Run slider1 which controls value of epsilon
plotterfcn(vars)

function plotterfcn(vars)
% get data
Points=vars.Points;
Epochs=vars.Epochs;
Cells=vars.Cells;
Image=vars.Image;
PIVresults=vars.PIVresults;
Op=vars.Op;

if ~isempty(PIVresults)
    ShowOriginalGrid=1;
else
    ShowOriginalGrid=0;
end

% get number
i=round(get(vars.slider1_handle,'Value'));
ii=Epochs.Index(i);
if isfield(Image,'TimeCorrIndex')
    iimage=Image.TimeCorrIndex(i);
else
    iimage=ii;
end
% Plots the image
% original vectors



% plot image
imshow(Image.Im{iimage},Image.RI)
if Op.UpwardYAxis
    set(gca,'YDir','normal')
else
    set(gca,'YDir','reverse') % because PIV results are in reverse y direction
end

hold on
% plot grid
if ShowOriginalGrid
    plot(PIVresults.x{ii}(PIVresults.mask{ii}==1),PIVresults.y{ii}(PIVresults.mask{ii}==1),'.k')
end

% plus one gives updated grid
vertices=[Points.x{ii+1}(Cells.IndexingVertices) Points.y{ii+1}(Cells.IndexingVertices)];
patch('vertices',vertices, 'faces',Cells.Connectivity,'FaceVertexCData',NaN(length(Cells.Connectivity),1),'FaceColor','flat','EdgeColor',ones(3,1)*0.5)
hold off
%title(strcat({'deformed grid, epoch '},{num2str(ii)}));
if isdatetime(Epochs.Time(ii))
            timestr=['grid:' datestr(Epochs.Time(ii)) ' image: ' datestr(Image.Dates(iimage))];
else
            timestr=strcat({' time: '},{num2str(Epochs.Time(ii),'%4.1f\n')},{' '},{Param.TimeUnit});
end
title(timestr)
% axes
axis equal
xlim([Points.minx Points.maxx]);
ylim([Points.miny Points.maxy]);





