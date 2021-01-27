function ViewPIVresults(PIVresults,Epochs,Param,Op,PIVcorrected)
%ViewPIVresults shows the displacements as vectors. It is possible to
%scroll through time.
%
% ViewPIVresults(PIVresults,Epochs,Param,Op)
%       PIVresults, structure with coordinates and displacement
%       data.
%       Epochs, time structure.
%       Op, structure with  options. See SetDefaults
%       Param, structure with parameters
%       Param.VectorScale=number magnifies the vectors in the plot
%  ViewPIVresults(PIVresults,Epochs,Param,Op,PIVcorrected)
%       If outlier detection has been performed, corrected PIV
%       displacements can be input with PIVcorrected
%
%
% STRAINMAP
%
% programmed with MATLAB version 2018b
% by Taco Broerse, 2020
% for more information on the theory:
% https://doi.org/10.31223/X5FS3H
%
if nargin == 4
    Outliers=0;
    nRowPlot=1;
elseif nargin == 5
    Outliers=1;
    nRowPlot=2;
else
    error('incorrect number of outputs')
end

% select displacements and apply mask
[PIV.u,PIV.v]=SelectDisplacements(PIVresults,Epochs,Op);

% image indices
EpochsIndex=Epochs.Index;
nepochs=length(EpochsIndex);

if nepochs < 2
    warning('too few epochs ,skipping')
    return
end


if Outliers
    PIV.correctedu=PIVcorrected.u;
    PIV.correctedv=PIVcorrected.v;
    PIV.Outlier=PIVcorrected.Outlier;
end



% coordinates
PIV.x=PIVresults.x;
PIV.y=PIVresults.y;


% steps in coordinates
dx=median(abs(diff(unique(PIVresults.x{EpochsIndex(1)}))));
dy=median(abs(diff(unique(PIVresults.y{EpochsIndex(1)}))));

factoredge=5;
PIV.minx = min(PIV.x{EpochsIndex(1)}(:))-factoredge*dx;
PIV.maxx = max(PIV.x{EpochsIndex(1)}(:))+factoredge*dx;

PIV.miny = min(PIV.y{EpochsIndex(1)}(:))-factoredge*dy;
PIV.maxy = max(PIV.y{EpochsIndex(1)}(:))+factoredge*dy;

% quiver scale
PIV.scale=Param.VectorScale;

% y-axis direction
PIV.UpwardYAxis= Op.UpwardYAxis;

PIV.dstep = sqrt(dx^2+dy^2);
% make figure
fig=figure;
set(fig,'Position',[0.1 0.1 900 900])
set(fig,'Name','Image','Toolbar','figure',...
    'NumberTitle','off')

if Outliers
    pos1=[0.1 0.55 0.8 0.35];
else
    pos1=[0.1 0.1 .8 .8];
end
    axSub(1)=subplot('Position',pos1);
    hold on

% sliders for sliding through time
slider1_handle=uicontrol(fig,'Style','slider','Max',nepochs,'Min',1,...
    'Value',1,'SliderStep',[1/(nepochs-1) 10/(nepochs-1)],...
    'Units','normalized','Position',[.02 .02 .4 .05]);
uicontrol(fig,'Style','text','Units','normalized','Position',[.02 0 .14 .04],...
    'String','Choose epoch','FontSize',14);



% Set up callbacks
vars=struct('slider1_handle',slider1_handle,'PIV',PIV,'EpochsIndex',EpochsIndex,'axSub',axSub);
set(slider1_handle,'Callback',{@slider1_callback,vars});
% call plot function
plotterfcn(vars)



if Outliers
    % create second (u) and third (v) plot to show time series

    pos2=[0.1 0.1 0.4 0.35];
    axSub(2)=subplot('Position',pos2);box on
    pos3=[0.55 0.1 0.4 0.35];
    axSub(3)=subplot('Position',pos3);box on


    
    % button for starting to select a point in the upper panel

    button_handle = uicontrol(fig,'style','push',...
        'units','normalized',...
        'position',[.5 .02 .2 .05],...
        'fontsize',14,...
        'string','select point');
    
    % Set up callbacks
    vars2=struct('button_handle',button_handle,'PIV',PIV,'EpochsIndex',EpochsIndex,'axSub',axSub,'Param',Param);
    set(button_handle,'Callback',{@button_callback,vars2});
    
    
    
end


% End of main file

% Callback subfunctions to support UI actions
    function slider1_callback(~,~,vars)
        % Run slider1
        plotterfcn(vars)
    end
    function plotterfcn(vars)
        % get data
        PIV=vars.PIV;
        Epochs=vars.EpochsIndex;
        
        
        
        if isfield(PIV,'correctedu')
            Outliers=1;
            Color1 = 'r';
            Color2 = 'k';
        else
            Outliers=0;
            Color1 = 'k';
        end
        
        % get number
        i=round(get(vars.slider1_handle,'Value'));
        ii=Epochs(i);
        % Plots the image
        % original vectors
        subplot(vars.axSub(1))
        cla(vars.axSub(1))
        quiver(PIV.x{ii},PIV.y{ii},PIV.u{ii}*PIV.scale,PIV.v{ii}*PIV.scale,0,'Color',Color1)
        if Outliers
            hold on
            % corrected vectors
            quiver(PIV.x{ii},PIV.y{ii},PIV.correctedu{ii}*PIV.scale,PIV.correctedv{ii}*PIV.scale,0,'Color',Color2)
            plot(PIV.x{ii}(find(squeeze(PIV.Outlier(ii,:,:)))),PIV.y{ii}(find(squeeze(PIV.Outlier(ii,:,:)))),'ro')
            hold off
             title(strcat({'displacements epoch '},{num2str(ii)},{' scaled by '},{num2str(PIV.scale)},{' outliers in red'}),'FontSize',Param.TitleFontSize,'FontWeight','Normal');

        else
             title(strcat({'displacements epoch '},{num2str(ii)},{' scaled by '},{num2str(PIV.scale)}),'FontSize',Param.TitleFontSize,'FontWeight','Normal');

        end
    
               % axes
        if PIV.UpwardYAxis==0
            set(gca,'YDir','reverse') % because PIV results are in reverse y direction
        end
        
        axis equal
        xlim([PIV.minx PIV.maxx]);
        ylim([PIV.miny PIV.maxy]);
    end

    function button_callback(~,~,vars2)
        % Run button
        plottertimeseriesfcn(vars2)
    end

    function plottertimeseriesfcn(vars)
        % function for selecting data in upper panel
        % time series will be shown in lower panel
         % get data
        PIV=vars.PIV;
        Epochs=vars.EpochsIndex;
        Param=vars.Param;
        % get colors
        [cmap]=makecolormap('qualitative2');
        hPoint=[];
        while 0<1
            subplot(vars.axSub(1))
            % allow for clicking multiple times
            disp('click on cell to inspect')
            disp('type 1 to zoom out')
            disp('type 2 to zoom in')
         
           % somehow ginput does not take keyboard input for the first try
                [x,y,b] = ginput(1);
          
                if isempty(b)
                    break; % get out of while loop
                elseif b==49 % (key 1)
                    % zoom out
                    ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
                    axis([x-width/2 x+width/2 y-height/2 y+height/2]);
                    zoom(1/2);
                elseif b==50 % key 2
                    % zoom in
                    ax = axis; width=ax(2)-ax(1); height=ax(4)-ax(3);
                    axis([x-width/2 x+width/2 y-height/2 y+height/2]);
                    zoom(2);
                else
                    % area clicked, proceed
                    % get corresponding point from coordinates
                    alldistances=sqrt((PIV.x{1}(:)-x).^2+(PIV.y{1}(:)-y).^2);
                    % find minimum distance
                    [mindist,indexpoint]=min(alldistances);
                    % get axis limits
                   
                    xmin=vars.axSub(1).XLim(1);xmax=vars.axSub(1).XLim(2);
                    ymin=vars.axSub(1).YLim(1);ymax=vars.axSub(1).YLim(2);
                    
                    % stop if there has been clicked outside axes
                    if x < xmin || x > xmax || y < ymin || y > ymax
                     
                        break
                    end
                    
                    % delete old point from plot
                    if ~isempty(hPoint)
                        delete(hPoint);
                    end
                    
                    % plot selected point
                    hold on
                    hPoint=plot(PIV.x{1}(indexpoint),PIV.y{1}(indexpoint),'o','MarkerSize',5);
                    
                    if isempty(indexpoint)
                        disp('did not find a corresponding cell, try again')
                    else
                        
                        % create vectors at single locations
                        clear uvec vvec
                        
                        for iTime=EpochsIndex
                            uvec(iTime)=PIV.u{iTime}(indexpoint);
                            vvec(iTime)=PIV.v{iTime}(indexpoint);
                            uveccor(iTime)=PIV.correctedu{iTime}(indexpoint);
                            vveccor(iTime)=PIV.correctedv{iTime}(indexpoint);
                        end
                        
                        % run FindReplaceOutliers again
                        [OutlierVecu,uinterp,umaddifffullvec,uInterpolated,uMovMedianThreshold]=FindReplaceOutliers(EpochsIndex,uvec,Param);
                        [OutlierVecv,vinterp,vmaddifffullvec,vInterpolated,vMovMedianThreshold]=FindReplaceOutliers(EpochsIndex,vvec,Param);
                       
                        OutlierVec = OutlierVecu & OutlierVecv;
                        % make separate plots for u and v
                        % make time series figure
                        
                        % plot u (displacement x direction)
                        subplot(vars.axSub(2))
                        cla(vars.axSub(2))
                        
                        hold on; box on
                        plot(EpochsIndex,uvec,':','Color',cmap(1,:),'LineWidth',2)
                        plot(EpochsIndex,uveccor,'-','Color',cmap(2,:),'LineWidth',1)
                        
                        % plot median absolute differences in moving window
                        % and threshold
                        plot(EpochsIndex,umaddifffullvec,'-','Color',cmap(3,:),'LineWidth',1)
                        plot(EpochsIndex,uMovMedianThreshold,'-','Color',cmap(4,:),'LineWidth',1)
                        legendstr={'u','u interpolated',strcat('MAD moving window:',num2str(Param.WidthWindow)),strcat('threshold, with scale:',num2str(Param.ThresholdFactor))};
                        
                        if ~isempty(find(OutlierVec))
                            plot(EpochsIndex(find(OutlierVec)),umaddifffullvec(find(OutlierVec)),'ro')
                            legendstr{5}='outlier';
                        end
                        if ~isempty(find(OutlierVecu))
                            plot(EpochsIndex(find(OutlierVecu)),umaddifffullvec(find(OutlierVecu)),'r.')
                            if isempty(find(OutlierVec))
                                legendstr{5}='outlier only on u';
                            else
                                legendstr{6}='outlier only on u';
                            end
                        end
                        

                        
                        legend(legendstr)
                        legend('boxoff')
                        legend('Location','SouthOutside','FontSize',Param.LabelFontSize)
                        xlabel('epoch')
                        ylabel('displacement')
                        [ix,iy]=ind2sub(size(PIV.u{iTime})',indexpoint);
                        title(['u timeseries for point ix ' num2str(ix) ' iy ' num2str(iy)],'FontSize',Param.TitleFontSize,'FontWeight','Normal')
                        axis tight
                        grid on
                        
                        
                        % plot v (displacement x direction)
                        subplot(vars.axSub(3))
                        cla(vars.axSub(3))
                        
                        hold on; box on
                        plot(EpochsIndex,vvec,':','Color',cmap(1,:),'LineWidth',2)
                        plot(EpochsIndex,vveccor,'-','Color',cmap(2,:),'LineWidth',1)
                        
                        % plot median absolute differences in moving window
                        % and threshold
                        plot(EpochsIndex,vmaddifffullvec,'-','Color',cmap(3,:),'LineWidth',1)
                        plot(EpochsIndex,vMovMedianThreshold,'-','Color',cmap(4,:),'LineWidth',1)
                        
                        legendstr={'v','v interpolated',strcat('MAD moving window:',num2str(Param.WidthWindow)),strcat('threshold, with scale:',num2str(Param.ThresholdFactor))};
                        
                        if ~isempty(find(OutlierVec))
                            plot(EpochsIndex(find(OutlierVec)),vmaddifffullvec(find(OutlierVec)),'ro')
                            legendstr{5}='outlier';
                        end
                        if ~isempty(find(OutlierVecv))
                            plot(EpochsIndex(find(OutlierVecv)),vmaddifffullvec(find(OutlierVecv)),'r.')
                            if isempty(find(OutlierVec))
                                legendstr{5}='outlier only on v';
                            else
                                legendstr{6}='outlier only on v';
                            end
                        end
                        
                        legend(legendstr)
                        legend('boxoff')
                        legend('Location','SouthOutside','FontSize',Param.LabelFontSize)
                        xlabel('epoch')
                        ylabel('displacement')
                        axis tight
                        grid on
                        
                        title(['v timeseries for point ix ' num2str(ix) ' iy ' num2str(iy)],'FontSize',Param.TitleFontSize,'FontWeight','Normal')
                        % axis limits
                        ymin2=vars.axSub(2).YLim(1);ymax2=vars.axSub(2).YLim(2);
                        ymin3=vars.axSub(3).YLim(1);ymax3=vars.axSub(3).YLim(2);
                        
                        ymax=max(abs([ymin2 ymax2 ymin3 ymax3]));
                        % use this to set the y limits
                        set(vars.axSub(2),'Ylim',[-ymax ymax]);
                        set(vars.axSub(3),'Ylim',[-ymax ymax]);
                        
                    end
                  
                end
         
  
        end
    end

end
