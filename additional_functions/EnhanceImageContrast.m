function [Image]=EnhanceImageContrast(Image,Param,Op)
% enhance the color of images by stretching the colorspace such that it
% more effectively fills the colorspace
% use a percentile of the color distribution to determine stretching
% 100 % results in the original colormap
% 75 % results in that half of the distribution is mapped to the full color
% range

if Op.EnhanceContrastImage
    if Param.ColorRangePerc <= 50
        error('Param.ColorRangePerc has to be larger than 50%')
        
    end
    if Param.ColorRangePerc > 100
        error('Param.ColorRangePerc has to be smaller or equal than 100')
        
    end
    nPixels=numel(Image.Scaled{1}(:,:,1));
    
    nImages=Image.nImages;
    % loop on images to determine mean distribution
    for i=1:nImages
        j=Image.EpochIndex(i);
        % take image
        I=Image.Scaled{j};
        % red
        [Rcounts,RbinLocations] = imhist(I(:,:,1));
        % green
        [Gcounts,GbinLocations] = imhist(I(:,:,2));
        % blue
        [Bcounts,BbinLocations] = imhist(I(:,:,3));
        % find normalisation value
        normval=sum(Rcounts);
        % normval = max([Rcounts ; Gcounts; Bcounts])
        % normalise and save
        Rcountnorm{i}=Rcounts/normval;
        Gcountnorm{i}=Gcounts/normval;
        Bcountnorm{i}=Bcounts/normval;
        
    end
    
    binLocations=RbinLocations;
    % now sum all normalised counts
    Rcountsum=Rcountnorm{1};
    Gcountsum=Gcountnorm{1};
    Bcountsum=Bcountnorm{1};
    for i=2:nImages
        Rcountsum=Rcountsum+Rcountnorm{i};
        Gcountsum=Gcountsum+Rcountnorm{i};
        Bcountsum=Bcountsum+Rcountnorm{i};
    end
    Rcountsum=Rcountsum/nImages;
    Gcountsum=Gcountsum/nImages;
    Bcountsum=Bcountsum/nImages;
    % now integrate sums
    nbins=length(binLocations);
    Rintegral=zeros(nbins,1);
    Gintegral=zeros(nbins,1);
    Bintegral=zeros(nbins,1);
    for j=1:nbins
        Rintegral(j)=sum(Rcountsum(1:j));
        Gintegral(j)=sum(Gcountsum(1:j));
        Bintegral(j)=sum(Bcountsum(1:j));
    end
    % find percentile range
    %  Rlowbound=prctile(Rintegral,100-Param.ColorRangePerc);
    %  Rhighbound=prctile(Rintegral,Param.ColorRangePerc);
    % Glowbound=prctile(Gintegral,100-Param.ColorRangePerc);
    % Ghighbound=prctile(Gintegral,Param.ColorRangePerc);
    %  Blowbound=prctile(Bintegral,100-Param.ColorRangePerc);
    %  Bhighbound=prctile(Bintegral,Param.ColorRangePerc);
    lowbound=(100-Param.ColorRangePerc)/100;
    highbound=(Param.ColorRangePerc)/100;
    
    [~,iRlow] = min(abs(Rintegral-lowbound));
    [~,iRhigh] = min(abs(Rintegral-highbound));
    
    [~,iGlow] = min(abs(Gintegral-lowbound));
    [~,iGhigh] = min(abs(Gintegral-highbound));
    
    [~,iBlow] = min(abs(Bintegral-lowbound));
    [~,iBhigh] = min(abs(Bintegral-highbound));
    
    fig=figure;
    subplot(3,1,1);plot(binLocations,Rcountsum,'r')
    
    
    hold on
    ylims=ylim;
    xlim([0 binLocations(end)])
    plot(iRlow*ones(2,1),ylims,':r')
    plot(iRhigh*ones(2,1),ylims,':r')
    text(10,0.9*ylims(2),strcat('bounds for :',num2str(Param.ColorRangePerc),' percentile'));
    
    xlabel('pixel intensity value')
    ylabel('fraction per bin')
    yyaxis right
    plot(binLocations,Rintegral,'-')
    ylabel('integrated fraction')
    ylim([0 1])
    title('red')
    
    subplot(3,1,2);plot(binLocations,Gcountsum,'g')
    hold on
    plot(iGlow*ones(2,1),ylims,':g')
    plot(iGhigh*ones(2,1),ylims,':g')
    title('green')
    xlabel('pixel intensity value')
    ylabel('fraction per bin')
    xlim([0 binLocations(end)])
    yyaxis right
    
    plot(binLocations,Gintegral,'-')
    ylim([0 1])
    ylabel('integrated fraction')
    
    subplot(3,1,3);plot(binLocations,Bcountsum,'b')
    hold on
    plot(iBlow*ones(2,1),ylims,':b')
    plot(iBhigh*ones(2,1),ylims,':b')
    title('blue')
    xlabel('pixel intensity value')
    ylabel('fraction per bin')
    xlim([0 binLocations(end)])
    xlim([0 binLocations(end)])
    yyaxis right
    plot(binLocations,Bintegral,'-')
    ylim([0 1])
    ylabel('integrated fraction')
    sgtitle('average color intensity for all images')
    
    if Op.SaveFigures
        if ~isfolder(Param.SaveDir)
            mkdir(Param.SaveDir)
            disp(strcat('making folder:',Param.SaveDir))
        end
        
        SaveFigName=strcat(Param.SaveDir,'/','color_intensity_mapping');
        savefig(fig,SaveFigName)
        if Op.SavePng
            % save as png
            print(fig,SaveFigName,'-dpng',Param.FigureResolution)
        end
    end
    %% now apply bounds to all images
    
    for i=1:Image.nImages
        j=Image.EpochIndex(i);
        I=Image.Scaled{j};
        % set up target histogram
        hgramR=[iRlow iRhigh]/nbins;
        hgramG=[iGlow iGhigh]/nbins;
        hgramB=[iBlow iBhigh]/nbins;
        hgram=[hgramR ; hgramG ; hgramB]';
        Image.EnhancedImage{j}=imadjust(I,hgram);
        
    end
    
    % show resulting image
    fig=figure;
    imshowpair(I,Image.EnhancedImage{Image.EpochIndex(nImages)},'montage')
    title(strcat('original (left) vs. enhanced image (right), color intensity percentile:',num2str(Param.ColorRangePerc) ))
    
    SaveFigName=strcat(Param.SaveDir,'/','change_in_image_color_intensity');
    savefig(fig,SaveFigName)
    if Op.SavePng
        % save as png
        print(fig,SaveFigName,'-dpng',Param.FigureResolution)
    end
else
    disp('skip image contrast enhancement')
    
end
end