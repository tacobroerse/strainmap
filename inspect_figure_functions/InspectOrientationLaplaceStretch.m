function InspectOrientationLaplaceStretch(Cells,Points,Epochs,Op,ix,iy)
% check the results for the Laplace Stretch under different orientations

% number of orientations
nrefs=200;

showevery=nrefs/20;
reforientation=linspace(0,1-1/(nrefs),nrefs)*pi;

ntimes=length(Cells.F);


% color map
Cmap=makecolormap('redyellowgreen',nrefs);
Cmap=makescientificcolormap('romaO',nrefs);
% rotation matrix

RotMat=@(theta) [cos(theta) -sin(theta) ;...
    sin(theta)  cos(theta)];


% loop on different reference frame orientations

Frot=cell(nrefs,ntimes);
U=cell(nrefs,ntimes);
Uinv=cell(nrefs,ntimes);
IncrU=cell(nrefs,ntimes);
R=cell(nrefs,ntimes);
a=zeros(nrefs,ntimes);
b=zeros(nrefs,ntimes);
gamma=zeros(nrefs,ntimes);
theta=zeros(nrefs,ntimes);


aincr=[];
bincr=[];
gammaincr=[];

for iref=1:nrefs
    for itime=1:ntimes
        
        % rotate to a different reference frame orientation
        refrot=reforientation(iref);
        Frot{iref,itime}=RotMat(refrot)*Cells.F{itime}(:,:,iy,ix)*RotMat(refrot)';
        
        % now do a QR decomposition in this reference orientation
        
        [U{iref,itime},R{iref,itime},a(iref,itime),b(iref,itime),gamma(iref,itime),theta(iref,itime),Uinv{iref,itime}]=QRdecompositionUpTriang(Frot{iref,itime});
        
        
        % make incremental distortion
        if itime>1
            % U(i)=Uincr(i)*U(i-1)
            IncrU{iref,itime-1}=U{iref,itime}*Uinv{iref,itime-1};
            aincr(iref,itime-1)=IncrU{iref,itime-1}(1,1);
            bincr(iref,itime-1)=IncrU{iref,itime-1}(2,2);
            gammaincr(iref,itime-1)=IncrU{iref,itime-1}(1,2)/aincr(iref,itime-1);
            
            thetaincr(iref,itime-1)=theta(iref,itime)-theta(iref,itime-1);
        end
        %             % first derivative
        %              aincr(iref,:)=diff(a(iref,:));
        %              bincr(iref,:)=diff(b(iref,:));
        %              gammaincr(iref,:)=diff(gamma(iref,:));
        
        
    end
    % second derivative
    aincr2(iref,:)=diff(aincr(iref,:));
    bincr2(iref,:)=diff(bincr(iref,:));
    gammaincr2(iref,:)=diff(gammaincr(iref,:));
    thetaincr2(iref,:)=diff(thetaincr(iref,:));
    
    
    % rms
    rmsdistortionchange(iref)=sqrt(mean((aincr2(iref,:)).^2+(bincr2(iref,:)).^2+(gammaincr2(iref,:)).^2));
    
    rmsrotationchange(iref)=sqrt(mean((thetaincr2(iref,:)*180/pi).^2));
    
    rmsrotationchange(iref)=sqrt(mean((thetaincr(iref,:)*180/pi).^2));
   
    
end

% find minimum value
[minrmsdist,imindist]=min(rmsdistortionchange);

[minrmsrot,iminrot]=min(rmsrotationchange);

%% make figures

titlestr=['point x:',num2str(ix),' y:',num2str(iy)];

LineWidth=1.5;
MarkerSize=6;
% figure cumulative distortion
figure
n=2;
m=2;
legendstr=[];
kref=0;
for iref=1:showevery:nrefs
    kref=kref+1;
    legendstr{kref}=[num2str(reforientation(iref)/pi) '\pi'];
end

subplot(n,m,1)
% plot extension 1 (a)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes],a(iref,:),'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
% plot optimal value

plot([1:ntimes],a(imindist,:),'Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
plot([1:ntimes],a(iminrot,:),':','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)

title('extension a')
legend(legendstr,'Location','Best')
subplot(n,m,2)
% plot shear (gamma)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes],gamma(iref,:),'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
% plot optimal value

plot([1:ntimes],gamma(imindist,:),'Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
plot([1:ntimes],gamma(iminrot,:),':','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)
title('shear \gamma')
subplot(n,m,3)
% plot extension 2 (b)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes],b(iref,:),'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
title('extension b')
% plot optimal value

plot([1:ntimes],b(imindist,:),'Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
plot([1:ntimes],b(iminrot,:),':','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)

subplot(n,m,4)
% plot rotation (theta)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes],theta(iref,:)*180/pi,'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
% plot optimal value

plot([1:ntimes],theta(imindist,:)*180/pi,'Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
plot([1:ntimes],theta(iminrot,:)*180/pi,':','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)
title('rotation angle \theta [deg]')

sgtitle(titlestr)

%% figure incremental distortion
figure
n=3;
m=2;
legendstr=[];
kref=0;
for iref=1:showevery:nrefs
    kref=kref+1;
    legendstr{kref}=[num2str(reforientation(iref)/pi) '\pi'];
end

subplot(n,m,1)

% plot extension 1 (a)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes-1],aincr(iref,:),'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
% plot optimal value

plot([1:ntimes-1],aincr(imindist,:),'Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
plot([1:ntimes-1],aincr(iminrot,:),':','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)

title('incremental extension a')
legend(legendstr,'Location','Best')
subplot(n,m,2)
% plot shear (gamma)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes-1],gammaincr(iref,:),'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
title('incrementalshear \gamma')
% plot optimal value

plot([1:ntimes-1],gammaincr(imindist,:),'Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
plot([1:ntimes-1],gammaincr(iminrot,:),':','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)

subplot(n,m,3)
% plot extension 2 (b)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes-1],bincr(iref,:),'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
title('incremental extension b')
% plot optimal value

plot([1:ntimes-1],bincr(imindist,:),'Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
plot([1:ntimes-1],bincr(iminrot,:)',':','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)
subplot(n,m,4)
% % plot rotation (theta)
%
hold on
grid on
% % loop on reference frames
for iref=1:showevery:nrefs
    
    
    
    plot([1:ntimes-1],thetaincr(iref,:)*180/pi,'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
title('incremental rotation angle [deg]')
% plot optimal value

plot([1:ntimes-1],thetaincr(imindist,:)*180/pi,'Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
plot([1:ntimes-1],thetaincr(iminrot,:)*180/pi,':','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)
subplot(n,m,5)
% % plot rotation (theta)
%
hold on
grid on
% % loop on reference frames
for iref=1:showevery:nrefs
    distortionchange=abs(aincr(iref,:)-1)+abs(bincr(iref,:)-1)+abs(gammaincr(iref,:));
    
    % distortionchange=(1/3*((aincr(iref,:)-1).^2+(bincr(iref,:)-1).^2+(gammaincr(iref,:)).^2));
    plot([1:ntimes-1],distortionchange,'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
title('distortion change')
% plot optimal value
%plot([1:ntimes-1],distortionchange,'Color',Cmap(iref,:),'LineWidth',LineWidth*2)
%plot([1:ntimes-1],distortionchange,'Color',Cmap(iref,:),'LineWidth',LineWidth*2)

sgtitle(titlestr)

%% figure second derivative distortion parameters
figure
n=3;
m=2;
legendstr=[];
for iref=1:nrefs
    legendstr{iref}=[num2str(reforientation(iref)/pi) '\pi'];
end

subplot(n,m,1)
% plot extension 1 (a)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes-2],aincr2(iref,:),'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
title('change in incremental extension a')

subplot(n,m,2)
% plot shear (gamma)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes-2],gammaincr2(iref,:),'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
title('change in incremental shear \gamma')
subplot(n,m,3)
% plot extension 2 (b)

hold on
grid on
% loop on reference frames
for iref=1:showevery:nrefs
    plot([1:ntimes-2],bincr2(iref,:),'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
title('change in incremental extension b')

subplot(n,m,4)
% % plot rotation (theta)
%
hold on
grid on
% % loop on reference frames
for iref=1:showevery:nrefs
    distortionchange=abs(aincr2(iref,:))+abs(bincr2(iref,:))+abs(gammaincr2(iref,:));
    
    % distortionchange=(1/3*((aincr(iref,:)-1).^2+(bincr(iref,:)-1).^2+(gammaincr(iref,:)).^2));
    plot([1:ntimes-2],distortionchange,'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
title('change in distortion ')

subplot(n,m,5)
hold on
grid on
for iref=1:nrefs
    
    
    % distortionchange=(1/3*((aincr(iref,:)-1).^2+(bincr(iref,:)-1).^2+(gammaincr(iref,:)).^2));
    plot(reforientation(iref)*180/pi,rmsdistortionchange(iref),'.','MarkerSize',MarkerSize,'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
%legend(legendstr,'Location','Best')
title('rms change in incremental laplace stretch')
plot(reforientation(imindist)*180/pi,rmsdistortionchange(imindist),'o','Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
plot(reforientation(iminrot)*180/pi,rmsdistortionchange(iminrot),'x','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)
xlabel('orientation [deg]')
xlim([0 180])
subplot(n,m,6)
hold on
grid on
for iref=1:nrefs
    
    
    % distortionchange=(1/3*((aincr(iref,:)-1).^2+(bincr(iref,:)-1).^2+(gammaincr(iref,:)).^2));
    plot(reforientation(iref)*180/pi,rmsrotationchange(iref),'.','MarkerSize',MarkerSize,'Color',Cmap(iref,:),'LineWidth',LineWidth)
end
plot(reforientation(iminrot)*180/pi,rmsrotationchange(iminrot),'o','Color',Cmap(iminrot,:),'LineWidth',LineWidth*2)
plot(reforientation(imindist)*180/pi,rmsrotationchange(imindist),'x','Color',Cmap(imindist,:),'LineWidth',LineWidth*2)
xlim([0 180])
xlabel('orientation [deg]')
title('rms change in incremental rotation [deg]')
%   subplot(n,m,6)
%  hold on
%  grid on
%  for iref=1:nrefs
%
%     distortion=(aincr(iref,:)).^2+(bincr(iref,:)).^2+(gammaincr(iref,:)).^2;
%    % distortionchange=(1/3*((aincr(iref,:)-1).^2+(bincr(iref,:)-1).^2+(gammaincr(iref,:)).^2));
%      plot(reforientation(iref)/pi,sqrt(mean(distortion)),'.','MarkerSize',10,'Color',Cmap(iref,:),'LineWidth',LineWidth)
% end
%  title('rms incremental laplace stretch')
sgtitle(titlestr)

end

