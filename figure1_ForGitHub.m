%% FIGURE 1 - Papier V1 - Basic info on the data
clear all
close all
catArray=['IJXX';'IJA2';'IJB2'];    % Cat Name
xSignalArray=[510,498,510];         % Size, x
ySignalArray=[330,330,348];         % Size, y

TF_values_original=[0.3800 0.8800 1.2400 1.5300 2.0000 2.7800 3.2600 3.8500 5.3700];
SF_values_original=[0.15 0.29 0.5 0.9 1.65];

TF_tested=1:9;
SF_tested=1:4;

nTF=length(TF_tested);
nSF=length(SF_tested);

TF_values=TF_values_original(TF_tested);
SF_values=SF_values_original(SF_tested);

catId=1; 

xSignal=xSignalArray(catId);
ySignal=ySignalArray(catId);


pathDir='~/Dropbox/Papier_TF_V1/figure_code/hannah_data/';
pathDirCat=[pathDir catArray(catId,:) '/'];

xWin=60:299;
yWin=50:227;

offset=0;
load([pathDirCat 'cache.mat']);
cache=cache(yWin,xWin);

%% Load Response Curves

disp(pathDirCat);
signal = zeros([9 5 200 ySignal*xSignal]);

% Load data 
for sf=SF_tested
    dataSFname=['dataSF' sprintf('%g',sf)];
    pathDirCatSF =[pathDirCat dataSFname '/'];
    load([pathDirCatSF dataSFname '_1-9_SFTF-smoothed.mat'],'signalS');
    signal(:,sf,:,:) = signalS;
end


%% Fig 1B, TF-amplitude map
Signal_Shaped = reshape(signal,9,5,200,ySignal,xSignal);

fig = figure('Position',[0,1000,1000,800]);
hold on 

for sf=1:length(SF_values)

    Signal_Shaped_Window(:,sf,:,:,:) = Signal_Shaped(:,sf,:,yWin,xWin);
    Signal_Mod_Pi_Window(:,sf,:,:,:) = Signal_Shaped_Window(:,sf,1:100,:,:)...
        +Signal_Shaped_Window(:,sf,101:200,:,:);
    Signal_Mean = squeeze(sum(Signal_Mod_Pi_Window,2)/length(SF_values));

    % Max response per TF, normalized pixel by pixel.
    AmplitudeMatrix=zeros(9,length(SF_values),length(yWin),length(xWin));
    MeanAmplitudeMatrix = zeros(9,length(yWin),length(xWin)); % Mean over SFs
    for i=1:length(yWin) 
        for j=1:length(xWin)
            % Max value on TF of the max on OR (i.e., global max)
            MaximalResponse(sf)=max(max(squeeze(Signal_Mod_Pi_Window(:,sf,:,i,j)))); 
            MeanMaximalResponse=max(max(Signal_Mean(:,:,i,j))); 
            % Min value on TF of the max on OR. 
            MinimalResponse(sf)=min(max(squeeze(Signal_Mod_Pi_Window(:,sf,:,i,j)),[],2));
            MeanMinimalResponse=min(max(Signal_Mean(:,:,i,j),[],2));
            for TF=TF_tested
                AmplitudeMatrix(TF,sf,i,j)=(max(Signal_Mod_Pi_Window(TF,sf,:,i,j))-MinimalResponse(sf))/(MaximalResponse(sf)-MinimalResponse(sf));
                MeanAmplitudeMatrix(TF,i,j)=(max(Signal_Mean(TF,:,i,j))-MeanMinimalResponse)/(MeanMaximalResponse-MeanMinimalResponse);
            end
        end
    end
    
    for TF=TF_tested
        ax(TF,sf) = subplot(5,9,TF+9*(sf-1));
        imagesc(squeeze(AmplitudeMatrix(TF,sf,:,:)));
        hold on
        caxis([0,1]);
        colormap(parula)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);  
        set(gca,'xaxisLocation','top');
        if sf == 1
           xlabel(ax(TF,1),sprintf('%.2f',TF_values(TF)),...
               'Interpreter','latex','FontSize',14);
        end
    end
    ylabel(ax(1,sf),sprintf('%.2f',SF_values(sf))...
        ,'Interpreter','latex','FontSize',14);
end

% Plotting average over SFs
for TF=TF_tested
    subplot(5,9,TF+36)
    imagesc(squeeze(MeanAmplitudeMatrix(TF,:,:)));
    hold on
    caxis([0,1]);
    colormap(parula)
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);  
    if TF==1
        ylabel('Average over SFs','Interpreter','latex','FontSize',14);
    end
end
%saveas(gca, 'Fig1B.fig');
c = colorbar;
c.Position = [0.93 0.168 0.022 0.7];
print -depsc epsFig1B

%% Fig 1D, TF-amplitude map 
figure; 
s = 1; % Control Parameter 
amp_TF = zeros(length(TF_tested),1);
for tf=TF_tested
    for i=1:s:length(yWin)
        for j=1:s:length(xWin)
            iy = min(i+s-1,length(yWin));
            ix = min(j+s-1,length(xWin));
            amp_spaced = (MeanMaximalResponse-MeanMinimalResponse)*MeanAmplitudeMatrix(tf,i:iy,j:ix);
            amp_TF(tf,end+1) = mean(amp_spaced(:));
        end
    end
end
amp_TF(:,1)=[];
errorbar(TF_tested,mean(amp_TF,2),std(amp_TF,[],2)/sqrt(sum(cache(:))));
grid on
set(gca,'XScale','log')
set(gca,'YTickLabel',[]);
xlabel('TF (Hz)','Interpreter','latex','FontSize',14);
ylabel('Mean Amplitude','Interpreter','latex','FontSize',14);
print -depsc epsFig1D

%% Fig 1C, On window shaped signal
% For chosen SF
[~,TF_Preferred]=max(AmplitudeMatrix(:,4,:,:),[],1);
%M_s = min(max(squeeze(Signal_Mod_Pi_Window(1,4,:,:))),1000);
%m_s = min(squeeze(Signal_Mod_Pi_Window(1,4,:,:)));
% Mean 
%[~,TF_Preferred]= max(MeanAmplitudeMatrix,[],1);
M_s = min(max(squeeze(Signal_Mean(1,:,:))),1000);
m_s = min(squeeze(Signal_Mean(1,:,:)));

amp = reshape(tanh((M_s-m_s)./(max(M_s(:))/0.7)),length(yWin),length(xWin));

N_pa = 7;
N_sh = 100;
[data intensityMap] = intensitymap(squeeze(TF_Preferred),amp,hsv(N_pa),N_sh);

figure;
imagesc(data')
colormap(intensityMap);
c = colorbar('location','southoutside');
caxis([0 N_pa*N_sh])
c.TickLabels = TF_values;    
c.TickLabelInterpreter = 'latex';
c.TickDirection = 'out';
set(gca,'FontSize',10,'LineWidth',1);
xlabel('Preferred TF (Hz)','Interpreter','latex','FontSize',14);
c.Ticks = [50 150 250 350 450 550 650];
c.TickLength = 0;
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
axis equal tight
%%
% Hand-picked cortical points
hold on
X= -60+[160 165 170 175];Y=[88 95 100 105]; 
s = scatter(X,Y,'filled','white','o');
s.LineWidth = 1.6;
s.MarkerEdgeColor = 'k';
print -depsc epsFig1C

%% Fig 1E
figure;
[a,b] = histc(TF_values(TF_Preferred(:)),TF_values);
bar(a/sum(a));
xticklabels(TF_values);
ylim auto
xlabel('Preferred TF (Hz)','Interpreter','latex','FontSize',14);
ylabel('Pixel Proportion','Interpreter','latex','FontSize',14);
print -depsc epsFig1E

%% Fig 1F
fig = figure; 
OR=[35 40 50 55]; % Arbitrary orientations 
for k=1:length(X)
    subplot(4,1,k);
    plot(TF_values,squeeze(mean(Signal_Shaped_Window(:,:,OR(k),Y(k),X(k)),2)),'.-','color',[0, 0.5, 0]); % Orientation picked arbitrarily
    hold on 
    plot(TF_values,squeeze(mean(Signal_Shaped_Window(:,:,OR(k)+100,Y(k),X(k)),2)),'.-b'); 
    set(gca, ...
      'Box'         , 'off'     , ...
      'TickLength'  , [0.02 0]   , ...
      'TickDir'     , 'in',...
      'XMinorTick'  , 'off'   , ...
      'YMinorTick'  , 'off'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'YTick'       , [], ...
      'XTick'       , [0 2 4], ...
      'LineWidth'   , 1         );
    ylim([5 30]);
    xlim([0 4]);
    yticklabels({'','',''});
    if k == 2
        ylabel('Response Amplitude','POSITION',[-.15 7],'Interpreter','latex','FontSize',14);
    end
    if k == 1
        legend('Chosen Direction', 'Chosen Direction $+ \pi$','Interpreter','latex');
    end
end
xlabel('TF (Hz)','Interpreter','latex','FontSize',14);
print -depsc epsFig1F

%% Fig 1A - Calcul de la carte d'orientation
OR_Maps=zeros(length(TF_tested),length(SF_tested),length(yWin),length(xWin));

for sf=SF_tested
    for TF=TF_tested
        for i=1:length(yWin)
            for j=1:length(xWin)
                if cache(i,j) 
                    S=Signal_Mod_Pi_Window(TF,sf,:,i,j);
                    fftSignal= fft(squeeze(S-mean(S)));
                    phase=angle(fftSignal(2));
                    if phase<0
                        phase=phase+2*pi;
                    end
                    OR_Maps(TF,sf,i,j)=phase/2;eli
                end
            end
        end
    end
end

figure;
for sf=SF_tested
    for TF=TF_tested
        ax(TF,sf) = subplot(length(SF_tested),9,TF+9*(sf-1));
        imagesc(squeeze(OR_Maps(TF,sf,:,:))');    
        colormap(hsv)
        set(gca,'YTickLabel',[]);
        set(gca,'XTickLabel',[]);
        set(gca,'xaxisLocation','top');
        if sf==1
            xlabel(sprintf('%.2f',TF_values(TF)),'Interpreter','latex','FontSize',14);
        end
    end
    ylabel(ax(1,sf),sprintf('%.2f',SF_values(sf)),'Interpreter','latex','FontSize',14);
end
print -depsc epsFig1A

%% 
Pix1x=50;
Pix1y=70;
Pix2x=20;
Pix2y=30;

figure();plot(Signal_Shaped_Window(:,2,20,Pix1x,Pix1y));
hold on
plot(Signal_Shaped_Window(:,2,20,Pix2x,Pix2y));

figure();plot(Signal_Shaped_Window(2,:,20,Pix1x,Pix1y));
hold on
plot(Signal_Shaped_Window(2,:,20,Pix2x,Pix2y));

figure();plot(squeeze(Signal_Shaped_Window(2,2,:,Pix1x,Pix1y)));
hold on
plot(squeeze(Signal_Shaped_Window(2,2,:,Pix2x,Pix2y)));
%%
figure()
kfig=0;
M=max(max(squeeze(Signal_Shaped_Window(:,:,:,Pix2x,Pix2y))),max(squeeze(Signal_Shaped_Window(:,:,:,Pix2x,Pix2y))));
for TF=1:9
    for SF=1:4
        kfig=kfig+1;
        subplot(9,4,kfig) 
        hold on
        plot(squeeze(Signal_Shaped_Window(TF,SF,:,Pix2x,Pix2y))/M);
        plot(squeeze(Signal_Shaped_Window(TF,SF,:,Pix1x,Pix1y))/M);
        end
end
%%
% Needs blocks 1,2,3 - Done with Isaiah and Jonathan on 5/7
figure();
ax1=subplot(3,1,1);
imagesc(squeeze(Signal_Shaped_Window(2,2,50,:,:)))
ax2=subplot(3,1,2);
imagesc(squeeze(Signal_Shaped_Window(2,2,100,:,:)))
ax2=subplot(3,1,3);
imagesc(squeeze(Signal_Shaped_Window(6,2,50,:,:)))
ax3=subplot(3,1,3);
imagesc(squeeze(Signal_Shaped_Window(6,2,50,:,:)))
linkaxes([ax1,ax2,ax3],'xy')

%%
Ref_SF=2;
Ref_TF=2;
Ref_OR=100;

MM_TF=squeeze( Signal_Shaped_Window(:     , Ref_SF ,Ref_OR ,:,:));
MM_DIR=squeeze(Signal_Shaped_Window(Ref_TF, Ref_SF ,  :    ,:,:));
MM_SF=squeeze( Signal_Shaped_Window(Ref_TF,  :     ,Ref_OR ,:,:));

MM_TF_Normalized=(MM_TF-min(MM_TF(:)))/(max(MM_TF(:))-min(MM_TF(:)));
MM_SF_Normalized=(MM_SF-min(MM_SF(:)))/(max(MM_SF(:))-min(MM_SF(:)));
MM_DIR_Normalized=(MM_DIR-min(MM_DIR(:)))/(max(MM_DIR(:))-min(MM_DIR(:)));

STD_TF =squeeze(std(MM_TF,[],1));
STD_DIR=squeeze(std(MM_DIR,[],1));
STD_SF =squeeze(std(MM_SF,[],1));

STD_TF_Normalized=squeeze(std(MM_TF_Normalized,[],1));
STD_DIR_Normalized=squeeze(std(MM_DIR_Normalized,[],1));
STD_SF_Normalized =squeeze(std(MM_SF_Normalized,[],1));

% Mean_STDs_TFs=squeeze(mean(squeeze(STD_TF),1));
% Mean_STDs_DIRss=squeeze(mean(squeeze(STD_DIR),1));
% Mean_STDs_SF=squeeze(mean(squeeze(STD_SF),1));

figure;pcolor(STD_TF );shading flat;colorbar()
figure;pcolor(STD_DIR);shading flat;colorbar()
figure;pcolor(STD_SF );shading flat;colorbar()


figure()
axes()
hold on

boxchart(repmat(1,size(STD_TF(:),1),1),STD_TF(:))
boxchart(repmat(3,size(STD_TF(:),1),1),STD_DIR(:))
boxchart(repmat(2,size(STD_TF(:),1),1),STD_SF(:))


figure()
axes()
hold on

boxchart(repmat(1,size(STD_TF_Normalized(:),1),1),STD_TF_Normalized(:))
boxchart(repmat(3,size(STD_TF(:),1),1),STD_DIR_Normalized(:))
boxchart(repmat(2,size(STD_TF(:),1),1),STD_SF_Normalized(:))

%%
 p = randperm(length(STD_TF(:)));
 p = p(1:100)
[h,p_tf_dir]=ttest2(STD_TF(p)  ,STD_DIR(p))
[h,p_tf_sf] =ttest2(STD_TF(p)  , STD_SF(p))
[h,p_sf_dir]=ttest2(STD_SF(p)  ,STD_DIR(p))

[h,p_tf_dir_norm]=ttest2(STD_TF_Normalized(p)  ,STD_DIR_Normalized(p))
[h,p_tf_sf_norm] =ttest2(STD_TF_Normalized(p)  , STD_SF_Normalized(p))
[h,p_sf_dir_norm]=ttest2(STD_SF_Normalized(p)  ,STD_DIR_Normalized(p))

%%
[h,p_tf_dir]=ztest(STD_TF(:)-STD_DIR(:),0,std(STD_TF(:)-STD_DIR(:)))
[h,p_tf_sf] =ztest(STD_TF(:)-STD_SF(:),0,std(STD_TF(:)-STD_SF(:)))
[h,p_sf_dir]=ztest(STD_SF(:)-STD_DIR(:),0,std(STD_SF(:)-STD_DIR(:)))

%%
cmin=min(min(min(squeeze(Signal_Shaped_Window(:,2,25,:,:)))));
cmax=max(max(max(squeeze(Signal_Shaped_Window(:,2,25,:,:)))));

figure;

subplot(5,3,1)
pcolor(squeeze(Signal_Shaped_Window(1,2,25,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[])
set(gca,'YTick',[])
axis image
shading interp
subplot(5,3,4)
pcolor(squeeze(Signal_Shaped_Window(3,2,25,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[])
set(gca,'YTick',[])
axis image
shading interp
subplot(5,3,7)
pcolor(squeeze(Signal_Shaped_Window(5,2,25,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,10)
pcolor(squeeze(Signal_Shaped_Window(7,2,25,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,13)
pcolor(squeeze(Signal_Shaped_Window(9,2,25,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp

maxima=zeros(5,1)
subplot(5,3,2)
pcolor(squeeze(Signal_Shaped_Window(1,2,25,:,:)))
maxima(1)=max(max(squeeze(Signal_Shaped_Window(1,2,25,:,:))))
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,5)
pcolor(squeeze(Signal_Shaped_Window(3,2,25,:,:)))
maxima(2)=max(max(squeeze(Signal_Shaped_Window(3,2,25,:,:))))
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,8)
pcolor(squeeze(Signal_Shaped_Window(5,2,25,:,:)))
maxima(3)=max(max(squeeze(Signal_Shaped_Window(5,2,25,:,:))))
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,11)
pcolor(squeeze(Signal_Shaped_Window(7,2,25,:,:)))
maxima(4)=max(max(squeeze(Signal_Shaped_Window(7,2,25,:,:))))
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,14)
pcolor(squeeze(Signal_Shaped_Window(9,2,25,:,:)))
maxima(5)=max(max(squeeze(Signal_Shaped_Window(9,2,25,:,:))))
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp

cmin=min(min(min(squeeze(Signal_Shaped_Window(5,2,:,:,:)))));
cmax=max(max(max(squeeze(Signal_Shaped_Window(5,2,:,:,:)))));

maxima2=zeros(5,1)
subplot(5,3,3)
pcolor(squeeze(Signal_Shaped_Window(5,2,10,:,:)))
maxima2(1)=max(max(squeeze(Signal_Shaped_Window(5,2,10,:,:))))

% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,6)
pcolor(squeeze(Signal_Shaped_Window(5,2,60,:,:)))
maxima2(2)=max(max(squeeze(Signal_Shaped_Window(5,2,60,:,:))))
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,9)
pcolor(squeeze(Signal_Shaped_Window(5,2,110,:,:)))
maxima2(3)=max(max(squeeze(Signal_Shaped_Window(5,2,110,:,:))))
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,12)
pcolor(squeeze(Signal_Shaped_Window(5,2,160,:,:)))
maxima2(4)=max(max(squeeze(Signal_Shaped_Window(5,2,160,:,:))))
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,3,15)
pcolor(squeeze(Signal_Shaped_Window(5,2,200,:,:)))
maxima2(5)=max(max(squeeze(Signal_Shaped_Window(5,2,200,:,:))))
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp

figure;
subplot(5,1,1)
pcolor(squeeze(Signal_Shaped_Window(5,2,10,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,1,2)
pcolor(squeeze(Signal_Shaped_Window(5,2,60,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,1,3)
pcolor(squeeze(Signal_Shaped_Window(5,2,110,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,1,4)
pcolor(squeeze(Signal_Shaped_Window(5,2,160,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(5,1,5)
pcolor(squeeze(Signal_Shaped_Window(5,2,200,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
%%
maxima3=zeros(5,1)

figure;
% subplot(4,1,1)
figure;
pcolor(squeeze(Signal_Shaped_Window(2,1,25,:,:)))
maxima3(1)=max(max(squeeze(Signal_Shaped_Window(2,1,25,:,:))))
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
% subplot(4,1,2)
figure;
pcolor(squeeze(Signal_Shaped_Window(2,2,25,:,:)))
maxima3(2)=max(max(squeeze(Signal_Shaped_Window(2,2,25,:,:))))
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
% subplot(4,1,3)
figure;
pcolor(squeeze(Signal_Shaped_Window(2,3,25,:,:)))
maxima3(3)=max(max(squeeze(Signal_Shaped_Window(2,3,25,:,:))))
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
% subplot(4,1,4)
figure;
pcolor(squeeze(Signal_Shaped_Window(2,4,25,:,:)))
maxima3(4)=max(max(squeeze(Signal_Shaped_Window(2,4,25,:,:))))
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp


cmin=min(min(min(squeeze(Signal_Shaped_Window(2,:,25,:,:)))));
cmax=max(max(max(squeeze(Signal_Shaped_Window(2,:,25,:,:)))));
%%
figure;
subplot(4,1,1)
pcolor(squeeze(Signal_Shaped_Window(2,1,25,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(4,1,2)
pcolor(squeeze(Signal_Shaped_Window(2,2,25,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(4,1,3)
pcolor(squeeze(Signal_Shaped_Window(2,3,25,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
subplot(4,1,4)
pcolor(squeeze(Signal_Shaped_Window(2,4,25,:,:)))
caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp

%%
figure;pcolor(squeeze(Signal_Shaped_Window(2,4,25,:,:)))
% caxis([cmin cmax])
% set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
hold on;scatter(X,Y,'filled','white','o');

% subplot(5,1,5)
% pcolor(squeeze(Signal_Shaped_Window(2,5,25,:,:)))
% axis image
% shading interp
%%
X= -50+[160 165 170 175];Y=[88 95 100 105]; 

TF=1
OR=zeros(1,4)
for k=1:4
signa=squeeze(mean(Signal_Shaped_Window(TF,:,:,Y(k),X(k)),2));
figure;plot(signa)
[~,OR_pref]=max(signa);
OR(k)=OR_pref;
end

%%
fig = figure; 
% OR=[74 67 23 22]; % Arbitrary orientations 
for k=1:length(X)
    subplot(4,1,k);
    plot(TF_values,squeeze(mean(Signal_Shaped_Window(:,:,OR(k),Y(k),X(k)),2)),'.-','color',[0, 0.5, 0]); % Orientation picked arbitrarily
    hold on 
    plot(TF_values,squeeze(mean(Signal_Shaped_Window(:,:,OR(k)+100,Y(k),X(k)),2)),'.-b'); 
    plot(TF_values,squeeze(mean(Signal_Shaped_Window(:,:,OR(k)+20,Y(k),X(k)),2)),'.-','color',[0.9290 0.6940 0.1250]); 
    set(gca, ...
      'Box'         , 'off'     , ...
      'TickLength'  , [0.02 0]   , ...
      'TickDir'     , 'in',...
      'XMinorTick'  , 'off'   , ...
      'YMinorTick'  , 'off'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'YTick'       , [], ...
      'XTick'       , [0 2 4], ...
      'LineWidth'   , 1         );
%     ylim([11 520]);
%     xlim([0 6]);
    yticklabels({'','',''});
    if k == 2
        ylabel('Response Amplitude','POSITION',[-.15 7],'Interpreter','latex','FontSize',14);
    end
    if k == 1
        legend('Preferred Direction', 'Preferred Direction $+ \pi$', 'Preferred Direction $+ \frac{\pi}{2}$','Interpreter','latex');
    end
end
xlabel('TF (Hz)','Interpreter','latex','FontSize',14);
print -depsc epsFig1F

%%
figure;

subplot(5,1,1)
M1=squeeze(Signal_Shaped_Window(1,2,50,:,:));
M1=(M1-min(M1(:)))/(max(M1(:))-min(M1(:)));
imagesc(M1)
colormap(bluewhitered)
caxis([-0.1,0.1])
% caxis([cmin cmax])
set(gca,'XTick',[])
set(gca,'YTick',[])
axis image
shading interp
% subplot(5,1,2)

figure;
M2=squeeze(Signal_Shaped_Window(3,2,50,:,:));
M2=(M2-min(M2(:)))/(max(M2(:))-min(M2(:)));
pcolor(M2-M1)
caxis([-0.1,0.1])
colormap(bluewhitered)
% caxis([cmin cmax])
set(gca,'XTick',[])
set(gca,'YTick',[])
axis image
shading interp


% subplot(5,1,3)
figure()
M2=squeeze(Signal_Shaped_Window(5,2,50,:,:));
M2=(M2-min(M2(:)))/(max(M2(:))-min(M2(:)));
pcolor(M2-M1)
caxis([-0.1,0.1])
colormap(bluewhitered)
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp


% subplot(5,1,4)
figure;
M2=squeeze(Signal_Shaped_Window(7,2,50,:,:));
M2=(M2-min(M2(:)))/(max(M2(:))-min(M2(:)));
pcolor(M2-M1)
caxis([-0.1,0.1])
colormap(bluewhitered)
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp
% subplot(5,1,5)


figure
M2=squeeze(Signal_Shaped_Window(9,2,50,:,:));
M2=(M2-min(M2(:)))/(max(M2(:))-min(M2(:)));
pcolor(M2-M1)
caxis([-0.1,0.1])
colormap(bluewhitered)
% caxis([cmin cmax])
set(gca,'XTick',[]); set(gca,'YTick',[]); axis image;
shading interp

%%

figure;
OR=25
M1=squeeze(Signal_Shaped_Window(1,2,OR,:,:));
M1=(M1-min(M1(:)))/(max(M1(:))-min(M1(:)));

for i=2:9
%     subplot(4,2,i-1)
figure;
M2=squeeze(Signal_Shaped_Window(i,2,OR,:,:));
M2=(M2-min(M2(:)))/(max(M2(:))-min(M2(:)));
pcolor(M2-M1)
caxis([-0.1,0.1])
colormap(bluewhitered)
% caxis([cmin cmax])
set(gca,'XTick',[])
set(gca,'YTick',[])
axis image
shading interp
end