%% File with fracture workflow combined
%Hannah's combination of Charlie's different files to make a single set of
%codeblocks to track pinwheels, fractures, calculate fracture stats, and
%plot the fractures for the paper.

%% INITIALIZE CONSTANTS (Always run this before running any other section)
clear all
close all
global cur_tf;
global path1;
global path2;
max_tf = 8; %created a maximum TF variable to use throughout program
cur_tf = 1; %created a current TF variable to use throughout program
max_sf = 4; %created a maximum SF variable to use throughout program
cur_sf = 1; %created a current SF variable to use throughout program
catId = 1; %Select your cat!

%root_path= '/Users/hannahgermaine/Desktop/dataSF/'; %location of dataSF folder
root_path = '../figure_code/hannah_data/';
catArray = ['IJXX';'IJA2';'IJB2'];
xSignalArray = [510,498,510];  %510, 498, 510
ySignalArray = [330,330,348];  %330, 330, 348
xWinStartArray = [60,50,210];
yWinStartArray = [50,100,90];
xWinEndArray = [299,250,330];
yWinEndArray = [227,250,220];
cat = catArray(catId,:);
xSignal = xSignalArray(catId);
ySignal = ySignalArray(catId);
path1 = string(strcat(root_path,cat,'/dataSF',string(cur_sf),'/')); %data folder 1
path2 = string(strcat(root_path,cat,'/')); %data folder 2
TFvalues = [0.3800 0.8800 1.2400 1.5300 2.0000 2.7800 3.2600 3.8500 5.3700];
SFvalues=[0.15 0.29 0.5 0.9 1.65];

TF_tested=1:max_tf;
SF_tested=1:max_sf;
DIR_tested=round(10:10:200);
DIR_values = DIR_tested/200*2*pi;
all_DIR_values = (1:200)/200*2*pi;
nTF=length(TF_tested);
nSF=length(SF_tested);
TF_values=TFvalues(TF_tested);
SF_values=SFvalues(SF_tested);
Signalpx = xSignal*ySignal;
xWin=xWinStartArray(catId):xWinEndArray(catId);
yWin=yWinStartArray(catId):yWinEndArray(catId);
Windowpx = length(yWin)*length(xWin);

m_px = 12.3e-6; %meters per pixel
x_um = m_px*length(xWin);
y_um = m_px*length(yWin);

%% PINWHEEL position tracking

%%To track properly, set the cur_tf variable above to 1. From there, the 
%%program automatically runs through the different TF values, and allows
%%you to track the PW locations, save a file, and then continue at the next
%%TF until you've completed through max_tf.
global pwTrack;
global pw_loc;
global ORdenoised;
global DIRall;
screenw = 1000;
screenh = 2000;

cat_data_dir = ['~/Dropbox/Papier_TF_V1/figure_code/hannah_data/', cat, '/'];

load(strcat(cat_data_dir,'dataSF',string(cur_sf),'/DIRall_test.mat'));

pwTrack = struct;
pwTrack(1).coord = NaN*ones(max_tf,2);
%%

pwTrack = trackPW(DIRall, pwTrack, cur_tf, max_tf, xSignal, ySignal)

%%
% increase cur_tf until = max_tf and keep re-running the above block to
%finish filling in.
if cur_tf < max_tf
    cur_tf = cur_tf + 1
else
    disp('Move on to next code block')
end

%%
%Clear first NaN placeholder
if pwTrack(1).coord(1,1) == NaN
    pwTrack(1) = [];
end
%Number all pinwheels
for i = 1:length(pwTrack)
    pwTrack(i).pwID = i;
end    
%Save
name = strcat('pwTrack_tf',string(max_tf),'.mat');
save(strcat(path1,'/',name), 'pwTrack'); 
%% TRACK FRACTURES (see PW TRACKING first ^ )

tf = cur_tf;
name = strcat('pwTrack_tf',string(max_tf),'.mat');
load(strcat(path1,'/',name), 'pwTrack');
global pwTrack; %#ok<REDEFGG> %in case PW Tracking not run first
global fracopen;
global fracclose;
global pw_loc;
global DIRall;
global cur_tf;
cat_data_dir = ['/Users/hannahgermaine/Documents/PhD/Past_Research/new_code/', cat, '/'];
load(strcat(cat_data_dir,'dataSF',string(cur_sf),'/DIRall_test.mat'));

fracopen=struct();
fracclose=struct();
%% 

trackFrac(cat,pwTrack,path1,DIRall,cur_tf,max_tf,xSignal,ySignal,fracopen,fracclose)

%%
%Save fractures
save(strcat(path1,'/fracopen_tf',sprintf('%.1g',cur_tf),'.mat'),'fracopen');
save(strcat(path1,'/fracclose_tf',sprintf('%.1g',cur_tf),'.mat'),'fracclose');
%When tracking, if you have a corner, because the derivative there is not
%smooth, try to track two points as close as possible to the corner on
%either side so the derivative is found at each point.
%%
% increase cur_tf until = max_tf and keep re-running the above block to
%finish filling in.
if cur_tf < max_tf
    cur_tf = cur_tf + 1
else
    disp('Move on to next code block')
end

%% Remove empty rows

for cur_tf = 1:max_tf
    name = strcat('pwTrack_tf',string(max_tf),'.mat');
    load(strcat(path1,'fracopen_tf',sprintf('%.1g',cur_tf),'.mat'),'fracopen');
    load(strcat(path1,'fracclose_tf',sprintf('%.1g',cur_tf),'.mat'),'fracclose');
    if isempty(fracclose(1).path)
        fracclose(1) = [];
    end
    if isempty(fracopen(1).path)
        fracopen(1) = [];
    end
    save(strcat(path1,'/fracopen_tf',sprintf('%.1g',cur_tf),'.mat'),'fracopen');
    save(strcat(path1,'/fracclose_tf',sprintf('%.1g',cur_tf),'.mat'),'fracclose');
end    

%% CREATE fS by JOINING ALL FRACTURE TF DATA

global pwTrack;
name = strcat('pwTrack_tf',string(max_tf),'.mat');
load(strcat(path1,name), 'pwTrack');

%Store in an array the pinwheel locations, the pinwheel number, and
%the temporal frequency
pw=zeros(1,4);
for i=1:length(pwTrack)
    for tf=1:length(pwTrack(i).coord)
        pw(end+1,1:2)=pwTrack(i).coord(tf,:); %x and y coordinates in first 2 columns
        pw(end,3)=i; %number of pinwheel in 3rd column
        pw(end,4)=tf; %temporal frequency of pinwheel location in 4th column
    end
end
pw(1,:)=[]; %remove first empty row


%Additional cache on the boundary of the image
A=zeros(330,510);
for i=1:size(A,1)
    A(i,:)=i;
end
cache=A>30;
[y, x] = find(cache);
a=[y,x];


%Create a structure to store fractures by temporal frequency
fS=struct();
pwPair=zeros(1,5);
for tf=1:max_tf
    load(strcat(path1,'fracopen_tf',sprintf('%.1g',tf),'.mat'),'fracopen');
    load(strcat(path1,'fracclose_tf',sprintf('%.1g',tf),'.mat'),'fracclose');
    fS(tf).frac=struct();
    for i=1:length(fracopen)
        fS(tf).frac(end+1).type='open';
        fS(tf).frac(end).nbPWs=0;
        pwCoord=[fracopen(i).pw.coord]; %concatenate all pinwheel coordinates into an array - len 4 means 2 pw, len 2 means 1 pw
        sFrac=size(fracopen(i).path,1); %how many points are in the traced path
        if length(pwCoord)>2 %if there are two pinwheels in the fracture
            if ~isempty(find(x==round(pwCoord(1)) & y==round(pwCoord(2)))) %Check if the first pinwheel is in the cache
                idPwTrack1=min(find(round(pw(:,1),4)==round(pwCoord(1),4))); %#ok<*MXFND> %Find the first id of the pinwheel
            end
            if ~isempty(find(x==round(pwCoord(3)) & y==round(pwCoord(4)))) %Check if the second pinwheel is in the cache
                idPwTrack2=min(find(round(pw(:,1),4)==round(pwCoord(3),4))); %Find the first id of the pinwheel
            end
            if exist('idPwTrack1','var') == 1 && exist('idPwTrack2','var') == 1 %if both are in the cache
                if ~isempty(idPwTrack1) && ~isempty(idPwTrack2)
                    pwPair(end+1,1)=pw(idPwTrack1,3); %number of the first pinwheel
                    pwPair(end,2)=pw(idPwTrack2,3); %number of the second pinwheel
                    pwPair(end,3)=pw(idPwTrack1,4); %first marked temporal frequency of the first pinwheel
                    pwPair(end,4)=pw(idPwTrack2,4); %first marked temporal frequency of the second pinwheel
                    pwPair(end,5)=tf; %current temporal frequency
                end
                distPWfrac=(pwCoord([1 3])-fracopen(i).path(1,1)).^2+(pwCoord([2 4])-fracopen(i).path(1,2)).^2;
                idPw=find(distPWfrac==min(distPWfrac),1); %find the closest pinwheel to the start of the fracture
                idPw2=mod(idPw,2)+1; %Calculate the id of the other pinwheel
                if ~isempty(find(x==round(pwCoord(1)) & y==round(pwCoord(2)))) %if the first pinwheel is in the cache    
                    fS(tf).frac(end).path(1,:)=pwCoord([2*idPw-1 2*idPw]); %store location of first pinwheel as beginning of path
                    fS(tf).frac(end).pw(1,:)=pwCoord([2*idPw-1 2*idPw]); %store location of first pinwheel
                end
                %Organize in fS the fractures
                if size(fracopen(i).path(:,:),1) > 1
                    for l=1:length(fracopen(i).path(:,:))
                        %x and y below are from the cache. We are checking if it is
                        %already in the cache or not.
                        if ~isempty(find(round(fracopen(i).path(l,1))==x & round(fracopen(i).path(l,2))==y)) %if in the cache
                            if ~isfield(fS(tf).frac(end),'path') %if the field of 'path' doesn't exist
                                fS(tf).frac(end).path(1,:)=fracopen(i).path(l,:); %set the first value of the path to match fracture path start
                            else %if the field of 'path' does exist
                                fS(tf).frac(end).path(end+1,:)=fracopen(i).path(l,:); %add to the field the values of the path
                            end
                        %The next portion appears unnecessary: the path is already
                        %marked open earlier and if a point is not in the cache, it
                        %simply isn't added to the path. This would overwrite a
                        %path with an empty one if any point is outside the cache,
                        %rather than saving only those in the cache.
        %                elseif isfield(fS(tf).frac(end),'path')
        %                     if ~isempty(fS(tf).frac(end).path) %if the path is not empty, and you go outside the cache, the fracture is open and you end
        %                         fS(tf).frac(end).type='open';
        %                         fS(tf).frac(end).path=[];
        %                     end
                        end
                    end
                end
                if ~isempty(find(x==round(pwCoord(3)) & y==round(pwCoord(4)))) %if the second pinwheel is also in the cache
                    fS(tf).frac(end).path(end+1,:)=pwCoord([2*idPw2-1 2*idPw2]); %add to the path the last point as that of the pinwheel
                    fS(tf).frac(end).pw(2,:)=pwCoord([2*idPw2-1 2*idPw2]); %add to the pinwheel collection the point of the second pinwheel
                end
                fS(tf).frac(end).distPWs=sqrt((pwCoord(1)-pwCoord(3)).^2+(pwCoord(2)-pwCoord(4)).^2); %calculate the pythagorean distance between pinwheels and store
                fS(tf).frac(end).nbPWs=2; %store the number of pinwheels as 2
            end
            clear idPwTrack1 idPwTrack2
        else %if we only have 1 pw with the fracture
            if ~isempty(find(x==round(pwCoord(1)) & y==round(pwCoord(2)))) %if the pw is in the cache
                fS(tf).frac(end).path(1,:)=pwCoord([1 2]); %store location of pw as first in path
                fS(tf).frac(end).pw(1,:)=pwCoord([1 2]); %store pw location
            end
            for l=1:length(fracopen(i).path(:,:))
                %if a point of the path is in the cache, add it to the path
                %in fS
                if ~isempty(find(round(fracopen(i).path(l,1))==x & round(fracopen(i).path(l,2))==y))
                    if ~isfield(fS(tf).frac(end),'path') %if the field doesn't exist, because the pinwheel was not in the cache, add the first point of the path
                        fS(tf).frac(end).path(1,:)=fracopen(i).path(l,:);
                    else %append the point to the path
                        fS(tf).frac(end).path(end+1,:)=fracopen(i).path(l,:);
                    end
                %The next portion appears unnecessary: the path is already
                %marked open earlier and if a point is not in the cache, it
                %simply isn't added to the path. This would overwrite a
                %path with an empty one if any point is outside the cache,
                %rather than saving only those in the cache.
%                 elseif isfield(fS(tf).frac(end),'path')
%                     if ~isempty(fS(tf).frac(end).path)
%                         fS(tf).frac(end).type='open';
%                         fS(tf).frac(end).path=[];
%                     end
                end
            end
            fS(tf).frac(end).nbPWs=1;
        end
    end

    %0 PW in the closed fracture
    for i=1:length(fracclose)
        fS(tf).frac(end+1).type='closed';
        fS(tf).frac(end).nbPWs=0;
        try
            for l=1:length(fracclose(i).path(:,:))
                if ~isempty(find(round(fracclose(i).path(l,1))==x & round(fracclose(i).path(l,2))==y)) %check if in the cache
                    if ~isfield(fS(tf).frac(end),'path') %if there is no path field
                        fS(tf).frac(end).path(1,:) = fracclose(i).path(l,:);
                    else %if path exists, add to it
                        fS(tf).frac(end).path(end+1,:)=fracclose(i).path(l,:);
                    end
                end
            end
        catch
            disp('No closed fractures. Moving on.')
        end
    end
    %clear any empty first rows
    if isempty(fS(tf).frac(1).type)
        fS(tf).frac(1) = [];
    end
end

clear A y x a i l sFrac tf idPw idPw2 pwCoord

save(strcat(path2,'fS_sf',string(cur_sf),'.mat'),'fS')
save(strcat(path2,'pwPair_sf',string(cur_sf),'.mat'),'pwPair')

%% Combine fS_sf# into fS_all
%create fS_all
fS_all = struct;
for i = 1:max_sf
    load(strcat(path2,'fS_sf',string(i),'.mat'));
    fS_all(i).sf = i;
    fS_all(i).tf = struct;
    for j = 1:length(fS)
        fS_all(i).tf(j).frac = fS(j).frac;
    end
end
save(strcat(path2,'fS_all'),'fS_all');

%% CHECK THAT EVERYTHING IS WELL DRAWN
global DIRall;
load(strcat(path1,'DIRall_test.mat'), 'DIRall');
load(strcat(path2,'fS_sf',string(cur_sf),'.mat'),'fS')

%close all
figure('Position',[1,1,600,300]);
tf_colors = hsv(8);

for tf=1:max_tf
    dim = ceil(sqrt(max_tf)); %comment out if looking at overlay of all tf
    subplot(dim,dim,tf) %comment out for overlay
    hold on;
    im = angle(exp(1i*reshape(squeeze(DIRall(tf,:)),ySignal,xSignal)));
    imagesc(im); %comment out for overlay
    colormap(hsv)
    for i=1:length(fS(tf).frac)
        if isfield(fS(tf).frac,'path') %Hannah Edit
            if ~isempty(fS(tf).frac(i).path)
                x=fS(tf).frac(i).path(:,1);
                y=fS(tf).frac(i).path(:,2);
                plot3(x,y,tf*ones(size(x)),'Color','k','linewidth',3); %switch 'k' to tf_colors(tf,:) if overlay
            end
            %comment out 'if statement' and title if overlay
            if ~isempty(fS(tf).frac(i).pw)
                xpw=fS(tf).frac(i).pw(:,1);
                ypw=fS(tf).frac(i).pw(:,2);
                scatter(xpw,ypw,50,'w','filled');
            end
            title(['tf = ' sprintf('%g',tf)]);
        else %Hannah Edit
            disp('No path stored.')
        end
    end
    xlim([1 xSignal]);
    ylim([1 ySignal]);
end

clear tf i x y xpw ypw dim
%% COUNT PAIRS

load(strcat(path2,'pwPair_sf',string(cur_sf),'.mat'),'pwPair')
name = strcat('pwTrack_tf',string(max_tf),'.mat');
load(strcat(path1,name), 'pwTrack');

pair_stats = struct;

nPw=length(pwTrack);
pair_stats.nPw = nPw;
partner=zeros(nPw,max_tf);
for i=1:nPw
    for tf=1:max_tf
        pId1=find(pwPair(:,1)==i & pwPair(:,5)==tf);
        pId2=find(pwPair(:,2)==i & pwPair(:,5)==tf);
        if length(pId1)<2 && length(pId2)<2
            if ~isempty(pId1) && isempty(pId2)
                partner(i,tf)=pwPair(pId1,2);
            elseif ~isempty(pId2) && isempty(pId1)
                partner(i,tf)=pwPair(pId2,1);
            end
        end
    end
end
pair_stats.partner = partner;

pwCache=ones(size(pair_stats.partner(:,1)));
for tfi=2:max_tf-2 %2:7
    pwCache = pwCache & (pair_stats.partner(:,tfi) ~= 0);
end
pCount=zeros(1,nPw); %count by pinwheel of #changes
ptfCount=zeros(1,max_tf); %count by tf of #changes
for j=1:sum(pwCache) %nPw
    A=find(pwCache);
    i=A(j);
    for tf=1:max_tf-1
        if partner(i,tf+1)~=partner(i,tf)
            pCount(i)=pCount(i)+1;
            ptfCount(tf)=ptfCount(tf)+1;
        end
    end
      num_change = length(unique(nonzeros(pair_stats.partner(j,:))));
      ind_change = find(pair_stats.partner(j,:));
      if ~isempty(ind_change)
         if ind_change(1) ~= 0
            cur_val = pair_stats.partner(j,ind_change(1));
            for k = 2:length(ind_change)
                if cur_val ~= pair_stats.partner(j,ind_change(k))
                    pCount(1,j) = pCount(1,j) + 1;
                    ptfCount(1,ind_change(k)) = ptfCount(1,ind_change(k)) + 1;
                    cur_val = pair_stats.partner(j,ind_change(k));
                end
            end
         end
      end
end
pair_stats.pCount = pCount;
pair_stats.ptfCount = ptfCount;
    

nPw2=0; %total number of paired pinwheels (instead of single)
for i=1:length(partner)
    if sum(partner(i,:)>0)
        nPw2=nPw2+1;
    end
end
pair_stats.nPw2 = sum(pwCache); %nPw2;

figure;
plot(TFvalues(1:max_tf),ptfCount(1:max_tf)/nPw2,'-o')
ylim([0 .2]);
xlim([0.3 4.5]);
title('ratio of pinwheel changing pairs = f(tf)');
set(gca, 'XScale', 'log');

save(strcat(path2,'pair_stats_sf',string(cur_sf),'.mat'),'pair_stats')

clear nPw i tf nPw2 pId1 pId2

%% PLOT AVG LENGTH AND CURVATURES

%This section creates statistics about the fractures and stores them

%If fS values are already saved on file somewhere
load(strcat(path2,'fS_sf',string(cur_sf),'.mat'));

%Variables:
min_size = 1;

%First find the slope at each point.
for tf=1:max_tf
    for i=1:length(fS(tf).frac)
        fracPath=fS(tf).frac(i).path;
        %create necessary fields
        fS(tf).frac(i).orPath = [];
        fS(tf).frac(i).normPath = [];
        fS(tf).frac(i).normPathFromS = [];
        if ~isempty(fracPath) && size(fS(tf).frac(i).path,1)>1 %if the path exists and is longer than 1 point
            for j=2:length(fracPath)
                fS(tf).frac(i).normPathFromS(j)=0;
                    %Next take the arctangent so that, instead of having a 2D path, you get a
                    %list of angles.
                    
                    %Hannah edit: this if statement ensures we don't look
                    %at double clicks
                    if (abs(round(fracPath(j,1)-fracPath(j-1,1))) > 0) && (abs(round(fracPath(j,2)-fracPath(j-1,2))) > 0)
                        fS(tf).frac(i).orPath(j-1)=atan((fracPath(j,1)-fracPath(j-1,1))/...
                                                        (fracPath(j,2)-fracPath(j-1,2))...
                                                        ); 
                        %Next find the distance between points and add to find the total distance.
                        fS(tf).frac(i).normPath(j-1)=sqrt((fracPath(j,2)-fracPath(j-1,2))^2+...
                                                        (fracPath(j,1)-fracPath(j-1,1))^2);
                        %The accumulation of the distance between points so you
                        %have total distance for the path.
                        fS(tf).frac(i).normPathFromS(j)=fS(tf).frac(i).normPathFromS(j-1)+...
                                                fS(tf).frac(i).normPath(j-1);
                    end
            end
        end
        g=fspecial('gaussian',70,10); %returns a rotationally symmetric Gaussian lowpass filter of size hsize with standard deviation sigma. Not recommended. Use imgaussfilt or imgaussfilt3 instead.
        dx = 1000*conv(g(8,:),[-1 1],'same');
        dx=dx(1:end-1); %remove last value of dx
        if length(fS(tf).frac(i).orPath)>2
           [Y,Ty]=resample(fS(tf).frac(i).orPath,fS(tf).frac(i).normPath,1); %resample at rate of (length of segment)*original rate
           fS(tf).frac(i).reSampledOrPath=Y;
           fS(tf).frac(i).reSampledNorm=Ty;
        else
           fS(tf).frac(i).reSampledOrPath=0;
           fS(tf).frac(i).reSampledNorm=1;
        end
        %Now you have resampled interpolated paths that are well spaced

        if isnan(sum(fS(tf).frac(i).reSampledOrPath))
           fS(tf).frac(i).cuPath=0;
        else
           fS(tf).frac(i).cuPath=conv(fS(tf).frac(i).reSampledOrPath,dx,'same');%*avgStep;
        end
        fS(tf).frac(i).curv=sum(fS(tf).frac(i).cuPath);
        fS(tf).frac(i).curvAbs=mean(abs(fS(tf).frac(i).cuPath));
    end
end

tfCount=zeros(1,max_tf);
tfCountbis=zeros(1,max_tf);
tfClosedCount=zeros(1,max_tf);
sizeFtf=zeros(1,max_tf);sizeFtf2=zeros(1,max_tf);sizeFVartf=zeros(1,max_tf);
sztodist=zeros(1,max_tf);sztodist2=zeros(1,max_tf);
distPWsavg = zeros(1,max_tf); distPWsVartf=zeros(1,max_tf);
curvFtf=zeros(1,max_tf);curvFtf2=zeros(1,max_tf);curvFVartf=zeros(1,max_tf);
curv2Ftf=zeros(1,max_tf);curv2Ftf2=zeros(1,max_tf);curv2FVartf=zeros(1,max_tf);

for tf=1:max_tf
    tfCount(tf)=sum(cellfun(@length,{fS(tf).frac.path})>0); %count how many have unempty paths
    tfCountbis(tf)=sum([fS(tf).frac.nbPWs]>1); %count how many have more than one pinwheel
    distPWsavg(tf) = sum([fS(tf).frac.distPWs])/length([fS(tf).frac.distPWs]); %average pw distance
    fS(tf).size = zeros(1,length(fS(tf).frac));
    for i=1:length(fS(tf).frac) %Hannah edit
        fSize=sum(fS(tf).frac(i).normPath);
        fS(tf).size(1,i) = fSize;
        if cellfun(@length,{fS(tf).frac(i).path})>min_size %look at only those with a path greater than a set minimum size
            if fSize>1 %to avoid noise and only work with large fractures
                sizeFtf(tf)=sizeFtf(tf)+fSize/tfCount(tf);
                sizeFtf2(tf)=sizeFtf2(tf)+fSize^2/tfCount(tf);
                if ~isempty(fS(tf).frac(i).distPWs) && fS(tf).frac(i).distPWs~=0 %look at those with two pinwheels who have a distance between them greater than 0
                    %calculate the variance in distance between pinwheels here
                    distPWsVartf(tf) = distPWsVartf(tf) + (fS(tf).frac(i).distPWs - distPWsavg(tf))^2/tfCountbis(tf);
                    sztodist(tf)=sztodist(tf)+ fSize/(fS(tf).frac(i).distPWs)/tfCountbis(tf); %+ pathlength/pwdistance/number-of-2pw-fractures
                    sztodist2(tf)=sztodist2(tf)+ (fSize/fS(tf).frac(i).distPWs)^2/tfCountbis(tf); %+ (pathlength/pwdistance)^2/number-of-2pw-fractures
                end
                %The curvFtf is the expectation of the curvature, and the curvFtf2 is
                %the absolute value of the curvature which should not be zero.
                %Then the square is the variance.
                curvFtf(tf)=curvFtf(tf)+fS(tf).frac(i).curv/fSize/tfCount(tf); %in order to give a weight to the curvature by the length of the fracture we divide
                curvFtf2(tf)=curvFtf2(tf)+(fS(tf).frac(i).curv/fSize)^2/tfCount(tf);
                %The absolute value is what we want (why?)
                curv2Ftf(tf)=curv2Ftf(tf)+fS(tf).frac(i).curvAbs/tfCount(tf);
                curv2Ftf2(tf)=curv2Ftf2(tf)+(fS(tf).frac(i).curvAbs)^2/tfCount(tf);

                if strcmp(fS(tf).frac(i).type,'closed')
                    tfClosedCount(tf)=tfClosedCount(tf)+1;
                end
            end
        end
    end
    sizeFVartf(tf)= (sizeFtf2(tf) - sizeFtf(tf).^2)*tfCount(tf)/(tfCount(tf)-1);
   %Modify to above changes to distance variable: distPWsVartf(tf)= (sztodist2(tf)- sztodist(tf).^2)*tfCountbis(tf)/(tfCountbis(tf)-1)/sizeFtf(tf)^2;
    curvFVartf(tf)= (curvFtf2(tf) - curvFtf(tf).^2)*tfCount(tf)/(tfCount(tf)-1);
    curv2FVartf(tf)= (curv2Ftf2(tf) - curv2Ftf(tf).^2)*tfCount(tf)/(tfCount(tf)-1);
end
curv2Err=sqrt(curv2FVartf(1:max_tf)./tfCount(1:max_tf)); %sqrt(error of sqrd curvature / num tf)
sizeErr=sqrt(sizeFVartf(1:max_tf)./tfCount(1:max_tf));

%Store all statistics
fracture_stats = struct;
fracture_stats.tfCount = tfCount;
fracture_stats.tfCountbis = tfCountbis;
fracture_stats.tfClosedCount = tfClosedCount;
fracture_stats.sizeErr = sizeErr;
fracture_stats.sizeFtf = sizeFtf;
fracture_stats.sizeFtf2 = sizeFtf2;
fracture_stats.sizeFVartf = sizeFVartf;
fracture_stats.sztodist = sztodist;
fracture_stats.sztodist2 = sztodist2;
fracture_stats.distPWsavg = distPWsavg;
fracture_stats.distPWsVartf = distPWsVartf;
fracture_stats.curvFtf = curvFtf;
fracture_stats.curvFtf2 = curvFtf2;
fracture_stats.curvFVartf = curvFVartf;
fracture_stats.curv2Ftf = curv2Ftf;
fracture_stats.curv2Ftf2 = curv2Ftf2;
fracture_stats.curv2FVartf = curv2FVartf;
fracture_stats.curv2Err = curv2Err;

curvGraph=figure('Position',[0,450,500,450]);
%subplot 1: 
subplot(2,2,1);
hold on
plot((TFvalues(1:max_tf)),sizeFtf(1:max_tf),'-ob');
errorbar((TFvalues(1:max_tf)),sizeFtf(1:max_tf),sizeErr,'ob');
hold on 
title('Average Fracture Size');
yticks([80 100 120 140]);
yticklabels([80 100 120 140].*0.0123);
xlim( [0.3    4.5]);
set(gca, 'XScale', 'log');
%subplot 2: 
subplot(2,2,2);
hold on
plot(log(TFvalues(1:max_tf)),tfCount(1:max_tf),'-^','LineWidth',3)
plot(log(TFvalues(1:max_tf)),tfClosedCount(1:max_tf),'-o','LineWidth',3)
title('Number of Fractures per TF');
legend('Total TF','Closed TF')
xlim( [-1.1    1.5]);
%subplot 3:
subplot(2,2,3);
hold on
plot((TFvalues(1:max_tf)),curv2Ftf(1:max_tf),'-o') %absolute value of curvature (theory that certain TF curve more)
errorbar((TFvalues(1:max_tf)),curv2Ftf(1:max_tf),curv2Err,'ob')
title('Absolute Curvature / Length');
    xlim( [0.3    4.5]);
set(gca, 'XScale', 'log');
%subplot 4:
subplot(2,2,4);
hold on
totSize=sizeFtf.*tfCount;totSizeErr=sqrt(sizeFVartf(1:max_tf).*tfCount(1:max_tf));
plot(log(TFvalues(1:max_tf)),log(totSize(1:max_tf)),'-o');
plot(log(TFvalues(1:max_tf)),log(totSize(1:max_tf)-totSizeErr),'v');
plot(log(TFvalues(1:max_tf)),log(totSize(1:max_tf)+totSizeErr),'^');
title(['<sizefrac> * #frac = f(TF)']);
xlim( [-1.1    1.5]);

save(strcat(path2,'fS_sf',string(cur_sf),'.mat'),'fS')
saveas(curvGraph,strcat(path1,'curvature_Stats.png'))
save(strcat(path1,'fracture_stats_sf',string(cur_sf),'.mat'),'fracture_stats');

clear tf i fracPath j g dx Y Ty fSize tfClosedCount tfCountbis tfCount
clear curv2Err curv2Ftf curvFtf curvFtf distPWsavg sztodist sizeFtf 
                     
%% Visualize Manual Fracture Changes - Full Map
%Plots all fractures on a black background. Each TF is a color of the
%rainbow (red for tf = 1, and onwards), and each SF has its own subplot.

%Load Data
load(strcat(path2,'fS_all.mat'))

%Create Variables
tf_colors = hsv(9);
grey = [0.7, 0.7, 0.7];
screenw = 1440;
screenh = 900;
plot_width = 1.5;

f = figure(1000);
f.Position = [1,1,screenw,screenh];
for sf = SF_tested
    ha(sf) = subplot(2,2,sf);
    for tf = TF_tested
        hold on
        for i = 1:length(fS_all(sf).tf(tf).frac)
            path = fS_all(sf).tf(tf).frac(i).path;
            if length(path) > 0
                color = tf_colors(tf,:);
                plot(path(:,1),path(:,2),'Color', color,'LineWidth',plot_width)
                if ~isempty(fS_all(sf).tf(tf).frac(i).pw)
                    for pw = 1:fS_all(sf).tf(tf).frac(i).nbPWs
                        scatter(fS_all(sf).tf(tf).frac(i).pw(pw,1),fS_all(sf).tf(tf).frac(i).pw(pw,2),[],grey,'filled')
                    end
                end
            end
        end
    end
    set(gca,'Color','k','XColor', 'none','YColor','none')
    title(strcat('SF = ',string(sf)))
end
linkaxes(ha)

clear sf tf i ha path

%% Visualize Manual Fracture Changes - Speed Colors

%Import existing fractures
load(strcat(path2,'fS_all.mat'))

%Calculate speed values
speed_vals = [];
match_tf = [];
match_sf = [];
for tf_i = TF_tested
    for sf_i = SF_tested
        speed_vals = [speed_vals, TFvalues(tf_i)/SFvalues(sf_i)];
        match_tf = [match_tf, tf_i];
        match_sf = [match_sf, sf_i];
    end
end
[sort_speed, sort_i] = sort(speed_vals);
sort_tf = match_tf(sort_i);
sort_sf = match_sf(sort_i);


%Plot fractures over blacked out figure
speed_colors = jet(length(speed_vals));
figure()
hold on
for s_i = 1:length(sort_speed)
    try
        speed_val = sort_speed(s_i);
        tf_i = sort_tf(s_i);
        sf_i = sort_sf(s_i);
        for i = 1:length(fS_all(sf_i).tf(tf_i).frac)
            frac_path = fS_all(sf_i).tf(tf_i).frac(i).path;
            try
                plot(frac_path(:,1),frac_path(:,2),'Color', speed_colors(s_i,:),'LineWidth',3)
                if ~isempty(fS_all(sf_i).tf(tf_i).frac(i).pw)
                    frac_pw = fS_all(sf_i).tf(tf_i).frac(i).pw;
                    scatter(frac_pw(:,1),frac_pw(:,2),'w','filled')
                end
            catch
                disp('Bad Frac.')
            end
        end
    catch
        disp(strcat('Skipping SF ',string(sf_i),' TF',string(tf_i)))
    end
end
%Add scalebar plot
mm_px = (1/1000)/m_px;
plot([1,mm_px],[10,10],'Color','white','LineWidth',3)

set(gca,'visible','on','position',[0,0,1,1],'color','k')
set(gcf,'Units','Pixels','Position',[0,0,2*xSignal,2*ySignal],'PaperUnits','points','PaperSize',[xSignal,ySignal])
iptsetpref('ImshowBorder','tight')
saveas(gcf,strcat(path2,'speed_frac_combined_pw.png'))
print(gcf,'-vector','-dsvg',strcat(path2,'speed_frac_combined_pw.svg'))

%Plot fractures over blacked out figure
speed_colors = jet(length(speed_vals));
figure()
hold on
for s_i = 1:length(sort_speed)
    try
        speed_val = sort_speed(s_i);
        tf_i = sort_tf(s_i);
        sf_i = sort_sf(s_i);
        for i = 1:length(fS_all(sf_i).tf(tf_i).frac)
            frac_path = fS_all(sf_i).tf(tf_i).frac(i).path;
            try
                plot(frac_path(:,1),frac_path(:,2),'Color', speed_colors(s_i,:),'LineWidth',3)
            catch
                disp('Bad Frac.')
            end
        end
    catch
        disp(strcat('Skipping SF ',string(sf_i),' TF',string(tf_i)))
    end
end
plot([1,mm_px],[10,10],'Color','white','LineWidth',3)

set(gca,'visible','on','position',[0,0,1,1],'color','k')
set(gcf,'Units','Pixels','Position',[0,0,2*xSignal,2*ySignal],'PaperUnits','points','PaperSize',[xSignal,ySignal])
iptsetpref('ImshowBorder','tight')
saveas(gcf,strcat(path2,'speed_frac_combined_nopw.png'))
print(gcf,'-vector','-dsvg',strcat(path2,'speed_frac_combined_nopw.svg'))

%% Import DIR Maps and use Gradient to Auto Find Fractures

frac_all = zeros(nSF,nTF,ySignal,xSignal);
for sf_i = SF_tested
    %import DIR map as DIRall
    load(strcat(path2,'dataSF',string(sf_i),'/DIRall_test.mat'));
    for tf_i = TF_tested
        DIR_pref = angle(exp(1i*reshape(squeeze(DIRall(tf_i,:)),ySignal,xSignal)));
        %calculate DIR gradient for frac
        [Gmag_DIR,Gdir_DIR] = imgradient(DIR_pref ,"prewitt");
        %Binarize
        DIR_frac_wide = Gmag_DIR > pi;
        SE = strel('disk',1);
        DIR_frac_close = imclose(DIR_frac_wide,SE);
        DIR_frac_open = imopen(DIR_frac_close,SE);
        DIR_skel = bwskel(DIR_frac_open);
        DIR_frac = DIR_skel;
        frac_all(sf_i,tf_i,:,:) = DIR_frac;
    end
end

%% Plot fractures with speed coloring

mkdir(strcat(path2,'frac_scatter/'))

%Get possible speed values
speed_vals = [];
for tf_i = TF_values
    for sf_i = SF_values
        speed_vals = [speed_vals, tf_i/sf_i];
    end
end
speed_vals = unique(sort(speed_vals));
speed_colors = jet(length(speed_vals));
%Reorg fracs by speed
frac_speed = zeros(length(speed_vals),ySignal,xSignal);
for sf_i = SF_tested
    for tf_i = TF_tested
        speed_val = TF_values(tf_i)/SF_values(sf_i);
        speed_ind = find(speed_vals == speed_val,1);
        frac_speed(speed_ind,:,:) = squeeze(frac_all(sf_i,tf_i,:,:));
    end
end

%Create RGB image of fractures
A = zeros(ySignal,xSignal,3); %black canvas
for speed_i = 1:length(speed_vals)
    frac_i = squeeze(frac_speed(speed_i,:,:));
    [frac_row,frac_col] = find(squeeze(frac_speed(speed_i,:,:)));
    for i_pair = 1:length(frac_row)
        A(frac_row(i_pair),frac_col(i_pair),1) = speed_colors(speed_i,1)*255;
        A(frac_row(i_pair),frac_col(i_pair),2) = speed_colors(speed_i,2)*255;
        A(frac_row(i_pair),frac_col(i_pair),3) = speed_colors(speed_i,3)*255;
    end
end
figure();
imagesc(A)
set(gcf,'Units','Pixels','Position',[0,0,xSignal,ySignal],'PaperUnits','points','PaperSize',[xSignal,ySignal])
iptsetpref('ImshowBorder','tight')
saveas(gcf,strcat(path2,'frac_scatter/','auto_fractures.svg'))
saveas(gcf,strcat(path2,'frac_scatter/','auto_fractures.png'))