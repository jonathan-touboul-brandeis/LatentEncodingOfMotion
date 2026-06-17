%% FIGURE 3 - Papier V1- Basic info on the data
clear all
close all
catArray=['IJXX';'IJA2';'IJB2'];
xSignalArray=[510,498,510];  %510, 498, 510
ySignalArray=[330,330,348];  %330, 330, 348

TFvalues=[0.3800 0.8800 1.2400 1.5300 2.0000 2.7800 3.2600 3.8500 5.3700];
SFvalues=[0.15 0.29 0.5 0.9 1.65];

TF_tested=1:9;
SF_tested=1:4;

nTF=length(TF_tested);
nSF=length(SF_tested);

%TF_values=TF_values_original(TF_tested);
%SF_values=SF_values_original(SF_tested);

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
ONLY_NECESSARY=1;

%% Load Direction Maps

disp(pathDirCat)
DIRmap = zeros(4,9,ySignal*xSignal);
for sf=SF_tested
    dataSFname=['dataSF' sprintf('%g',sf)];
    pathDirCatSF =[pathDirCat dataSFname '/'];
    load([pathDirCatSF 'DIRall_test.mat'],'DIRall');
    DIRmap(sf,:,:) = DIRall;
end
DIRmap_Shaped = reshape(DIRmap,4,9,ySignal,xSignal);

DIRmap_Shaped_Win=DIRmap_Shaped(:,:,yWin,xWin);
%% Fig 3A- Similarity between direction maps with 2D Gaussian Fit (to be cleaned up)
dirmapdist = zeros(length(SF_tested),length(TF_tested),length(SF_tested),length(TF_tested));
N=20;
dsf_grid_number = 2*length(SF_tested);
dtf_grid_number_list = 17;% Obtained via parameter search on the (10:2:18) space
dtf_cutoff_list = [2]; % Obtained via parameter search on the (0.7:0.2:3) space

N_bstrp=100;
p=0.4;
coeffs=zeros(1,N_bstrp);
rhos=zeros(1,N_bstrp);
PriebeAlpha=zeros(1,N_bstrp);
NN=size(DIRmap_Shaped_Win(1,1,:),3);
N=ceil(p*NN);
for bootstrap=1:N_bstrp
    U=randperm(NN);
    U=U(1:N);
    % U=1:NN;
    % for dtf_grid_number_i = 1:length(dtf_grid_number_list)
    dtf_grid_number_i =1;
    % for cutoff_i = 1:length(dtf_cutoff_list)
    cutoff_i=1;
        dtf_grid_number = dtf_grid_number_list(dtf_grid_number_i);
        dtf_cutoff = dtf_cutoff_list(cutoff_i);
        
        dsf_grid=linspace(-2.1,2.1,dsf_grid_number);
        dtf_grid=linspace(-2.6,2.6,dtf_grid_number);
        dsf = dsf_grid(2:end);
        dtf = dtf_grid(2:end);

        trace_diff = zeros(dsf_grid_number-1,dtf_grid_number-1);
        trace_diff_S = zeros(dsf_grid_number-1,dtf_grid_number-1);
        coeff_diff = zeros(dsf_grid_number-1,dtf_grid_number-1);
        Similarity_index=NaN*ones(45,dsf_grid_number-1,dtf_grid_number-1);
        AllSims=zeros(40,dsf_grid_number-1,dtf_grid_number-1);
        
        for tf = TF_tested
            for TF = TF_tested
                for sf = SF_tested
                    for SF = SF_tested
                        % Log-likelihood ratio
                        dtf_idx = find(dtf - log(TFvalues(tf)/TFvalues(TF))>0,1);
                        dsf_idx = find(dsf - log(SFvalues(sf)/SFvalues(SF))>0,1);
                        dirmapdist(sf,tf,SF,TF) = mean(abs(circ_dist(DIRmap_Shaped_Win(sf,tf,U),DIRmap_Shaped_Win(SF,TF,U))));
                        trace_diff(dsf_idx,dtf_idx) = trace_diff(dsf_idx,dtf_idx) + (1-dirmapdist(sf,tf,SF,TF));                         
                        trace_diff_S(dsf_idx,dtf_idx) = trace_diff_S(dsf_idx,dtf_idx) + dirmapdist(sf,tf,SF,TF).^2;
                        coeff_diff(dsf_idx,dtf_idx) = coeff_diff(dsf_idx,dtf_idx) + 1; 
                        Similarity_index(coeff_diff(dsf_idx,dtf_idx),dsf_idx,dtf_idx)=mean(mean(abs(circ_dist(squeeze(DIRmap_Shaped_Win(sf,tf,:,:)),squeeze(DIRmap_Shaped_Win(SF,TF,:,:))))>pi/2));
                        AllSims(coeff_diff(dsf_idx,dtf_idx),dsf_idx,dtf_idx)=1-dirmapdist(sf,tf,SF,TF);
                    end
                end
            end
        end
%%
        trace_diff = trace_diff./coeff_diff;
        % trace_diff=1-squeeze(nanmean(Similarity_index));
        M_trace_diff = trace_diff;
        
        dtf_p = (dtf_grid(2:end) + dtf_grid(1:end-1))/2;
        dsf_p = (dsf_grid(2:end) + dsf_grid(1:end-1))/2;
        
        [X,Y] = meshgrid(dtf_p,dsf_p);
        M_trace_diff_copy = M_trace_diff;
        M_trace_diff(abs(X)>dtf_cutoff | abs(Y)>dtf_cutoff) = NaN;

        % figure;hold off
        % imagesc(dsf_p,dtf_p,M_trace_diff_copy');hold on
        % colormap(parula)
        % caxis([0 1])
        % colorbar;

        % min_M_trace = nanmin(M_trace_diff(:));
        % sum_M_trace = nansum(M_trace_diff(:)-min_M_trace);
        % p_M_trace = (M_trace_diff(:) - min_M_trace)/sum_M_trace;
        % 
        % p_M_trace(isnan(p_M_trace)) = 0;
        % mean_M_trace = ([X(:),Y(:)]' .* p_M_trace(:))';
        % cov_M_trace = ([X(:),Y(:)] - mean_M_trace)'.* repmat(p_M_trace(:),1,2)' * ([X(:),Y(:)] - mean_M_trace);

        % n_grid = 20;
        % [X2,Y2] = meshgrid(linspace(min(dtf_p),max(dtf_p),n_grid),linspace(min(dsf_p),max(dsf_p),n_grid));
        %%
        [params, Z_fit, Xi, axes_info] = fit_gaussian_to_matrix(M_trace_diff_copy, dtf_p, dsf_p,0);
        fprintf('Xi = %.3f, angle = %.1f deg, slope=%.2f \n', Xi, axes_info.angle_deg,axes_info.slope);
        axes_info.Sigma
        %%
        % for k=1:N_bstrp
        %     NewMat=zeros(size(M_trace_diff_copy));
        %     for i=1:size(NewMat,1)
        %         for j=1:size(NewMat,2)
        %             N=ceil(p*coeff_diff(i,j));
        %             U=randperm(coeff_diff(i,j));
        %             NewMat(i,j)=mean(AllSims(U(1:N),i,j));
        %         end
        %     end
        %     [params, Z_fit, Xi, axes_info] = fit_gaussian_to_matrix(NewMat, dtf_p, dsf_p,0);
        %     fprintf('Xi = %.3f, angle = %.1f deg\n', Xi, axes_info.angle_deg);
            coeffs(bootstrap)=Xi;
            rhos(bootstrap)=axes_info.rho;
            PriebeAlpha(bootstrap)=axes_info.priebe_alpha;
        end
        p_val0=mean(rhos<0)
        p_val1=mean(rhos>1);
        figure;hist(rhos);
        
        % pdf_M_trace = mvnpdf([X2(:),Y2(:)],mean_M_trace,cov_M_trace);
        % fit_M_trace = pdf_M_trace * sum_M_trace+ min_M_trace;
        % fit_M_trace = reshape(fit_M_trace,size(X2,1),size(X2,2));

        % 
        % 
        % for i=1:2
        %     plot([0,((D(i,i))*v(1,i))],[0,((D(i,i))*v(2,i))],'-o','linewidth',1.5)
        % end
        % 
        % plot([-3,3],[-3,3],'-.ok','linewidth',.75)
        % contour(X2,Y2,fit_M_trace,2,'r--','LineWidth',1)
        % xlabel('$\Delta$SF','Interpreter','latex','FontSize',14);
        % ylabel('$\Delta$TF','Interpreter','latex','FontSize',14);
        % set(gca, 'YDir','normal')
        % title('Similarity between direction maps','Interpreter','latex','FontSize',14);
        % figure;pcolor(1-squeeze(nanmean(Similarity_index)));shading flat
        % figure;pcolor(squeeze(nanstd(Similarity_index)));shading flat; title('Standard deviation')
        % figure;pcolor(squeeze(nanstd(Similarity_index))./(1-squeeze(nanmean(Similarity_index))));shading flat;title('coefficient of variation')
        
%     end      
% 
% end

% print -depsc epsFig3A

%%
size(dirmapdist)
n_conds=nSF*nTF;
AllDistances=zeros(n_conds);

i=0;
LogDifferenceTF=zeros(n_conds);
LogDifferenceSF=zeros(n_conds);
for tf = TF_tested
    for sf = SF_tested
        i=i+1;
        j=0;
        for TF = TF_tested
            for SF = SF_tested
                j=j+1;
                AllDistances(i,j)=dirmapdist(sf,tf,SF,TF);
                LogDifferenceTF(i,j)=abs(log(TFvalues(TF)/TFvalues(tf)));
                LogDifferenceSF(i,j)=abs(log(SFvalues(SF)/SFvalues(sf)));
            end
        end
    end
end
% figure;
% imagesc(AllDistances)
% shading flat
% hold on
% 
% x=0:nSF:n_conds;
% [X,Y] = meshgrid(x+0.5,x+0.5);
% plot(X,Y,'k-','Linewidth',2)
% hold on
% plot(Y,X,'k-','Linewidth',2);
% 
% figure;
% imagesc(LogDifferenceTF)
% shading flat
% hold on
% 
% x=0:nSF:n_conds;
% [X,Y] = meshgrid(x+0.5,x+0.5);
% plot(X,Y,'k-','Linewidth',2)
% hold on
% plot(Y,X,'k-','Linewidth',2);
% 
% figure;
% imagesc(LogDifferenceSF)
% shading flat
% hold on
% 
% x=0:nSF:n_conds;
% [X,Y] = meshgrid(x+0.5,x+0.5);
% plot(X,Y,'k-','Linewidth',2)
% hold on
% plot(Y,X,'k-','Linewidth',2);
%%
[params, Z_fit, Xi, axes_info] = fit_gaussian_to_matrix(M_trace_diff_copy, dtf_p, dsf_p,1);
fprintf('Xi = %.3f, angle = %.1f deg, slope=%.2f \n', Xi, axes_info.angle_deg,axes_info.slope);
axes_info.Sigma
axes_info.priebe_alpha

%%
N_speeds=8;
Speeds_min=min(TFvalues)/max(SFvalues);
Speeds_max=max(TFvalues)/min(SFvalues);
Speeds_binned=linspace(log(Speeds_min),log(Speeds_max/2.5),N_speeds);
SpeedDistances=zeros(N_speeds);
SpeedNumbers=zeros(N_speeds);
AllDists=NaN*ones(100,N_speeds,N_speeds);
NumEltsTaken=zeros(N_speeds);
Near=[];
Mids=[];
Fars=[];
for tf = TF_tested
    for sf = SF_tested
        i=i+1;
        j=0;
        speed=tf/sf;
        speed1=find(log(speed)>Speeds_binned,1,'last');
        for TF = TF_tested
            for SF = SF_tested
                if abs(sf-SF)>0 || abs(tf-TF)>0
                    speed=TF/SF;
                    speed2=find(log(speed)>Speeds_binned,1,'last');
                    SpeedNumbers(speed1,speed2) = SpeedNumbers(speed1,speed2) + 1;
                    SpeedDistances(speed1,speed2) = SpeedDistances(speed1,speed2) + dirmapdist(sf,tf,SF,TF);
                    NumEltsTaken(speed1,speed2)=NumEltsTaken(speed1,speed2)+1;
                    AllDists(NumEltsTaken(speed1,speed2),speed1,speed2)=dirmapdist(sf,tf,SF,TF);
                    if abs(speed1-speed2)<2
                        Near(end+1)=dirmapdist(sf,tf,SF,TF);
                    elseif abs(speed1-speed2)<4
                        Mids(end+1)=dirmapdist(sf,tf,SF,TF);
                    else
                        Fars(end+1)=dirmapdist(sf,tf,SF,TF);
                    end
                end
            end
        end
    end
end
figure;pcolor(Speeds_binned,Speeds_binned,SpeedDistances./SpeedNumbers)
hold on;plot(Speeds_binned,Speeds_binned,'w','LineWidth',2)
figure;pcolor(SpeedNumbers)

figure;
% plot([1,2,3],[mean(Near),mean(Mids),mean(Fars)],'square')
% hold on
errorbar([1,2,3],[mean(Near),mean(Mids),mean(Fars)],...
    [std(Near)/sqrt(length(Near)),std(Mids)/sqrt(length(Mids)),std(Fars)/sqrt(length(Fars))])
xlim([0.5,3.5])
ylim([0.5,1.5])


function [params, Z_fit, Xi, axes_info] = fit_gaussian_to_matrix(Z_data, tf_vec, sf_vec,Plots)
% Fit a 2D Gaussian to a TF x SF response matrix
%
% INPUTS:
%   Z_data : matrix of values (e.g. firing rates), size [numel(sf_vec) x numel(tf_vec)]
%   tf_vec : temporal frequency axis values
%   sf_vec : spatial frequency axis values
%
% OUTPUTS:
%   params     : [mu_tf, mu_sf, sigma_tf, sigma_sf, rho, A]
%   Z_fit      : fitted gaussian matrix
%   Xi         : alignment of principal axis with TF axis (1 = pure speed tuning)
%   axes_info  : struct with eigenvalues, eigenvectors, principal axis angle

[TF, SF] = meshgrid(tf_vec, sf_vec);

TF_flat = TF(:);
SF_flat = SF(:);
Z_flat  = Z_data(:);

% 2D Gaussian
gauss = @(p, x) p(6) * exp(- ( ...
    ((x(:,1)-p(1))/p(3)).^2 ...
  - 2*p(5)*((x(:,1)-p(1))/p(3)).*((x(:,2)-p(2))/p(4)) ...
  + ((x(:,2)-p(2))/p(4)).^2 ) / (2*(1-p(5)^2)));

gauss_priebe = @(p, x) p(5) * exp( ...
    - (log2(x(:,2)) - log2(p(2))).^2 / (2*p(4)^2) ...
    - (log2(x(:,1)) - log2(p(1)) - p(6)*log2(x(:,2)/p(2))).^2 / (2*p(3)^2) );

xy = [TF_flat, SF_flat];

% Initial guess
[~, idx] = max(Z_flat);
p0 = [TF_flat(idx), SF_flat(idx), ...
      range(tf_vec)/4, range(sf_vec)/4, ...
      0, max(Z_flat)];

p02 = [TF_flat(idx), SF_flat(idx), ...
      range(tf_vec)/4, range(sf_vec)/4, ...
      0, max(Z_flat),0.3];

lb = [min(tf_vec), min(sf_vec), 0,   0,   -0.99, 0];
ub = [max(tf_vec), max(sf_vec), inf, inf,  0.99, inf];

opts = optimoptions('lsqcurvefit', 'Display', 'off');
params = lsqcurvefit(gauss, p0, xy, Z_flat, lb, ub, opts);

params2 = lsqcurvefit(gauss, p02, xy, Z_flat, lb, ub, opts);


Z_fit = reshape(gauss(params, xy), size(Z_data));

Z_fit2 = reshape(gauss_priebe (params2, xy), size(Z_data));


% --- Covariance matrix from fitted params ---
sigma_tf = params(3);
sigma_sf = params(4);
rho      = params(5);

PriebeAlpha=params2(6);

Sigma = [sigma_tf^2,              rho*sigma_tf*sigma_sf;
         rho*sigma_tf*sigma_sf,   sigma_sf^2];

% --- Eigendecomposition -> principal axes ---
[V, D] = eig(Sigma);
eigenvalues = diag(D);
[~, idx] = max(eigenvalues);
principal_axis = V(:, idx);

% Xi: how aligned is the principal axis with the TF axis (speed tuning)?
% Xi = 1 -> axis aligned with TF (pure speed tuning: TF/SF = const)
% Xi = 0 -> axis aligned with SF
Xi = 1-abs(principal_axis(1));
% abs(principal_axis(2)/principal_axis(1));

% Angle of principal axis (degrees, relative to TF axis)
angle_deg = atan2d(principal_axis(2), principal_axis(1));

axes_info.eigenvectors  = V;
axes_info.eigenvalues   = eigenvalues;
axes_info.principal_axis = principal_axis;
axes_info.angle_deg     = angle_deg;
axes_info.slope         = 1/tan(angle_deg);
axes_info.Sigma         = Sigma;
axes_info.rho           = rho;
axes_info.priebe_alpha  = PriebeAlpha;

if Plots==1
figure;

% Heatmap of data and fit side by side
subplot(1,2,1)
imagesc(tf_vec, sf_vec, Z_data)
hold on
axis xy
xlabel('Temporal Frequency (Hz)')
ylabel('Spatial Frequency (cpd)')
title('Data')
colorbar

subplot(1,2,2)
imagesc(tf_vec, sf_vec, Z_fit)
axis xy
hold on

% Plot principal axes through the mean
mu_tf = params(1);
mu_sf = params(2);

% Scale eigenvectors by sqrt(eigenvalue) for visualization
for i = 2:2
    vec = V(:,i);% * sqrt(eigenvalues(i));
    % Draw axis as line through mean (both directions)
    quiver(mu_tf, mu_sf,  vec(2),  vec(1), 2, 'w', 'LineWidth', 2, 'MaxHeadSize', 0.5)
    quiver(mu_tf, mu_sf, -vec(2), -vec(1), 2, 'w', 'LineWidth', 2, 'MaxHeadSize', 0.5)
    subplot(1,2,1)
    quiver(mu_tf, mu_sf,  vec(2),  vec(1), 2, 'w', 'LineWidth', 2, 'MaxHeadSize', 0.5)
    quiver(mu_tf, mu_sf, -vec(2), -vec(1), 2, 'w', 'LineWidth', 2, 'MaxHeadSize', 0.5)


end

% Mark the mean
plot(mu_tf, mu_sf, 'w+', 'MarkerSize', 10, 'LineWidth', 2)

xlabel('Temporal Frequency (Hz)')
ylabel('Spatial Frequency (cpd)')
title(sprintf('Fit  |  Xi = %.2f  |  angle = %.1f°', Xi, axes_info.angle_deg))
colorbar

figure;
imagesc(tf_vec, sf_vec, Z_data)
hold on
plot([mu_tf-2.5,mu_tf+2.5], [mu_sf-2.5*PriebeAlpha,mu_sf+2.5*PriebeAlpha], 'w', 'LineWidth', 2)
plot([mu_tf-2.5,mu_tf+2.5], [mu_sf-2.5*1,mu_sf+2.5*1], 'k--', 'LineWidth', 2)
plot([mu_tf,mu_tf], [mu_sf-2.5,mu_sf+2.5], 'r--', 'LineWidth', 2)

axis xy
xlabel('Temporal Frequency (Hz)')
ylabel('Spatial Frequency (cpd)')
title('Data')
colorbar
end
end