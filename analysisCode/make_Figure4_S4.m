%% Figure 4AB
% A: Object motion vs. Surface motion from an example neuron (210623_n25) 
%    ==> At 3 speeds
%    ==> Direction Tuning curve
% B: Population histograms of the correlation between tuning curves 
%    for Object (left column)/Surface (right column) motion 
%    versus those of other motion types and speeds (rows)

%%% cellData_NPX_ObjSurf.mat

%%% 49 stimulus conditions
% 1-8  : Grating Object Fast
% 9-16 : Grating Object Middle
% 17-24: Grating Object Slow
% 25-32: Grating Surface Fast
% 33-40: Grating Surface Midle
% 41-48: Grating Surface Slow
% 49   : No Stim

%% Figure 4C
% C: Direction tuning curves for two example neurons (columns) measured 
%    with translating chromatic and achromatic bars (rows) 
%    at five luminance contrasts (see legend). 
%    Direction selectivity is evident for chromatic bars and noise patches 
%    (bottom row) but not achromatic bars.

%%% cellData_NPX_EqLum.mat
%
% 5 luminance levels from Blue
% 5 luminance levels from Green
% 5 luminance levels from Red
% 5 luminance levels from Gray
% Noise patch


%% load data file
load ../dataFiles/cellData_NPX_ObjSurf.mat; 


%% plot Figure 4A
% Object motion vs. Surface motion from an example neuron (210623_n25) 
%    ==> At 3 speeds
%    ==> Direction Tuning curve

figure; 
set(gcf,'position',[100 100 300 400]); 
sgtitle('Figure 4A'); 

for n=1:3 % for three speeds
    baseline = nanmean(cellData_NPX_ObjSurf(25).respMtx(:,49)); 

    %%%%%%%%%%%%%%%
    % object motion
    subplot(3,2,n*2-1); 
    dataNow = cellData_NPX_ObjSurf(25).dir_tuning(:,n);
    idx_start = (n-1)*8 + 1; 
    idx_end = n*8; 
    steNow = nanstd(cellData_NPX_ObjSurf(25).respMtx(:,idx_start:idx_end)) / ...
             sqrt(size(cellData_NPX_ObjSurf(25).respMtx,1)); 
    if size(steNow,1)==1, steNow = steNow'; end

    plot(0:45:360,[dataNow; dataNow(1)],'k.-','LineWidth',1);   hold on;
    errorbar(0:45:360,[dataNow; dataNow(1)],[steNow; steNow(1)],'k');       
    plot([0 360],[baseline baseline],'k:', 'LineWidth',1); 

    set(gca,'ylim',[0,31],'xlim',[-10, 370],...
        'xtick',0:90:360,'box','off','TickDir','out'); 
    title(['DI: ',num2str(cellData_NPX_ObjSurf(25).DI_base(n))]); 

    if n==3
        xlabel('Direction of motion (deg)');
        ylabel('Response (spk/s)');
    end

    %%%%%%%%%%%%%%%%
    % surface motion
    subplot(3,2,n*2);     
    dataNow = cellData_NPX_ObjSurf(25).dir_tuning(:,3+n); 
    idx_start = (n+2)*8 + 1; 
    idx_end = (n+3)*8; 
    steNow = nanstd(cellData_NPX_ObjSurf(25).respMtx(:,idx_start:idx_end)) / ...
             sqrt(size(cellData_NPX_ObjSurf(25).respMtx,1)); 
    if size(steNow,1)==1, steNow = steNow'; end    

    plot(0:45:360,[dataNow; dataNow(1)],'k.-','LineWidth',1);   hold on; 
    errorbar(0:45:360,[dataNow; dataNow(1)],[steNow; steNow(1)],'k');       
    plot([0 360],[baseline baseline],'k:', 'LineWidth',1); 

    set(gca,'ylim',[0,31],'xlim',[-10, 370],...
        'xtick',0:90:360,'box','off','TickDir','out'); 
    title(['DI: ',num2str(cellData_NPX_ObjSurf(25).DI_base(3+n))]); 
end
clearvars -except cellData_NPX_*; 

%% plot Figure 4B / Figure S4A
%  Population histograms of the correlation between tuning curves 
%  for Object (left column)/Surface (right column) motion 
%  versus those of other motion types and speeds (rows)

corrMtx = zeros([6,6,length(cellData_NPX_ObjSurf)]); 
pMtx = zeros([6,6,length(cellData_NPX_ObjSurf)]); 
for n = 1:length(cellData_NPX_ObjSurf)
    [r_mtxNow, p_mtxNow] = corrcoef(cellData_NPX_ObjSurf(n).dir_tuning); 
    corrMtx(:,:,n) = r_mtxNow; 
    pMtx(:,:,n) = p_mtxNow; 
end
med_corrMtx = median(corrMtx,3)

figure; 
set(gcf,'position',[100,100,500,300]); 
sgtitle('Figure 4B / S4A'); 

for n=1:36
    rowNum = ceil(n/6); 
    colNum = rem(n-1,6) + 1; 
    %if rowNum>=colNum
    subplot(6,6,(rowNum-1)*6 + colNum); 
    histogram(corrMtx(rowNum,colNum,:),-1:0.1:1,...
              'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;
    signi = find((pMtx(rowNum,colNum,:)<0.05) & (corrMtx(rowNum,colNum,:)>0)); 
    histogram(corrMtx(rowNum,colNum,signi),-1:0.1:1,...
              'FaceColor',[0,0,0],'EdgeColor',[1,1,1],'FaceAlpha',1); 
    yMax = max(histcounts(corrMtx(rowNum,colNum,:),-1:0.1:1));        
    if rowNum ~= colNum
        plot(med_corrMtx(rowNum,colNum),yMax,'kv','MarkerFaceColor','k'); 
    end
    if rowNum<6
        set(gca,'XTickLabel',[]);             
    end
    set(gca,'Ylim',[0,yMax],'YTick',[0, yMax],'Box','off','TickDir','out'); 
    %end
end

%% Figure S4B
%  from Figure S4A. 
%  Correlation histograms are divided into three groups
%  "within Object", "within Surface", "Object and Surface"

figure; 
set(gcf,'position',[100,100,150,300]); 
sgtitle('Figure S4B'); 

subplot(3,1,1); 
corr_withinObj = [reshape(corrMtx(2,1,:),[58,1]) reshape(pMtx(2,1,:),[58,1]); ...
                  reshape(corrMtx(3,1,:),[58,1]) reshape(pMtx(3,1,:),[58,1]); ...    
                  reshape(corrMtx(3,2,:),[58,1]) reshape(pMtx(3,2,:),[58,1])]; 
histogram(corr_withinObj(:,1),-1:0.1:1,...
          'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;    
signi = find((corr_withinObj(:,2)<0.05) & (corr_withinObj(:,1)>0)); 
histogram(corr_withinObj(signi,1),-1:0.1:1,...
          'FaceColor',[0,0,0],'EdgeColor',[1,1,1],'FaceAlpha',1);    
yMax = max(histcounts(corr_withinObj(:,1),-1:0.1:1));    
plot(median(corr_withinObj(:,1)),yMax,'kv','MarkerFaceColor','k'); 
set(gca,'XTickLabel',[],'Ylim',[0,yMax],'YTick',[0, yMax],'Box','off','TickDir','out'); 
title('within Object'); 

subplot(3,1,2); 
corr_withinSurf = [reshape(corrMtx(5,4,:),[58,1]) reshape(pMtx(5,4,:),[58,1]); ...
                   reshape(corrMtx(6,4,:),[58,1]) reshape(pMtx(6,4,:),[58,1]); ...    
                   reshape(corrMtx(6,5,:),[58,1]) reshape(pMtx(6,5,:),[58,1])]; 
histogram(corr_withinSurf(:,1),-1:0.1:1,...
          'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;    
signi = find((corr_withinSurf(:,2)<0.05) & (corr_withinSurf(:,1)>0)); 
histogram(corr_withinSurf(signi,1),-1:0.1:1,...
          'FaceColor',[0,0,0],'EdgeColor',[1,1,1],'FaceAlpha',1);    
yMax = max(histcounts(corr_withinSurf(:,1),-1:0.1:1));  
plot(median(corr_withinSurf(:,1)),yMax,'kv','MarkerFaceColor','k'); 
set(gca,'XTickLabel',[],'Ylim',[0,yMax],'YTick',[0, yMax],'Box','off','TickDir','out'); 
title('within Surface'); 

subplot(3,1,3); 
corr_between = [reshape(corrMtx(4,1,:),[58,1]) reshape(pMtx(4,1,:),[58,1]); ...
                reshape(corrMtx(4,2,:),[58,1]) reshape(pMtx(4,2,:),[58,1]); ...    
                reshape(corrMtx(4,3,:),[58,1]) reshape(pMtx(4,3,:),[58,1]); ...
                reshape(corrMtx(5,1,:),[58,1]) reshape(pMtx(5,1,:),[58,1]); ...
                reshape(corrMtx(5,2,:),[58,1]) reshape(pMtx(5,2,:),[58,1]); ...    
                reshape(corrMtx(5,3,:),[58,1]) reshape(pMtx(5,3,:),[58,1]); ...
                reshape(corrMtx(6,1,:),[58,1]) reshape(pMtx(6,1,:),[58,1]); ...
                reshape(corrMtx(6,2,:),[58,1]) reshape(pMtx(6,2,:),[58,1]); ...    
                reshape(corrMtx(6,3,:),[58,1]) reshape(pMtx(6,3,:),[58,1])]; 
histogram(corr_between(:,1),-1:0.1:1,...
          'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;    
signi = find((corr_between(:,2)<0.05) & (corr_between(:,1)>0)); 
histogram(corr_between(signi,1),-1:0.1:1,...
          'FaceColor',[0,0,0],'EdgeColor',[1,1,1],'FaceAlpha',1);    
yMax = max(histcounts(corr_between(:,1),-1:0.1:1));    
plot(median(corr_between(:,1)),yMax,'kv','MarkerFaceColor','k'); 
set(gca,'Ylim',[0,yMax],'YTick',[0, yMax],'Box','off','TickDir','out'); 
title('Object and Surface'); 
xlabel('Correlation coefficient (r)'); 
ylabel('Number of neurons'); 

clearvars -except cellData_NPX_ObjSurf; 

%% Figure 4C
% C: Direction tuning curves. 5 stimulus types. two example neurons

load ../dataFiles/cellData_NPX_EqLum.mat; 

figure; 
set(gcf,'position',[100,100,400,500]); 
sgtitle('Figure 4C'); 

gCells = [27 + 38,42 + 38];  % 38 cells in exp_id z200811
yMax = [50, 40];             % yMax for n27, n42

grad_now = [50:40:240]; 
grays = 50:30:170;

for m=1:length(gCells)
    n = gCells(m); 

    baseline = nanmean(cellData_NPX_EqLum(n).baseline); 

    subplot(5, 2, m); 
    for k = 1:5
        numRepeat = size(eval(['cellData_NPX_EqLum(n).blue_lum',num2str(k)]),1); 
        mean_dir = nanmean(eval(['cellData_NPX_EqLum(n).blue_lum',num2str(k)])); 
        ste_dir = nanstd(eval(['cellData_NPX_EqLum(n).blue_lum',num2str(k)])) / sqrt(numRepeat); 
        errorbar(0:45:360, mean_dir([1:8,1]), ste_dir([1:8,1]),...
                 'Color',[25 35 grad_now(k)]/255,'LineWidth',1);   hold on; 
    end
    plot([0,360],[baseline baseline],'k:'); 
    set(gca,'Box','off','Xlim',[-10,360],'Ylim',[0,yMax(m)],'TickDir','out',...
        'XTick',0:90:360,'YTick',[0,yMax(m)],'YTickLabel',[0,yMax(m)]);    
    title(['unit# ',num2str(n-38)]); 
    
    subplot(5, 2, m+2); 
    for k = 1:5
        numRepeat = size(eval(['cellData_NPX_EqLum(n).green_lum',num2str(k)]),1); 
        mean_dir = nanmean(eval(['cellData_NPX_EqLum(n).green_lum',num2str(k)])); 
        ste_dir = nanstd(eval(['cellData_NPX_EqLum(n).green_lum',num2str(k)])) / sqrt(numRepeat); 
        errorbar(0:45:360, mean_dir([1:8,1]), ste_dir([1:8,1]),...
                 'Color',[25 grad_now(k) 35]/255,'LineWidth',1);   hold on; 
    end
    plot([0,360],[baseline baseline],'k:'); 
    set(gca,'Box','off','Xlim',[-10,360],'Ylim',[0,yMax(m)],'TickDir','out',...
        'XTick',0:90:360,'YTick',[0,yMax(m)],'YTickLabel',[0,yMax(m)]);    

    subplot(5, 2, m+4); 
    for k = 1:5
        numRepeat = size(eval(['cellData_NPX_EqLum(n).red_lum',num2str(k)]),1); 
        mean_dir = nanmean(eval(['cellData_NPX_EqLum(n).red_lum',num2str(k)])); 
        ste_dir = nanstd(eval(['cellData_NPX_EqLum(n).red_lum',num2str(k)])) / sqrt(numRepeat); 
        errorbar(0:45:360, mean_dir([1:8,1]), ste_dir([1:8,1]),...
                 'Color',[grad_now(k) 25 35]/255,'LineWidth',1);   hold on; 
    end
    plot([0,360],[baseline baseline],'k:'); 
    set(gca,'Box','off','Xlim',[-10,360],'Ylim',[0,yMax(m)],'TickDir','out',...
        'XTick',0:90:360,'YTick',[0,yMax(m)],'YTickLabel',[0,yMax(m)]);    

    subplot(5, 2, m+6); 
    for k = 1:5
        numRepeat = size(eval(['cellData_NPX_EqLum(n).gray_lum',num2str(k)]),1); 
        mean_dir = nanmean(eval(['cellData_NPX_EqLum(n).gray_lum',num2str(k)])); 
        ste_dir = nanstd(eval(['cellData_NPX_EqLum(n).gray_lum',num2str(k)])) / sqrt(numRepeat); 
        errorbar(0:45:360, mean_dir([1:8,1]), ste_dir([1:8,1]),...
                 'Color',[grays(m) grays(m) grays(m)]/255,'LineWidth',1);   hold on; 
    end
    plot([0,360],[baseline baseline],'k:'); 
    set(gca,'Box','off','Xlim',[-10,360],'Ylim',[0,yMax(m)],'TickDir','out',...
        'XTick',0:90:360,'YTick',[0,yMax(m)],'YTickLabel',[0,yMax(m)]);    

    subplot(5, 2, m+8);     
    numRepeat = size(cellData_NPX_EqLum(n).noise,1); 
    mean_dir = nanmean(cellData_NPX_EqLum(n).noise); 
    ste_dir = nanstd(cellData_NPX_EqLum(n).noise) / sqrt(numRepeat); 
    errorbar(0:45:360, mean_dir([1:8,1]), ste_dir([1:8,1]),...
             'Color',[grays(m) grays(m) grays(m)]/255,'LineWidth',1);   hold on; 
    plot([0,360],[baseline baseline],'k:'); 
    set(gca,'Box','off','Xlim',[-10,360],'Ylim',[0,yMax(m)],'TickDir','out',...
        'XTick',0:90:360,'YTick',[0,yMax(m)],'YTickLabel',[0,yMax(m)]);    

    if m==1
        xlabel('Direction of motion (deg)');
        ylabel('Response (spk/s)');
    end
end    

clearvars -except cellData_NPX_*; 