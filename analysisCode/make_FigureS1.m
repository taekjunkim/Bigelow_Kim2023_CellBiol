%% Code to make_FigureS1


%% Figure S1A-C
%
% A. raster plot: DI (Local) vs. DI (LRM noise)
% B. raster plot: DI (Local) vs. DI (LRM sinusoid)
% C. corr btw LRM-sameDrift and LRM-antiDrift 
%             as a function of Drift resposne
%             as a function of Global resposne

load ../dataFiles/cellData_sua.mat; 

DI_base_mtx = []; 
corr_same_anti = []; 
LRM_Local_resp = [];    % maximum responses from LRM-sinusoid, Local

condNames = [{'LRM-noise'},{'LRM-sinusoid'},{'Local'},...
             {'LRM-sinu-Local-same'},{'LRM-sinu-Local-opp'}]; 

for n=1:length(cellData_sua)
    DI_base_mtx = [DI_base_mtx; cellData_sua(n).DI_base];

    [r,p] = corr(cellData_sua(n).dir_tuning(:,4),cellData_sua(n).dir_tuning(:,5));         
    corr_same_anti = [corr_same_anti; r p n]; 

    dir_tuning = cellData_sua(n).dir_tuning; 
    LRM_Local_resp = [LRM_Local_resp; nanmax(dir_tuning(:,[2,3]),[],1)]; 
end

figure; 
set(gcf,'Position',[100 100 800 400]);
sgtitle('Figure S1A-C');                 

%%% Figure S1A
subplot(2,3,2); 
plot(DI_base_mtx(:,3),DI_base_mtx(:,1),'ko','MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',5);   hold on;
[r,p] = corr(DI_base_mtx(:,3),DI_base_mtx(:,1)); 
plot([0.5,0.5],[0,2],'r:'); 
plot([0,2],[0.5,0.5],'r:'); 
title(['r: ',num2str(round(r,3)),', p:',num2str(round(p,3))]); 
xlabel('DI (Drift)'); 
ylabel('DI (LRM noise)'); 
set(gca,'box','off','TickDir','out','ylim',[0,2],'xlim',[0,2]); 

%%% Figure S1B
subplot(2,3,3); 
plot(DI_base_mtx(:,3),DI_base_mtx(:,2),'ko','MarkerFaceColor',[0.5,0.5,0.5],'MarkerSize',5);   hold on; 
[r,p] = corr(DI_base_mtx(:,3),DI_base_mtx(:,2)); 
plot([0.5,0.5],[0,2],'r:'); 
plot([0,2],[0.5,0.5],'r:'); 
title(['r: ',num2str(round(r,3)),', p:',num2str(round(p,3))]); 
xlabel('DI (Drift)'); 
ylabel('DI (LRM sinusoid)'); 
set(gca,'box','off','TickDir','out','ylim',[0,2],'xlim',[0,2]); 


%%% Figure S1C
subplot(2,3,4); 
histogram(corr_same_anti(:,1),-0.5:0.05:1,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on
signi = find(corr_same_anti(:,2)<0.05); 
histogram(corr_same_anti(signi,1),-0.5:0.05:1,'FaceColor',[0,0,0],'EdgeColor',[1,1,1],'FaceAlpha',1);
plot(median(corr_same_anti(:,1)),9,'k>','MarkerFaceColor','k')
ylabel('# Units'); 
title([num2str(length(signi)),'/115'])
set(gca,'box','off','TickDir','out','Xlim',[-0.51,1.01]); 
set(gca,'view',[-90 90],'XAxisLocation','top'); 

subplot(2,3,5); 
plot(log10(LRM_Local_resp(:,1)),corr_same_anti(:,1),'ko','MarkerFaceColor',[0.5,0.5,0.5]);   hold on;
plot(log10(LRM_Local_resp(signi,1)),corr_same_anti(signi,1),'ko','MarkerFaceColor','k'); 
set(gca,'box','off','TickDir','out','XTick',[0,1,2],'XTickLabel',[0,10,100],'Xlim',[0,2.1]); 
xlabel('Global Response (spikes/s)'); 
ylabel('r btw same, anti'); 
subplot(2,3,6); 
plot(log10(LRM_Local_resp(:,2)),corr_same_anti(:,1),'ko','MarkerFaceColor',[0.5,0.5,0.5]);   hold on; 
plot(log10(LRM_Local_resp(signi,2)),corr_same_anti(signi,1),'ko','MarkerFaceColor','k'); 
set(gca,'box','off','TickDir','out','XTick',[0,1,2],'XTickLabel',[0,10,100],'Xlim',[0,2.1]); 
xlabel('Local Response (spikes/s)'); 
ylabel('r btw same, anti'); 


%% Figure S1D. DI histogram for monkeys - M1 (Left), M2 (Right)
% LRM-noise
% LRM-sinusoid
% Local

DI_base_mtx = []; 
anova_mtx = []; 
monkey_num = []; 

condNames = [{'LRM-noise'},{'LRM-sinusoid'},{'Local'},...
             {'LRM-sinu-Local-same'},{'LRM-sinu-Local-opp'}]; 

for n=1:length(cellData_sua)
    DI_base_mtx = [DI_base_mtx; cellData_sua(n).DI_base];

    anova_mtx = [anova_mtx; cellData_sua(n).anova(2,:)]; 
    if cellData_sua(n).file_id(1)=='z'
        monkey_num = [monkey_num; 1];
    else
        monkey_num = [monkey_num; 2];
    end
end

figure; 
set(gcf,'Position',[100 1000 400 400]); 
sgtitle('Figure S1D');              

for n=1:3
    monkey1 = find(monkey_num(:)==1); 
    monkey2 = find(monkey_num(:)==2); 
    
    % monkey1
    subplot(3,2,(n-1)*2+1);
    histogram(DI_base_mtx(monkey1,n),'BinEdges',0:0.1:2,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;

    sig_anova = find(anova_mtx(monkey1,n)<0.05); 
    sig_DI_base = find(anova_mtx(monkey1,n)<0.05 & DI_base_mtx(monkey1,n)>=0.5);     
    nonsig_anova = find(anova_mtx(monkey1,n)>=0.05); 

    histogram(DI_base_mtx(monkey1(sig_anova),n),'BinEdges',0:0.1:2,'FaceColor',[0,0,0],'EdgeColor',[1,1,1],'FaceAlpha',1);  
    plot(median(DI_base_mtx(monkey1,n)), 13, 'kv', 'markerfacecolor', 'k');    
    disp(['sig. median = ',num2str(round(median(DI_base_mtx(monkey1(sig_anova),n)),2))]); 
    disp(['all median = ',num2str(round(median(DI_base_mtx(monkey1,n)),2))]);     
    title(['#sig: ',num2str(length(sig_anova)), ...
           '(',num2str(length(sig_DI_base)),'), #non-sig: ',num2str(length(nonsig_anova))]); 
    xlabel('Direction Index (DI)');     
    ylabel('# Neurons')   
    set(gca,'box','off','TickDir','out','Ylim',[0 15],'Xlim',[0 2]); 

    % monkey2    
    subplot(3,2,n*2);
    histogram(DI_base_mtx(monkey2,n),'BinEdges',0:0.1:2,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;

    sig_anova = find(anova_mtx(monkey2,n)<0.05); 
    sig_DI_base = find(anova_mtx(monkey2,n)<0.05 & DI_base_mtx(monkey2,n)>=0.5);     
    nonsig_anova = find(anova_mtx(monkey2,n)>=0.05); 

    histogram(DI_base_mtx(monkey2(sig_anova),n),'BinEdges',0:0.1:2,'FaceColor',[0,0,0],'EdgeColor',[1,1,1],'FaceAlpha',1);  
    plot(median(DI_base_mtx(monkey2,n)), 13, 'kv', 'markerfacecolor', 'k');    
    disp(['sig. median = ',num2str(round(median(DI_base_mtx(monkey2(sig_anova),n)),2))]); 
    disp(['all median = ',num2str(round(median(DI_base_mtx(monkey2,n)),2))]); 
    title(['#sig: ',num2str(length(sig_anova)), ...
           '(',num2str(length(sig_DI_base)),'), #non-sig: ',num2str(length(nonsig_anova))]); 
    xlabel('Direction Index (DI)');     
    ylabel('# Neurons')   
    set(gca,'box','off','TickDir','out','Ylim',[0 15],'Xlim',[0 2]);     
end    

clearvars -except cellData_sua; 