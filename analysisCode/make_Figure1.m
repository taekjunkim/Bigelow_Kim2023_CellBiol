%% Code to make_Figure1

% direction tuning curves from 5 neurons, 5 stim types
% Neurons: P190520, Z181009, P190410, P190429, P190516
% Stim Types: LRM-noise, LRM-sinusoid, Local, 
%             LRM-sinusoid-Local-same, 
%             LRM-sinusoid-Local-opp

% DI_base histogram for 5 stim types


%% direction tuning curves
% direction tuning curves from 5 neurons, 5 stim types
% Neurons: P190520, Z181009, P190410, P190429, P190516
% Stim Types: LRM-noise, LRM-sinusoid, Local, 
%             LRM-sinusoid-Local-same, LRM-sinusoid-Local-opp

load ../dataFiles/cellData_sua.mat; 
% P190520 (89), Z181009 (45), P190410 (65), P190429 (80), P190516 (86)

figure; 
set(gcf,'Position',[100 100 1000 700]);
sgtitle('Figure 1B'); 

cNums = [89, 45, 65, 80, 86]; 
for n=1:length(cNums)

    cNum = cNums(n); 
    respMtx = cellData_sua(cNum).respMtx; 

    % no stim response 
    baseline = nanmean(respMtx(:,end)); 

    % stim condition: five stim types
    condNames = [{'LRM-noise'},{'LRM-sinusoid'},{'Local'},...
                 {'LRM-sinu-Local-same'},{'LRM-sinu-Local-opp'}]; 

    % file_ID
    file_ID = cellData_sua(cNum).file_id; 

    for m=1:5   % 5 stim types
        
        if (n<4) & (m>=4)
            continue;
        end

        subplot(5,5,m+(n-1)*5);   
        condNum = (m-1)*8+1:m*8;

        dataNow = nanmean(respMtx(:,condNum),1); 
        yMax = max(nanmean(respMtx(:,1:40),1)); 
        yMin = min(nanmean(respMtx(:,1:40),1)); 
        steAll = nanstd(respMtx(:,:),1)/sqrt(size(respMtx,1));
        steNow = nanstd(respMtx(:,condNum),1)/sqrt(size(respMtx,1)); 

        plot(0:45:360,[dataNow dataNow(1)],'.-','LineWidth',1);   hold on; 
        errorbar(0:45:360,[dataNow dataNow(1)],[steNow steNow(1)]);   
        if baseline>yMin-max(steAll)
            plot([0 360],[baseline baseline],'--', 'LineWidth',1); 
        end
        set(gca,'ylim',[yMin-max(steAll) yMax+max(steAll)],...
            'xlim',[-10, 370],'xtick',0:90:360,'box','off','TickDir','out'); 
        if m==1
            title([file_ID,': ',condNames{m}]); 
        else
            title([condNames{m}]); 
        end
    end
end
clearvars -except cellData_sua; 


%% make DI histograms
load ../dataFiles/cellData_sua.mat; 

DI_base_mtx = []; 
anova_mtx = []; 

condNames = [{'LRM-noise'},{'LRM-sinusoid'},{'Local'},...
             {'LRM-sinu-Local-same'},{'LRM-sinu-Local-opp'}]; 

for n=1:length(cellData_sua)
    DI_base_mtx = [DI_base_mtx; cellData_sua(n).DI_base n];
    anova_mtx = [anova_mtx; cellData_sua(n).anova(2,:) n]; 
end

figure; 
set(gcf,'Position',[100 100 200 500]); 
sgtitle('Figure 1C'); 

for n=1:3
    sig_anova = find(anova_mtx(:,n)<0.05); 
    sig_DI_base = find(anova_mtx(:,n)<0.05 & DI_base_mtx(:,n)>=0.5);     
    nonsig_anova = find(anova_mtx(:,n)>=0.05); 

    subplot(3,1,n); 
    histogram(DI_base_mtx(:,n),'BinEdges',0:0.1:2.0,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;
    histogram(DI_base_mtx(sig_anova,n),'BinEdges',0:0.1:2.0,'FaceColor',[0,0,0],'EdgeColor',[1,1,1],'FaceAlpha',1);  
    plot(median(DI_base_mtx(:,n)), 21, 'kv', 'markerfacecolor', 'k');    
    disp(['median of all DI: ',num2str(round(median(DI_base_mtx(:,n)),2))]); 
    title(['#sig: ',num2str(length(sig_anova)), ...
           '(',num2str(length(sig_DI_base)),'), #non-sig: ',num2str(length(nonsig_anova))]); 
    xlabel('Direction Index (DI)');     
    ylabel('# Neurons')   
    set(gca,'box','off','TickDir','out','Ylim',[0 25],'Xlim',[0 2.0]); 
    
end    

clearvars -except cellData_sua; 

