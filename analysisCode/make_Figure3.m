%% Code to make_Figure3

% C. direction tuning curves from 5 neurons, 5 stim types
%    Neurons: 31+45, 31+35, 31+34, 31+19, 31+17, 31+16
%    Stim Types: RF/3, RF/6, RF/12, RF/36, Drift (RF/160)

% D. DI_base histogram for 5 stim types
%    Stim Types: RF/3, RF/6, RF/12, RF/36, Drift (RF/160)


load ../dataFiles/cellData_NPX_dX.mat; 

%% C. direction tuning curves
% direction tuning curves from 5 neurons, 5 stim types
% Neurons: 31+45, 31+35, 31+34, 31+19, 31+17, 31+16
% Stim Types: RF/3, RF/6, RF/12, RF/36, Drift (RF/160)

figure; 
set(gcf,'Position',[100 100 1000 700]);
sgtitle('Figure 3C'); 

cNums = [76, 66, 65, 50, 48, 47]; 
for n=1:length(cNums)

    cNum = cNums(n); 

    % no stim response 
    baseline = mean(cellData_NPX_dX(cNum).respMtx(:,41)); 

    % dir_tuning_mean
    dir_mean = reshape(mean(cellData_NPX_dX(cNum).respMtx(:,1:40)),[8,5]); 

    % dir_tuning_ste
    dir_ste = reshape(std(cellData_NPX_dX(cNum).respMtx(:,1:40))/sqrt(20),[8,5]); 

    % stim conditions: five stim types
    condNames = [{'RF/3'},{'RF/6'},{'RF/12'},...
                 {'RF/36'},{'Local (RF/160)'}]; 

    for m=1:5   % 5 stim types
        subplot(6,5,m+(n-1)*5);   

        dataNow = dir_mean(:,m); 
        yMax = max(max(dir_mean)); 
        yMin = min(min(dir_mean));       

        steNow = dir_ste(:,m); 

        plot(0:45:360,[dataNow; dataNow(1)],'.-','LineWidth',1);   hold on; 
        errorbar(0:45:360,[dataNow; dataNow(1)],[steNow; steNow(1)]);   
        if baseline>yMin-max(max(dir_ste))
            plot([0 360],[baseline baseline],'--', 'LineWidth',1); 
        end
        set(gca,'ylim',[yMin-max(max(dir_ste)) yMax+max(max(dir_ste))],...
            'xlim',[-10, 370],'xtick',0:90:360,'box','off','TickDir','out'); 
        title([num2str(cNum),': ',condNames{m}]); 
    end
end



%% D. make DI histograms
DI_base_mtx = []; 
anova_mtx = []; 
for n=1:length(cellData_NPX_dX)
    DI_base_mtx = [DI_base_mtx; cellData_NPX_dX(n).DI_base];
    anova_mtx = [anova_mtx; cellData_NPX_dX(n).anova(2,:)]; 
end

figure; 
set(gcf,'Position',[100 100 1000 140]);
sgtitle('Figure 3D'); 

for n=1:5
    sig_anova = find(anova_mtx(:,n)<0.05); 
    sig_DI_base = find(anova_mtx(:,n)<0.05 & DI_base_mtx(:,n)>=0.5);     
    nonsig_anova = find(anova_mtx(:,n)>=0.05); 

    subplot(1,5,n);
    histogram(DI_base_mtx(:,n),0:0.1:2.0,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1],'FaceAlpha',1);   hold on;
    histogram(DI_base_mtx(sig_anova,n),0:0.1:2.0,'FaceColor',[0,0,0],'EdgeColor',[1,1,1],'FaceAlpha',1);  
    title(['#sig: ',num2str(length(sig_anova)), ...
           '(',num2str(length(sig_DI_base)),'), #non-sig: ',num2str(length(nonsig_anova))]); 
    plot(median(DI_base_mtx(:,n)), 14, 'kv', 'markerfacecolor', 'k');   
    disp(['median of all DI: ',num2str(round(median(DI_base_mtx(:,n)),2))]);    
    xlabel('Direction Index (DI)');     
    ylabel('# Neurons');    
    set(gca,'box','off','TickDir','out','Ylim',[0,15],'Xlim',[0,2],'XTick',0:0.5:2); 
end

clearvars -except cellData_NPX_dX; 