%% Code to make_Figure2

% A. DI vs dX
% LRMnoise_DI vs. dX (RF/3)

% B. Drift dominance histogram
% to facilitate direct comparison with MT data
% compute global DI by averaging across all local motion conditions that 
% were paired with each global motion and vice versa for local motion DI
%
% Local DI (avg. of same/opp LRM) - LRM DI (avg. of same/opp Local)

load ../dataFiles/cellData_sua.mat; 


%% A. DI vs dX
% LRMnoise_DI vs. dX (RF/3)

DI_dX_mtx = []; 
% col1: dX
% col2: LRMnoise_DI_base; 
% col3: LRMsinusoid_DI_base; 
% col4: LRMnoise_ANOVA_p
% col5: LRMsinusoid_ANOVA_p

condNames = [{'LRM-noise'},{'LRM-sinusoid'},{'Local'},...
             {'LRM-sinu-Local-same'},{'LRM-sinu-Local-opp'}]; 
pix_per_deg = 37; 

DI_dX_mtx = []; 
for n=1:length(cellData_sua)
    RF_diameter = sqrt(cellData_sua(n).rf_pix(1)^2 + cellData_sua(n).rf_pix(2)^2)*0.625 + pix_per_deg; 
    dX = RF_diameter / 3; 

    DI_dX_mtx = [DI_dX_mtx; 
                 dX cellData_sua(n).DI_base(1) cellData_sua(n).anova(2,1)]; 
end

figure; 
set(gcf,'Position',[100 100 250 350]);
sgtitle('Figure 2AB'); 

subplot(4,1,[1 2]); 
plot(DI_dX_mtx(:,1)/pix_per_deg, DI_dX_mtx(:,2),'ko','MarkerSize',6);   hold on; 
sigDI1 = find(DI_dX_mtx(:,3)<0.05); 
sigDI1B = find((DI_dX_mtx(:,3)<0.05) & (DI_dX_mtx(:,2)>=0.5)); 
plot(DI_dX_mtx(sigDI1,1)/pix_per_deg, DI_dX_mtx(sigDI1,2),'ko','MarkerFaceColor','k','MarkerSize',6);   
plot([1,1],[0,1.5],'r:'); 
xlabel('dX (step size in degree)'); 
ylabel('Direction index (DI)'); 
title(['LRM noise: ',num2str(length(sigDI1)),'(',num2str(length(sigDI1B)),')']); 
set(gca,'box','off','TickDir','out','ylim',[0,1.5]); 



%% B. Drift dominance histogram
% to facilitate direct comparison with MT data
% compute global DI by averaging across all local motion conditions that 
% were paired with each global motion and vice versa for local motion DI

% compute avg. LRM DI_base, ANOVA avg.LRM 
for n=1:length(cellData_sua)

    respMtx = cellData_sua(n).respMtx; 

    % no stim response 
    baseline = nanmean(respMtx(:,end)); 

    % stim condition: five stim types
    condNames = [{'LRM-noise'},{'LRM-sinusoid'},{'Local'},...
                 {'LRM-sinu-Local-same'},{'LRM-sinu-Local-opp'}]; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LRM DI by averaging same/opp Local
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % for avg.LRM
    respMtx_for_avg = [respMtx(:,25:32); respMtx(:,33:40)]; 
   
    % Compute DI-base for avg.LRM
    dir_tuning = nanmean(respMtx_for_avg,1); 

    abs_modulation = abs(dir_tuning - baseline); 
    pref_id = find(abs_modulation(:)==max(abs_modulation));         
    mod_depth = []; 
    for p=1:length(pref_id) % when multiple points were found
        if pref_id(p) <= 4
            non_pref_id(p) = pref_id(p) + 4; 
        else
            non_pref_id(p) = pref_id(p) - 4; 
        end

        mod_depth(p) = abs(dir_tuning(pref_id(p)) - dir_tuning(non_pref_id(p))); 
    end
    pref_id = pref_id(find(mod_depth(:)==max(mod_depth))); 
    non_pref_id = non_pref_id(find(mod_depth(:)==max(mod_depth)));         

    numerator = abs(dir_tuning(pref_id) - dir_tuning(non_pref_id)); 
    denominator = abs(dir_tuning(pref_id)-baseline);         
    if(size(numerator/denominator)> 1)
        DI_base = numerator(1)/denominator(1);
    else
        DI_base = numerator / denominator; 
    end        
   
    % replace NaNs with appropriate values
    if isnan(DI_base)
        DI_base = 0; 
    end
    cellData_sua(n).DI_base_avgLRM = DI_base;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Local DI by averaging same/opp LRM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

    % for avg.Local
    respMtx_for_avgLocal = [respMtx(:,25:32); respMtx(:,[37 38 39 40 33 34 35 36])]; 
    
    % Compute DI-base for avg.Local
    dir_tuning = nanmean(respMtx_for_avgLocal,1); 

    abs_modulation = abs(dir_tuning - baseline); 
    pref_id = find(abs_modulation(:)==max(abs_modulation));         
    mod_depth = []; 
    for p=1:length(pref_id) % when multiple points were found
        if pref_id(p) <= 4
            non_pref_id(p) = pref_id(p) + 4; 
        else
            non_pref_id(p) = pref_id(p) - 4; 
        end

        mod_depth(p) = abs(dir_tuning(pref_id(p)) - dir_tuning(non_pref_id(p))); 
    end
    pref_id = pref_id(find(mod_depth(:)==max(mod_depth))); 
    non_pref_id = non_pref_id(find(mod_depth(:)==max(mod_depth)));         

    numerator = abs(dir_tuning(pref_id) - dir_tuning(non_pref_id)); 
    denominator = abs(dir_tuning(pref_id)-baseline);         
    if(size(numerator/denominator)> 1)
        DI_base = numerator(1)/denominator(1);
    else
        DI_base = numerator / denominator; 
    end        
   
    % replace NaNs with appropriate values
    if isnan(DI_base)
        DI_base = 0; 
    end
    cellData_sua(n).DI_base_avgLocal = DI_base;    

end
    

% accumulate values
avgLRM_DI_base = []; 
avgLocal_DI_base = []; 

for n=1:length(cellData_sua)
    avgLRM_DI_base = [avgLRM_DI_base; cellData_sua(n).DI_base_avgLRM];
    avgLocal_DI_base = [avgLocal_DI_base; cellData_sua(n).DI_base_avgLocal];         
end

% local Dominance
subplot(3,1,3); 
histogram(avgLocal_DI_base - avgLRM_DI_base, -2:0.1:2, 'FaceColor',[0.5,0.5,0.5],'EdgeColor',[0,0,0]); 
hold on; 
plot(median(avgLocal_DI_base - avgLRM_DI_base),24,'kv','MarkerFaceColor','k'); 
set(gca,'box','off','TickDir','out','Ylim',[0,25]); 
xlabel('Local dominance'); 
ylabel('Number of units'); 
title('V4 data'); 

clearvars -except cellData_sua; 
