
clear;
LoadMyColor;
LoadDefaultSettings;
% set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);
set(0,'DefaultAxesFontSize',13,'DefaultTextFontSize',13);
% set(0,'DefaultAxesFontSize',15,'DefaultTextFontSize',15);

%% Set working directory

base_path = fullfile('..'); % If your current working directory is 'demo', use this. Otherwise, set your own path.
addpath(fullfile(base_path, 'functions', 'utils'))

% If your current working directory is 'demo', use this. Otherwise, set your own path that includes 'data'.
cd('..'); 

%% setting

% ----- vs behavior -----

% !!! for shattering analysis !!!
% bhv_sfx = [];
bhv_sfx = ' (G)';   % G trials
% bhv_sfx = ' (H)';   % H trials
% bhv_sfx = ' (G Early)';   % G trials, half early
% bhv_sfx = ' (G Late)';   % G trials, half late
% bhv_sfx = ' (G uL-uH)';   % G trials, uLow - uHigh
% bhv_sfx = ' (G uL-uE)';   % G trials, uLate - uEarly
% bhv_sfx = ' (G L-E)';   % G trials, late - early

% !!! for PCA !!!
% bhv_sfx = [];

% BEHAV = {'ChoSwitch1n2(GoalSwitch)', 'ChoConsist1n2', 'ChoOptBin1n2'}; bhv_horz = 1;
BEHAV = {'ChoSwitch1n2(GoalSwitch)', 'ChoConsist1n2', 'ChoOptBinNaN1n2'}; bhv_horz = 1;
% % BEHAV = {'ChoiceSwitch1n2(GoalSwitch)', 'ChoConsist1n2', 'ChoOptBin1n2'}; bhv_horz = 3;

% BEHAV = {'ChoOptBin1', 'ChoiceSwitch1(GoalSwitch)', 'ChoConsist1'}; bhv_horz = 3;
% BEHAV = {'ChoOptBin2', 'ChoiceSwitch2(GoalSwitch)', 'ChoConsist2'}; bhv_horz = 3;

% BEHAV = {'ChoOptBin1n2', 'ChoiceSwitch1n2(GoalSwitch)', 'ChoConsist1n2', 'PMB'};

% % % BEHAV = {'ChoOptBin1', 'ChoOptBin2', 'ChoOptBin1n2', ...
% % %     'ChoiceSwitch1(GoalSwitch)', 'ChoiceSwitch2(GoalSwitch)', 'ChoiceSwitch1n2(GoalSwitch)', ...
% % %     'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', ...
% % %     'R', 'PMB'}; bhv_horz = 3;
% % 
% % % BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'R' ...
% % %     'ChoiceSwitch1(GoalSwitch)', 'ChoiceSwitch2(GoalSwitch)', 'ChoiceSwitch1n2(GoalSwitch)', 'PMB', ...
% % %     'ChoiceSwitch1(GoalStay)', 'ChoiceSwitch2(GoalStay)', 'ChoiceSwitch1n2(GoalStay)', 'lr', ...
% % %     'ChoiceSwitch1(GoalSw-St)', 'ChoiceSwitch2(GoalSw-St)', 'ChoiceSwitch1n2(GoalSw-St)', 'invTemp'};

% ----- Goal-directed behaviors + important RL parameters -----
% BEHAV = {'ChoOptBin1', 'ChoOptBin2', 'ChoOptBin1n2', 'R' ...
%     'ChoiceSwitch1(GoalSwitch)', 'ChoiceSwitch2(GoalSwitch)', 'ChoiceSwitch1n2(GoalSwitch)', 'PMB', ...
%     'ChoiceSwitch1(GoalStay)', 'ChoiceSwitch2(GoalStay)', 'ChoiceSwitch1n2(GoalStay)', 'lr', ...
%     'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'invTemp'}; bhv_horz = 4;

% % ----- UC effect on choice behavior -----
% BEHAV = {'UncCond--A1(anovan)', 'UncCond--A2 state2(anovan)', ... 
%             'UncCond--A2 state3(anovan)', 'UncCond--A2 state4(anovan)', 'UncCond--A2 state5(anovan)'}; bhv_horz = 5;

% ----- RL parameters & variables (per event) -----
% BEHAV = {'ThrSPE', 'lrRPE', 'MB2MF', 'MF2MB', 'invTemp', 'lr', ...
%     'SPE2', 'SPE3', 'SPE', 'RPE2', 'RPE3', 'RPE', ...
%     'relMB2', 'relMB3', 'relMB', 'relMF2', 'relMF3', 'relMF', ...
%     'relMAX2', 'relMAX3', 'relMAX', 'relDiff2', 'relDiff3', 'relDiff'}; bhv_horz = 6;

% % ----- RL parameters & variables (event average) -----
% BEHAV = {'ThrSPE', 'lrRPE', 'MB2MF', 'MF2MB', 'invTemp', 'lr', ...
%     'SPE', 'RPE', 'relMB', 'relMF', 'relMAX', 'relDiff'}; bhv_horz = 6;


% BEHAV = {'ChoOpt1', 'ChoOpt2', 'ChoOpt1n2', 'R' ...
%     'ChoiceSwitch1(GoalSwitch)', 'ChoiceSwitch2(GoalSwitch)', 'ChoiceSwitch1n2(GoalSwitch)', 'PMB', ...
%     'ChoiceSwitch1(GoalStay)', 'ChoiceSwitch2(GoalStay)', 'ChoiceSwitch1n2(GoalStay)', 'lr', ...
%     'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'invTemp'}; bhv_horz = 4;

% BEHAV = {'ChoOptBin1', 'ChoOptBin2', 'ChoOptBin1n2', 'R' ...
%     'ChoiceSwitch1(GoalSwitch)', 'ChoiceSwitch2(GoalSwitch)', 'ChoiceSwitch1n2(GoalSwitch)', 'PMB', ...
%     'ChoConsist1', 'ChoConsist2', 'ChoConsist1n2', 'invTemp'};

% ----- vs behavior -----
% NEUR = {'Goal(6,7,8) LOROV_Separability dimension'}; BEHAV = NEUR;

% NEUR = {'Goal(6,7,8) LOROV_Separability dimension', 'Goal(6,7,8) LOROV_CCGP across UncCond', ... 
%     'Goal(6,7,8) uL LOROV vs Goal(6,7,8) uH LOROV_Cosine similarity'}; BEHAV = NEUR; bhv_horz = 3;

% NEUR = {'UncCond G LOROV_Separability dimension'}; BEHAV = NEUR; bhv_horz = 1;


% ----- Supervised dimensionality -----
label_name = 'Goal(6,7,8)'; temp_color = [125 125 125]/255; % repeat 100 CCGP
% label_name = 'Goal(6,7,8) LOROV';
% label_name = 'Goal(6,7,8) LOROV HalfEarly';
% label_name = 'Goal(6,7,8) LOROV HalfLate';
% label_name = 'Goal(6,7,8) HalfLate - Goal(6,7,8) HalfEarly';
% label_name = 'S2 LOROV G';
% label_name = 'S3 LOROV G';
% label_name = 'S3(6,7,8) LOROV G';
% label_name = 'R LOROV G';
% label_name = 'R(10,20,40) LOROV G';
% label_name = 'Goal(6,7,8) x Unc LOROV';
% label_name = 'Goal(6,7,8) x Unc LOROV goal'; temp_color = [77 161 169]/255;
% label_name = 'Goal(6,7,8) x Unc LOROV unc'; temp_color = [226 16 103]/255;
% label_name = 'Goal(6,7,8) x Unc LOROV linear'; temp_color = [155 80 148]/255;
% label_name = 'Goal(6,7,8) x Unc LOROV nonlinear'; temp_color = [94 35 157]/255;                    
% label_name = 'Goal(6,7,8)xUC LOROV goal'; temp_color = [77 161 169]/255; % repeat 100 SD
% label_name = 'Goal(6,7,8)xUC LOROV unc'; temp_color = [226 16 103]/255; % repeat 100 SD
% label_name = 'Goal(6,7,8)xUC LOROV linear'; temp_color = [155 80 148]/255; % repeat 100 SD
% label_name = 'Goal(6,7,8)xUC LOROV nonlinear'; temp_color = [94 35 157]/255; % repeat 100 SD      
% label_name = 'Goal(6,7,8)xUC LOROV HalfEarly goal';
% label_name = 'Goal(6,7,8)xUC LOROV HalfEarly unc';
% label_name = 'Goal(6,7,8)xUC LOROV HalfEarly linear';
% label_name = 'Goal(6,7,8)xUC LOROV HalfEarly nonlinear';
% label_name = 'Goal(6,7,8)xUC LOROV HalfLate goal';
% label_name = 'Goal(6,7,8)xUC LOROV HalfLate unc';
% label_name = 'Goal(6,7,8)xUC LOROV HalfLate linear';
% label_name = 'Goal(6,7,8)xUC LOROV HalfLate nonlinear';
% label_name = 'Goal(6,7,8)xUC HalfLate goal - Goal(6,7,8)xUC HalfEarly goal';
% label_name = 'Goal(6,7,8)xUC HalfLate unc - Goal(6,7,8)xUC HalfEarly unc';
% label_name = 'Goal(6,7,8)xUC HalfLate linear - Goal(6,7,8)xUC HalfEarly linear';
% label_name = 'Goal(6,7,8)xUC HalfLate nonlinear - Goal(6,7,8)xUC HalfEarly nonlinear';
% label_name = 'UncCond G LOROV';
% label_name = 'UncCond H LOROV';
% label_name = 'UncCond at state 4 H';
% label_name = 'binSPE2 G LOROV';
% label_name = 'binSPE3 G LOROV';
% label_name = 'binRPE2 G LOROV';
% label_name = 'binRPE3 G LOROV';
% label_name = 'relComp2 G LOROV';
% label_name = 'relComp3 G LOROV';
% ----- Unsupervised dimensionality -----
% label_name = 'GoalCond1';
% label_name = 'GoalCond2';
% label_name = 'Goal(6,7,8)Gall';
% ----- Cosine similarity -----
% label_name = 'Goal(6,7,8) uL LOROV vs Goal(6,7,8) uH LOROV'; temp_color = [125 125 125]/255; % repeat 100 CS
% % label_name = 'Goal(6,7,8) uL vs Goal(6,7,8) uH'; 

ROI = {'vlPFC', 'dlPFC', 'OFC', 'ACC', 'preSMA', 'V1', 'HPC', 'vStr'}; LR_integration = 1;
% ROI = {'vlPFC', 'dlPFC', 'OFC', 'V1', 'HPC', 'vStr'}; LR_integration = 1;
% ROI = {'F3TL', 'F3TR', 'DLPFCL', 'DLPFCR', 'V1L', 'V1R', 'HIPPOL', 'HIPPOR', 'VentStrL', 'VentStrR', 'OFC'}; LR_integration = 0;
% ROI = {'vlPFC', 'dlPFC', 'OFC', 'V1'}; LR_integration = 1; % +DLPFC & R/L integrated OFC
N_ROI = length(ROI);

% EVENT = [1 2 4 5 7 8 9 10]; EVENT_NAME = {  'fix ', 'S1 ', 'A1 ', 'S2 ', 'A2 ', 'S3 ', 'fix` ', 'S1` '}; strict_pca = 0; % f2 = A1, f3 = A2
% EVENT = [2 4 5 7 8]; EVENT_NAME = {'S1 ', 'A1 ', 'S2 ', 'A2 ', 'S3 '}; strict_pca = 0; % f2 = A1, f3 = A2
% EVENT = 24578; EVENT_NAME = {  'S1 A1 S2 A2 S3' };
EVENT = 2457; EVENT_NAME = {  'S1 A1 S2 A2' };
% EVENT = 24; EVENT_NAME = {  'S1 A1' };
% EVENT = 57; EVENT_NAME = {  'S2 A2' };
% EVENT = 12345678; EVENT_NAME = {  'Entire events' }; strict_pca = 1;

% corr_behav_dir = 'corr_behav_original'; neural_measure = 'Separability dimension';
corr_behav_dir = 'corr_behav_original'; neural_measure = 'CCGP across UncCond';
% corr_behav_dir = 'corr_behav_original'; neural_measure = 'Cosine similarity';
% corr_behav_dir = 'corr_behav_original'; neural_measure = 'SD diff';
% corr_behav_dir = 'corr_behav_original'; neural_measure = 'CCGP diff';
% corr_behav_dir = 'corr_behav_original'; neural_measure = 'Generalization index across UncCond';
% corr_behav_dir = 'corr_behav_original'; neural_measure = 'Ratio(CCGP across UncCond, SD)';
% corr_behav_dir = 'corr_behav_ROIeffects'; neural_measure = 'Separability dimension';
% corr_behav_dir = 'corr_behav_dichotomy'; neural_measure = 'Separability dimension';
% corr_behav_dir = 'corr_behav_unsupervised'; neural_measure = 'PR'; norm_pca = 'basic'; shuffle_pca = [];
% corr_behav_dir = 'corr_behav_unsupervised'; neural_measure = 'ER'; norm_pca = 'basic'; shuffle_pca = [];
% corr_behav_dir = 'corr_behav_unsupervised'; neural_measure = 'nPC'; norm_pca = 'basic'; shuffle_pca = [];
% corr_behav_dir = 'corr_supervised_unsupervised'; neural_measure = 'PR'; norm_pca = 'basic'; shuffle_pca = [];
% corr_behav_dir = 'corr_supervised_unsupervised'; neural_measure = 'ER'; norm_pca = 'basic'; shuffle_pca = [];
% corr_behav_dir = 'corr_supervised_unsupervised'; neural_measure = 'nPC'; norm_pca = 'basic'; shuffle_pca = [];
% corr_behav_dir = 'corr_supervised_unsupervised'; neural_measure = 'pat_size'; norm_pca = 'basic'; shuffle_pca = [];
% corr_behav_dir = 'corr_supervised_unsupervised'; neural_measure = 'N_eigs'; norm_pca = 'basic'; shuffle_pca = [];

% fdr = 0; % no multiple testing correction
fdr = 'Benjamini-Hochberg';


if ~LR_integration
    colors = reshape(cat(1, colors, cellfun(@(x) (x .* 0.8), colors, 'UniformOutput', false)), 1, []);
end

%% loading

disp('BHV loading')

switch corr_behav_dir
    case 'corr_behav_original'
%         coef_cell = cell(length(ROI), length(BEHAV));
%         pval_cell = cell(length(ROI), length(BEHAV));
        coef_mat = nan(length(ROI), length(BEHAV), length(EVENT));
        pval_mat = nan(length(ROI), length(BEHAV), length(EVENT));
        
        for roii = 1:length(ROI)
            for bhvi = 1:length(BEHAV)
                for evi = 1:length(EVENT)
                    
                    exp_name = 'Lee2014';
                    % exp_name = 'Lee2014pp2';
                    file_name = ['data\ccgp_result\' corr_behav_dir '\' num2str(EVENT(evi)) '\' ...
                        exp_name '_' neural_measure '_' label_name '_' ROI{roii} '_' BEHAV{bhvi} bhv_sfx '.mat'];
                    
                    load(file_name, 'R')
                    load(file_name, 'p')
                    
                    coef_mat(roii, bhvi, evi) = R(1,2);
                    pval_mat(roii, bhvi, evi) = p;
                    % try
                    %     pval_mat(roii, bhvi, evi) = p(1,2);
                    % catch
                    %     pval_mat(roii, bhvi, evi) = p;
                    % end
                    
                end
                
            end
            fprintf('%d ', roii)
        end
        
    case 'corr_behav_ROIeffects'
        coef_cell = cell(length(BEHAV), 1);
        pval_cell = cell(length(BEHAV), 1);
        coef_mat = nan(length(EVENT), length(ROI));
        pval_mat = nan(length(EVENT), length(ROI));
        
        for bhvi = 1:length(BEHAV)
            for evi = 1:length(EVENT)
                
                exp_name = 'Lee2014';
                behavior_measure = BEHAV{bhvi};
                file_name = ['data\ccgp_result\' corr_behav_dir '\' num2str(EVENT(evi)) '\' ...
                    exp_name '_' neural_measure '_' label_name '_' behavior_measure bhv_sfx '.mat'];
                
                load(file_name, 'mdl')
                
                coef_mat(evi, :) = mdl.Coefficients.Estimate(2:end)';
                pval_mat(evi, :) = mdl.Coefficients.pValue(2:end)';
                
            end
            
            coef_cell{bhvi} = coef_mat;
            pval_cell{bhvi} = pval_mat;
            fprintf('%d ', bhvi)
        end
        
    case 'corr_behav_dichotomy'
        coef_cell = cell(length(ROI), length(BEHAV));
        pval_cell = cell(length(ROI), length(BEHAV));
        coef_mat = nan(length(EVENT), 4);
        pval_mat = nan(length(EVENT), 4);
%         coef_mat = nan(length(EVENT), 31);
%         pval_mat = nan(length(EVENT), 31);
        
        for roii = 1:length(ROI)
            for bhvi = 1:length(BEHAV)
                for evi = 1:length(EVENT)
                    
                    exp_name = 'Lee2014';
                    file_name = ['data\ccgp_result\' corr_behav_dir '\' num2str(EVENT(evi)) '\' ...
                        exp_name '_' neural_measure '_' label_name '_' ROI{roii} '_' BEHAV{bhvi} bhv_sfx '.mat'];
                    
                    load(file_name, 'mdl')
                    
                    coef_mat(evi, :) = mdl.Coefficients.Estimate(2:end)';
                    pval_mat(evi, :) = mdl.Coefficients.pValue(2:end)';
                    
                end
                
                coef_cell{roii, bhvi} = coef_mat;
                pval_cell{roii, bhvi} = pval_mat;
            end
            fprintf('%d ', roii)
        end
        
    case {'corr_behav_unsupervised', 'corr_supervised_unsupervised'}
%         coef_cell = cell(length(ROI), length(BEHAV));
%         pval_cell = cell(length(ROI), length(BEHAV));
        coef_mat = nan(length(ROI), length(BEHAV), length(EVENT));
        pval_mat = nan(length(ROI), length(BEHAV), length(EVENT));
        
        for roii = 1:length(ROI)
            for bhvi = 1:length(BEHAV)
                for evi = 1:length(EVENT)
                    
                    exp_name = 'Lee2014';                    
                    file_name = ['boldpat_dimension_result\seed2021\' corr_behav_dir '\' ...
                        ['strict' num2str(strict_pca) ' ' norm_pca shuffle_pca] '\' ...
                        num2str(EVENT(evi)) '\' ...
                        exp_name '_' neural_measure '_' label_name '_' ROI{roii} '_' BEHAV{bhvi} bhv_sfx '.mat'];
                    
                    load(file_name, 'R')
                    load(file_name, 'p')
                    
                    coef_mat(roii, bhvi, evi) = R(1,2);
                    pval_mat(roii, bhvi, evi) = p(1,2);
                    
                end
                
            end
            fprintf('%d ', roii)
        end
        
end
disp(' done')

%% visualization

close all

switch corr_behav_dir
    case {'corr_behav_original', 'corr_behav_unsupervised', 'corr_supervised_unsupervised'}
        
%         for roii = 1:length(ROI)
%             figure('Name', ROI{roii});
            %
            % figure('Position', [0 0 1920/8 1080-80]) % [left bottom width height] 
            % figure('Position', [0 1080/2-80 1920/12 1080/2]) % [left bottom width height] % dim2024 fig3c
            % figure('Position', [0 1080/2-80 1920/8 1080/2]) % [left bottom width height]
            figure('Position', [0 0 1920/8 1080]) % [left bottom width height] % dim2024 fig4b
            % figure('Position', [0 1080-1080/3*2-80 1920/8 1080/3*2]) % [left bottom width height] % dim2024 fig3d
%             figure('Position', [0 -40 1920/3 1080]) % [left bottom width height] 
            
            for bhvi = 1:length(BEHAV)
                
                c_mat = squeeze(coef_mat(:, bhvi, :))'; c_mat_sig = c_mat; % N_event x N_ROI
                p_mat = squeeze(pval_mat(:, bhvi, :))'; % (nEVENT, nROI)

                % For multiple testing correction across ROIs
                nEVENT = length(EVENT);
                for ei = 1:nEVENT
                    pVals = squeeze(p_mat(ei, :)); % (1, nROI)
                    switch fdr
                        case 'Benjamini-Hochberg'
                            [cpVals, cAlpha] = fdr_BH(pVals, .5);
                            p_mat(ei, :) = cpVals;
                        case 'Bonferroni'
                            [cpVals, cAlpha] = fwer_bonf(pVals, .5);
                            p_mat(ei, :) = cpVals;
                    end
                end

                sig_ind = (p_mat < 0.05);
                c_mat_sig(p_mat >= 0.05) = nan;
                
                subplot(ceil(length(BEHAV)/bhv_horz), bhv_horz, bhvi)
                %     plot(1:length(EVENT), coef_cell{bhvi}); hold on
%                     plot(1:length(EVENT), c_mat, '-ok'); hold off
                
                if length(EVENT) > 1                    
                    
                    % colororder(cell2mat(colors(1:N_ROI)'))
                    plot(1:length(EVENT), c_mat_sig, ... 
                        'LineStyle', '-', ... 
                        'Marker', '.', ...
                        'MarkerSize', 20); hold on
                    plot(1:length(EVENT), c_mat, ':'); hold off
                    
                    title({erase(label_name, ' LOROV'), [BEHAV{bhvi} bhv_sfx]}, 'Interpreter', 'none')
                    xlim([1 length(EVENT)])
                    xticks(1:length(EVENT));
                    xticklabels(EVENT_NAME)
                    % ylabel({['Corr. coef. (' neural_measure], ' vs behav.)'}, ...
                    %     'FontSize', 8)

                    box off
                    
                    if bhvi == length(BEHAV)
                        legend(ROI{:})
                    end
                    
                else
%                     bar([c_mat', c_mat_sig'])
                    %
                    sig_text = cell(1, N_ROI);
                    for roii = 1:N_ROI
                        if p_mat(roii)<0.05; sig_text{roii} = '*'; end
                        if p_mat(roii)<0.01; sig_text{roii} = '**'; end
                        if p_mat(roii)<0.001; sig_text{roii} = '***'; end
                    end
                    
                    b = bar(c_mat, 'FaceColor', 'flat', 'LineWidth', 1);
                    b.CData(isnan(c_mat_sig), :) = repmat([1 1 1], sum(isnan(c_mat_sig)), 1);
%                     b.CData(~isnan(c_mat_sig), :) = repmat([0 .8 .8], sum(~isnan(c_mat_sig)), 1);
                    b.CData(~isnan(c_mat_sig), :) = repmat(temp_color, sum(~isnan(c_mat_sig)), 1);
                    b.EdgeColor = temp_color;
                    b.LineWidth = 2;
                    hold on
                    text(find(~isnan(c_mat_sig)), ...
                        c_mat(~isnan(c_mat_sig)), ...
                        sig_text(~isnan(c_mat_sig)), 'FontSize', 25, ...
                        'HorizontalAlignment', 'center')
                    hold off
                    box off
%                     title({erase(label_name, ' LOROV'), [BEHAV{bhvi} bhv_sfx], EVENT_NAME{1}}, 'Interpreter', 'none') % commented out for COSYNE 2023 poster fig 3
                    xticks(1:length(ROI));
                    % xtickangle(0)
%                     xtickangle(45) % commented out for COSYNE 2023 poster fig 3
                    xticklabels(ROI)
%                     ylabel({['Corr. coef. (' neural_measure], ' vs behav.)'}, ...
%                         'FontSize', 8) % commented out for COSYNE 2023 poster fig 3
                    
                    % % COSYNE 2023 poster fig 3
                    % ylim([-.4 .8])
                    
                    % COSYNE 2023 poster fig 4
%                     ylim([-.5 .8])
                    
                end

            end

            % if fdr
            %     sgtitle(fdr)
            % end

%         end
    
    case 'corr_behav_ROIeffects'
        
        figure
        for bhvi = 1:length(BEHAV)
            
            c_mat = coef_cell{bhvi};
            p_mat = pval_cell{bhvi};
            
            sig_ind = (p_mat < 0.05);
            c_mat(p_mat >= 0.05) = nan;
            
            subplot(ceil(length(BEHAV)/bhv_horz), bhv_horz, bhvi)
            %     plot(1:length(EVENT), coef_cell{bhvi}); hold on
            %     plot(1:length(EVENT), c_mat, '-ok'); hold off
            plot(1:length(EVENT), c_mat, '-o');
            title([BEHAV{bhvi} bhv_sfx])
            xlim([1 8])
            xticks(1:8);
            xticklabels(EVENT_NAME)
            
            if bhvi == length(BEHAV)
                legend(ROI{:})
            end
        end
        
    case 'corr_behav_dichotomy'
        
        for roii = 1:length(ROI)
            figure('Name', ROI{roii});
            
            for bhvi = 1:length(BEHAV)
                
                c_mat = coef_cell{roii, bhvi};
                p_mat = pval_cell{roii, bhvi};
                
                sig_ind = (p_mat < 0.05);
%                 c_mat(p_mat >= 0.05) = nan;
                
                subplot(ceil(length(BEHAV)/bhv_horz), bhv_horz, bhvi)
                %     plot(1:length(EVENT), coef_cell{bhvi}); hold on
                %     plot(1:length(EVENT), c_mat, '-ok'); hold off
                
                if length(EVENT) > 1
                    plot(1:length(EVENT), c_mat, '-o');
                    title([BEHAV{bhvi} bhv_sfx])
                    xlim([1 8])
                    xticks(1:8);
                    xticklabels(EVENT_NAME)
                    
                    if bhvi == length(BEHAV)
                        legend('Goal', 'UC', 'NonlinearIT', 'LinearIT')
                    end
                    
                else
                    bar(c_mat)
                    title({[BEHAV{bhvi} bhv_sfx], EVENT_NAME{1}})
                    xticks(1:4);
                    xtickangle(45)
                    xticklabels({'Goal', 'UC', 'NonlinearIT', 'LinearIT'})
                    
                end
            end
        end
        
end
disp('done')

