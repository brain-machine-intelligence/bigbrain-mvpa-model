
clear all

%% Set working directory

base_path = fullfile('..'); % If your current working directory is 'demo', use this. Otherwise, set your own path.
addpath(fullfile(base_path, 'functions', 'utils'))

% If your current working directory is 'demo', use this. Otherwise, set your own path that includes 'data'.
cd('..'); 


%% taskVar profile (continuous)

% plot var2 time series
% observe var2 class distribution for each var1 condition
% subplot per var1 class
% legend per var2 class

% within-subject measures (mean var2 per var1 condition)
% across-subject statistical test (mean var2s ANOVA/paired t-test)

LoadDefaultSettings;
close all
warning off
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);

disp('start')

% Setting =================================================================

NameValue1 = {};
NameValue2 = {};

% ===== Experiment =====
Exp = 'Lee2014';  idsize = 20;  

% Choice optimality, specific goal condition, uncertainty low vs high
var1_name = 'UncCond'; N_var1_global = 2; var1_cond = []; var2_name = 'ChoOptBinNaN1n2'; var2_cond = 'G';

% Choice optimality, non-specific goal condition, uncertainty low vs high
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = []; var2_name = 'ChoOptBinNaN1n2'; var2_cond = 'H';

% Choice consistency, specific goal condition, uncertainty low vs high
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = []; var2_name = 'ChoConsist1n2'; var2_cond = 'G';

% Choice consistency, non-specific goal condition, uncertainty low vs high
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = []; var2_name = 'ChoConsist1n2'; var2_cond = 'H';


var1_unq_pre = []; % default setting - some var1s are special cases
categorize1 = [];
categorize2 = [];


% ===== visualization parameters =====
vis_mode = 'distribution_per_condition'; horz_slot = 11;

% --- individual result ---
% indiv_vis = 0;
indiv_vis = 1;
sig_ind_mark = 0;
% sig_ind_mark = 1;

% plot_mode = 'plot';
plot_mode = 'stem';
% plot_mode = 'movmean';
% horz_slot = 11;       % mainly for distribution
% horz_slot = 4;        % mainly for time series
% horz_slot = 3;
% horz_slot = 2;

% group result histogram
nbins = 8;
FaceAlpha = 0.4;
color_map = {'r', 'g', 'b', 'c', 'm', 'y'};

% Visualization ===========================================================

var2_means = nan * ones(idsize, N_var1_global);
var1_storage = cell(idsize, 1);
var2_storage = cell(idsize, 1);
id_storage = cell(idsize, 1);

% Q. Is var1 effect on var2 significant for each subject?: N_sig_sbj
% only for vis_mode = 'distribution_per_condition';
N_sig_sbj = nan(idsize, 1);
EL_storage = cell(idsize, 1);
GT_storage = cell(idsize, 1);
UC_storage = cell(idsize, 1);
GS_storage = cell(idsize, 1);
storages.EL_storage = EL_storage;
storages.GT_storage = GT_storage;
storages.UC_storage = UC_storage;
storages.GS_storage = GS_storage;

% figure
if indiv_vis
    figure('Position', [0 1080/2-80 1920 1080/2])
end

% ===== Individual result =====
for id = 1:idsize
    fprintf('%d ', id)
    
    % ===== Data load =====
    [var1, var2, var1_unq, N_var1] = load_var1var2(Exp, id, ... 
        var1_name, var2_name, ...
        var1_cond, var2_cond, ...
        NameValue1, NameValue2, ...
        'categorize1', categorize1, ...
        'categorize2', categorize2, ...
        'var1_unq_pre', var1_unq_pre);
        
    % =============== var2 mean per var1 class ================
    var1_class_names = cell(1, N_var1);
    
    if N_var1 ~= N_var1_global
        fprintf('var1 class lack: ')
        fprintf('(%d) ', N_var1)
    end
    
    for eli = 1:N_var1        
        var1_class_names{eli} = [var1_name var1_cond ' ' num2str(var1_unq(eli))];
        
        sub_var2 = var2(var1==var1_unq(eli));
        if ~isempty(sub_var2)
%             var2_means(id, eli) = mean(sub_var2);
            var2_means(id, eli) = nanmean(sub_var2);
        end
    end
    
    
    
    if id==idsize
        
        switch vis_mode
            case 'time_series' % time series per sbj
                % ===== legend =====
                legend(stem_per_var1, var1_class_names)
                title_text2 = [];
                
            case 'distribution_per_condition' % distribution per sbj
                title_text2 = [];

                switch size(N_sig_sbj, 2)
                    case 1
                        % fprintf('\n # of subject with significant var1 effect on var2: %d \n', ...
                        %     nansum(any(N_sig_sbj, 2)))
                        % title_text2 = ['Sig. sbj N=' sprintf('%d', nansum(any(N_sig_sbj, 2)))];
                    case 3
                        % fprintf('\n # of subject with significant var1 effect on var2: %d \n', ...
                        %     nansum(N_sig_sbj))
                        % title_text2 = ['Sig. sbj N=' ... 
                        %     sprintf('%d, ', nansum(N_sig_sbj(:, 1))), ...
                        %     sprintf('%d, ', nansum(N_sig_sbj(:, 2))), ...
                        %     sprintf('%d', nansum(N_sig_sbj(:, 3)))];
                    case 6
                        % fprintf('\n # of subject with significant var1 effect on var2: %d \n', ...
                        %     nansum(N_sig_sbj))
                        % title_text2 = ['Sig. sbj N=' ... 
                        %     sprintf('%d, ', nansum(N_sig_sbj(:, 1))), ...
                        %     sprintf('%d, ', nansum(N_sig_sbj(:, 2))), ...
                        %     sprintf('%d, ', nansum(N_sig_sbj(:, 3))), ...
                        %     sprintf('%d, ', nansum(N_sig_sbj(:, 4))), ...
                        %     sprintf('%d, ', nansum(N_sig_sbj(:, 5))), ...
                        %     sprintf('%d', nansum(N_sig_sbj(:, 6)))];
                end
        end
    end
    
end % id


% #########################################################################
% =================== Group result (statistical test) =====================
% #########################################################################


set(0,'DefaultAxesFontSize',20,'DefaultTextFontSize',12); titleFontSize = 14; % for paper figure
% set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
% set(0,'DefaultAxesFontSize',10,'DefaultTextFontSize',10);

switch N_var1_global

    case 2 % binary condition

        % across-subject paired t-test (var2 mean values for var1 class 1 vs 2)
        [~, p, ci, stats] = ttest(var2_means(:, 1), var2_means(:, 2));
        title_prefix = 'paired t-test p=';
        % figure
        % figure('Position', [0 1080/2-80 1920/8 1080/2]) % [left bottom width height], main figure 1 
        figure('Position', [0 1080/2-80 1920/6 1080/2]) % [left bottom width height] 
        % figure('Position', [0 1080/2-80 1920/8 1080/3]) % [left bottom width height] 
        
        
        try
            % MATLAB 2020
            
            % Boxchart color control
            if isempty(var2_cond)
                BoxFaceColor = [.4 .4 .4];
            elseif contains(var2_cond, 'G')
                BoxFaceColor = [.75 .05 .55];
            elseif contains(var2_cond, 'H')
                BoxFaceColor = [.45 .3 .8];
            else
                BoxFaceColor = [.4 .4 .4];
            end
            
            % Boxchart
            if any(ismember(var1_unq, 0))
                % if var1_unq is [0, 1] -> change into [1 2] for boxchart
                var1_unq2 = var1_unq+1;
            else
                var1_unq2 = var1_unq;
            end
            bc = boxchart(reshape(repmat(var1_unq2, idsize, 1), [], 1), ...
                reshape(var2_means, [], 1), ...
                'GroupByColor', reshape(repmat(var1_unq2, idsize, 1), [], 1), ...
                'BoxFaceAlpha', .4, ...
                'LineWidth', 4, ...
                'MarkerStyle', 'none'); 
            hold on;
            bc(1).BoxFaceColor = [1 1 1];
            bc(1).BoxEdgeColor = [0 0 0];
            bc(1).XData = 1.2 * ones(idsize, 1);
            bc(2).BoxFaceColor = [0 0 0];
            bc(2).BoxEdgeColor = [0 0 0];
            bc(2).XData = 1.8 * ones(idsize, 1);
            % bc = boxchart(reshape(repmat(var1_unq, idsize, 1), [], 1), ...
            %     reshape(var2_means, [], 1), ...
            %     'BoxFaceColor', BoxFaceColor, ... % G: [.85 .35 .45], H: [.25 .3 .9], Grey: [.4 .4 .4]
            %     'LineWidth', 4, ...
            %     'MarkerStyle', 'none'); hold off
            xlabel(var1_name)
            xticks(var1_unq2)
            xticklabels(var1_unq)

            % plot(var2_means', ... 
            %     'LineWidth', 2, ... 
            %     'Marker', 'o', ... 
            %     'MarkerSize', 7, ...
            %     'Color', [.85 .85 .85]); hold off

            % For Figure 2
%             plot(var2_means(isnan(N_sig_sbj), :)', 'LineWidth', 1, 'Marker', 'o', 'Color', [.9 .9 .9]); hold on;
%             plot(var2_means(~isnan(N_sig_sbj), :)', 'LineWidth', 1, 'Marker', 'o', 'Color', [.6 .6 .6]);
            plot([var1_unq2(1) + .2  var1_unq2(2) - .2], ...
                var2_means', ... 
                'LineWidth', 2, ... 
                'Marker', 'o', ... 
                'MarkerSize', 7, ...
                'Color', [.85 .85 .85]); 
            
            % Significance mark
%             if p < 0.05
%                 sig_text = sprintf('%.2e', p);
%             else
%                 sig_text = 'n.s.';
%             end
            if p < 0.0001
                sig_text = '****';
            elseif p < 0.001
                sig_text = '***';
            elseif p < 0.01
                sig_text = '**';
            elseif p < 0.05
                sig_text = '*';
            else
                sig_text = 'n.s.';
            end
            
            % Significance mark display
            if contains(var2_name, {'ChoOptBinNaN', 'ChoConsist'})
                
                % For opt/consist                
                line(var1_unq2, [1 + .04,  1 + .04], ...
                    'Color', 'k')
                text(1.5 - .055 * length(sig_text), 1 + .07, sig_text, ...
                    'FontSize', 20)
            
                if contains(var2_name, 'ChoOptBinNaN')
                    ylim([.3 1.1])
                end
                
                if contains(var2_name, 'ChoConsist')
                    ylim([.5 1.1])
                end
                
                xlim([.7 2.3])
                % xlim([.6 2.4])
                % xlim([.5 2.5])
                xticks(1:2); xticklabels({'Low', 'High'})

                % print statistic and p-value
                fprintf(['t(%d)=%.3f\t p=%.3e (' sig_text ')\n'], stats.df, stats.tstat, p)
                % fprintf('t(%d)=%.3f\n', stats.df, stats.tstat)
                % fprintf('p=%.3e\n', p)
                
            else
                
%                 set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
                
                % For normal cases
                line(var1_unq2, [max(var2_means(:)) + .04,  max(var2_means(:)) + .04], ...
                    'Color', 'k')
                text(1.5 - .04 * length(sig_text), max(var2_means(:)) + .06, sig_text, ...
                    'FontSize', 15)
                ylim([min(var2_means(:)) - .04, max(var2_means(:)) + .08])
                
            end

            hold off

        catch err
            disp(err)
            
            % MATLAB 2018
            plot(var2_means(isnan(N_sig_sbj), :)', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.9 .9 .9]); hold on;
            plot(var2_means(~isnan(N_sig_sbj), :)', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.3 .7 .3]);
            boxplot(var2_means, ... 
                'BoxStyle', 'outline', ...
                'Colors', 'k', ... [0.7 0.1 0.6]
                'Width', .5, ...
                'Notch', 'on', ...
                'Labels', var1_class_names); hold off            
        end
%         boxplot(var2_means, ... 
%             'Colors', 'k', ... [0.7 0.1 0.6]
%             'Width', 0.3, ...
%             'Notch', 'on', ...
%             'Labels', var1_class_names);
% %         boxplot(var2_means, 'Colors', [0.7 0.1 0.6], 'PlotStyle', 'compact');
% %         title({[var2_name var2_cond], [' per ' var1_name var1_cond]})
% %         title_text = [];
        title_text = [title_prefix, num2str(p)];
        % rldm 2022
        hold on
%         bar(nanmean(var2_means), 'BarWidth', 0.3);
%         errorbar(nanmean(var2_means, 1), nanstd(var2_means, 0 ,1)./sqrt(idsize), ...
%             'LineStyle', 'none', ...
%             'LineWidth', 1.2);
%         ylim([0.6 .9])
%         ylim([0.5 1])
%         ylim([.4 .7])
%         ylim([.5 .7])
%         ylim([.4 1])
        box off

        % =================================================================
    otherwise % N_var1_global > 2
        % =================================================================

%         % across-subject one-way ANOVA (var2 mean values per var1 class)
%         var2mean_vec = reshape(var2_means, [], 1);
%         var1class_vec = reshape(repmat(1:N_var1_global, idsize, 1), [], 1);
%         [p, tbl, stats] = anova1(var2mean_vec, var1class_vec, 'off'); % auto plot off
%         title_prefix = 'ANOVA p=';
        
        % 2023-04-30 
        % N-way ANOVA for repeated measure anova
        var1_total_vec = cat(1, var1_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
        var2_total_vec = cell2mat(var2_storage); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] mat 
        id_total_vec = cat(1, id_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
        [p, tbl, stats] = anovan(var2_total_vec, {var1_total_vec, id_total_vec}, ...
            'model', 'interaction', ...
            'varnames', {var1_name, 'Subject'}, ...
            'display', 'off');
        p = p(1); % main effect of var1
        title_prefix = 'N-way ANOVA p=';
        

        
        % #################################################################
        % ------------------------ Special cases ------------------------
        % #################################################################
        
        title_text = [];
        
        [title_text] = anovan_group(idsize, N_var1_global, ...
            var1_name, var1_cond, var2_name, var2_cond, ...
            var2_means, ...
            var2_total_vec, ...
            'storages', storages, ...
            'id_total_vec', id_total_vec);

        if isempty(title_text)
            title_text = [title_prefix num2str(p)];
        end
        
        % #################################################################
        % ------------------------ General cases ------------------------ 
        % #################################################################

        figure('Position', [0 1080/2-80 1920/4 1080/2]) 
%         subplot(1,3,1:2)
        %
        
        if indiv_vis
            if sig_ind_mark
                plot(var2_means(any(~isnan(N_sig_sbj), 2), :)', ... 
                    'LineWidth', 1.5, ... 
                    'Marker', 'o', ... 
                    'MarkerSize', 7, ...
                    'Color', [.3 .7 .3]);
            else
                plot(var1_unq, var2_means(any(~isnan(N_sig_sbj), 2), :)', ... 
                    'LineWidth', 1.5, ... 
                    'Marker', 'o', ... 
                    'MarkerSize', 7, ...
                    'Color', [.85 .85 .85]);
            end
            hold on;
        end
        
        
        % #################################################################
        % ------------------------ Visualization ------------------------ 
        % #################################################################

        try
            % MATLAB 2020, 2023
            boxchart(reshape(repmat(var1_unq, idsize, 1), [], 1), ...
                reshape(var2_means, [], 1), ...
                'BoxFaceColor', 'none', ...
                'BoxEdgeColor', 'k', ...
                'LineWidth', 4, ...
                'MarkerStyle', 'none'); hold off
            xticks(1:length(var1_unq));
            xticklabels(var1_class_names);

        catch
            % MATLAB 2018
            boxplot(var2_means, ...
                'BoxStyle', 'outline', ...
                'Colors', 'k', ... [0.7 0.1 0.6]
                'Width', .5, ...
                'Notch', 'on', ...
                'Labels', var1_class_names); hold off
        end
%         x = reshape(repmat(1:N_var1_global, idsize, 1), [], 1);
%         y = reshape(var2_means, [], 1);
%         beeswarm(x, y, 'corral_style', 'random', 'overlay_style', 'ci')
        %
        box off
end

title(sprintf([var1_name var1_cond cell2mat(cellfun(@num2str, NameValue1, 'UniformOutput', false)) ... 
                ' %d ' var2_name var2_cond]), ... 
                'FontSize', 14)

% Title the group result figure
if iscell(title_text)
    title(cat(2, {[var2_name var2_cond ... 
        cell2mat(cellfun(@num2str, NameValue2, 'UniformOutput', false)) ... 
        ' per ' var1_name var1_cond]}, ...
        title_text, title_text2), 'FontSize', titleFontSize)
else
    title({[var2_name var2_cond ... 
        cell2mat(cellfun(@num2str, NameValue2, 'UniformOutput', false)) ... 
        ' per ' var1_name var1_cond], ...
        title_text, title_text2}, 'FontSize', titleFontSize)
end

% x-label the group result figure
xlabel([var1_name var1_cond])
% xticks(var1_unq2)
% xticklabels(var1_unq)

% ylim
if contains(var2_name, 'A1') || contains(var2_name, 'A2')
    ylim([1 2])
end



disp('end')




%% helper function: N-way anova for an individual participant

function [N_sig_sbj, storages] = anovan_indiv(id, var1_name, var1_cond, var2_name, var2_cond, ...
    g, x, N_sig_sbj, storages, varargin)

% for 'taskVar profile (continuous)' section

% Default input parameters
options = struct('indiv_vis', [], ...
    'var1_unq', []); 
% Read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(anovan_indiv) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(anovan_indiv) %s is not a recognized parameter name', pair{1})
    end
end

indiv_vis = options.indiv_vis;
var1_unq = options.var1_unq;

if strcmp(var1_name, 'blkCond')

    % Two-way ANOVA
    GT = cell(1, length(g));
    GT(ismember(g, [1, 2])) = repmat({'Spec'}, 1, sum(ismember(g, [1, 2])));
    GT(ismember(g, [3, 4])) = repmat({'Flex'}, 1, sum(ismember(g, [3, 4])));
    UC = cell(1, length(g));
    UC(ismember(g, [1, 4])) = repmat({'Low'}, 1, sum(ismember(g, [1, 4])));
    UC(ismember(g, [2, 3])) = repmat({'High'}, 1, sum(ismember(g, [2, 3])));

    storages.GT_storage{id} = GT';
    storages.UC_storage{id} = UC';
    % storages.GT_storage = GT_storage;
    % storages.UC_storage = UC_storage;

    [p,~,stats, ~] = anovan(x, {GT, UC}, ...
        'model', 'interaction', ...
        'varnames', {'GT', 'UC'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GT p=%0.2e \n(%g)', p(1), stats.coeffs(2))]; % coeff: G effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC p=%0.2e \n(%g)', p(2), stats.coeffs(5))]; % coeff: UC=High effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GT*UC p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'GT x UC')

    % Two-way ANOVA
    GT = cell(1, length(g));
    GT(ismember(g, [1, 2])) = repmat({'Spec'}, 1, sum(ismember(g, [1, 2])));
    GT(ismember(g, [3, 4])) = repmat({'Flex'}, 1, sum(ismember(g, [3, 4])));
    UC = cell(1, length(g));
    UC(ismember(g, [1, 3])) = repmat({'Low'}, 1, sum(ismember(g, [1, 3])));
    UC(ismember(g, [2, 4])) = repmat({'High'}, 1, sum(ismember(g, [2, 4])));

    storages.GT_storage{id} = GT';
    storages.UC_storage{id} = UC';
    % storages.GT_storage = GT_storage;
    % storages.UC_storage = UC_storage;

    [p,~,stats, ~] = anovan(x, {GT, UC}, ...
        'model', 'interaction', ...
        'varnames', {'GT', 'UC'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GT p=%0.2e \n(%g)', p(1), stats.coeffs(2))]; % coeff: G effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC p=%0.2e \n(%g)', p(2), stats.coeffs(5))]; % coeff: UC=High effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GT*UC p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'EL x GT x UC')

    % Three-way ANOVA
    EL = cell(1, length(g));
    EL(ismember(g, 1:4)) = repmat({'E'}, 1, sum(ismember(g, 1:4)));
    EL(ismember(g, 5:8)) = repmat({'L'}, 1, sum(ismember(g, 5:8)));
    GT = cell(1, length(g));
    GT(ismember(g, [1 2 5 6])) = repmat({'G'}, 1, sum(ismember(g, [1 2 5 6])));
    GT(ismember(g, [3 4 7 8])) = repmat({'H'}, 1, sum(ismember(g, [3 4 7 8])));
    UC = cell(1, length(g));
    UC(ismember(g, [1 3 5 7])) = repmat({'L'}, 1, sum(ismember(g, [1 3 5 7])));
    UC(ismember(g, [2 4 6 8])) = repmat({'H'}, 1, sum(ismember(g, [2 4 6 8])));

    storages.EL_storage{id} = EL';
    storages.GT_storage{id} = GT';
    storages.UC_storage{id} = UC';
    % storages.EL_storage = EL_storage;
    % storages.GT_storage = GT_storage;
    % storages.UC_storage = UC_storage;

    [p,~,stats, ~] = anovan(x, {EL, GT, UC}, ...
        'model', 'interaction', ...
        'varnames', {'EL', 'GT', 'UC'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('EL p=%0.2e \n(%g)', p(1), stats.coeffs(3))]; % coeff: Late effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GT p=%0.2e \n(%g)', p(2), stats.coeffs(4))]; % coeff: G effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC p=%0.2e \n(%g)', p(3), stats.coeffs(7))]; % coeff: UC=High effect
        N_sig_sbj(id, 3) = 1;
    end

    if p(4) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('EL*GT p=%0.2e', p(4))];
        N_sig_sbj(id, 4) = 1;
    end

    if p(5) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('EL*UC p=%0.2e', p(5))];
        N_sig_sbj(id, 5) = 1;
    end

    if p(6) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GT*UC p=%0.2e', p(6))];
        N_sig_sbj(id, 6) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'uLuH EL')

    % Two-way ANOVA
    %                 UC = cell(1, length(g));
    %                 UC(ismember(g, [1, 2])) = repmat({'Low'}, 1, sum(ismember(g, [1, 2])));
    %                 UC(ismember(g, [3, 4])) = repmat({'High'}, 1, sum(ismember(g, [3, 4])));
    %                 EL = cell(1, length(g));
    %                 EL(ismember(g, [1, 3])) = repmat({'Early'}, 1, sum(ismember(g, [1, 3])));
    %                 EL(ismember(g, [2, 4])) = repmat({'Late'}, 1, sum(ismember(g, [2, 4])));
    UC = cell(1, length(g));
    UC(ismember(g, [1, 3])) = repmat({'Low'}, 1, sum(ismember(g, [1, 3])));
    UC(ismember(g, [2, 4])) = repmat({'High'}, 1, sum(ismember(g, [2, 4])));
    EL = cell(1, length(g));
    EL(ismember(g, [1, 2])) = repmat({'Early'}, 1, sum(ismember(g, [1, 2])));
    EL(ismember(g, [3, 4])) = repmat({'Late'}, 1, sum(ismember(g, [3, 4])));

    storages.EL_storage{id} = EL';
    storages.UC_storage{id} = UC';
    % storages.EL_storage = EL_storage;
    % storages.UC_storage = UC_storage;

    [p,~,stats] = anovan(x, {UC, EL}, ...
        'model', 'interaction', ...
        'varnames', {'UC', 'EL'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC p=%0.2e \n(%g)', p(1), stats.coeffs(3))]; % coeff: UC=High effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('EL p=%0.2e \n(%g)', p(2), stats.coeffs(5))]; % coeff: EL=Late effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC*EL p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'uLuH ucEL')

    % Two-way ANOVA
    UC = cell(1, length(g));
    UC(ismember(g, [1, 2])) = repmat({'Low'}, 1, sum(ismember(g, [1, 2])));
    UC(ismember(g, [3, 4])) = repmat({'High'}, 1, sum(ismember(g, [3, 4])));
    EL = cell(1, length(g));
    EL(ismember(g, [1, 3])) = repmat({'Early'}, 1, sum(ismember(g, [1, 3])));
    EL(ismember(g, [2, 4])) = repmat({'Late'}, 1, sum(ismember(g, [2, 4])));

    [p,~,stats] = anovan(x, {UC, EL}, ...
        'model', 'interaction', ...
        'varnames', {'UC', 'ucEL'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC p=%0.2e \n(%g)', p(1), stats.coeffs(3))]; % coeff: UC=High effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('ucEL p=%0.2e \n(%g)', p(2), stats.coeffs(5))]; % coeff: ucEL=Late effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC*ucEL p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'GoalxUC')

    % Two-way ANOVA
    UC = cell(1, length(g));
    UC(ismember(g, [1 2 3 4])) = repmat({'Low'}, 1, sum(ismember(g, [1 2 3 4])));
    UC(ismember(g, [5 6 7 8])) = repmat({'High'}, 1, sum(ismember(g, [5 6 7 8])));
    G678H = cell(1, length(g));
    G678H(ismember(g, [1 4])) = repmat({'6'}, 1, sum(ismember(g, [1 4])));
    G678H(ismember(g, [2 5])) = repmat({'7'}, 1, sum(ismember(g, [2 5])));
    G678H(ismember(g, [3 6])) = repmat({'8'}, 1, sum(ismember(g, [3 6])));
    G678H(ismember(g, [7 8])) = repmat({'H'}, 1, sum(ismember(g, [7 8])));

    [p,~,stats] = anovan(x, {UC, G678H}, ...
        'model', 'interaction', ...
        'varnames', {'UC', 'G678H'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC p=%0.2e \n(%g)', p(1), stats.coeffs(3))]; % coeff: UC=High effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('G678H p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(4:7))]; % coeff: G678 effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC*G678H p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'Goal(6,7,8)xUC')

    % Two-way ANOVA
    UC = cell(1, length(g));
    UC(ismember(g, [1 2 3])) = repmat({'Low'}, 1, sum(ismember(g, [1 2 3])));
    UC(ismember(g, [4 5 6])) = repmat({'High'}, 1, sum(ismember(g, [4 5 6])));
    G678 = cell(1, length(g));
    G678(ismember(g, [1 4])) = repmat({'6'}, 1, sum(ismember(g, [1 4])));
    G678(ismember(g, [2 5])) = repmat({'7'}, 1, sum(ismember(g, [2 5])));
    G678(ismember(g, [3 6])) = repmat({'8'}, 1, sum(ismember(g, [3 6])));

    [p,~,stats] = anovan(x, {UC, G678}, ...
        'model', 'interaction', ...
        'varnames', {'UC', 'G678'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC p=%0.2e \n(%g)', p(1), stats.coeffs(3))]; % coeff: UC=High effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('G678 p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(4:6))]; % coeff: G678 effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC*G678 p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'uLuH Goal(6,7,8)Switch')

    %                 var1(UC == 1 & GS == 1) = 1; % uL, switch
    %                 var1(UC == 1 & GS == 0) = 2; % uL, stay
    %                 var1(UC == 2 & GS == 1) = 3; % uH, switch
    %                 var1(UC == 2 & GS == 0) = 4; % uH, stay
    %

    % Two-way ANOVA
    UC = cell(1, length(g));
    UC(ismember(g, [1, 2])) = repmat({'Low'}, 1, sum(ismember(g, [1, 2])));
    UC(ismember(g, [3, 4])) = repmat({'High'}, 1, sum(ismember(g, [3, 4])));
    GS = cell(1, length(g));
    GS(ismember(g, [1, 3])) = repmat({'Switch'}, 1, sum(ismember(g, [1, 3])));
    GS(ismember(g, [2, 4])) = repmat({'Stay'}, 1, sum(ismember(g, [2, 4])));

    storages.UC_storage{id, 1} = UC';
    storages.GS_storage{id, 1} = GS';

    [p,~,stats] = anovan(x, {UC, GS}, ...
        'model', 'interaction', ...
        'varnames', {'UC', 'GS'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC p=%0.2e \n(%g)', p(1), stats.coeffs(3))]; % coeff: UC=High effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GS p=%0.2e \n(%g)', p(2), stats.coeffs(4))]; % coeff: GS=Switch effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC*GS p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'Session x GT')

    % Two-way ANOVA
    Ss = cell(1, length(g));
    Ss(ismember(g, [1 6])) = repmat({'1'}, 1, sum(ismember(g, [1 6])));
    Ss(ismember(g, [2 7])) = repmat({'2'}, 1, sum(ismember(g, [2 7])));
    Ss(ismember(g, [3 8])) = repmat({'3'}, 1, sum(ismember(g, [3 8])));
    Ss(ismember(g, [4 9])) = repmat({'4'}, 1, sum(ismember(g, [4 9])));
    Ss(ismember(g, [5 10])) = repmat({'5'}, 1, sum(ismember(g, [5 10])));
    GT = cell(1, length(g));
    GT(ismember(g, 1:5)) = repmat({'1'}, 1, sum(ismember(g, 1:5)));
    GT(ismember(g, 6:10)) = repmat({'2'}, 1, sum(ismember(g, 6:10)));

    [p,~,stats] = anovan(x, {Ss, GT}, ...
        'model', 'interaction', ...
        'varnames', {'Ss', 'GT'}, ...
        'display', 'off');

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('Ss p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(2:6))]; % coeff: Session effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GT p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(7))]; % coeff: specific goal (GT==1) effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('Ss*GT p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'Session x UC')

    % Two-way ANOVA
    Ss = cell(1, length(g));
    Ss(ismember(g, [1 6])) = repmat({'1'}, 1, sum(ismember(g, [1 6])));
    Ss(ismember(g, [2 7])) = repmat({'2'}, 1, sum(ismember(g, [2 7])));
    Ss(ismember(g, [3 8])) = repmat({'3'}, 1, sum(ismember(g, [3 8])));
    Ss(ismember(g, [4 9])) = repmat({'4'}, 1, sum(ismember(g, [4 9])));
    Ss(ismember(g, [5 10])) = repmat({'5'}, 1, sum(ismember(g, [5 10])));
    UC = cell(1, length(g));
    UC(ismember(g, 1:5)) = repmat({'1'}, 1, sum(ismember(g, 1:5)));
    UC(ismember(g, 6:10)) = repmat({'2'}, 1, sum(ismember(g, 6:10)));

    [p,~,stats] = anovan(x, {Ss, UC}, ...
        'model', 'interaction', ...
        'varnames', {'Ss', 'UC'}, ...
        'display', 'off');

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('Ss p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(2:6))]; % coeff: Session effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('UC p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(8))]; % coeff: high uncertainty (UC==2) effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('Ss*UC p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'Session x Goal(6,7,8)')

    % Two-way ANOVA
    Ss = cell(1, length(g));
    Ss(ismember(g, [1 6 11])) = repmat({'1'}, 1, sum(ismember(g, [1 6 11])));
    Ss(ismember(g, [2 7 12])) = repmat({'2'}, 1, sum(ismember(g, [2 7 12])));
    Ss(ismember(g, [3 8 13])) = repmat({'3'}, 1, sum(ismember(g, [3 8 13])));
    Ss(ismember(g, [4 9 14])) = repmat({'4'}, 1, sum(ismember(g, [4 9 14])));
    Ss(ismember(g, [5 10 15])) = repmat({'5'}, 1, sum(ismember(g, [5 10 15])));
    G678 = cell(1, length(g));
    G678(ismember(g, 1:5)) = repmat({'6'}, 1, sum(ismember(g, 1:5)));
    G678(ismember(g, 6:10)) = repmat({'7'}, 1, sum(ismember(g, 6:10)));
    G678(ismember(g, 11:15)) = repmat({'8'}, 1, sum(ismember(g, 11:15)));

    [p,~,stats] = anovan(x, {Ss, G678}, ...
        'model', 'interaction', ...
        'varnames', {'Ss', 'G678'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('Ss p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(2:6))]; % coeff: Session effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('G678 p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(7:9))]; % coeff: G678 effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('Ss*G678 p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'Session x GS')

    % Two-way ANOVA
    Ss = cell(1, length(g));
    Ss(ismember(g, [1 6])) = repmat({'1'}, 1, sum(ismember(g, [1 6])));
    Ss(ismember(g, [2 7])) = repmat({'2'}, 1, sum(ismember(g, [2 7])));
    Ss(ismember(g, [3 8])) = repmat({'3'}, 1, sum(ismember(g, [3 8])));
    Ss(ismember(g, [4 9])) = repmat({'4'}, 1, sum(ismember(g, [4 9])));
    Ss(ismember(g, [5 10])) = repmat({'5'}, 1, sum(ismember(g, [5 10])));
    GS = cell(1, length(g));
    GS(ismember(g, 1:5)) = repmat({'0'}, 1, sum(ismember(g, 1:5)));
    GS(ismember(g, 6:10)) = repmat({'1'}, 1, sum(ismember(g, 6:10)));

    [p,~,stats] = anovan(x, {Ss, GS}, ...
        'model', 'interaction', ...
        'varnames', {'Ss', 'GS'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('Ss p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(2:6))]; % coeff: Session effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GS p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(8))]; % coeff: goal switch (GS==1) effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('Ss*GS p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

elseif strcmp(var1_name, 'EL x GS')

    % Two-way ANOVA
    EL = cell(1, length(g));
    EL(ismember(g, [1 3])) = repmat({'1'}, 1, sum(ismember(g, [1 3])));
    EL(ismember(g, [2 4])) = repmat({'2'}, 1, sum(ismember(g, [2 4])));
    GS = cell(1, length(g));
    GS(ismember(g, 1:2)) = repmat({'0'}, 1, sum(ismember(g, 1:2)));
    GS(ismember(g, 3:4)) = repmat({'1'}, 1, sum(ismember(g, 3:4)));

    [p,~,stats] = anovan(x, {EL, GS}, ...
        'model', 'interaction', ...
        'varnames', {'EL', 'GS'}, ...
        'display', 'off');
    %                 disp([var2_name var2_cond ' 2-way ANOVA: ' var1_name var1_cond ' effects'])
    %                 disp(tbl)

    if indiv_vis; boxplot(x, g); end

    sbj_title_text = {};

    if p(1) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('EL p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(3))]; % coeff: Learning (EL==2) effect
        N_sig_sbj(id, 1) = 1;
    end

    if p(2) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('GS p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(5))]; % coeff: goal switch (GS==1) effect
        N_sig_sbj(id, 2) = 1;
    end

    if p(3) < 0.05
        sbj_title_text = [sbj_title_text, ...
            sprintf('EL*GS p=%0.2e', p(3))];
        N_sig_sbj(id, 3) = 1;
    end

    title(sbj_title_text);

else

    % within-subject one-way ANOVA
    [p, ~, stats] = anova1(x, g, 'off'); % auto plot off

    boxplot(x, g, 'Positions', var1_unq);
    %                 xticks(1:length(var1_unq));
    %                 xticklabels(var1_unq);
    if p < 0.05
        title(num2str(p))
        N_sig_sbj(id) = 1;
    end
    xlabel([var1_name var1_cond])
    ylabel([var2_name var2_cond])

end

end % function end

%% helper function: N-way anova for group analysis

function [title_text] = anovan_group(idsize, N_var1_global, ...
    var1_name, var1_cond, var2_name, var2_cond, ... 
    var2_means, var2_total_vec, ...
    varargin)

% Default input parameters
options = struct('storages', [], ...
    'id_total_vec', [], ...
    'id_vec', []);
% Read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(anovan_group) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(anovan_group) %s is not a recognized parameter name', pair{1})
    end
end

except_ids = [];

if ~isempty(options.storages)
    storages = options.storages;
end
id_total_vec = options.id_total_vec;
% id_vec = options.id_vec;

var2mean_vec = reshape(var2_means, [], 1); % for visualization

if strcmp(var1_name, 'blkCond') || strcmp(var1_name, 'GT x UC')

    %             % two-way ANOVA
    %             GT = reshape(repmat({'Spec', 'Spec', 'Flex', 'Flex'}, idsize, 1), [], 1);
    %             UC = reshape(repmat({'Low', 'High', 'High', 'Low'}, idsize, 1), [], 1);
    %             id_vec = cellfun(@num2str, num2cell(repmat(1:idsize, 1, N_var1_global)'), ...
    %                 'UniformOutput', false);
    %             [p,tbl,stats,terms] = anovan(var2mean_vec, {GT, UC, id_vec}, ...
    %                 'model', 'interaction', ...
    %                 'varnames', {'GT', 'UC', 'sbj'}, ...
    %                 'display', 'off');

    if isfield(storages, 'GT_storage'); GT_storage = storages.GT_storage; end
    if isfield(storages, 'UC_storage'); UC_storage = storages.UC_storage; end
    GT_total_vec = cat(1, GT_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
    UC_total_vec = cat(1, UC_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
    [p,tbl,stats,~] = anovan(var2_total_vec, ...
        {GT_total_vec, UC_total_vec, id_total_vec}, ...
        'model', 'interaction', ...
        'varnames', {'GT', 'UC', 'sbj'}, ...
        'display', 'off');

    disp([var2_name var2_cond ' 3-way ANOVA: ' var1_name var1_cond ' effects'])
    disp(tbl)
    fprintf(strjoin(stats.coeffnames(2:5), ', '));
    disp(stats.coeffs(2:5)');

    p = p([1, 2, 4]); % GT, UC, GT * UC

    title_text = {};

    if p(1) < 0.05
        title_text = [title_text, ...
            sprintf('GT p=%0.2e (%.3g)', p(1), stats.coeffs(2))]; % coeff: specific goal (GT==Spec) effect
    end

    if p(2) < 0.05
        title_text = [title_text, ...
            sprintf('UC p=%0.2e (%.3g)', p(2), stats.coeffs(5))]; % coeff: high uncertainty (UC==2) effect
    end

    if p(3) < 0.05
        title_text = [title_text, ...
            sprintf('GT*UC p=%0.2e', p(3))];
    end

elseif strcmp(var1_name, 'EL x GT x UC')

    if isfield(storages, 'EL_storage'); EL_storage = storages.EL_storage; end
    if isfield(storages, 'GT_storage'); GT_storage = storages.GT_storage; end
    if isfield(storages, 'UC_storage'); UC_storage = storages.UC_storage; end
    EL_total_vec = cat(1, EL_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
    GT_total_vec = cat(1, GT_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
    UC_total_vec = cat(1, UC_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
    [p,tbl,stats,~] = anovan(var2_total_vec, ...
        {EL_total_vec, GT_total_vec, UC_total_vec, id_total_vec}, ...
        'model', 'interaction', ...
        'varnames', {'EL', 'GT', 'UC', 'sbj'}, ...
        'display', 'off');
    disp([var2_name var2_cond ' n-way ANOVA: ' var1_name var1_cond ' effects'])
    disp(tbl)
    fprintf(strjoin(stats.coeffnames(2:7), ', '));
    disp(stats.coeffs(2:7)');

    p = p([1:3 5 6 8]); % EL, GT, UC, EL*GT, EL*UC, GT*UC

    title_text = {};
    title_text3 = [];

    if p(1) < 0.05
        title_text = [title_text, ...
            sprintf('EL p=%0.2e (%.3g)', p(1), stats.coeffs(3))]; % coeff: EL==L effect
    end

    if p(2) < 0.05
        title_text = [title_text, ...
            sprintf('GT p=%0.2e (%.3g)', p(2), stats.coeffs(4))]; % coeff: G trial (GT==1) effect
    end

    if p(3) < 0.05
        title_text = [title_text, ...
            sprintf('UC p=%0.2e (%.3g)', p(3), stats.coeffs(7))]; % coeff: UC==H effect
    end

    if p(4) < 0.05
        title_text3 = [title_text3, ...
            sprintf('EL*GT p=%0.2e', p(4))];
    end

    if p(5) < 0.05
        title_text3 = [title_text3, ...
            sprintf(', EL*UC p=%0.2e', p(5))];
    end

    if p(6) < 0.05
        title_text3 = [title_text3, ...
            sprintf(', GT*UC p=%0.2e', p(6))];
    end

    if ~isempty(title_text3)
        title_text = [title_text, title_text3];
    end


elseif strcmp(var1_name, 'uLuH EL')

    % Three-way ANOVA
    if isfield(storages, 'EL_storage'); EL_storage = storages.EL_storage; end
    if isfield(storages, 'UC_storage'); UC_storage = storages.UC_storage; end
    EL_total_vec = cat(1, EL_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
    UC_total_vec = cat(1, UC_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
    [p,tbl,stats,~] = anovan(var2_total_vec, ...
        {UC_total_vec, EL_total_vec, id_total_vec}, ...
        'model', 'interaction', ...
        'varnames', {'UC', 'EL', 'sbj'}, ...
        'display', 'off');
    %             UC = reshape(repmat({'Low', 'High'}, 2 * idsize, 1), [], 1);
    %             EL = reshape(repmat({'Early', 'Late', 'Early', 'Late'}, idsize, 1), [], 1);
    %             id_vec = cellfun(@num2str, num2cell(repmat(1:idsize, 1, N_var1_global)'), ...
    %                 'UniformOutput', false);
    %             [p,tbl,stats,terms] = anovan(var2mean_vec, {UC, EL, id_vec}, ...
    %                 'model', 'interaction', ...
    %                 'varnames', {'UC', 'EL', 'sbj'}, ...
    %                 'display', 'off');
    disp([var2_name var2_cond ' 3-way ANOVA: ' var1_name var1_cond ' effects'])
    disp(tbl)
    fprintf(strjoin(stats.coeffnames(2:5), ', '));
    disp(stats.coeffs(2:5)');

    p = p([1, 2, 4]); % Session, UC, Session * UC
    %             title_text = [sprintf('UC p=%0.2e ', p(1)), ...
    %                 sprintf('EL p=%0.2e ', p(2)), ...
    %                 sprintf('interaction p=%0.2e', p(3))];

    title_text = {};

    if p(1) < 0.05
        title_text = [title_text, ...
            sprintf('UC p=%0.2e (%.3g)', p(1), stats.coeffs(3))]; % coeff: UC==2 effect
    end

    if p(2) < 0.05
        title_text = [title_text, ...
            sprintf('EL p=%0.2e (%.3g)', p(2), stats.coeffs(5))]; % coeff: EL==2 effect
    end

    if p(3) < 0.05
        title_text = [title_text, ...
            sprintf('UC*EL p=%0.2e', p(3))];
    end

elseif strcmp(var1_name, 'uLuH ucEL')

    % Three-way ANOVA
    UC = reshape(repmat({'Low', 'High'}, 2 * idsize, 1), [], 1);
    EL = reshape(repmat({'Early', 'Late', 'Early', 'Late'}, idsize, 1), [], 1);
    id_vec = cellfun(@num2str, num2cell(repmat(1:idsize, 1, N_var1_global)'), ...
        'UniformOutput', false);
    [p,tbl,stats,~] = anovan(var2mean_vec, {UC, EL, id_vec}, ...
        'model', 'interaction', ...
        'varnames', {'UC', 'ucEL', 'sbj'}, ...
        'display', 'off');
    disp([var2_name var2_cond ' 3-way ANOVA: ' var1_name var1_cond ' effects'])
    disp(tbl)
    fprintf(strjoin(stats.coeffnames(2:5), ', '));
    disp(stats.coeffs(2:5)');

    p = p([1, 2, 4]); % Session, UC, Session * UC
    title_text = [sprintf('UC p=%0.2e ', p(1)), ...
        sprintf('ucEL p=%0.2e ', p(2)), ...
        sprintf('interaction p=%0.2e', p(3))];

    %             % repeated ANOVA
    %             UC = {'Low', 'Low', 'High', 'High'}';
    %             EL = {'Early', 'Late', 'Early', 'Late'}';
    %             UC_x_EL = table(UC, EL);
    % %             id_vec = (1:idsize)';
    %             id_vec = cellfun(@num2str, num2cell((1:idsize)'), 'UniformOutput', false);
    % %             baseline = ones(idsize, 1);
    %             var2_LE = var2_means(:, 1);
    %             var2_LL = var2_means(:, 2);
    %             var2_HE = var2_means(:, 3);
    %             var2_HL = var2_means(:, 4);
    %             var2_tbl = table(id_vec, var2_LE, var2_LL, var2_HE, var2_HL);
    % %             var2_tbl = table(baseline, var2_LE, var2_LL, var2_HE, var2_HL);
    %             rm = fitrm(var2_tbl, 'var2_LE-var2_HL ~ id_vec', 'WithinDesign', UC_x_EL);
    % %             rm = fitrm(var2_tbl, 'var2_LE-var2_HL ~ baseline', 'WithinDesign', UC_x_EL);
    %             rm.Coefficients
    %             ranova_tbl = ranova(rm, 'WithinModel', 'UC+EL+UC:EL');
    %
    %             % significance check
    %             disp([var2_name var2_cond ' repeated ANOVA: ' var1_name var1_cond ' effects'])
    %             pVal_tbl = ranova_tbl(:, 5);
    %             for tbl_row = 1:size(pVal_tbl, 1)
    %                 if ~isempty(pVal_tbl{tbl_row, 1}) && pVal_tbl{tbl_row, 1} < 0.05
    %                     disp(ranova_tbl(tbl_row, 1:5))
    %                 end
    %             end
    %

elseif strcmp(var1_name, 'Goal(6,7,8)xUC')
    UC = reshape(repmat({'Low', 'Low', 'Low', 'High', 'High', 'High'}, idsize, 1), [], 1);
    Goal678 = reshape(repmat({'6', '7', '8', '6', '7', '8'}, idsize, 1), [], 1);
    id_vec = cellfun(@num2str, num2cell(repmat(1:idsize, 1, N_var1_global)'), ...
        'UniformOutput', false);
    [p,tbl,stats,~] = anovan(var2mean_vec, {UC, Goal678, id_vec}, ...
        'model', 'interaction', ...
        'varnames', {'UC', 'Goal678', 'sbj'}, ...
        'display', 'off');
    disp(' ')
    disp([var2_name var2_cond ' N-way ANOVA: ' var1_name var1_cond ' effects'])
    disp(tbl)
    fprintf(strjoin(stats.coeffnames(2:6), ', '));
    disp(stats.coeffs(2:6)');

    p = p([1, 2, 4]); % UC, Goal678, Interaction
    title_text = [sprintf('UC p=%0.2e ', p(1)), ...
        sprintf('G678 p=%0.2e ', p(2)), ...
        sprintf('UC*678 p=%0.2e', p(3))];

elseif strcmp(var1_name, 'uLuH Goal(6,7,8)Switch')
    
    % % N-way ANOVA (trial-level) 2023-12-09
    % if isfield(storages, 'UC_storage'); UC_storage = storages.UC_storage; end
    % if isfield(storages, 'GS_storage'); GS_storage = storages.GS_storage; end
    % UC_total_vec = cat(1, UC_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
    % GS_total_vec = cat(1, GS_storage{:}); % [N_sbj x 1] cell -> [(N_sbj*N_trial) x 1] cell
    % [p,tbl,stats,~] = anovan(var2_total_vec, ...
    %     {UC_total_vec, GS_total_vec, id_total_vec}, ...
    %     'model', 'interaction', ...
    %     'varnames', {'UC', 'GS', 'sbj'}, ...
    %     'display', 'off');
    % 
    % disp([var2_name var2_cond ' 3-way ANOVA: ' var1_name var1_cond ' effects'])
    % disp(tbl)
    % fprintf(strjoin(stats.coeffnames(2:5), ', '));
    % disp(stats.coeffs(2:5)');
    % 
    % p = p([1, 2, 4]); % UC, GS, UC * GS
    % 
    % title_text = {};
    % 
    % if p(1) < 0.05
    %     title_text = [title_text, ...
    %         sprintf('UC p=%0.2e (%.3g)', p(1), stats.coeffs(3))]; % coeff: high uncertainty (UC==High) effect
    % end
    % 
    % if p(2) < 0.05
    %     title_text = [title_text, ...
    %         sprintf('GS p=%0.2e (%.3g)', p(2), stats.coeffs(4))]; % coeff: goal switch (GS==Switch) effect
    % end
    % 
    % if p(3) < 0.05
    %     title_text = [title_text, ...
    %         sprintf('UC*GS p=%0.2e', p(3))];
    % end

    % two-way ANOVA (trial-averaged)
    UC = reshape(repmat({'Low', 'High'}, 2 * idsize, 1), [], 1);
    GS = reshape(repmat({'Switch', 'Stay', 'Switch', 'Stay'}, idsize, 1), [], 1);
    id_vec = cellfun(@num2str, num2cell(repmat(1:idsize, 1, N_var1_global)'), ...
        'UniformOutput', false);
    [p,tbl,stats,~] = anovan(var2mean_vec, {UC, GS, id_vec}, ...
        'model', 'interaction', ...
        'varnames', {'UC', 'GS', 'sbj'}, ...
        'display', 'off');
    p = p([1, 2, 4]); % UC, GS, UC * GS
    title_text = [sprintf('UC p=%0.2e ', p(1)), ...
        sprintf('GS p=%0.2e ', p(2)), ...
        sprintf('Interaction p=%0.2e', p(3))];
    disp(' ')
    disp([var2_name var2_cond ' N-way ANOVA: ' var1_name var1_cond ' effects'])
    disp(tbl)
    fprintf(strjoin(stats.coeffnames(2:5), ', '));
    disp(stats.coeffs(2:5)');
    for cfi = 2:5
        disp([stats.coeffnames(cfi) ' ' stats.coeffs(cfi)])
    end
    for cfi = 6+idsize-length(except_ids):6+idsize-length(except_ids)+3
        disp([stats.coeffnames(cfi) ' ' stats.coeffs(cfi)])
    end

    %             % repeated ANOVA
    %             UC = {'Low', 'Low', 'High', 'High'}';
    %             GS = {'Switch', 'Stay', 'Switch', 'Stay'}';
    %             UC_x_GS = table(UC, GS);
    % %             id_vec = (1:idsize)';
    %             id_vec = cellfun(@num2str, num2cell((1:idsize)'), 'UniformOutput', false);
    % %             baseline = ones(idsize, 1);
    %             var2_LW = var2_means(:, 1);
    %             var2_LT = var2_means(:, 2);
    %             var2_HW = var2_means(:, 3);
    %             var2_HT = var2_means(:, 4);
    %             var2_tbl = table(id_vec, var2_LW, var2_LT, var2_HW, var2_HT);
    % %             var2_tbl = table(baseline, var2_LE, var2_LL, var2_HE, var2_HL);
    %             rm = fitrm(var2_tbl, 'var2_LW-var2_HT ~ id_vec', 'WithinDesign', UC_x_GS);
    % %             rm = fitrm(var2_tbl, 'var2_LE-var2_HL ~ baseline', 'WithinDesign', UC_x_EL);
    %             rm.Coefficients
    %             ranova_tbl = ranova(rm, 'WithinModel', 'UC+GS+UC:GS');
    %
    %             % significance check
    %             disp([var2_name ' repeated ANOVA: ' var1_name ' effects'])
    %             pVal_tbl = ranova_tbl(:, 5);
    %             for tbl_row = 1:size(pVal_tbl, 1)
    %                 if ~isempty(pVal_tbl{tbl_row, 1}) && pVal_tbl{tbl_row, 1} < 0.05
    %                     disp(ranova_tbl(tbl_row, 1:5))
    %                 end
    %             end

elseif strcmp(var1_name, 'Session x GT')

    % 3-way ANOVA: session, GT, subject effect
    Session_vec = reshape(repmat({'1', '2', '3', '4', '5', '1', '2', '3', '4', '5'}, idsize, 1), [], 1); % column
    GT_vec = reshape(repmat({'G', 'G', 'G', 'G', 'G', 'H', 'H', 'H', 'H', 'H'}, idsize, 1), [], 1); % column
    id_vec = cellfun(@num2str, num2cell(repmat(1:idsize, 1, N_var1_global)'), ...
        'UniformOutput', false);

    [p,tbl,stats,~] = anovan(var2mean_vec, {Session_vec, GT_vec, id_vec}, ...
        'model', 'interaction', ...
        'varnames', {'Session', 'GT', 'ID'}, ...
        'display', 'off');
    disp(tbl)
    p = p([1, 2, 4]); % Session, GT, Session * GT

    % Significance
    title_text = {};

    if p(1) < 0.05
        title_text = [title_text, ...
            sprintf('Ss p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(2:6))]; % coeff: Session effect
    end

    if p(2) < 0.05
        title_text = [title_text, ...
            sprintf('GT p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(7:8))]; % coeff: GT effect
    end

    if p(3) < 0.05
        title_text = [title_text, ...
            sprintf('Ss*GT p=%0.2e', p(3))];
    end

    % Visualization
    figure('Position', [0 1080/2-80 1920/4 1080/2])

    %             % Session 1~5, GT low
    %             plot(var2_means(:, 1:5)', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.9 .6 .6]); hold on;
    %             % Session 1~5, GT high
    %             plot(var2_means(:, 6:10)', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.6 .6 .9]);

    boxplot(var2_means(:, 1:5), ... % Session 1~5, GT low
        'Colors', [.9 .4 .4], ...
        'Width', 0.3, ...
        'Notch', 'on'); hold on;
    boxplot(var2_means(:, 6:10), ... % Session 1~5, GT high
        'Colors', [.4 .4 .9], ...
        'Width', 0.3, ...
        'Notch', 'on');
    %             legend({'GT low', '1', '11', 'OL', 'GT high', '2', '22', 'OL'})
    %             legend({'GT high'})
    hold off;
    box off

    title(cat(2, {[var2_name var2_cond ' per ' var1_name var1_cond]}, ...
        title_text, title_text2))

    % Gap visualization
    figure('Position', [0 1080/2-80 1920/4 1080/2])

    % Absolute gap

    % Additional ANOVA: session effect on GT gap
    var2mean_vec2 = reshape(abs(var2_means(:, 1:5) - var2_means(:, 6:10)), [], 1); % N_sbj x N_Ss
    [p2,~,~,~] = anovan(var2mean_vec2, {Session_vec(1:end/2), id_vec(1:end/2)}, ...
        'model', 'linear', ...
        'varnames', {'Session', 'ID'}, ...
        'display', 'off');
    if p2(1) < 0.05
        title_text3 = sprintf('Ss efct on gap p=%0.2e', p2(1));
    else
        title_text3 = [];
    end

    plot(abs(var2_means(:, 1:5)' - var2_means(:, 6:10)'), 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.9 .9 .9]); hold on
    boxplot(abs(var2_means(:, 1:5) - var2_means(:, 6:10)), ... % Session 1~5, GT G - H
        'Colors', [.9 .5 .9], ...
        'Width', 0.3, ...
        'Notch', 'on'); hold off;
    box off

    title(cat(2, {[var2_name var2_cond '|GT G - H| per Session']}, ...
        title_text3))


elseif strcmp(var1_name, 'Session x UC')

    % 3-way ANOVA: session, UC, subject effect
    Session_vec = reshape(repmat({'1', '2', '3', '4', '5', '1', '2', '3', '4', '5'}, idsize, 1), [], 1); % column
    UC_vec = reshape(repmat({'L', 'L', 'L', 'L', 'L', 'H', 'H', 'H', 'H', 'H'}, idsize, 1), [], 1); % column
    id_vec = cellfun(@num2str, num2cell(repmat(1:idsize, 1, N_var1_global)'), ...
        'UniformOutput', false);

    [p,tbl,stats,~] = anovan(var2mean_vec, {Session_vec, UC_vec, id_vec}, ...
        'model', 'interaction', ...
        'varnames', {'Session', 'UC', 'ID'}, ...
        'display', 'off');
    disp(tbl)
    p = p([1, 2, 4]); % Session, UC, Session * UC

    % Significance
    title_text = {};

    if p(1) < 0.05
        title_text = [title_text, ...
            sprintf('Ss p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(2:6))]; % coeff: Session effect
    end

    if p(2) < 0.05
        title_text = [title_text, ...
            sprintf('UC p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(7:8))]; % coeff: UC effect
    end

    if p(3) < 0.05
        title_text = [title_text, ...
            sprintf('Ss*UC p=%0.2e', p(3))];
    end

    % Visualization
    figure('Position', [0 1080/2-80 1920/4 1080/2])

    %             % Session 1~5, UC low
    %             plot(var2_means(:, 1:5)', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.9 .6 .6]); hold on;
    %             % Session 1~5, UC high
    %             plot(var2_means(:, 6:10)', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.6 .6 .9]);

    boxplot(var2_means(:, 1:5), ... % Session 1~5, UC low
        'Colors', [.9 .4 .4], ...
        'Width', 0.3, ...
        'Notch', 'on'); hold on;
    boxplot(var2_means(:, 6:10), ... % Session 1~5, UC high
        'Colors', [.4 .4 .9], ...
        'Width', 0.3, ...
        'Notch', 'on');
    %             legend({'UC low', '1', '11', 'OL', 'UC high', '2', '22', 'OL'})
    %             legend({'UC high'})
    hold off;
    box off

    title(cat(2, {[var2_name var2_cond ' per ' var1_name var1_cond]}, ...
        title_text, title_text2))

    % Gap visualization
    figure('Position', [0 1080/2-80 1920/4 1080/2])

    %             % Normal gap

    %             % Additional ANOVA: session effect on UC gap
    %             var2mean_vec2 = reshape(var2_means(:, 1:5) - var2_means(:, 6:10), [], 1); % N_sbj x N_Ss
    %             [p2,tbl2,stats2,terms2] = anovan(var2mean_vec2, {Session_vec(1:end/2), id_vec(1:end/2)}, ...
    %                 'model', 'linear', ...
    %                 'varnames', {'Session', 'ID'}, ...
    %                 'display', 'off');
    %             title_text3 = sprintf('Ss efct on gap p=%0.2e', p2(1));
    %             plot(var2_means(:, 1:5)' - var2_means(:, 6:10)', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.9 .9 .9]); hold on
    %             boxplot(var2_means(:, 1:5) - var2_means(:, 6:10), ... % Session 1~5, UC low - high
    %                 'Colors', [.9 .5 .9], ...
    %                 'Width', 0.3, ...
    %                 'Notch', 'on'); hold off;
    %             box off
    %
    %             title(cat(2, {[var2_name var2_cond '(UC low - high) per Session']}, ...
    %             title_text3))


    % Absolute gap

    % Additional ANOVA: session effect on UC gap
    var2mean_vec2 = reshape(abs(var2_means(:, 1:5) - var2_means(:, 6:10)), [], 1); % N_sbj x N_Ss
    [p2,~,~,~] = anovan(var2mean_vec2, {Session_vec(1:end/2), id_vec(1:end/2)}, ...
        'model', 'linear', ...
        'varnames', {'Session', 'ID'}, ...
        'display', 'off');
    if p2(1) < 0.05
        title_text3 = sprintf('Ss efct on gap p=%0.2e', p2(1));
    else
        title_text3 = [];
    end

    plot(abs(var2_means(:, 1:5)' - var2_means(:, 6:10)'), 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.9 .9 .9]); hold on
    boxplot(abs(var2_means(:, 1:5) - var2_means(:, 6:10)), ... % Session 1~5, UC low - high
        'Colors', [.9 .5 .9], ...
        'Width', 0.3, ...
        'Notch', 'on'); hold off;
    box off

    title(cat(2, {[var2_name var2_cond '|UC low - high| per Session']}, ...
        title_text3))

    %
    %
    %             % repeated ANOVA
    %             Session = {'1', '2', '3', '4', '5'}; Session_vec = reshape(repmat(Session, 1, 2), [], 1); % column
    %             UC = {'Low', 'High'}; UC_vec = reshape(repmat(UC, 5, 1), [], 1); % column
    %             id_vec = cellfun(@num2str, num2cell((1:idsize)'), 'UniformOutput', false); % within design
    %             y_name_vec = cell(1, idsize);
    %             for yi = 1:idsize
    %                 y_name_vec{yi} = ['y' num2str(yi)];
    %             end
    %             var2_mean_cell = num2cell(var2_means'); % [(N_UC * N_session) x N_sbj]
    %
    %             vartype_cell = cat(2, {'cellstr'}, {'cellstr'}, repmat({'double'}, 1, idsize)); % [1 x (N_sbj+2)]
    %             varname_cell = cat(2, {'Session'}, {'UC'}, y_name_vec); % [1 x (N_sbj+2)]
    %
    %             var2_tbl = table('Size', [2*5 idsize+2], 'VariableTypes', vartype_cell, 'VariableNames', varname_cell);
    %             var2_tbl(:,:) = cat(2, Session_vec, UC_vec, var2_mean_cell);
    %
    %             rm = fitrm(var2_tbl, sprintf('y1-y%d ~ Session+UC+Session:UC', idsize), ...
    %                 'WithinDesign', id_vec);
    %             rm.Coefficients
    %             ranova_tbl = ranova(rm, 'WithinModel', 'id_vec');
    %
    %             % significance check
    %             disp([var2_name var2_cond ' repeated ANOVA: ' var1_name var1_cond ' effects'])
    %             pVal_tbl = ranova_tbl(:, 5);
    %             for tbl_row = 1:size(pVal_tbl, 1)
    %                 if ~isempty(pVal_tbl{tbl_row, 1}) && pVal_tbl{tbl_row, 1} < 0.05
    %                     disp(ranova_tbl(tbl_row, 1:5))
    %                 end
    %             end
    %

elseif strcmp(var1_name, 'Session x Goal(6,7,8)')

    % 3-way ANOVA: session, G678, subject effect
    Session_vec = reshape(repmat({'1', '2', '3', '4', '5', '1', '2', '3', '4', '5', '1', '2', '3', '4', '5'}, idsize, 1), [], 1);
    G678_vec = reshape(repmat({'6', '6', '6', '6', '6', '7', '7', '7', '7', '7', '8', '8', '8', '8', '8'}, idsize, 1), [], 1);

    [p,tbl,stats,~] = anovan(var2mean_vec, {Session_vec, G678_vec, id_vec}, ...
        'model', 'interaction', ...
        'varnames', {'Session', 'G678', 'ID'}, ...
        'display', 'off');
    disp(tbl)
    p = p([1, 2, 4]); % Session, G678, Session * G678
    %
    %             title_text = [sprintf('Session p=%0.2e ', p(1)), ...
    %                 sprintf('G678 p=%0.2e ', p(2)), ...
    %                 sprintf('interaction p=%0.2e', p(3))];
    %
    title_text = {};

    if p(1) < 0.05
        title_text = [title_text, ...
            sprintf('Ss p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(2:6))]; % coeff: Session effect
    end

    if p(2) < 0.05
        title_text = [title_text, ...
            sprintf('G678 p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(7:9))]; % coeff: G678 effect
    end

    if p(3) < 0.05
        title_text = [title_text, ...
            sprintf('Ss*G678 p=%0.2e', p(3))];
    end

elseif strcmp(var1_name, 'Session x GS')

    % 3-way ANOVA: session, G678, subject effect
    Session_vec = reshape(repmat({'1', '2', '3', '4', '5', '1', '2', '3', '4', '5'}, idsize, 1), [], 1);
    GS_vec = reshape(repmat({'0', '0', '0', '0', '0', '1', '1', '1', '1', '1'}, idsize, 1), [], 1);

    [p,tbl,stats,~] = anovan(var2mean_vec, {Session_vec, GS_vec, id_vec}, ...
        'model', 'interaction', ...
        'varnames', {'Session', 'GS', 'ID'}, ...
        'display', 'off');
    disp(tbl)
    p = p([1, 2, 4]); % Session, G678, Session * G678

    title_text = {};

    if p(1) < 0.05
        title_text = [title_text, ...
            sprintf('Ss p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(2:6))]; % coeff: Session effect
    end

    if p(2) < 0.05
        title_text = [title_text, ...
            sprintf('GS p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(8))]; % coeff: goal switch (GS==1) effect
    end

    if p(3) < 0.05
        title_text = [title_text, ...
            sprintf('Ss*GS p=%0.2e', p(3))];
    end

elseif strcmp(var1_name, 'EL x GS')

    % 3-way ANOVA: session, G678, subject effect
    EL_vec = reshape(repmat({'1', '2', '1', '2'}, idsize, 1), [], 1);
    GS_vec = reshape(repmat({'0', '0', '1', '1'}, idsize, 1), [], 1);

    [p,tbl,stats,~] = anovan(var2mean_vec, {EL_vec, GS_vec, id_vec}, ...
        'model', 'interaction', ...
        'varnames', {'Session', 'GS', 'ID'}, ...
        'display', 'off');
    disp(tbl)
    p = p([1, 2, 4]); % Session, G678, Session * G678

    title_text = {};

    if p(1) < 0.05
        title_text = [title_text, ...
            sprintf('EL p=%0.2e', p(1)), ...
            sprintf('(%.3g)', stats.coeffs(3))]; % coeff: learning (EL==2) effect
    end

    if p(2) < 0.05
        title_text = [title_text, ...
            sprintf('GS p=%0.2e', p(2)), ...
            sprintf('(%.3g)', stats.coeffs(5))]; % coeff: goal switch (GS==1) effect
    end

    if p(3) < 0.05
        title_text = [title_text, ...
            sprintf('EL*GS p=%0.2e', p(3))];
    end

end

end % function end

%% helper function: load var1 and var2

function [var1, var2, var1_unq, N_var1] = load_var1var2(Exp, id, ... 
    var1_name, var2_name, ... 
    var1_cond, var2_cond, ...
    NameValue1, NameValue2, ... 
    varargin)

% for 'taskVar profile (continuous)' section

% Default input parameters
options = struct('categorize1', [], ...
    'categorize2', [], ...
    'var1_unq_pre', []); 
% Read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(load_var1var2) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(load_var1var2) %s is not a recognized parameter name', pair{1})
    end
end

categorize1 = options.categorize1;
categorize2 = options.categorize2;
var1_unq_pre = options.var1_unq_pre;

if strcmp(var1_name(end-1:end), ' G')
    var1_name_new = erase(var1_name, ' G');
else
    var1_name_new = var1_name;
end
switch var1_name_new
    case 'Half'
        var1 = arbMBMF_load_var(Exp, 'Session', id, var1_cond);
        var1(1:floor(length(var1)/2)) = 1;
        var1(ceil(length(var1)/2):end) = 2;
    case 'UncCond G'
        var1 = arbMBMF_load_var(Exp, 'UncCond', id, 'G');
    case 'uLuH EL'
        UC = arbMBMF_load_var(Exp, 'UncCond', id, var1_cond);
        EL = arbMBMF_load_var(Exp, 'Half', id, var1_cond); % 1, 2
        var1 = nan(size(UC));
        var1(UC == 1 & EL == 1) = 1; % uL, early
        var1(UC == 2 & EL == 1) = 2; % uH, early
        var1(UC == 1 & EL == 2) = 3; % uL, late
        var1(UC == 2 & EL == 2) = 4; % uH, late
    case 'uLuH ucEL'
        UC = arbMBMF_load_var(Exp, 'UncCond', id, var1_cond);
        ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var1_cond);
        var1 = nan(size(UC));
        var1(UC == 1 & ucEL == 0) = 1; % uL, early
        var1(UC == 1 & ucEL == 1) = 2; % uL, late
        var1(UC == 2 & ucEL == 0) = 3; % uH, early
        var1(UC == 2 & ucEL == 1) = 4; % uH, late
    case 'uLuH Goal(6,7,8)Switch'
        UC = arbMBMF_load_var(Exp, 'UncCond', id, var1_cond);
        GS = arbMBMF_load_var(Exp, 'Goal(6,7,8)Switch', id, var1_cond);
        var1 = nan(size(UC));
        var1(UC == 1 & GS == 1) = 1; % uL, switch
        var1(UC == 1 & GS == 0) = 2; % uL, stay
        var1(UC == 2 & GS == 1) = 3; % uH, switch
        var1(UC == 2 & GS == 0) = 4; % uH, stay
    case 'uLL vs uHE' % uLL: 1, uHE: 2
        UC = arbMBMF_load_var(Exp, 'UncCond', id, var1_cond);
        ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var1_cond);
        var1 = nan(size(UC));
        var1(UC == 1 & ucEL == 1) = 1; % uL, late
        var1(UC == 2 & ucEL == 0) = 2; % uH, early
    case 'uHL vs uLE'
        UC = arbMBMF_load_var(Exp, 'UncCond', id, var1_cond);
        ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var1_cond);
        var1 = nan(size(UC));
        var1(UC == 2 & ucEL == 1) = 1; % uH, late
        var1(UC == 1 & ucEL == 0) = 2; % uL, early
    otherwise
        if isempty(NameValue1)
            var1 = arbMBMF_load_var(Exp, var1_name_new, id, var1_cond);
        else
            var1 = arbMBMF_load_var(Exp, var1_name_new, id, var1_cond, NameValue1{:});
            % disp('test')
        end
end
%     if strcmp(var1_name(end-1:end), ' G')
%         GT = arbMBMF_load_var(Exp, 'GoalCond', id, []);
%         var1(~ismember(GT, 1)) = nan;
%         disp('Var1 limitation to G trials')
%     end
if ~isempty(categorize1)   % continuous to categorical
    [var1, cutoff] = regr2class(var1, categorize1, 'Cutoff', 'median', 'Verbose', 1);
    %         [var1, cutoff] = regr2class(var1, categorize1, 'Cutoff', 'mean', 'Verbose', 1);
    %         [var1, cutoff] = regr2class(var1, categorize1);
    %         [~, N_el] = unq_elms(var1);
    %         disp([cutoff.type ' ' num2str(N_el')])
end

% Classes
if isempty(var1_unq_pre) % default setting - var1 from data

    var1_unq = unq_elms(var1);

else % special cases - predefined var1 classes

    if any(~ismember(unq_elms(var1), var1_unq))
        error('(arbMBMF_taskVar_profile) exceptions to predefined var1_unq')
    end

end

N_var1 = length(var1_unq);

switch var2_name
    case 'ChoCon1'
        var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
        var2 = var2(1,:);
    case 'ChoCon2'
        var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
        var2 = var2(2,:);
    otherwise
        if isempty(NameValue2)
            var2 = arbMBMF_load_var(Exp, var2_name, id, var2_cond);
        else
            var2 = arbMBMF_load_var(Exp, var2_name, id, var2_cond, NameValue2{:});
            % disp('test')
        end
end
if ~isempty(categorize2)   % continuous to categorical
    [var2, cutoff] = regr2class(var2, categorize2);
    [~, N_el] = unq_elms(var2);
    disp([cutoff.type ' ' num2str(N_el')])
end
%     [~, N_el] = unq_elms(var2);
%     disp([num2str(N_el')])

%     if strcmp(var2_name, 'ChoOpt1n2') || strcmp(var2_name, 'ChoConsist1n2') ...
%             || strcmp(var2_name, 'ChoCon ignoring context') || strcmp(var2_name, 'ChoCon ignoring context G')
%         var2 = reshape(var2, 1, []);
%         var1 = reshape([var1; var1], 1, []);
%     end

if size(var2, 1) == 2
    var2 = reshape(var2, 1, []);
    var1 = reshape([var1; var1], 1, []);
end

end % function end
