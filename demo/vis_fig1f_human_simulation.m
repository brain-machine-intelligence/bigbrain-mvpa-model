
clear all

%% Set working directory

base_path = fullfile('..'); % If your current working directory is 'demo', use this. Otherwise, set your own path.
addpath(fullfile(base_path, 'functions', 'utils'))

% If your current working directory is 'demo', use this. Otherwise, set your own path that includes 'data'.
cd('..'); 


%% Load human data

clear
LoadDefaultSettings;

% ============================= Data loading ============================

% experiment info
sbjID = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,23,24}; % full 20 subjects
validIdx = 1:20;
nSbj = length(validIdx);
sbjNames = {'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
    'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};  % Lee2014
nSessions = 5*ones(1,length(sbjNames)); nSessions(11)=4; nSessions(16)=3; nSessions(20)=2; nSessions(21)=4;     % Lee2014

% load human data (session combined)
humanRoot = 'data';
humanDir = 'human_behavior';
humanName = 'SBJ_structure_each_exp_BESTcollectionJan29.mat';
load([humanRoot '/' humanDir '/' humanName], 'SBJ'); % 22 subjects data

% make the cell(sbj) of behavior matrices
humanBhvs = cell(nSbj, 1);
humanBlks = cell(nSbj, 1);
humanTensor = cell(nSbj, 1);
for i = 1:nSbj   % subject loop
    humanBhvs{i} = cell2mat(SBJ{validIdx(i)}.HIST_behavior_info')'; % [18 x nTrial] double
    humanBlks{i} = cell2mat(SBJ{validIdx(i)}.HIST_block_condition); % [2 x nTrial] double
    
    hs1 = humanBhvs{i}(4, :);
    ha1 = humanBhvs{i}(7, :);
    hs2 = humanBhvs{i}(5, :);
    ha2 = humanBhvs{i}(8, :);
    hGoal = humanBhvs{i}(18, :);
    hBlk = humanBlks{i}(2, :);

    hGT = (hBlk > 2.5) + 1; % G(1): 1 & 2, H(2): 3 & 4
    hUC = ismember(hBlk, [2 3]) + 1; % 1:Low, 2:High

    % ******************** choice optimality ********************
    [opt, opt1, opt2, optnan, optnan1, optnan2] = optimality_temp(hs1, ha1, hs2, ha2, hGoal, hBlk);

    % ******************** choice versatility ********************
    [flex, flex1, flex2] = versatility_temp(ha1, hs2, ha2, hGoal);

    % ******************** choice consistency ********************
    [stab, stab1, stab2] = consistency_temp(ha1, hs2, ha2, hGoal);

    % sanity check
    tempOpt = arbMBMF_load_var('Lee2014', 'ChoOpt1n2', validIdx(i), []);                 tempOpt1 = tempOpt(1, :); tempOpt2 = tempOpt(2, :); tempOpt = mean(tempOpt, 1, 'omitnan');
    tempOptNan = arbMBMF_load_var('Lee2014', 'ChoOptBinNaN1n2', validIdx(i), []);        tempOptNan1 = tempOptNan(1, :); tempOptNan2 = tempOptNan(2, :); tempOptNan = mean(tempOptNan, 1, 'omitnan');
    tempFlex = arbMBMF_load_var('Lee2014', 'ChoSwitch1n2(GoalSwitch)', validIdx(i), []); tempFlex1 = tempFlex(1, :); tempFlex2 = tempFlex(2, :); tempFlex = mean(tempFlex, 1, 'omitnan');
    tempStab = arbMBMF_load_var('Lee2014', 'ChoConsist1n2', validIdx(i), []);            tempStab1 = tempStab(1, :); tempStab2 = tempStab(2, :); tempStab = mean(tempStab, 1, 'omitnan');
    if any(opt ~= tempOpt) || any(opt1 ~= tempOpt1) || any(opt2 ~= tempOpt2); error('Optimality error'); end
    if any(optnan(~isnan(optnan)) ~= tempOptNan(~isnan(tempOptNan))) || any(optnan1(~isnan(optnan1)) ~= tempOptNan1(~isnan(tempOptNan1))) || any(optnan2(~isnan(optnan2)) ~= tempOptNan2(~isnan(tempOptNan2))); error('OptimalityNaN error'); end
    if any(flex(~isnan(flex)) ~= tempFlex(~isnan(tempFlex))) || any(flex1(~isnan(flex1)) ~= tempFlex1(~isnan(tempFlex1))) || any(flex2(~isnan(flex2)) ~= tempFlex2(~isnan(tempFlex2))); error('Flexibility error'); end
    if any(stab(~isnan(stab)) ~= tempStab(~isnan(tempStab))) || any(stab1(~isnan(stab1)) ~= tempStab1(~isnan(tempStab1))) || any(stab2(~isnan(stab2)) ~= tempStab2(~isnan(tempStab2))); error('Stability error'); end

    % storing
    humanTensor{i} = [hs1; hs2; ha1; ha2; ... % 1~4
        hGT; hUC; hGoal; ... % 5~7
        optnan; flex; stab; ...  % 8~10
        optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
        opt; opt1; opt2; ... % 17~19
        hBlk]; % 20
end
% col1 - block #, col2 - trial # (in each block), col3 - block condition - 1: G(with low uncertainty), 2: G'(with high uncertainty), 3:H(with high uncertainty), 4:H'(with low uncertainty), 
% col4 - S1, col5 - S2, col6 - S3, col7 - A1 (action in S1), col8 - A2 (action in S2), 
% col9 - RT(A1), col10 - RT(A2), 
% col11 - onset (S1) from the trial start, col12 - onset (S2) from the trial start, col13 - onset (S3) from the trial start, 
% col14 - onset (A1) from the trial start, col15 - onset (A2) from the trial start, 
% col16 - reward amount (0/10/20/40) at S3, col17 - total amount (total in the current session), 
% col18 - goal state=outcome state. 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state,  
disp('Human loading done')

%% Load simulation data (S1, S2, A1, A2, Qcho1, Qcho2 Qunc1, Qunc2, Model1, Model2)

% simulRoot = pwd;
simulRoot = 'data';
simulDirs = {'mb_virtual_episode1000', 'mf_virtual_episode1000'}; nSimul = 1000; nParam = 6;
% simulDirs = {'mb_virtual_episode100', 'mf_virtual_episode100'}; nSimul = 100; nParam = 6;
simulName = 'SBJ';
nFeat = 29; % # of features
nMdl = length(simulDirs);
epiCell = cell(nSbj, nMdl);
paramMat = nan(nMdl, nSbj, nSimul, nParam);
disp('Virtual simulation result loading')
for sim = 1:nMdl
    fprintf(simulDirs{sim})

    for i = 1:nSbj

        % sanity check - block condition
        if any(humanBhvs{i}(3, :) ~= humanBlks{i}(2, :)); error('Inconsistent block condition in SBJ!'); end

        nTrial = size(humanBhvs{i}, 2); % *** nTrial from the human data, for validation ***
        epiTensor = nan(nSimul, nFeat, nTrial);
        load([simulRoot '/' simulDirs{sim} '/' simulName '_' sbjNames{sbjID{validIdx(i)}}], 'SBJtot') % a single simulation from one subject experiment
        for j = 1:nSimul

            % sanity check - block condition (human vs virtual)
            tempBlk = cell2mat(SBJtot{j}{1}.HIST_block_condition); tempBlk = tempBlk(2, :);
            if any(tempBlk ~= humanBlks{i}(2, :)); error('Inconsistent block condition: human vs virtual'); end
            tempGT = (tempBlk > 2.5) + 1; % G(1): 1 & 2, H(2): 3 & 4
            tempUC = ismember(tempBlk, [2 3]) + 1; % 1:Low, 2:High

            % sanity check - goal condition (human vs virtual)
            tempBhv = cell2mat(SBJtot{j}{1}.HIST_behavior_info')';
            tempGoal = tempBhv(18, :);
            if any(humanBhvs{i}(18, :) ~= tempGoal); error('Inconsistent goal condition: human vs virtual'); end

            % parameter saving
            if isfield(SBJtot{j}{1}.model_BayesArb, 'param')
                paramMat(sim, i, j, :) = SBJtot{j}{1}.model_BayesArb.param(1:nParam);
            end

            % episode info.: states, actions, chosen value, unchosen value, active model (fwd vs sarsa)
            tempEpisode = cell2mat(SBJtot{j}{1}.model_BayesArb.episode(:)); % [N_session x 1] struct
            ss = [tempEpisode.HIST_ArbState]; ss1 = ss(1:2:end); ss2 = ss(2:2:end); % simulation states, [nTrial x 1] double
            sa = [tempEpisode.HIST_ArbAction]; sa1 = sa(1:2:end); sa2 = sa(2:2:end); % simulation actions
            sqc = [tempEpisode.HIST_ArbQcho]; sqc1 = sqc(1:2:end); sqc2 = sqc(2:2:end); % simulation Q(chosen)
            squ = [tempEpisode.HIST_ArbQunc]; squ1 = squ(1:2:end); squ2 = squ(2:2:end); % simulation Q(unchosen)
            sm = [tempEpisode.HIST_ArbMode]; sm1 = sm(1:2:end); sm2 = sm(2:2:end); % simulation active model
            if nFeat > 29
                nll = [tempEpisode.HIST_NegLogLik]; nll1 = nll(1:2:end); nll2 = nll(2:2:end); % human choice negative log likelihood
                lik = exp(-nll); lik1 = lik(1:2:end); lik2 = lik(2:2:end); % human choice likelihood
            end

            % ******************** choice optimality ********************
            [opt, opt1, opt2, optnan, optnan1, optnan2] = optimality_temp(ss1, sa1, ss2, sa2, tempGoal, tempBlk);

            % ******************** choice versatility ********************
            [flex, flex1, flex2] = versatility_temp(sa1, ss2, sa2, tempGoal);

            % ******************** choice consistency ********************
            [stab, stab1, stab2] = consistency_temp(sa1, ss2, sa2, tempGoal);

            % storing the values

            if nFeat > 29
                epiTensor(j, :, :) = [ss1; ss2; sa1; sa2; ... % 1~4
                    tempGT; tempUC; tempGoal; ... % 5~7
                    optnan; flex; stab; ...  % 8~10
                    optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
                    opt; opt1; opt2; ... % 17~19
                    tempBlk; ... % 20
                    sqc1; sqc2; squ1; squ2; ... % 21~24
                    abs(sqc1 - squ1); abs(sqc2 - squ2); mean([abs(sqc1 - squ1); abs(sqc2 - squ2)], 1); ... % 25~27
                    sm1; sm2; ... % 28~29
                    mean([nll1; nll2], 1); nll1; nll2; ... % 30~32
                    mean([lik1; lik2], 1); lik1; lik2]; % 33~35
            else

                epiTensor(j, :, :) = [ss1; ss2; sa1; sa2; ... % 1~4
                    tempGT; tempUC; tempGoal; ... % 5~7
                    optnan; flex; stab; ...  % 8~10
                    optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
                    opt; opt1; opt2; ... % 17~19
                    tempBlk; ... % 20
                    sqc1; sqc2; squ1; squ2; ... % 21~24
                    abs(sqc1 - squ1); abs(sqc2 - squ2); mean([abs(sqc1 - squ1); abs(sqc2 - squ2)], 1); ... % 25~27
                    sm1; sm2]; % 28~29
            end

            if nTrial ~= size(epiTensor(j, :, :), 3); error('Different nTrial for human vs simulation'); end
        end
        epiCell{i, sim} = epiTensor;
        fprintf(' %d', i)
    end

    disp(' ')
end

disp(' loading done')



%% Visualization (scatter plot)

close all

disp(' ')
disp('Start')
disp(' ')

% humanTensor{i} = [hs1; hs2; ha1; ha2; ... % 1~4
%     hGT; hUC; hGoal; ... % 5~7
%     optnan; flex; stab; ...  % 8~10
%     optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
%     opt; opt1; opt2]; % 17~19

% epiCell{N_sbj}(nSimul, nFeat, nTrial)

% epiCell{i} = epiTensor;
% epiTensor(j, :, :) = [ss1; ss2; sa1; sa2; ... % 1~4
%             tempGT; tempUC; tempGoal; ... % 5~7
%             optnan; flex; stab; ...  % 8~10
%             optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
%             opt; opt1; opt2; ... % 17~19
%             sqc1; sqc2; squ1; squ2; ... % 20~23
%             sm1; sm2]; % 24~25

subjRestVar = []; subjRestFeat = []; subjFeatRange = [];
% subjRestVar = 'Optimality'; subjRestFeat = 8; subjFeatRange = [.7614 .9932];

% trialRestVar = []; trialRestFeat = []; trialFeatSet = [];
trialRestVar = 'GoalCond'; trialRestFeat = 5; trialFeatSet = 1;
% trialRestVar = 'GoalCond'; trialRestFeat = 5; trialFeatSet = 2;
% trialRestVar = 'Uncertainty'; trialRestFeat = 6; trialFeatSet = 1;
% trialRestVar = 'Uncertainty'; trialRestFeat = 6; trialFeatSet = 2;
% trialRestVar = 'Block'; trialRestFeat = 20; trialFeatSet = 1;
% trialRestVar = 'Block'; trialRestFeat = 20; trialFeatSet = 2;
% trialRestVar = 'Block'; trialRestFeat = 20; trialFeatSet = 3;
% trialRestVar = 'Block'; trialRestFeat = 20; trialFeatSet = 4;

% Flexibility vs performance
name1 = 'Flexibility'; name2 = 'Performance'; feat1 = [13 14]; feat2 = [11 12]; name11 = 'Flexibility (choice versatility)'; name22 = 'Performance (choice optimality)';

% Stability vs performance
% name1 = 'Stability'; name2 = 'Performance'; feat1 = [15 16]; feat2 = [11 12]; name11 = 'Stability (choice consistency)'; name22 = 'Performance (choice optimality)';

% Flexibility vs stability
% name1 = 'Flexibility'; name2 = 'Stability'; feat1 = [13 14]; feat2 = [15 16]; name11 = 'Flexibility (choice versatility)'; name22 = 'Stability (choice consistency)';

figure('Position', [0 400-80 400 400]) % [left bottom width height]

% scatter virtual simulation results
for sim = 1:nMdl

    var1 = nan(nSbj, nSimul);
    var2 = nan(nSbj, nSimul);
    restVar = ones(nSbj, nSimul);

    for i = 1:nSbj
        
        var1_temp = epiCell{i, sim}(:, feat1, :); % (nSimul, nFeat, nTrial)
        var2_temp = epiCell{i, sim}(:, feat2, :); % (nSimul, nFeat, nTrial)
        if ~isempty(trialRestVar)
            trialRestIdx = ismember(epiCell{i, sim}(:, trialRestFeat, :), trialFeatSet); % [nSimul x 1 x nTrial]
            var1_temp(repmat(~trialRestIdx, 1, length(feat1), 1)) = nan;
            var2_temp(repmat(~trialRestIdx, 1, length(feat2), 1)) = nan;
        end

        var1(i, :) = mean(var1_temp, [2 3], 'omitnan');
        var2(i, :) = mean(var2_temp, [2 3], 'omitnan');
        if ~isempty(subjRestVar)
            restVar(i, :) = (mean(epiCell{i, sim}(:, subjRestFeat, :), 3, 'omitnan') > subjFeatRange(1)) ...
                & (mean(epiCell{i, sim}(:, subjRestFeat, :), 3, 'omitnan') < subjFeatRange(2));
        end
    end

    var1Vec = var1(logical(restVar(:)));
    var2Vec = var2(logical(restVar(:)));
    [simulColor, sbjColor] = meshgrid(1:nSimul, 1:nSbj); % x and y grid
    simulColor = simulColor(logical(restVar(:)));
    sbjColor = sbjColor(logical(restVar(:)));

    % % sanity check
    % p1 = anovan(var1Vec, {simulColor, sbjColor})
    % p2 = anovan(var2Vec, {simulColor, sbjColor})
    
     [R, p, stats] = scatter_corr(var1Vec, var2Vec, name1, name2, name11, name22, ...
        'MarkerSize', 5, ...
        'LinearTrend', false, ...
        'ColorWeight', 2.25 - .75 * sim, ...
        'MarkerFaceAlpha', .1);

     % To print statistic and p value
     p = p(1,2);
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
     disp(simulDirs{sim})
     fprintf(['r=%.3f, t(%d)=%.3f\t p=%.3e (' sig_text ')\n'], R(1,2), stats.df, stats.tstat, p)

    % scatter_corr(var1Vec, var2Vec, name1, name2, name11, name22, ...
    %     'MarkerSize', 5, ...
    %     'LinearTrend', false, ...
    %     'ColorScalars', colors(7 + sim, :), ...
    %     'MarkerFaceAlpha', .1);
    hold on
end


% scatter human results
hvar1 = nan(nSbj, 1);
hvar2 = nan(nSbj, 1);

for i = 1:nSbj

    hvar1_temp = humanTensor{i}(feat1, :); % [N_feat x nTrial]
    hvar2_temp = humanTensor{i}(feat2, :); % [N_feat x nTrial]
    if ~isempty(trialRestVar)
        trialRestIdx = ismember(humanTensor{i}(trialRestFeat, :), trialFeatSet); % [1 x nTrial]
        hvar1_temp(repmat(~trialRestIdx, length(feat1), 1)) = nan;
        hvar2_temp(repmat(~trialRestIdx, length(feat2), 1)) = nan;
    end
    
    hvar1(i) = mean(hvar1_temp, [1 2], 'omitnan');
    hvar2(i) = mean(hvar2_temp, [1 2], 'omitnan');

end

[R, p, stats] = scatter_corr(hvar1, hvar2, name1, name2, name11, name22, ...
        'ExperimentName', ['Human, MB, MF ' trialRestVar ' ' num2str(trialFeatSet)], ...
        'MarkerSize', 15, ...
        'LinearTrend', false, ...
        'ColorScalars', colors(1, :));
hold off;

% To print statistic and p value
p = p(1,2);
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
disp('Human')
fprintf(['r=%.3f, t(%d)=%.3f\t p=%.3e (' sig_text ')\n'], R(1,2), stats.df, stats.tstat, p)

disp(' ')
disp('End')
disp(' ')



%% Helper functions

% ******************** choice optimality ********************
function [opt, opt1, opt2, optnan, optnan1, optnan2] = optimality_temp(s1, a1, s2, a2, goals, blocks)

nTrial = length(s1);

% sanity check
if length(a1) ~= nTrial; error('S1 length ~= A1 length'); end
if length(s2) ~= nTrial; error('S1 length ~= S2 length'); end
if length(a2) ~= nTrial; error('S1 length ~= A2 length'); end
if length(goals) ~= nTrial; error('S1 length ~= Goals length'); end
if length(blocks) ~= nTrial; error('S1 length ~= Blocks length'); end

opt1 = nan(1, nTrial); optnan1 = nan(1, nTrial);
opt2 = nan(1, nTrial); optnan2 = nan(1, nTrial);
for k = 1:nTrial

    opt1(k) = opt_normalization_1st2ndstage_exnan(s1(k), a1(k), goals(k), blocks(k));
    opt2(k) = opt_normalization_1st2ndstage_exnan(s2(k), a2(k), goals(k), blocks(k));

    if opt1(k)==.5
        optnan1(k)=NaN;
    elseif opt1(k)>.5
        optnan1(k)=1;
    else
        optnan1(k)=0;
    end

    if opt2(k)==.5
        optnan2(k)=NaN;
    elseif opt2(k)>.5
        optnan2(k)=1;
    else
        optnan2(k)=0;
    end
end
opt = mean([opt1; opt2], 1, 'omitnan');
optnan = mean([optnan1; optnan2], 1, 'omitnan');
end

% ******************** choice versatility ********************
function [flex, flex1, flex2] = versatility_temp(a1, s2, a2, goals)

% ===== stage 1 =====
gs1 = [1 diff(goals)~=0]; % goal switch

flex1 = 1 - [0 diff(a1)==0]; % A1 choice switch
flex1(:, ~ismember(gs1, 1)) = nan; % for A1 choice switch at goal switch

% ===== stage 2 =====
flex2 = nan(1, length(a2));
gs2 = nan(1, length(a2));
S2_set = [2 3 4 5];
for si = 1:length(S2_set)
    s2_idx = (s2 == S2_set(si));
    if any(s2_idx)
        a2_s2 = a2(s2_idx);
        flex2(s2_idx) = 1 - [0 diff(a2_s2)==0]; % A2 switch at S2
        g_s2 = goals(s2_idx);
        gs2(s2_idx) = 1 - [0 diff(g_s2)==0]; % Goal switch at S2
    end
end
if any(isnan(flex2)); error('stage 2 versatility error'); end % sanity check
flex2(:, ~ismember(gs2, 1)) = NaN; % for A2 switch at goal switch
flex = mean([flex1; flex2], 1, 'omitnan');
end

% ******************** choice consistency ********************
function [stab, stab1, stab2] = consistency_temp(a1, s2, a2, goals)

nTrial = length(a1);

% sanity check
if length(s2) ~= nTrial; error('A1 length ~= S2 length'); end
if length(a2) ~= nTrial; error('A1 length ~= A2 length'); end
if length(goals) ~= nTrial; error('A1 length ~= Goals length'); end

stab1 = nan(1, nTrial);
stab2 = nan(1, nTrial);
Goal_set = [6 7 8 -1];
for gi = 1:length(Goal_set)

    % ===== stage 1 =====
    g_idx = (goals == Goal_set(gi)); % goal index

    stab1(g_idx) = [0 diff(a1(g_idx))==0]; % stage 1 choice consistency

    % ===== stage 2 =====
    S2_temp = s2(g_idx);
    A2_temp = a2(g_idx);
    ChoCon2_temp = NaN(1,length(A2_temp));
    S2_set = [2 3 4 5];
    for si = 1:length(S2_set)
        s2_idx = (S2_temp == S2_set(si));
        if any(s2_idx)
            a2_temp = A2_temp(s2_idx);
            ChoCon2_temp(s2_idx) = [0 diff(a2_temp)==0];
        end
    end
    stab2(g_idx) = ChoCon2_temp;
end
if any(isnan(stab2)); error('stage 2 consistency error'); end % sanity check
stab = mean([stab1; stab2], 1, 'omitnan');
end

% ******************** visualization by scatter plot ********************
function [R, p, stats, mdl, tbl] = scatter_corr(var1Vec, var2Vec, name1, name2, name11, name22, varargin)

% Default input parameters
options = struct('ExperimentName', [], ...
    'MarkerSize', 3, ... 
    'MarkerFaceAlpha', 1, ...
    'ColorScalars', [77 161 169]/255, ...
    'ColorWeight', 1, ...
    'LinearTrend', true);
% Read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(scatter_corr) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(scatter_corr) %s is not a recognized parameter name', pair{1})
    end
end

scatter(var1Vec, var2Vec, ...
    options.MarkerSize, ...
    options.ColorWeight * options.ColorScalars, ...    
    'filled', ...
    'MarkerFaceAlpha', options.MarkerFaceAlpha)
hold on;

% significant correlation check
[R, p] = corrcoef(var1Vec, var2Vec);
X = var1Vec; Y = var2Vec; % they must be columns
ttl_sfx = [];
if p(2,1) < .05 && options.LinearTrend
    ttl_sfx = [' r=' num2str(R(2,1),'%2.2f') ', p=' num2str(p(2,1), '%2.4f')];

    tbl = table(X, Y); mdl = fitlm(tbl,'Y~X');
    crr_fit = plotAdded(mdl,'X', 'Marker', 'none');
    crr_fit(2).Color = 'k';
    crr_fit(3).Color = 'k';
    legend('off')

    title({options.ExperimentName, [name1 ' vs ' name2], ...
        [' r=' sprintf('%.3f', R(2,1)) ', p=' sprintf('%.2e', p(2,1))]}, ...
        'Interpreter', 'none')
elseif p(2,1) < .05
    % disp({options.ExperimentName, [name1 ' vs ' name2], ...
    %     [' r=' sprintf('%.3f', R(2,1)) ', p=' sprintf('%.2e', p(2,1))]})
    title({options.ExperimentName, [name1 ' vs ' name2]}, ...
        'Interpreter', 'none')
else
    title({options.ExperimentName, [name1 ' vs ' name2]}, ...
        'Interpreter', 'none')
end
hold off
xlabel(name11)
ylabel(name22)

% To print statistic and p-value
% Extract the correlation coefficient
r = R(1,2);  % off-diagonal element
n = length(var1Vec);  % or y, assuming x and y are the same length

% Compute t-statistic
tstat = r * sqrt((n - 2) / (1 - r^2));

% Degrees of freedom
df = n - 2;

stats.tstat = tstat;
stats.df = df;
stats.p = 2 * (1 - tcdf(abs(tstat), df));

end


% ******************** vis & compare two group stats ********************

function [] = compare2(col1, col2, nSbj, varargin)

% Default input parameters
options = struct('figTitle', [], ...
    'figXlabel', [], ...
    'figXticklabels', []);
% Read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(scatter_corr) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(scatter_corr) %s is not a recognized parameter name', pair{1})
    end
end

[~, p] = ttest(col1, col2);
visMat = [col1, col2]; % (nSbj, 2)
xTicks = 1:2;
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

figure
bc = boxchart(reshape(repmat(xTicks, nSbj, 1), [], 1), ...
    reshape(visMat, [], 1), ...
    'GroupByColor', reshape(repmat(xTicks, nSbj, 1), [], 1), ...
    'BoxFaceAlpha', .4, ...
    'LineWidth', 4, ...
    'MarkerStyle', 'none');
hold on;
bc(1).BoxFaceColor = [1 1 1];
bc(1).BoxEdgeColor = [0 0 0];
bc(1).XData = 1.2 * ones(nSbj, 1);
bc(2).BoxFaceColor = [0 0 0];
bc(2).BoxEdgeColor = [0 0 0];
bc(2).XData = 1.8 * ones(nSbj, 1);

title(options.figTitle, 'Interpreter', 'none')
xlabel(options.figXlabel)
xticks(xTicks)
xticklabels(options.figXticklabels)

ax = gca;
ax.TickLabelInterpreter = 'none';
plot([xTicks(1) + .2  xTicks(2) - .2], ...
    visMat', ...
    'LineWidth', 2, ...
    'Marker', 'o', ...
    'MarkerSize', 7, ...
    'Color', [.85 .85 .85]);
line(xTicks, [max(visMat(:)) + .04,  max(visMat(:)) + .04], ...
    'Color', 'k')
text(1.5 - .02 * length(sig_text), max(visMat(:)) + .06, sig_text, ...
    'FontSize', 15)
ylim([min(visMat(:)) - .04, max(visMat(:)) + .08])
hold off

end
