
%% Comparison for validation: virtual simulation episode vs human behavior data

clear 
LoadDefaultSettings;

%% Data loading

% experiment info
sbjID = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24};
N_sbj = length(sbjID);
sbjNames = {'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
    'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};  % Lee2014
nSessions = 5*ones(1,length(sbjNames)); nSessions(11)=4; nSessions(16)=3; nSessions(20)=2; nSessions(21)=4;     % Lee2014

% load human data (session combined)
humanRoot = '\\143.248.30.94\bmlsamba\ydsung/A_Research';
humanDir = 'Dim_control_PFC_metaRL/M3_2014subs_ori';
humanName = 'SBJ_structure_each_exp_BESTcollectionJan29.mat';
load([humanRoot '/' humanDir '/' humanName], 'SBJ'); % 22 subjects data

% make the cell(sbj) of behavior matrices
humanBhvs = cell(length(SBJ), 1);
humanBlks = cell(length(SBJ), 1);
for i = 1:N_sbj   % subject loop
    humanBhvs{i} = cell2mat(SBJ{i}.HIST_behavior_info')'; % [18 x N_trial] double
    humanBlks{i} = cell2mat(SBJ{i}.HIST_block_condition); % [2 x N_trial] double
end
% col1 - block #, col2 - trial # (in each block), col3 - block condition - 1: G(with low uncertainty), 2: G'(with high uncertainty), 3:H(with high uncertainty), 4:H'(with low uncertainty), 
% col4 - S1, col5 - S2, col6 - S3, col7 - A1 (action in S1), col8 - A2 (action in S2), 
% col9 - RT(A1), col10 - RT(A2), 
% col11 - onset (S1) from the trial start, col12 - onset (S2) from the trial start, col13 - onset (S3) from the trial start, 
% col14 - onset (A1) from the trial start, col15 - onset (A2) from the trial start, 
% col16 - reward amount (0/10/20/40) at S3, col17 - total amount (total in the current session), 
% col18 - goal state=outcome state. 6: 40reward state, 7: 20reward state, 8: 10reward state, -1: universal state,  

% load virtual simulation data (S1, S2, A1, A2, Qcho1, Qcho2 Qunc1, Qunc2, Model1, Model2)
simulRoot = 'D:\[소스코드예제]\[소스코드] arbitration model\neural-dim_stab-flex_simul';
% simulDir = 'arb_virtual_episode100'; N_simul = 100;
% simulDir = 'arb_virtual_episode1000'; N_simul = 1000;
% simulDir = 'mb_virtual_episode100'; N_simul = 100;
simulDir = 'mb_virtual_episode1000'; N_simul = 1000;
% simulDir = 'mf_virtual_episode1000'; N_simul = 1000;
simulName = 'SBJ';
N_feat = 25; % # of features
epiCell = cell(N_sbj, 1);
disp('Virtual simulation result loading')
for i = 1:N_sbj

    % sanity check - block condition
    if any(humanBhvs{i}(3, :) ~= humanBlks{i}(2, :)); error('Inconsistent block condition in SBJ!'); end
    
    N_trial = size(humanBhvs{i}, 2); % *** N_trial from the human data, for validation ***
    epiTensor = nan(N_simul, N_feat, N_trial);
    load([simulRoot '/' simulDir '/' simulName '_' sbjNames{sbjID{i}}], 'SBJtot') % a single simulation from one subject experiment
    for j = 1:N_simul
        
        % sanity check - block condition (human vs virtual)
        tempBlk = cell2mat(SBJtot{j}{1}.HIST_block_condition); tempBlk = tempBlk(2, :);
        if any(tempBlk ~= humanBlks{i}(2, :)); error('Inconsistent block condition: human vs virtual'); end
        tempGT = (tempBlk > 2.5) + 1; % G(1): 1 & 2, H(2): 3 & 4
        tempUC = ismember(tempBlk, [2 3]) + 1; % 1:Low, 2:High

        % sanity check - goal condition (human vs virtual)        
        tempBhv = cell2mat(SBJtot{j}{1}.HIST_behavior_info')';
        tempGoal = tempBhv(18, :);
        if any(humanBhvs{i}(18, :) ~= tempGoal); error('Inconsistent goal condition: human vs virtual'); end
        
        % episode info.: states, actions, chosen value, unchosen value, active model (fwd vs sarsa)
        tempEpisode = cell2mat(SBJtot{j}{1}.model_BayesArb.episode(:)); % [N_session x 1] struct
        ss = [tempEpisode.HIST_ArbState]; ss1 = ss(1:2:end); ss2 = ss(2:2:end); % simulation states, [N_trial x 1] double
        sa = [tempEpisode.HIST_ArbAction]; sa1 = sa(1:2:end); sa2 = sa(2:2:end); % simulation actions
        sqc = [tempEpisode.HIST_ArbQcho]; sqc1 = sqc(1:2:end); sqc2 = sqc(2:2:end); % simulation Q(chosen)
        squ = [tempEpisode.HIST_ArbQunc]; squ1 = squ(1:2:end); squ2 = squ(2:2:end); % simulation Q(unchosen)
        sm = [tempEpisode.HIST_ArbMode]; sm1 = sm(1:2:end); sm2 = sm(2:2:end); % simulation active model
        
        % ******************** choice optimality ********************
        opt1 = nan(1, N_trial); optnan1 = nan(1, N_trial);
        opt2 = nan(1, N_trial); optnan2 = nan(1, N_trial);
        for k = 1:N_trial
            try
                opt1(k) = opt_normalization_1st2ndstage_exnan(ss1(k), sa1(k), tempGoal(k), tempBlk(k));
            catch
                opt1(k) = opt_normalization_1st2ndstage_exnan(ss1(k), sa1(k), tempGoal(k), tempBlk(k));
            end
            opt2(k) = opt_normalization_1st2ndstage_exnan(ss2(k), sa2(k), tempGoal(k), tempBlk(k));

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

        % ******************** choice versatility ********************
        gs1 = [1 diff(tempGoal)~=0]; % goal switch

        flex1 = 1 - [0 diff(sa1)==0]; % A1 choice switch
        flex1(:, ~ismember(gs1, 1)) = nan; % for A1 choice switch at goal switch

        % ===== stage 2 choice versatility =====
        flex2 = nan(1, length(sa2));
        gs2 = nan(1, length(sa2));
        S2_set = [2 3 4 5];
        for si = 1:length(S2_set)
            s2_idx = (ss2 == S2_set(si));
            if any(s2_idx)
                a2_s2 = sa2(s2_idx);
                flex2(s2_idx) = 1 - [0 diff(a2_s2)==0]; % A2 switch at S2
                g_s2 = tempGoal(s2_idx);
                gs2(s2_idx) = 1 - [0 diff(g_s2)==0]; % Goal switch at S2
            end
        end
        if any(isnan(flex2)); error('stage 2 versatility error'); end % sanity check
        flex2(:, ~ismember(gs2, 1)) = NaN; % for A2 switch at goal switch
        flex = mean([flex1; flex2], 1, 'omitnan');

        % ******************** choice consistency ********************
        stab1 = nan(1, N_trial);
        stab2 = nan(1, N_trial);
        Goal_set = [6 7 8 -1];
        for gi = 1:length(Goal_set)

            g_idx = (tempGoal == Goal_set(gi)); % goal index
            
            stab1(g_idx) = [0 diff(sa1(g_idx))==0]; % stage 1 choice consistency 

            % ===== stage 2 choice consistency =====
            S2_temp = ss2(g_idx);
            A2_temp = sa2(g_idx);
            ChoCon2_temp = NaN(1,length(A2_temp));
            S2_set = [2 3 4 5];
            for si = 1:length(S2_set)
                s2_idx = (S2_temp == S2_set(si));
                if any(s2_idx)
                    a2 = A2_temp(s2_idx);
                    ChoCon2_temp(s2_idx) = [0 diff(a2)==0];
                end
            end
            stab2(g_idx) = ChoCon2_temp;
        end
        if any(isnan(stab2)); error('stage 2 consistency error'); end % sanity check
        stab = mean([stab1; stab2], 1, 'omitnan');
        
        % storing the values
        epiTensor(j, :, :) = [ss1; ss2; sa1; sa2; ... % 1~4
            tempGT; tempUC; tempGoal; ... % 5~7
            optnan; flex; stab; ...  % 8~10
            optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
            opt; opt1; opt2; ... % 17~19
            sqc1; sqc2; squ1; squ2; ... % 20~23
            sm1; sm2]; % 24~25
        if N_trial ~= size(epiTensor(j, :, :), 3); error('Different N_trial for human vs simulation'); end
    end
    epiCell{i} = epiTensor;
    fprintf(' %d', i)
end
disp(' loading done')



%% Simple comparison with human data

close all

% N_trial: identical

% S2: different
figure
for i = 1:N_sbj

    dispMatS2 = [humanBhvs{i}(5, :); squeeze(epiCell{i}(:, 2, :))];
    if ~all(ismember([2 3 4 5], dispMatS2(:))); error('Different S2'); end

    % [N_trial x N_sbj]
    subplot(1, N_sbj, i)
    imagesc(dispMatS2)
    title([num2str(i) 'S2'])
    axis off

end

disp('First row: human data, below rows: virtual simulation')

% A1, A2: different
figure
for i = 1:N_sbj

    dispMatA1 = [humanBhvs{i}(7, :); squeeze(epiCell{i}(:, 3, :))];
    if ~all(ismember([1 2], dispMatA1(:))); error('Different A1'); end

    % [N_trial x N_sbj]
    subplot(1, N_sbj, i)
    imagesc(dispMatA1)
    title([num2str(i) 'A1'])
    axis off

end

figure
for i = 1:N_sbj

    dispMatA2 = [humanBhvs{i}(8, :); squeeze(epiCell{i}(:, 4, :))];
    if ~all(ismember([1 2], dispMatA2(:))); error('Different A2'); end

    % [N_trial x N_sbj]
    subplot(1, N_sbj, i)
    imagesc(dispMatA2)
    title([num2str(i) 'A2'])
    axis off

end

disp('done')

%% Conditional distribution of a categorical variable

% Goal: goal-dependency of S2 and A2 check - goal-directedness of virtual agents

close all

% epiCell{i} = epiTensor;
% epiTensor(j, :, :) = [ss1; ss2; sa1; sa2; ... % 1~4
%             tempGT; tempUC; tempGoal; ... % 5~7
%             optnan; flex; stab; ...  % 8~10
%             optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
%             opt; opt1; opt2; ... % 17~19
%             sqc1; sqc2; squ1; squ2; ... % 20~23
%             sm1; sm2]; % 24~25

% ========== Target variable's dependency on condition variable ===========

% condName = 'Goal'; targName = 'S2'; condIndex = 7; targIndex = 2;  condClass = [-1 6 7 8]; targClass = [2 3 4 5];
condName = 'Goal'; targName = 'A1'; condIndex = 7; targIndex = 3;  condClass = [-1 6 7 8]; targClass = [1 2];

N_cond = length(condClass);
N_targ = length(targClass);
arraySbjSimulCondTarg = nan(N_sbj, N_simul, N_cond, N_targ); % storing category frequency

for i = 1:N_sbj

    condMat = epiCell{i}(:, condIndex, :);
    targMat = epiCell{i}(:, targIndex, :); % [N_simul x 1 x N_trial]
    N_trial = size(condMat, 3);

    % sanity check
    if any(~ismember(targMat, targClass)); error('target error'); end
    if any(~ismember(condMat, condClass)); error('condition error'); end

    for ci = 1:N_cond % condition loop

        c_idx = (condMat == condClass(ci)); % [N_simul x 1 x N_trial]

        for ti = 1:N_targ % target loop

            arraySbjSimulCondTarg(i, :, ci, ti) = ...
                sum((targMat == targClass(ti)) & c_idx, 3) / N_trial;
            % [N_simul x 1 x N_trial] -> [N_simul x 1]

        end
    end
end

% Subject plot (simulations are averaged)
figure('Position', [0 1080/2-80 1920 1080/2]) % [left bottom width height] 
for ci = 1:N_cond
    subplot(1, N_cond, ci)
    boxchart(reshape(mean(arraySbjSimulCondTarg(:, :, ci, :), 2), [], N_targ), ...
        'JitterOutliers','on', ...
        'MarkerStyle','.');
    xticks(categorical(1:N_targ)); xticklabels(targClass)
    title({simulDir, [condName ' ' num2str(condClass(ci)) ' ' targName ' distribution']}, ...
        'Interpreter', 'none')
end

% Simulation plot (subjects are averaged)
figure('Position', [0 1080/2-80 1920 1080/2]) % [left bottom width height] 
for ci = 1:N_cond
    subplot(1, N_cond, ci)
    boxchart(reshape(mean(arraySbjSimulCondTarg(:, :, ci, :), 1), [], N_targ), ...
        'JitterOutliers','on', ...
        'MarkerStyle','.');
    xticks(categorical(1:N_targ)); xticklabels(targClass)
    title({simulDir, [condName ' ' num2str(condClass(ci)) ' ' targName ' distribution']}, ...
        'Interpreter', 'none')
end

% Subject x simulations plot
figure('Position', [0 1080/2-80 1920 1080/2]) % [left bottom width height] 
for ci = 1:N_cond
    subplot(1, N_cond, ci)
    boxchart(reshape(arraySbjSimulCondTarg(:, :, ci, :), [], N_targ), ...
        'JitterOutliers','on', ...
        'MarkerStyle','.');
    xticks(categorical(1:N_targ)); xticklabels(targClass)
    title({simulDir, [condName ' ' num2str(condClass(ci)) ' ' targName ' distribution']}, ...
        'Interpreter', 'none')
end



disp('done')

%% FlexStabOpt correlation check

close all

% epiCell{N_sbj}(N_simul, N_feature, N_trial)

% epiCell{i} = epiTensor;
% epiTensor(j, :, :) = [ss1; ss2; sa1; sa2; ... % 1~4
%             tempGT; tempUC; tempGoal; ... % 5~7
%             optnan; flex; stab; ...  % 8~10
%             optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
%             opt; opt1; opt2; ... % 17~19
%             sqc1; sqc2; squ1; squ2; ... % 20~23
%             sm1; sm2]; % 24~25

% restVarName = []; featR = []; featRange = [];
subjRestVar = 'Optimality'; subjRestFeat = 8; subjFeatRange = [.7614 .9932];

% name1 = 'Flexibility'; name2 = 'Optimality'; feat1 = 9; feat2 = 8;
% % name1 = 'Stability'; name2 = 'Optimality'; feat1 = 10; feat2 = 8;
% % name1 = 'Flexibility'; name2 = 'Stability'; feat1 = 9; feat2 = 10;

name1 = 'Flexibility'; name2 = 'Performance'; feat1 = 9; feat2 = 8; name11 = 'Flexibility (choice versatility)'; name22 = 'Performance (choice optimality)';
% name1 = 'Stability'; name2 = 'Performance'; feat1 = 10; feat2 = 8; name11 = 'Stability (choice consistency)'; name22 = 'Performance (choice optimality)';
% name1 = 'Flexibility'; name2 = 'Stability'; feat1 = 9; feat2 = 10; name11 = 'Flexibility (choice versatility)'; name22 = 'Stability (choice consistency)';

var1 = nan(N_sbj, N_simul);
var2 = nan(N_sbj, N_simul);
restVar = ones(N_sbj, N_simul);

for i = 1:N_sbj
    var1(i, :) = mean(epiCell{i}(:, feat1, :), 3, 'omitnan');
    var2(i, :) = mean(epiCell{i}(:, feat2, :), 3, 'omitnan');
    if ~isempty(subjRestVar)
        restVar(i, :) = (mean(epiCell{i}(:, subjRestFeat, :), 3, 'omitnan') > subjFeatRange(1)) ... 
            & (mean(epiCell{i}(:, subjRestFeat, :), 3, 'omitnan') < subjFeatRange(2));
    end
end

var1Vec = var1(logical(restVar(:)));
var2Vec = var2(logical(restVar(:)));
[simulColor, sbjColor] = meshgrid(1:N_simul, 1:N_sbj); % x and y grid
simulColor = simulColor(logical(restVar(:)));
sbjColor = sbjColor(logical(restVar(:)));

% % sanity check
% p1 = anovan(var1Vec, {simulColor, sbjColor})
% p2 = anovan(var2Vec, {simulColor, sbjColor})

scatter_corr(var1Vec, var2Vec, name1, name2, name11, name22, ...
    'ExperimentName', simulDir);
scatter_corr(var1Vec, var2Vec, name1, name2, name11, name22, ...
    'ExperimentName', simulDir, ...
    'ColorScalars', sbjColor);
scatter_corr(var1Vec, var2Vec, name1, name2, name11, name22, ...
    'ExperimentName', simulDir, ...
    'ColorScalars', simulColor, ...
    'MarkerSize', 3);

disp('done')

%% Overlay FlexStabOpt plot (human, MB, MF)

clear 
LoadDefaultSettings;
% LoadMyColor;

% ============================= Data loading ============================

% experiment info
sbjID = {2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23,24}; % full 22 subjects
validIdx = [1:18, 21:22];
N_sbj = length(validIdx);
sbjNames = {'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
    'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};  % Lee2014
nSessions = 5*ones(1,length(sbjNames)); nSessions(11)=4; nSessions(16)=3; nSessions(20)=2; nSessions(21)=4;     % Lee2014

% load human data (session combined)
humanRoot = '\\143.248.73.45\bmlsamba\ydsung/A_Research';
humanDir = 'Dim_control_PFC_metaRL/M3_2014subs_ori';
humanName = 'SBJ_structure_each_exp_BESTcollectionJan29.mat';
load([humanRoot '/' humanDir '/' humanName], 'SBJ'); % 22 subjects data

% make the cell(sbj) of behavior matrices
humanBhvs = cell(N_sbj, 1);
humanBlks = cell(N_sbj, 1);
humanTensor = cell(N_sbj, 1);
for i = 1:N_sbj   % subject loop
    humanBhvs{i} = cell2mat(SBJ{validIdx(i)}.HIST_behavior_info')'; % [18 x N_trial] double
    humanBlks{i} = cell2mat(SBJ{validIdx(i)}.HIST_block_condition); % [2 x N_trial] double
    
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

% load virtual simulation data (S1, S2, A1, A2, Qcho1, Qcho2 Qunc1, Qunc2, Model1, Model2)
% simulRoot = 'D:\[소스코드예제]\[소스코드] arbitration model\neural-dim_stab-flex_simul';
simulRoot = 'D:\neural-dim-github-v2\functions\analysis_simulation';
% simulDirs = {'mb_virtual_episode1000', 'mf_virtual_episode1000'}; N_simul = 1000;
simulDirs = {'mb_human_episode100'}; N_simul = 100;
simulName = 'SBJ';
N_feat = 29; % # of features
epiCells = cell(N_sbj, 2);
disp('Virtual simulation result loading')
for sim = 1:2
    for i = 1:N_sbj

        % sanity check - block condition
        if any(humanBhvs{i}(3, :) ~= humanBlks{i}(2, :)); error('Inconsistent block condition in SBJ!'); end

        N_trial = size(humanBhvs{i}, 2); % *** N_trial from the human data, for validation ***
        epiTensor = nan(N_simul, N_feat, N_trial);
        load([simulRoot '/' simulDirs{sim} '/' simulName '_' sbjNames{sbjID{validIdx(i)}}], 'SBJtot') % a single simulation from one subject experiment
        for j = 1:N_simul

            % sanity check - block condition (human vs virtual)
            tempBlk = cell2mat(SBJtot{j}{1}.HIST_block_condition); tempBlk = tempBlk(2, :);
            if any(tempBlk ~= humanBlks{i}(2, :)); error('Inconsistent block condition: human vs virtual'); end
            tempGT = (tempBlk > 2.5) + 1; % G(1): 1 & 2, H(2): 3 & 4
            tempUC = ismember(tempBlk, [2 3]) + 1; % 1:Low, 2:High

            % sanity check - goal condition (human vs virtual)
            tempBhv = cell2mat(SBJtot{j}{1}.HIST_behavior_info')';
            tempGoal = tempBhv(18, :);
            if any(humanBhvs{i}(18, :) ~= tempGoal); error('Inconsistent goal condition: human vs virtual'); end

            % episode info.: states, actions, chosen value, unchosen value, active model (fwd vs sarsa)
            tempEpisode = cell2mat(SBJtot{j}{1}.model_BayesArb.episode(:)); % [N_session x 1] struct
            ss = [tempEpisode.HIST_ArbState]; ss1 = ss(1:2:end); ss2 = ss(2:2:end); % simulation states, [N_trial x 1] double
            sa = [tempEpisode.HIST_ArbAction]; sa1 = sa(1:2:end); sa2 = sa(2:2:end); % simulation actions
            sqc = [tempEpisode.HIST_ArbQcho]; sqc1 = sqc(1:2:end); sqc2 = sqc(2:2:end); % simulation Q(chosen)
            squ = [tempEpisode.HIST_ArbQunc]; squ1 = squ(1:2:end); squ2 = squ(2:2:end); % simulation Q(unchosen)
            sm = [tempEpisode.HIST_ArbMode]; sm1 = sm(1:2:end); sm2 = sm(2:2:end); % simulation active model

            % ******************** choice optimality ********************
            [opt, opt1, opt2, optnan, optnan1, optnan2] = optimality_temp(ss1, sa1, ss2, sa2, tempGoal, tempBlk);

            % ******************** choice versatility ********************
            [flex, flex1, flex2] = versatility_temp(sa1, ss2, sa2, tempGoal);

            % ******************** choice consistency ********************
            [stab, stab1, stab2] = consistency_temp(sa1, ss2, sa2, tempGoal);

            % storing the values
            epiTensor(j, :, :) = [ss1; ss2; sa1; sa2; ... % 1~4
                tempGT; tempUC; tempGoal; ... % 5~7
                optnan; flex; stab; ...  % 8~10
                optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
                opt; opt1; opt2; ... % 17~19
                tempBlk; ... % 20
                sqc1; sqc2; squ1; squ2; ... % 21~24
                abs(sqc1 - squ1); abs(sqc2 - squ2); mean([abs(sqc1 - squ1); abs(sqc2 - squ2)], 1); ... % 25~27
                sm1; sm2]; % 28~29
            if N_trial ~= size(epiTensor(j, :, :), 3); error('Different N_trial for human vs simulation'); end
        end
        epiCells{i, sim} = epiTensor;
        fprintf(' %d', i)
    end
end
disp(' loading done')

%% ============================= Visualization =============================
close all

% humanTensor{i} = [hs1; hs2; ha1; ha2; ... % 1~4
%     hGT; hUC; hGoal; ... % 5~7
%     optnan; flex; stab; ...  % 8~10
%     optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
%     opt; opt1; opt2]; % 17~19

% epiCell{N_sbj}(N_simul, N_feature, N_trial)

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

% name1 = 'Flexibility'; name2 = 'Performance'; feat1 = [13 14]; feat2 = [11 12]; name11 = 'Flexibility (choice versatility)'; name22 = 'Performance (choice optimality)';
% name1 = 'Stability'; name2 = 'Performance'; feat1 = [15 16]; feat2 = [11 12]; name11 = 'Stability (choice consistency)'; name22 = 'Performance (choice optimality)';
% name1 = 'Flexibility'; name2 = 'Stability'; feat1 = 9; feat2 = 10; name11 = 'Flexibility (choice versatility)'; name22 = 'Stability (choice consistency)';
name1 = 'Flexibility'; name2 = 'Stability'; feat1 = [13 14]; feat2 = [15 16]; name11 = 'Flexibility (choice versatility)'; name22 = 'Stability (choice consistency)';

figure('Position', [0 400-80 400 400]) % [left bottom width height]

% scatter virtual simulation results
for sim = 1:2

    var1 = nan(N_sbj, N_simul);
    var2 = nan(N_sbj, N_simul);
    restVar = ones(N_sbj, N_simul);

    for i = 1:N_sbj
        
        var1_temp = epiCells{i, sim}(:, feat1, :); % [N_simul x N_feat x N_trial]
        var2_temp = epiCells{i, sim}(:, feat2, :); % [N_simul x N_feat x N_trial]
        if ~isempty(trialRestVar)
            trialRestIdx = ismember(epiCells{i, sim}(:, trialRestFeat, :), trialFeatSet); % [N_simul x 1 x N_trial]
            var1_temp(repmat(~trialRestIdx, 1, length(feat1), 1)) = nan;
            var2_temp(repmat(~trialRestIdx, 1, length(feat2), 1)) = nan;
        end

        var1(i, :) = mean(var1_temp, [2 3], 'omitnan');
        var2(i, :) = mean(var2_temp, [2 3], 'omitnan');
        if ~isempty(subjRestVar)
            restVar(i, :) = (mean(epiCells{i, sim}(:, subjRestFeat, :), 3, 'omitnan') > subjFeatRange(1)) ...
                & (mean(epiCells{i, sim}(:, subjRestFeat, :), 3, 'omitnan') < subjFeatRange(2));
        end
    end

    var1Vec = var1(logical(restVar(:)));
    var2Vec = var2(logical(restVar(:)));
    [simulColor, sbjColor] = meshgrid(1:N_simul, 1:N_sbj); % x and y grid
    simulColor = simulColor(logical(restVar(:)));
    sbjColor = sbjColor(logical(restVar(:)));

    % % sanity check
    % p1 = anovan(var1Vec, {simulColor, sbjColor})
    % p2 = anovan(var2Vec, {simulColor, sbjColor})
    
    scatter_corr(var1Vec, var2Vec, name1, name2, name11, name22, ...
        'MarkerSize', 5, ...
        'LinearTrend', false, ...
        'ColorWeight', 2.25 - .75 * sim, ...
        'MarkerFaceAlpha', .1);
    % scatter_corr(var1Vec, var2Vec, name1, name2, name11, name22, ...
    %     'MarkerSize', 5, ...
    %     'LinearTrend', false, ...
    %     'ColorScalars', colors(7 + sim, :), ...
    %     'MarkerFaceAlpha', .1);
    hold on
end


% scatter human results
hvar1 = nan(N_sbj, 1);
hvar2 = nan(N_sbj, 1);

for i = 1:N_sbj

    hvar1_temp = humanTensor{i}(feat1, :); % [N_feat x N_trial]
    hvar2_temp = humanTensor{i}(feat2, :); % [N_feat x N_trial]
    if ~isempty(trialRestVar)
        trialRestIdx = ismember(humanTensor{i}(trialRestFeat, :), trialFeatSet); % [1 x N_trial]
        hvar1_temp(repmat(~trialRestIdx, length(feat1), 1)) = nan;
        hvar2_temp(repmat(~trialRestIdx, length(feat2), 1)) = nan;
    end
    
    hvar1(i) = mean(hvar1_temp, [1 2], 'omitnan');
    hvar2(i) = mean(hvar2_temp, [1 2], 'omitnan');

end

scatter_corr(hvar1, hvar2, name1, name2, name11, name22, ...
        'ExperimentName', ['Human, MB, MF ' trialRestVar ' ' num2str(trialFeatSet)], ...
        'MarkerSize', 15, ...
        'LinearTrend', false, ...
        'ColorScalars', colors(1, :));
hold off;

disp('done')

%% Conditional distribution of a continuous variable

% Goal: context-dependency of choice optimality and consistency

close all

% epiCell{i} = epiTensor;
% epiTensor(j, :, :) = [ss1; ss2; sa1; sa2; ... % 1~4
%                 tempGT; tempUC; tempGoal; ... % 5~7
%                 optnan; flex; stab; ...  % 8~10
%                 optnan1; optnan2; flex1; flex2; stab1; stab2; ... % 11~16
%                 opt; opt1; opt2; ... % 17~19
%                 tempBlk; ... % 20
%                 sqc1; sqc2; squ1; squ2; ... % 21~24
%                 abs(sqc1 - squ1); abs(sqc2 - squ2); mean([abs(sqc1 - squ1); abs(sqc2 - squ2)], 1) ... % 25~27
%                 sm1; sm2]; % 28~29

% ========== Target variable's dependency on condition variable ===========

% condName = 'Uncertainty'; targName = 'Value difference'; restName = 'G'; condIndex = 6; targIndex = 27; restIndex = 5; condClass = [1 2]; restClass = 1;
condName = 'Uncertainty'; targName = 'Value difference'; restName = 'H'; condIndex = 6; targIndex = 27; restIndex = 5; condClass = [1 2]; restClass = 2;
% condName = 'Uncertainty'; targName = 'Optimality'; restName = 'G'; condIndex = 6; targIndex = 8; restIndex = 5; condClass = [1 2]; restClass = 1;
% condName = 'Uncertainty'; targName = 'Optimality'; restName = 'H'; condIndex = 6; targIndex = 8; restIndex = 5; condClass = [1 2]; restClass = 2;
% condName = 'Goal condition'; targName = 'Optimality'; restName = 'H'; condIndex = 6; targIndex = 8; restIndex = 5; condClass = [1 2]; restClass = 2;

% for sim = 1:2
for sim = 1

    N_cond = length(condClass);
    arraySbjSimulCond = nan(N_sbj, N_simul, N_cond); % storing conditional mean value
    
    for i = 1:N_sbj

        condMat = epiCells{i, sim}(:, condIndex, :);
        targMat = epiCells{i, sim}(:, targIndex, :);
        restMat = epiCells{i, sim}(:, restIndex, :); % [N_simul x 1 x N_trial]
        N_trial = size(condMat, 3);
        
        % sanity check
        if any(~ismember(condMat, condClass)); error('condition error'); end

        for ci = 1:N_cond % condition loop

            c_idx = (condMat == condClass(ci)); % [N_simul x 1 x N_trial]
            targMatTemp = targMat;

            targMatTemp(~ismember(restMat, restClass) | ~c_idx) = nan;
            arraySbjSimulCond(i, :, ci) = mean(targMatTemp, 3, 'omitnan');
            % [N_simul x 1 x N_trial] -> [N_simul x 1]
            
        end
    end

    % % Subject plot (simulations are averaged)
    % figure('Position', [0 1080/2-80 1920 1080/2]) % [left bottom width height]
    % for ci = 1:N_cond
    %     subplot(1, N_cond, ci)
    %     boxchart(reshape(mean(arraySbjSimulCond(:, :, ci), 2), [], N_targ), ...
    %         'JitterOutliers','on', ...
    %         'MarkerStyle','.');
    %     xticks(categorical(1:N_targ)); xticklabels(targClass)
    %     title({simulDir, [condName ' ' num2str(condClass(ci)) ' ' targName ' distribution']}, ...
    %         'Interpreter', 'none')
    % end
    %
    % % Simulation plot (subjects are averaged)
    % figure('Position', [0 1080/2-80 1920 1080/2]) % [left bottom width height]
    % for ci = 1:N_cond
    %     subplot(1, N_cond, ci)
    %     boxchart(reshape(mean(arraySbjSimulCond(:, :, ci), 1), [], N_targ), ...
    %         'JitterOutliers','on', ...
    %         'MarkerStyle','.');
    %     xticks(categorical(1:N_targ)); xticklabels(targClass)
    %     title({simulDir, [condName ' ' num2str(condClass(ci)) ' ' targName ' distribution']}, ...
    %         'Interpreter', 'none')
    % end

    % across-subject paired t-test (var2 mean values for var1 class 1 vs 2)
    arrayReshape = reshape(arraySbjSimulCond, [], N_cond);
    [~, p] = ttest(arrayReshape(:, 1), arrayReshape(:, 2));

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

    % Subject x simulations plot
    figure('Position', [0 400-80 400 400]) % [left bottom width height]
    % figure('Position', [0 1080/2-80 1920 1080/2]) % [left bottom width height]
    % for ci = 1:N_cond
        % subplot(1, N_cond, ci)
        boxchart(reshape(arrayReshape, [], N_cond), ...
            'BoxFaceAlpha', .4, ...
            'LineWidth', 4, ...
            'MarkerStyle', 'none');
        xticks(categorical(1:N_cond)); xticklabels(condClass)
        hold on
        plot([condClass(1) + .2 condClass(2) - .2], ...
            reshape(arrayReshape, [], N_cond)', ...
            'LineWidth', .1, ...
            'Marker', 'o', ...
            'MarkerSize', 2, ...
            'Color', [.85 .85 .85])
        title({simulDirs{sim}, [condName ' ' targName ' ' restName ' distribution']}, ...
            'Interpreter', 'none')
        line(condClass, [1 + .04,  1 + .04], ...
            'Color', 'k')
        text(1.5 - .055 * length(sig_text), 1 + .07, sig_text, ...
            'FontSize', 20)
        % xlim([0.7, 2.3])
        % xlim([.6 2.4])
        % xlim([.5 2.5])
        % xticks(1:2); xticklabels({'Low', 'High'})
        hold off
    % end
end

disp('done')



%% Helper functions

% ******************** choice optimality ********************
function [opt, opt1, opt2, optnan, optnan1, optnan2] = optimality_temp(s1, a1, s2, a2, goals, blocks)

N_trial = length(s1);

% sanity check
if length(a1) ~= N_trial; error('S1 length ~= A1 length'); end
if length(s2) ~= N_trial; error('S1 length ~= S2 length'); end
if length(a2) ~= N_trial; error('S1 length ~= A2 length'); end
if length(goals) ~= N_trial; error('S1 length ~= Goals length'); end
if length(blocks) ~= N_trial; error('S1 length ~= Blocks length'); end

opt1 = nan(1, N_trial); optnan1 = nan(1, N_trial);
opt2 = nan(1, N_trial); optnan2 = nan(1, N_trial);
for k = 1:N_trial

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

N_trial = length(a1);

% sanity check
if length(s2) ~= N_trial; error('A1 length ~= S2 length'); end
if length(a2) ~= N_trial; error('A1 length ~= A2 length'); end
if length(goals) ~= N_trial; error('A1 length ~= Goals length'); end

stab1 = nan(1, N_trial);
stab2 = nan(1, N_trial);
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
function [R, p, mdl, tbl] = scatter_corr(var1Vec, var2Vec, name1, name2, name11, name22, varargin)

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
    disp({options.ExperimentName, [name1 ' vs ' name2], ...
        [' r=' sprintf('%.3f', R(2,1)) ', p=' sprintf('%.2e', p(2,1))]})
    title({options.ExperimentName, [name1 ' vs ' name2]}, ...
        'Interpreter', 'none')
else
    title({options.ExperimentName, [name1 ' vs ' name2]}, ...
        'Interpreter', 'none')
end
hold off
xlabel(name11)
ylabel(name22)

end
