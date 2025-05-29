
%% (Across sbj) taskVar corr

LoadDefaultSettings;
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);
close all

disp('start')

% Setting =================================================================

% ===== Experiment =====
% Exp = 'Lee2014';  idsize = 22;  except_ids = [];
Exp = 'Lee2014';  idsize = 22;  except_ids = [19 20];
% Exp = 'Heo2018';    idsize = 28;
% Exp = 'Kim2019';    idsize = 21;

% ===== taskVar names =====
% var1_name = 'R'; var1_cond = 'G';
% var1_name = 'R'; var1_cond = 'H';
% var1_name = 'R'; var1_cond = 'Goal6';
% var1_name = 'R'; var1_cond = 'Goal7';
% var1_name = 'R'; var1_cond = 'Goal8';
% var1_name = 'relMB'; var1_cond = 'H';
% var1_name = 'relMB'; var1_cond = 'Goal6';
% var1_name = 'relMB'; var1_cond = 'Goal7';
% var1_name = 'relMB'; var1_cond = 'Goal8';
% var1_name = 'SPE'; var1_cond = 'H';
% var1_name = 'SPE'; var1_cond = 'Goal6';
% var1_name = 'SPE'; var1_cond = 'Goal7';
% var1_name = 'SPE'; var1_cond = 'Goal8';
% var1_name = 'PMB'; var1_cond = 'H';
% var1_name = 'PMB'; var1_cond = 'Goal6';
% var1_name = 'PMB'; var1_cond = 'Goal7';
% var1_name = 'PMB'; var1_cond = 'Goal8';
% var1_name = 'ChoOpt2 uEarly';
% var1_name = 'PMB_G';
% var1_name = 'ChoOptG';
% var1_name = 'ChoOpt_G1';
% var1_name = 'ChoOpt_G2';
% var1_name = 'ChoOpt2 G uEarly';
% var1_name = 'ChoOptBin1'; var1_cond = 'G';
% var1_name = 'ChoOptBin1'; var1_cond = 'Goal(6,7,8)Switch';
% var1_name = 'ChoOptBin1'; var1_cond = 'Goal(6,7,8)Stay';
% var1_name = 'ChoOptBin2'; var1_cond = 'Goal(6,7,8)Switch';
% var1_name = 'ChoOptBin2'; var1_cond = 'Goal(6,7,8)Stay';
% var1_name = 'ChoOptBin G uEarly';
% var1_name = 'ChoOptBin1'; var1_cond = 'Goal6';
% var1_name = 'ChoOptBin1'; var1_cond = 'Goal7';
% var1_name = 'ChoOptBin1'; var1_cond = 'Goal8';
% var1_name = 'ChoOptBin1'; var1_cond = 'H';
% var1_name = 'ChoOptBin2'; var1_cond = 'Goal6';
% var1_name = 'ChoOptBin2'; var1_cond = 'Goal7';
% var1_name = 'ChoOptBin2'; var1_cond = 'Goal8';
% var1_name = 'ChoOptBin2'; var1_cond = 'H';
% var1_name = 'ChoOptBin1 G uEarly';
% var1_name = 'ChoOptBin2 G uEarly';
% var1_name = 'ChoiceSwitch1(GoalSwitch)'; var1_cond = 'G';
% var1_name = 'ChoiceSwitch2(GoalSwitch)'; var1_cond = 'G';
% var1_name = 'ChoiceSwitch1n2(GoalSwitch)'; var1_cond = 'G';
% var1_name = 'ChoiceSwitch1n2(GoalStay)'; var1_cond = 'G';
% var1_name = 'ChoiceSwitch1n2'; var1_cond = 'G';
% var1_name = 'ChoiceSwitch1n2'; var1_cond = 'Goal(6,7,8)Switch';
% var1_name = 'ChoiceSwitch1n2'; var1_cond = 'Goal(6,7,8)Stay';
% var1_name = 'ChoConsist1n2'; var1_cond = 'Goal(6,7,8)Switch';
% var1_name = 'ChoConsist1n2'; var1_cond = 'Goal(6,7,8)Stay';
var1_name = 'ChoConsist1n2'; var1_cond = 'G';
% var1_name = 'ChoConsist1'; var1_cond = 'Goal(6,7,8)Switch';
% var1_name = 'ChoConsist2'; var1_cond = 'Goal(6,7,8)Switch';
% var1_name = 'ChoConsist1'; var1_cond = 'Goal(6,7,8)Stay';
% var1_name = 'ChoConsist2'; var1_cond = 'Goal(6,7,8)Stay';
% var1_name = 'ChoCon1 per goal G'; var1_cond = [];

% var2_name = 'ucdiff:A1'; var2_cond = 'Goal6';
% var2_name = 'ucdiff:A2 state2'; var2_cond = 'Goal6';
% var2_name = 'ucdiff:A2 state3'; var2_cond = 'Goal6';
% var2_name = 'ucdiff:A2 state4'; var2_cond = 'Goal6';
% var2_name = 'ucdiff:A2 state5'; var2_cond = 'Goal6';
% var2_name = 'ucdiff:A1'; var2_cond = 'Goal7';
% var2_name = 'ucdiff:A2 state2'; var2_cond = 'Goal7';
% var2_name = 'ucdiff:A2 state3'; var2_cond = 'Goal7';
% var2_name = 'ucdiff:A2 state4'; var2_cond = 'Goal7';
% var2_name = 'ucdiff:A2 state5'; var2_cond = 'Goal7';
% var2_name = 'ucdiff:A1'; var2_cond = 'Goal8';
% var2_name = 'ucdiff:A2 state2'; var2_cond = 'Goal8';
% var2_name = 'ucdiff:A2 state3'; var2_cond = 'Goal8';
% var2_name = 'ucdiff:A2 state4'; var2_cond = 'Goal8';
% var2_name = 'ucdiff:A2 state5'; var2_cond = 'Goal8';
% var2_name = 'ucdiff:A1'; var2_cond = 'H';
% var2_name = 'ucdiff:A2 state2'; var2_cond = 'H';
% var2_name = 'ucdiff:A2 state3'; var2_cond = 'H';
% var2_name = 'ucdiff:A2 state4'; var2_cond = 'H';
% var2_name = 'ucdiff:A2 state5'; var2_cond = 'H';

% var2_name = 'UncCond--A1(anovan)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A2 state2(anovan)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A2 state3(anovan)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A2 state4(anovan)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A2 state5(anovan)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A1(anovan)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A2 state2(anovan)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A2 state3(anovan)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A2 state4(anovan)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A2 state5(anovan)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A1(anovan)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A2 state2(anovan)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A2 state3(anovan)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A2 state4(anovan)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A2 state5(anovan)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A1(anovan)'; var2_cond = 'H';
% var2_name = 'UncCond--A2 state2(anovan)'; var2_cond = 'H';
% var2_name = 'UncCond--A2 state3(anovan)'; var2_cond = 'H';
% var2_name = 'UncCond--A2 state4(anovan)'; var2_cond = 'H';
% var2_name = 'UncCond--A2 state5(anovan)'; var2_cond = 'H';

% var2_name = 'UncCond--A1(fitlm)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A2 state2(fitlm)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A2 state3(fitlm)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A2 state4(fitlm)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A2 state5(fitlm)'; var2_cond = 'Goal6';
% var2_name = 'UncCond--A1(fitlm)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A2 state2(fitlm)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A2 state3(fitlm)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A2 state4(fitlm)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A2 state5(fitlm)'; var2_cond = 'Goal7';
% var2_name = 'UncCond--A1(fitlm)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A2 state2(fitlm)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A2 state3(fitlm)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A2 state4(fitlm)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A2 state5(fitlm)'; var2_cond = 'Goal8';
% var2_name = 'UncCond--A1(fitlm)'; var2_cond = 'H';
% var2_name = 'UncCond--A2 state2(fitlm)'; var2_cond = 'H';
% var2_name = 'UncCond--A2 state3(fitlm)'; var2_cond = 'H';
% var2_name = 'UncCond--A2 state4(fitlm)'; var2_cond = 'H';
% var2_name = 'UncCond--A2 state5(fitlm)'; var2_cond = 'H';

% var2_name = 'ChoOptBin1n2'; var2_cond = 'Goal(6,7,8)Switch';
% var2_name = 'ChoOptBin1n2'; var2_cond = 'Goal(6,7,8)Stay';
% var2_name = 'ChoOptBin1n2'; var2_cond = 'G';
% var2_name = 'ChoOptBin1n2'; var2_cond = 'H';
% var2_name = 'ChoOptBin1';
% var2_name = 'ChoOptBin2';
% var2_name = 'ChoOptBin2'; var2_cond = 'Goal(6,7,8)Switch';
% var2_name = 'ChoCon2 uLate';
% var2_name = 'ChoOpt_G1';
% var2_name = 'ChoOpt_G2';
% var2_name = 'ChoiceSwitch1'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch2'; var2_cond = 'G';
%  var2_name = 'ChoiceSwitch1n2'; var2_cond = 'G';
 var2_name = 'ChoSwitch1n2(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch1n2'; var2_cond = 'Goal(6,7,8)Switch';
% var2_name = 'ChoiceSwitch1n2'; var2_cond = 'Goal(6,7,8)Stay';
% var2_name = 'ChoiceSwitch1(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoSwitch1(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch1(GoalStay)'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch2(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoSwitch2(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch2'; var2_cond = 'Goal(6,7,8)Stay';
% var2_name = 'ChoCon1 per goal G'; var2_cond = 'G';
% var2_name = 'ChoCon1 per goal G'; var2_cond = [];
% var2_name = 'ChoCon2 per goal G'; var2_cond = 'G';
% var2_name = 'ChoCon2 per goal G'; var2_cond = [];
% var2_name = 'ChoCon per goal G uLate';
% var2_name = 'ChoCon1 per goal G uLate';
% var2_name = 'ChoCon2 per goal G uLate';
% var2_name = 'ChoCon1 G';
% var2_name = 'ChoCon2 G';
% var2_name = 'ChoOptVarG';
% var2_name = 'ChoOptVarG1';
% var2_name = 'ChoOptVarG2';

% ===== visualization parameters =====

% Measure loading ===========================================================
var1_mean_vals = nan(1,idsize);
var2_mean_vals = nan(1,idsize);
var1_sig_vals = nan(1,idsize);
var2_sig_vals = nan(1,idsize);
for id = 1:idsize
    if ismember(id, except_ids); continue; end
    fprintf('%d ', id)
    
    % ===== Data load =====
    switch var1_name
        case 'ChoOpt2 uEarly'
            var1 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var1_cond); var1 = var1(2,:);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var1_cond); % 0:early, 1:late
            var1(ucEL~=0) = nan;
        case 'ChoOptBin1'
            var1 = arbMBMF_load_var(Exp, 'ChoOptBin1n2', id, var1_cond); var1 = var1(1,:);
        case 'ChoOptBin2'
            var1 = arbMBMF_load_var(Exp, 'ChoOptBin1n2', id, var1_cond); var1 = var1(2,:);
        case 'ChoOpt_G1'
            var1 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var1_cond); var1 = var1(1,:);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var1_cond);
            var1(gt~=1) = nan;
        case 'ChoOpt_G2'
            var1 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var1_cond); var1 = var1(2,:);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var1_cond);
            var1(gt~=1) = nan;
        case 'ChoOpt2 G uEarly'
            var1 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var1_cond); var1 = var1(2,:);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var1_cond);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var1_cond); % 0:early, 1:late
            var1(ucEL~=0 | gt~=1) = nan;
        case 'ChoOptBin G uEarly'
            var1 = arbMBMF_load_var(Exp, 'ChoOptBin1n2', id, var1_cond); var1 = mean(var1, 1);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var1_cond);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var1_cond); % 0:early, 1:late
            var1(ucEL~=0 | gt~=1) = nan;
        case 'ChoOptBin1 G uEarly'
            var1 = arbMBMF_load_var(Exp, 'ChoOptBin1n2', id, var1_cond); var1 = var1(1,:);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var1_cond);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var1_cond); % 0:early, 1:late
            var1(ucEL~=0 | gt~=1) = nan;
        case 'ChoOptBin2 G uEarly'
            var1 = arbMBMF_load_var(Exp, 'ChoOptBin1n2', id, var1_cond); var1 = var1(2,:);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var1_cond);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var1_cond); % 0:early, 1:late
            var1(ucEL~=0 | gt~=1) = nan;
        otherwise
            if contains(var1_name, '--')
                var1 = arbMBMF_glm_coeff([], [], ... 
                    'Experiment', Exp, ...
                    'SubjectID', id, ...
                    'Keyword', var1_name, ... 
                    'RestrictContext', var1_cond, ...
                    'Verbose', 1 ...
                    );
            else
                var1 = arbMBMF_load_var(Exp, var1_name, id, var1_cond);
            end
    end
        
    switch var2_name
        case 'ChoCon1'
            var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
            var2 = var2(1,:);
        case 'ChoCon2'
            var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
            var2 = var2(2,:);
        case 'ChoCon2 uLate'
            var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var2_cond); % 0:early, 1:late
            var2 = var2(2,:); var2(ucEL~=1) = nan;
        case 'ChoCon1 G'
            var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2 = var2(1,:);
            var2(GT~=1) = nan;
        case 'ChoCon2 G'
            var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2 = var2(2,:);
            var2(GT~=1) = nan;
        case 'ChoCon2 G uLate'
            var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var2_cond); % 0:early, 1:late
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2 = var2(2,:); var2(ucEL~=1 | GT~=1) = nan;
        case 'ChoCon1 per goal'
            var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
            var2 = var2(1,:);
        case 'ChoCon2 per goal'
            var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
            var2 = var2(2,:);
        case 'ChoCon1 per goal G'
            var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2 = var2(1,:);
            var2(GT~=1) = nan;
        case 'ChoCon2 per goal G'
            var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2 = var2(2,:);
            var2(GT~=1) = nan;
        case 'ChoCon per goal G uLate'
            var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var2_cond); % 0:early, 1:late
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2 = mean(var2, 1); var2(ucEL~=1 | GT~=1) = nan;
        case 'ChoCon1 per goal G uLate'
            var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var2_cond); % 0:early, 1:late
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2 = var2(1,:); var2(ucEL~=1 | GT~=1) = nan;
        case 'ChoCon2 per goal G uLate'
            var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
            ucEL = arbMBMF_load_var(Exp, 'UncPhaseEL', id, var2_cond); % 0:early, 1:late
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2 = var2(2,:); var2(ucEL~=1 | GT~=1) = nan;
        case 'StrictChoCon G'
            var2 = arbMBMF_load_var(Exp, 'StrictChoCon', id, var2_cond);
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2(GT~=1) = nan;
        case 'ChoOpt1'
            var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(1,:);
        case 'ChoOptBin1'
            var2 = arbMBMF_load_var(Exp, 'ChoOptBin1n2', id, var2_cond); var2 = var2(1,:);
        case 'ChoOpt2'
            var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(2,:);
        case 'ChoOptBin2'
            var2 = arbMBMF_load_var(Exp, 'ChoOptBin1n2', id, var2_cond); var2 = var2(2,:);
        case 'ChoOpt_G1'
            var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(1,:);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2(gt~=1) = nan;
        case 'ChoOpt_G2'
            var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(2,:);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2(gt~=1) = nan;
        case 'ChoOptVarG'
            opt = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); 
            opt = mean(opt, 1);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            opt(gt~=1) = nan;
            var2 = nanvar(opt);
        case 'ChoOptVarG1'
            opt = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); 
            opt = opt(1, :);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            opt(gt~=1) = nan;
            var2 = nanvar(opt);
        case 'ChoOptVarG2'
            opt = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); 
            opt = opt(2, :);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            opt(gt~=1) = nan;
            var2 = nanvar(opt);
        otherwise
            if contains(var2_name, '--')
                [var2, p_var2] = arbMBMF_glm_coeff([], [], ... 
                    'Experiment', Exp, ...
                    'SubjectID', id, ...
                    'Keyword', var2_name, ... 
                    'RestrictContext', var2_cond, ...
                    'Verbose', 1 ...
                    );
            elseif contains(var2_name, 'ucdiff')
                var2 = ucdiff(var2_name, ...
                    'Experiment', Exp, ...
                    'SubjectID', id, ...
                    'RestrictContext', var2_cond);
                p_var2 = [];
            else
                var2 = arbMBMF_load_var(Exp, var2_name, id, var2_cond);
                p_var2 = [];
            end
    end
    
    if strcmp(var2_name, 'ChoOpt1n2') || strcmp(var2_name, 'ChoConsist1n2') || size(var2, 1) == 2
        if size(var1, 1) == 2
            var1 = reshape(var1, 1, []);
            var2 = reshape(var2, 1, []);
        else
            var2 = reshape(var2, 1, []);
            var1 = reshape([var1; var1], 1, []);
        end
    end
    
    var1_mean_vals(id) = nanmean(var1);
    var2_mean_vals(id) = nanmean(var2);
    
    if p_var2 < .05
        var1_sig_vals(id) = nanmean(var1);
        var2_sig_vals(id) = nanmean(var2);
    elseif contains(var2_name, 'ucdiff') && (var2_mean_vals(id) == 0)
        var1_sig_vals(id) = nanmean(var1);
        var2_sig_vals(id) = nanmean(var2);
        var2_mean_vals(id) = nan;
    elseif isnan(p_var2)
        var2_mean_vals(id) = nan;
    end
            
end % id

% remove nan
if ~isempty(except_ids)    
    var1_mean_vals(except_ids) = [];
    var2_mean_vals(except_ids) = [];
end

% correlation
[R, p] = corrcoef([var1_mean_vals', var2_mean_vals'], 'Rows', 'pairwise');

% Visualization ===========================================================
figure
scatter(var1_mean_vals, var2_mean_vals)
if any(~isnan(var2_sig_vals))
    hold on
    scatter(var1_sig_vals, var2_sig_vals, 'filled') % visualize sig. sbj
    hold off
end

% Counting missing subjects
if any(isnan(var1_mean_vals)) || any(isnan(var2_mean_vals))
    
    if sum(isnan(var1_mean_vals)) > any(isnan(var2_mean_vals))
        NanIDs = isnan(var1_mean_vals);
    else
        NanIDs = isnan(var2_mean_vals);
    end
    ttl_sfx = ['Missing sbj: ' sprintf('%d ', find(NanIDs))];
    
elseif contains(var2_name, 'ucdiff') && any(var2_mean_vals == 0)
    
    % Counting no-UC-diff subjects
    ZeroIDs = var2_mean_vals==0;
    ttl_sfx = ['Zero-diff sbj: ' sprintf('%d ', find(ZeroIDs))];
    
else    
    ttl_sfx = [];
end

% Showing linear trend (for significant corr. only)
if p(1,2) < 0.05
    disp('test')
    hold on
    X = var1_mean_vals; Y = var2_mean_vals;
    tbl = table(X, Y);
    mdl = fitlm(tbl, 'Y~X');
    plotAdded(mdl); legend('off')
    hold off
else
    hold off
end
xlabel([var1_name var1_cond])
ylabel([var2_name var2_cond])
title({[var1_name var1_cond ' vs ' var2_name var2_cond], ... 
    ['r=' num2str(R(1,2)) ', p=' num2str(p(1,2))], ...
    ttl_sfx})

disp('end')

%% taskVar profile (continuous)
% plot var2 time series
% observe var2 class distribution for each var1 condition
% subplot per var1 class
% legend per var2 class

% within-subject measures (mean var2 per var1 condition)
% across-subject statistical test (mean var2s ANOVA/paired t-test)

LoadDefaultSettings;
close all
set(0,'DefaultAxesFontSize',12,'DefaultTextFontSize',12);

disp('start')

% Setting =================================================================

NameValue1 = {};
NameValue2 = {};

% ===== Experiment =====
% Exp = 'Lee2014';  idsize = 22;  except_ids = [];
Exp = 'Lee2014';  idsize = 22;  except_ids = [19 20];
% Exp = 'arb_virtual_episode100'; idsize = 22; except_ids = [19 20]; NameValue1 = {'SimulationNumber', 1}; NameValue2 = {'SimulationNumber', 1};
% Exp = 'arb_virtual_episode1'; idsize = 22; except_ids = [19 20]; NameValue1 = {'SimulationNumber', 1}; NameValue2 = {'SimulationNumber', 1};
% Exp = 'arb_virtual_episode1000'; idsize = 22; except_ids = [19 20]; NameValue1 = {'SimulationNumber', 100}; NameValue2 = {'SimulationNumber', 100};
% Exp = 'Heo2018';    idsize = 28;
% Exp = 'Kim2019';    idsize = 21;  except_ids = [];
% Exp = 'Kim2019 24';    idsize = 24;  except_ids = [];

% ===== condition variable (mainly categorical) =====
categorize1 = [];
% categorize1 = 2;
% categorize1 = 3;

var1_unq_pre = []; % default setting - some var1s are special cases

% var1_name = 'Half'; N_var1_global = 2; var1_cond = [];
% var1_name = 'EarlyLate'; N_var1_global = 2; var1_cond = [];
% var1_name = 'Phase3'; N_var1_global = 3; var1_cond = [];
% var1_name = 'EL x GS'; N_var1_global = 4; var1_cond = [];
% var1_name = 'EL x GT x UC'; N_var1_global = 8; var1_cond = [];
% var1_name = 'GT x UC'; N_var1_global = 4; var1_cond = [];

% var1_name = 'Session'; N_var1_global = 5; var1_cond = [];
% var1_name = 'Session'; N_var1_global = 6; 
% var1_name = 'blkCond'; N_var1_global = 4; var1_cond = [];
% var1_name = 'blkCond'; N_var1_global = 4; var1_cond = 'Goal(6,7,8)Switch';
% var1_name = 'blkCond'; N_var1_global = 4; var1_cond = 'Goal(6,7,8)Stay';
% var1_name = 'GoalCond'; N_var1_global = 2; var1_cond = []; 
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = [];
var1_name = 'UncCond'; N_var1_global = 2; var1_cond = 'G';
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = 'H';
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = 'Goal6';
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = 'Goal7';
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = 'Goal8';
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = 'H + Red';
% var1_name = 'UncCond'; N_var1_global = 2; var1_cond = 'uLate';
% var1_name = 'UncCond G'; N_var1_global = 2; var1_cond = [];
% var1_name = 'UncPhaseEL'; N_var1_global = 2; var1_cond = [];
% var1_name = 'UncPhaseEL2'; N_var1_global = 2; var1_cond = [];
% var1_name = 'UncPhaseEL2'; N_var1_global = 2; var1_cond = 'uL';
% var1_name = 'UncPhaseEL2'; N_var1_global = 2; var1_cond = 'uH';
% var1_name = 'UncPhase2'; N_var1_global = 2; var1_cond = [];
% var1_name = 'UncPhase'; N_var1_global = 3; var1_cond = [];
% var1_name = 'Session x GT'; N_var1_global = 10; var1_cond = [];
% var1_name = 'Session x UC'; N_var1_global = 10; var1_cond = [];
% var1_name = 'Session x UC'; N_var1_global = 10; var1_cond = 'G'; var1_unq_pre = 1:10;
% var1_name = 'CmplxCond'; N_var1_global = 2; 

% var1_name = 'GoaltypePhase'; N_var1_global = 3; 
% var1_name = 'GoaltypeSwitch'; N_var1_global = 2; 
% var1_name = 'UncPhase'; N_var1_global = 3; 
% var1_name = 'UncPhase G'; N_var1_global = 3; 
% var1_name = 'UncPhase2'; N_var1_global = 2; 
% var1_name = 'UncPhaseEL'; N_var1_global = 2; 
% var1_name = 'UncPhaseEL G'; N_var1_global = 2; 
% var1_name = 'uLuH EL'; N_var1_global = 4; var1_cond = [];
% var1_name = 'uLuH EL'; N_var1_global = 4; var1_cond = 'G';
% var1_name = 'uLuH EL'; N_var1_global = 4; var1_cond = 'H';
% var1_name = 'uLuH ucEL'; N_var1_global = 4; var1_cond = [];
% var1_name = 'uLuH ucEL'; N_var1_global = 4; var1_cond = 'G';
% var1_name = 'uLuH ucEL'; N_var1_global = 4; var1_cond = 'H';
% var1_name = 'uLL vs uHE'; N_var1_global = 2; 
% var1_name = 'uHL vs uLE'; N_var1_global = 2; 
% var1_name = 'UncSwitch'; N_var1_global = 2; 
% var1_name = 'GoalSwitch'; N_var1_global = 2; var1_cond = [];
% var1_name = 'GoalSwitch'; N_var1_global = 2; var1_cond = 'G';
% var1_name = 'Session x GS'; N_var1_global = 10; var1_cond = [];
% var1_name = 'Session x GS'; N_var1_global = 10; var1_cond = 'G';
% var1_name = 'Goal(6,7,8)Switch'; N_var1_global = 2; var1_cond = 'G';
% var1_name = 'Goal(6,7,8)Switch'; N_var1_global = 2; var1_cond = 'S1';
% var1_name = 'Goal(6,7,8)Switch'; N_var1_global = 2; var1_cond = 'S2';
% var1_name = 'Goal(6,7,8)Switch'; N_var1_global = 2; var1_cond = 'S3';
% var1_name = 'Goal(6,7,8)Switch'; N_var1_global = 2; var1_cond = 'S4';
% var1_name = 'Goal(6,7,8)Switch'; N_var1_global = 2; var1_cond = 'S5';
% var1_name = 'Goal(6,7,8)Switch'; N_var1_global = 2; var1_cond = 'uL';
% var1_name = 'Goal(6,7,8)Switch'; N_var1_global = 2; var1_cond = 'uH';
% var1_name = 'uLuH Goal(6,7,8)Switch'; N_var1_global = 4;  var1_cond = [];
% var1_name = 'uLuH Goal(6,7,8)Switch'; N_var1_global = 4;  var1_cond = 'G';

% var1_name = 'GoalPhase'; N_var1_global = 3;
% var1_name = 'Goal'; N_var1_global = 4; var1_cond = []; % Lee 2014
% var1_name = 'Goal'; N_var1_global = 3; % Kim 2019
% var1_name = 'Goal(6,7,8)'; N_var1_global = 3; var1_cond = [];
% var1_name = 'Goal(6,7,8)'; N_var1_global = 3; var1_cond = []; NameValue1 = {'SimulationNumber', 1};
% var1_name = 'Session x Goal(6,7,8)'; N_var1_global = 15; var1_cond = [];
% var1_name = 'GoalxUC'; N_var1_global = 8; var1_cond = []; 
% var1_name = 'Goal(6,7,8)xUC'; N_var1_global = 6; var1_cond = []; 
% var1_name = 'Goal(6,7,8)xUCSwitch'; N_var1_global = 2;
% var1_name = 'AchGoalS1RandA'; N_var1_global = 3;
% var1_name = 'AchGoalS1OptA'; N_var1_global = 3;
% var1_name = 'GoalValEntropy'; N_var1_global = 2;
% var1_name = 'S3';
% N_var1_global = 4;
% var1_name = 'R';
% N_var1_global = 4;

% var1_name = 'SPE'; N_var1_global = 3; var1_cond = []; 
% var1_name = 'RPE'; N_var1_global = 3; var1_cond = []; 
% var1_name = 'RPE3';
% N_var1_global = 4;

% var1_name = 'relMB_G'; N_var1_global = 2;  % caterorized
% var1_name = 'relMB'; N_var1_global = 3;
% var1_name = 'relMB2';
% N_var1_global = 2;  % caterorized
% var1_name = 'relMB3';
% N_var1_global = 4;
% var1_name = 'relMF'; N_var1_global = 3;
% var1_name = 'relMF2';
% N_var1_global = 4;
% var1_name = 'relMF3';
% N_var1_global = 4;
% var1_name = 'relMAX_G'; N_var1_gloabl = 2;
% var1_name = 'relMAX'; N_var1_global = 3;
% var1_name = 'relMAX2';
% N_var1_global = 4;
% var1_name = 'relMAX3';
% N_var1_global = 4;
% var1_name = 'relComp'; N_var1_global = 3;
% var1_name = 'relDiff'; N_var1_global = 2;
% var1_name = 'relDiff_G'; N_var1_global = 2;

% var1_name = 'PMB'; N_var1_global = 2; var1_cond = [];
% var1_name = 'PMB_G'; N_var1_global = 2;

% var1_name = 'ChoOpt'; N_var1_global = 2;
% var1_name = 'binChoOpt1n2'; N_var1_global = 2;
% var1_name = 'ChoOptBinNaN1n2'; N_var1_global = 2; var1_cond = 'G'; NameValue1 = {'MovingMean', [0 4], 'Mat2Row', 1};
% var1_name = 'ChoOptBinNaN1n2'; N_var1_global = 2; var1_cond = 'G'; NameValue1 = {'MovingMean', [4 0], 'Mat2Row', 1};
% var1_name = 'relComp2';
% N_var1_global = 2;

% var1_name = 'ChoConsist1n2'; N_var1_global = 2; var1_cond = 'G'; NameValue1 = {'MovingMean', [0 4], 'Mat2Row', 1};
% var1_name = 'ChoConsist1n2'; N_var1_global = 2; var1_cond = 'G'; NameValue1 = {'MovingMean', [4 0], 'Mat2Row', 1};

% var1_name = 'ChoConsist1'; N_var1_global = 2; var1_cond = 'H';
% var1_name = 'ChoConsist1'; N_var1_global = 2; var1_cond = 'Goal6';
% var1_name = 'ChoConsist1'; N_var1_global = 2; var1_cond = 'Goal7';
% var1_name = 'ChoConsist1'; N_var1_global = 2; var1_cond = 'Goal8';

% var1_name = 'ChoConsist2'; N_var1_global = 2; var1_cond = 'H state4';
% var1_name = 'ChoConsist2'; N_var1_global = 2; var1_cond = 'Goal6 state4';
% var1_name = 'ChoConsist2'; N_var1_global = 2; var1_cond = 'Goal7 state4';
% var1_name = 'ChoConsist2'; N_var1_global = 2; var1_cond = 'Goal8 state4';

% var1_name = 'ChoSwitch1n2(GoalSwitch)'; N_var1_global = 2; var1_cond = 'G'; NameValue1 = {'MovingMean', [0 4], 'Mat2Row', 1};
% var1_name = 'ChoSwitch1n2(GoalSwitch)'; N_var1_global = 2; var1_cond = 'G'; NameValue1 = {'MovingMean', [4 0], 'Mat2Row', 1};

% ===== target variable (time series) =====
categorize2 = [];
% categorize2 = 2;
% categorize2 = 3;

% var2_name = 'Session'; var2_cond = [];
% var2_name = 'Session x GT'; var2_cond = [];
% var2_name = 'Session x GS'; var2_cond = [];
% var2_name = 'Session x GS'; var2_cond = 'G';
% var2_name = 'EL x GS'; var2_cond = [];
% var2_name = 'GT x UC'; var2_cond = [];
% var2_name = 'EL x GT x UC'; var2_cond = [];

% var2_name = 'UncCond'; var2_cond = [];

% var2_name = 'blkSwitch';
% var2_name = 'blkPhase';

% var2_name = 'GoaltypeSwitch';
% var2_name = 'GoaltypePhase';

% var2_name = 'UncSwitch';
% var2_name = 'UncPhase';

% var2_name = 'CmplxSwitch';
% var2_name = 'CmplxPhase';

% var2_name = 'Goal'; var2_cond = [];
% var2_name = 'Goal(6,7,8)'; var2_cond = 'Goal(6,7,8)Switch';
% var2_name = 'Goal(6,7,8)'; var2_cond = 'Goal(6,7,8)Stay';
% var2_name = 'Session x Goal(6,7,8)'; var2_cond = [];
% var2_name = 'Goal(6,7,8)xUC'; var2_cond = [];
% var2_name = 'Goal(6,7,8)xUCd1'; var2_cond = [];
% var2_name = 'Goal(6,7,8)xUCd2'; var2_cond = [];
% var2_name = 'Goal(6,7,8)xUC'; var2_cond = 'uL';
% var2_name = 'Goal(6,7,8)xUC'; var2_cond = 'uH';
% var2_name = 'GoalSwitch'; var2_cond = 'G';
% var2_name = 'GoalPhase';
% var2_name = 'GoalValEntropy';

% var2_name = 'SPE'; var2_cond = [];
% var2_name = 'SPE'; var2_cond = 'G';
% var2_name = 'SPE'; var2_cond = 'H';
% var2_name = 'SPE2'; var2_cond = [];
% var2_name = 'SPE2'; var2_cond = 'G';
% var2_name = 'SPE3'; var2_cond = [];
% var2_name = 'SPE3'; var2_cond = 'G';

% var2_name = 'RPE';
% var2_name = 'absRPE'; var2_cond = [];
% var2_name = 'absRPE'; var2_cond = 'G';
% var2_name = 'absRPE'; var2_cond = 'H';

% var2_name = 'RPE2'; var2_cond = [];
% var2_name = 'RPE2'; var2_cond = 'G';
% var2_name = 'RPE3'; var2_cond = [];
% var2_name = 'RPE3'; var2_cond = 'G';
% var2_name = 'absRPE3';

% var2_name = 'relMB';
% var2_name = 'relMB_G';
% var2_name = 'relMF';
% var2_name = 'relMAX';
% var2_name = 'relMAX_G';
% var2_name = 'relDiff'; var2_cond = [];
% var2_name = 'relDiff2'; var2_cond = [];
% var2_name = 'relDiff2'; var2_cond = 'G';
% var2_name = 'relDiff3'; var2_cond = [];
% var2_name = 'relDiff3'; var2_cond = 'G';
% var2_name = 'relDiff_G';
% var2_name = 'relComp2'; var2_cond = [];
% var2_name = 'relComp2'; var2_cond = 'G';
% var2_name = 'relComp3'; var2_cond = [];
% var2_name = 'relComp3'; var2_cond = 'G';

% var2_name = 'relMB2';
% var2_name = 'relMB3';
% var2_name = 'relMF2';
% var2_name = 'relMF3';
% var2_name = 'relMAX2';
% var2_name = 'relMAX3';
% var2_name = 'relDiff2'; var2_cond = [];
% var2_name = 'relDiff3'; var2_cond = [];
% var2_name = 'relComp2';
% var2_name = 'relComp3';

% var2_name = 'PMB2'; var2_cond = [];   
% var2_name = 'PMB3'; var2_cond = [];   
% var2_name = 'PMB'; var2_cond = [];   
% var2_name = 'PMB'; var2_cond = 'G';  
% var2_name = 'PMB'; var2_cond = 'H';   
% var2_name = 'PMB_G';

% var2_name = 'ActiveModel'; var2_cond = [];   

% var2_name = 'Qarb'; var2_cond = [];
% var2_name = 'Qarb1'; var2_cond = [];
% var2_name = 'Qarb2'; var2_cond = [];

% var2_name = 'S2';
% var2_name = 'S3';

% var2_name = 'A1'; var2_cond = [];
% var2_name = 'A1'; var2_cond = 'G';
% var2_name = 'A1'; var2_cond = 'G';  NameValue2 = {'SimulationNumber', 1};
% var2_name = 'A1'; var2_cond = 'H';
% var2_name = 'A2'; var2_cond = [];
% var2_name = 'A2'; var2_cond = 'state2';
% var2_name = 'A2'; var2_cond = 'state3';
% var2_name = 'A2'; var2_cond = 'state4';
% var2_name = 'A2'; var2_cond = 'state5';
% var2_name = 'A2 at state 4'; var2_cond = 'H';
% var2_name = 'A2 at state 4'; var2_cond = 'H + Red';
% var2_name = 'A2 at state 4'; var2_cond = 'Goal6';
% var2_name = 'ChoiceSwitch1'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch1'; var2_cond = [];
% var2_name = 'ChoiceSwitch1(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch1(GoalSwitch)'; var2_cond = [];
% var2_name = 'ChoiceSwitch2'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch2'; var2_cond = [];
% var2_name = 'ChoiceSwitch2(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch2(GoalSwitch)'; var2_cond = [];
% var2_name = 'ChoiceSwitch1n2'; var2_cond = [];
% var2_name = 'ChoiceSwitch1n2'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch1n2'; var2_cond = 'H';
% var2_name = 'ChoiceSwitch1n2'; var2_cond = 'uL';
% var2_name = 'ChoiceSwitch1n2'; var2_cond = 'uH';
% var2_name = 'ChoiceSwitch1n2(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoiceSwitch1n2(GoalStay)'; var2_cond = 'G';

% var2_name = 'ChoSwitch1n2(GoalSwitch)'; var2_cond = [];
% var2_name = 'ChoSwitch1(GoalSwitch)'; var2_cond = [];
% var2_name = 'ChoSwitch2(GoalSwitch)'; var2_cond = [];
% var2_name = 'ChoSwitch1n2(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoSwitch1n2(GoalSwitch)'; var2_cond = 'G'; NameValue2 = {'MovingMean', [0 4], 'Mat2Row', 1};
% var2_name = 'ChoSwitch1n2(GoalSwitch)'; var2_cond = 'G'; NameValue2 = {'MovingMean', [4 0], 'Mat2Row', 1};
% var2_name = 'ChoSwitch1n2(GoalSwitch)'; var2_cond = 'H';
% var2_name = 'ChoSwitch1(GoalSwitch)'; var2_cond = 'G';
% var2_name = 'ChoSwitch2(GoalSwitch)'; var2_cond = 'G';

var2_name = 'RTA1A2'; var2_cond = [];
% var2_name = 'RTA1A2'; var2_cond = 'G'; NameValue2 = {'Mat2Row', 1};
% var2_name = 'RTA1A2'; var2_cond = 'H'; NameValue2 = {'Mat2Row', 1};
% var2_name = 'RTA1'; var2_cond = [];
% var2_name = 'RTA1'; var2_cond = 'G';
% var2_name = 'RTA1'; var2_cond = 'H';
% var2_name = 'RTA2'; var2_cond = [];
% var2_name = 'RTA2'; var2_cond = 'G';
% var2_name = 'RTA2'; var2_cond = 'H';

% var2_name = 'R'; var2_cond = [];
% var2_name = 'R'; var2_cond = 'G';
% var2_name = 'Hit';
% var2_name = 'ChoOpt';
% var2_name = 'ChoOptG';
% var2_name = 'ChoOptH';
% var2_name = 'ChoOpt_old';
% var2_name = 'ChoOpt1n2'; var2_cond = [];   % choice optimality in Kim2019
% var2_name = 'ChoOpt1n2'; var2_cond = 'G';   % choice optimality in Kim2019
% var2_name = 'ChoOpt1n2TransEst';  var2_cond = [];
% var2_name = 'ChoOptTransEst1';  var2_cond = []; variant = 'default';
% var2_name = 'ChoOptTransEst1';  var2_cond = []; variant = 'OptimalA2SoftmaxA1TransEstblkCond';
% var2_name = 'ChoOptTransEst2';  var2_cond = []; variant = 'default';
% var2_name = 'ChoOptTransEst2';  var2_cond = []; variant = 'OptimalA2SoftmaxA2TransEstblkCond';
% var2_name = 'ChoOptBin1n2'; var2_cond = []; % choice optimality in Kim2021
% var2_name = 'ChoOptBin1n2'; var2_cond = 'G'; % choice optimality in Kim2021
% var2_name = 'ChoOptBin1n2'; var2_cond = 'Early'; % choice optimality in Kim2021
% var2_name = 'ChoOptBin1n2'; var2_cond = 'G uL'; % choice optimality in Kim2021
% var2_name = 'ChoOptBin1n2'; var2_cond = 'G uH'; % choice optimality in Kim2021
% var2_name = 'ChoOptBinNew1n2'; var2_cond = [];
% var2_name = 'ChoOptBinNew1n2'; var2_cond = 'G';
% var2_name = 'ChoOpt1'; var2_cond = []; variant = 'default'; 
% var2_name = 'ChoOpt1'; var2_cond = 'G'; variant = 'default'; 
% var2_name = 'ChoOpt1'; var2_cond = 'H'; variant = 'default'; 
% var2_name = 'ChoOpt1'; var2_cond = 'Goal6'; variant = 'default'; 
% var2_name = 'ChoOpt1'; var2_cond = 'Goal7'; variant = 'default'; 
% var2_name = 'ChoOpt1'; var2_cond = 'Goal8'; variant = 'default'; 
% var2_name = 'ChoOpt1'; var2_cond = []; variant = 'OptimalA2SoftmaxA1blkCond';
% var2_name = 'ChoOptBin1'; var2_cond = [];
% var2_name = 'ChoOptBin1'; var2_cond = 'G';
% var2_name = 'ChoOptBin1'; var2_cond = 'G uL';
% var2_name = 'ChoOptBin1'; var2_cond = 'G uH';
% var2_name = 'ChoOptBin1'; var2_cond = 'H';
% var2_name = 'ChoOptBin1'; var2_cond = 'HalfEarly';
% var2_name = 'ChoOptBin1'; var2_cond = 'HalfLate';
% var2_name = 'ChoOptBin1'; var2_cond = 'Ss1';
% var2_name = 'ChoOptBin1'; var2_cond = 'Ss2';
% var2_name = 'ChoOptBinNew1'; var2_cond = []; 
% var2_name = 'ChoOptBinNew1'; var2_cond = 'G'; 
% var2_name = 'ChoOpt2'; var2_cond = []; variant = 'default';
% var2_name = 'ChoOpt2'; var2_cond = 'G'; variant = 'default';
% var2_name = 'ChoOpt2'; var2_cond = 'H'; variant = 'default';
% var2_name = 'ChoOpt2'; var2_cond = 'Goal6'; variant = 'default';
% var2_name = 'ChoOpt2'; var2_cond = 'Goal7'; variant = 'default';
% var2_name = 'ChoOpt2'; var2_cond = 'Goal8'; variant = 'default';
% var2_name = 'ChoOpt2'; var2_cond = []; variant = 'OptimalA2SoftmaxA2blkCond';
% var2_name = 'ChoOptBin2'; var2_cond = [];
% var2_name = 'ChoOptBin2'; var2_cond = 'G';
% var2_name = 'ChoOptBin2'; var2_cond = 'G uL';
% var2_name = 'ChoOptBin2'; var2_cond = 'G uH';
% var2_name = 'ChoOptBin2'; var2_cond = 'H';
% var2_name = 'ChoOptBin2'; var2_cond = 'Goal7';
% var2_name = 'ChoOptBin2'; var2_cond = 'HalfEarly';
% var2_name = 'ChoOptBin2'; var2_cond = 'HalfLate';
% var2_name = 'ChoOptBin2'; var2_cond = 'Ss1';
% var2_name = 'ChoOptBin2'; var2_cond = 'Ss2';
% var2_name = 'ChoOptBinNew2'; var2_cond = [];
% var2_name = 'ChoOptBinNew2'; var2_cond = 'G';
% var2_name = 'ChoOptBinNew2'; var2_cond = 'Goal7';
% var2_name = 'ChoOpt_G1';
% var2_name = 'ChoOpt_G2';
% var2_name = 'ChoOpt_cL1';
% var2_name = 'ChoOpt_cL2';
% var2_name = 'ChoOpt_cH1';
% var2_name = 'ChoOpt_cH2';
% var2_name = 'ChoOpt1n2G';    
% var2_name = 'ChoOpt1n2H';    
% var2_name = 'ChoOptA1LeftBias';
% var2_name = 'ChoOptA2LeftBias';
% var2_name = 'ChoOptMaxGoalA1LeftBias';
% var2_name = 'ChoOptMaxGoalA2LeftBias';
% var2_name = 'ChoOptAchGoalRanA1LeftBias';
% var2_name = 'ChoOptAchGoalRanA2LeftBias';
% var2_name = 'ChoOptAchGoalOptA1LeftBias';
% var2_name = 'ChoOptAchGoalOptA2LeftBias';
% var2_name = 'ChoOptBinNaN1n2'; var2_cond = [];
% var2_name = 'ChoOptBinNaN1n2'; var2_cond = []; NameValue2 = {'SimulationNumber', 1};
% var2_name = 'ChoOptBinNaN1n2'; var2_cond = 'G'; NameValue = [];
% var2_name = 'ChoOptBinNaN1n2'; var2_cond = 'G'; NameValue = {'MovingMean', 5};
% var2_name = 'ChoOptBinNaN1n2'; var2_cond = 'G'; NameValue2 = {'MovingMean', [0 4], 'Mat2Row', 1};
% var2_name = 'ChoOptBinNaN1n2'; var2_cond = 'G'; NameValue2 = {'MovingMean', [4 0], 'Mat2Row', 1};
% var2_name = 'ChoOptBinNaN1n2'; var2_cond = 'H';
% var2_name = 'ChoOptBinNaN1n2'; var2_cond = 'uL';
% var2_name = 'ChoOptBinNaN1n2'; var2_cond = 'uH';
% var2_name = 'ChoOptBinNaN1'; var2_cond = [];
% var2_name = 'ChoOptBinNaN1'; var2_cond = 'G';
% var2_name = 'ChoOptBinNaN2'; var2_cond = [];
% var2_name = 'ChoOptBinNaN2'; var2_cond = 'G';

% var2_name = 'ChoCon ignoring context';
% var2_name = 'ChoCon ignoring context G';
% var2_name = 'ChoCon1';
% var2_name = 'ChoCon2';
% var2_name = 'ChoCon1 G';
% var2_name = 'ChoCon2 G';

% var2_name = 'ChoConsist1n2'; var2_cond = [];
% var2_name = 'ChoConsist1n2'; var2_cond = 'Goal(6,7,8)Stay';
% var2_name = 'ChoConsist1n2'; var2_cond = 'G';
% var2_name = 'ChoConsist1n2'; var2_cond = 'G'; NameValue2 = {'MovingMean', [0 4], 'Mat2Row', 1};
% var2_name = 'ChoConsist1n2'; var2_cond = 'G'; NameValue2 = {'MovingMean', [4 0], 'Mat2Row', 1};
% var2_name = 'ChoConsist1n2'; var2_cond = 'H';
% var2_name = 'ChoConsist1n2'; var2_cond = 'uL';
% var2_name = 'ChoConsist1n2'; var2_cond = 'uH';
% var2_name = 'ChoCon1 per goal'; var2_cond = [];
% var2_name = 'ChoCon2 per goal'; var2_cond = [];
% var2_name = 'ChoCon1 per goal'; var2_cond = 'G';
% var2_name = 'ChoCon2 per goal'; var2_cond = 'G';
% var2_name = 'ChoCon2 per goal'; var2_cond = 'state4';
% var2_name = 'ChoCon2 per goal'; var2_cond = 'state2';
% var2_name = 'StrictChoCon';
% var2_name = 'ChoCon1 per goal G';
% var2_name = 'ChoCon2 per goal G';
% var2_name = 'StrictChoCon G';

% var2_name = 'ChoOpt1n2MaxGoal'; % ChoOpt policy: random A2 (default), ChoOpt likelihood: NormVal (default), DoubleInComplex: true (default)
% var2_name = 'ChoOpt1n2AchGoalOpt'; 
% var2_name = 'ChoOpt1AchGoalOpt'; 
% var2_name = 'ChoOpt1AchGoalRan'; 

% ===== visualization parameters =====
switch Exp
    case {'Lee2014', 'arb_virtual_episode1', 'arb_virtual_episode100', 'arb_virtual_episode1000'}
% % %         vis_mode = 'time_series';   MarkerSize = 2; horz_slot = 2;
% % %         vis_mode = 'time_series';   MarkerSize = 2; horz_slot = 8;
% % %         vis_mode = 'time_series';   MarkerSize = 3; horz_slot = 8;
% % %         vis_mode = 'distribution_per_condition'; horz_slot = 8;
        % vis_mode = 'time_series';   MarkerSize = 12; horz_slot = 2;
        % vis_mode = 'time_series';   MarkerSize = 12; horz_slot = 3;
        % vis_mode = 'time_series';   MarkerSize = 12; horz_slot = 4;
        % vis_mode = 'time_series';   MarkerSize = 12; horz_slot = 8;        
        vis_mode = 'distribution_per_condition'; horz_slot = 11;
        % vis_mode = 'distribution_per_condition'; horz_slot = 8;
    case 'Kim2019'
%         vis_mode = 'time_series';   MarkerSize = 2; horz_slot = 3;
%         vis_mode = 'time_series';   MarkerSize = 4; horz_slot = 3;
        vis_mode = 'distribution_per_condition'; horz_slot = 11;
end
% vis_mode = 'time_series';   MarkerSize = 2; horz_slot = 8;
% vis_mode = 'time_series';   MarkerSize = 2; horz_slot = 11;
% vis_mode = 'time_series';   MarkerSize = 3; horz_slot = 2;
% vis_mode = 'distribution_per_condition'; horz_slot = 11; 

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
if ismember(var1_name, {'blkCond', 'uLuH EL', 'uLuH ucEL', 'Goal(6,7,8)xUC', 'uLuH Goal(6,7,8)Switch', ... 
        'Session x GT', 'Session x UC', 'Session x Goal(6,7,8)', 'Session x GS', 'EL x GS', 'GT x UC'})
% if strcmp(var1_name, 'blkCond') || strcmp(var1_name, 'uLuH ucEL') || strcmp(var1_name, 'uLuH Goal(6,7,8)Switch')
    N_sig_sbj = nan(idsize, 3); % main effect 1, main 2, interaction effect
    EL_storage = cell(idsize, 1);
    GT_storage = cell(idsize, 1);
    UC_storage = cell(idsize, 1);
    GS_storage = cell(idsize, 1);
elseif ismember(var1_name, {'EL x GT x UC'})
    N_sig_sbj = nan(idsize, 6); % main effect 1~3, interaction effect 1~3
    EL_storage = cell(idsize, 1);
    GT_storage = cell(idsize, 1);
    UC_storage = cell(idsize, 1);
    GS_storage = cell(idsize, 1);
else
    N_sig_sbj = nan(idsize, 1);
    EL_storage = cell(idsize, 1);
    GT_storage = cell(idsize, 1);
    UC_storage = cell(idsize, 1);
    GS_storage = cell(idsize, 1);
end
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
    if ismember(id, except_ids); continue; end
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
    
    if strcmp(var1_name, 'blkCond')
%         var1_class_names = {'LowEarly', 'LowLate', 'HighEarly', 'HighLate'};
    elseif strcmp(var1_name, 'GT x UC')
        var1_class_names = {'G low', 'G high', 'H low', 'H high'};
    elseif strcmp(var1_name, 'EL x GT x UC')
        var1_class_names = {'EGL', 'EGH', 'EHL', 'EHH', 'LGL', 'LGH', 'LHL', 'LHH'};
    elseif strcmp(var1_name, 'uLuH EL')
%         var1_class_names = {'LowEarly', 'LowLate', 'HighEarly', 'HighLate'};
        var1_class_names = {'Early Low', 'Early High', 'Late Low', 'Late High'};
    elseif strcmp(var1_name, 'uLuH ucEL')
        var1_class_names = {'LowEarly', 'LowLate', 'HighEarly', 'HighLate'};
    elseif strcmp(var1_name, 'uLL vs uHE')
        var1_class_names = {'LowLate', 'HighEarly'};
    elseif strcmp(var1_name, 'uHL vs uLE')
        var1_class_names = {'HighLate', 'LowEarly'};
    elseif strcmp(var1_name, 'Goal(6,7,8)xUC')
        var1_class_names = {'G6 L', 'G7 L', 'G8 L', 'G6 H', 'G7 H', 'G8 H'};
    elseif strcmp(var1_name, 'Session') && strcmp(Exp, 'Kim2019')
        var1_class_names = {'1', '2', '3', '4', '5', '6'};
    elseif strcmp(var1_name, 'uLuH Goal(6,7,8)Switch')
        var1_class_names = {'LowSwitch', 'LowStay', 'HighSwitch', 'HighStay'};
    elseif strcmp(var1_name, 'Session x UC')
        var1_class_names = {'1L', '2L', '3L', '4L', '5L', '1H', '2H', '3H', '4H', '5H'};
    elseif strcmp(var1_name, 'Session x Goal(6,7,8)')
        var1_class_names = {'G6 1', 'G6 2', 'G6 3', 'G6 4', 'G6 5', 'G7 1', 'G7 2', 'G7 3', 'G7 4', 'G7 5', 'G8 1', 'G8 2', 'G8 3', 'G8 4', 'G8 5'};
    elseif strcmp(var1_name, 'EL x GS')
        var1_class_names = {'StEarly', 'StLate', 'SwEarly', 'SwLate'};
    elseif strcmp(var1_name, 'Session x GS')
        var1_class_names = {'St1', 'St2', 'St3', 'St4', 'St5', 'Sw1', 'Sw2', 'Sw3', 'Sw4', 'Sw5'};
    end
    
    % =============== Individual visualization ===============
    if indiv_vis
        subplot(ceil(idsize/horz_slot), horz_slot, id);
    end
    
    switch vis_mode
        case 'time_series' % time series per sbj
            stem_per_var1 = [];
            for eli = 1:N_var1
                temp_var2 = var2;
                temp_var2(var1~=var1_unq(eli)) = nan;
                
                switch plot_mode
                    case 'plot'
                        s = plot(1:length(temp_var2), temp_var2, ...
                            'LineWidth', 0.1);
                    case 'stem'
%                         s = stem(1:length(temp_var2), temp_var2, ...
%                             'LineStyle', 'none', 'MarkerFaceColor', 'auto', ...
%                             'MarkerSize', MarkerSize);
                        s = stem(1:length(temp_var2), temp_var2, ...
                            'LineStyle', 'none', ...
                            'Marker', ".", ...
                            'MarkerSize', MarkerSize);

                    case 'movmean'
                        s = plot(1:length(temp_var2), ...
                            movmean(temp_var2, 3), ...
                            'LineWidth', 0.1);
                        
                end
                
                hold on
%                 title([var2_name ' per ' var1_name])
                
                stem_per_var1 = [stem_per_var1, s];
                
            end % var1 el
            
            
            
            hold off
            ylabel([var2_name var2_cond]);
            
        case 'distribution_per_condition' % distribution per sbj
            
            x = []; g = [];
            for eli = 1:N_var1
                sub_var2 = var2(var1==var1_unq(eli));
                x = [x; sub_var2'];
                g = [g; var1_unq(eli)*ones(length(sub_var2),1)];
            end % var1 el
            
            var1_storage{id} = cellfun(@num2str, num2cell(g), 'UniformOutput', false);
            var2_storage{id} = x;
            id_storage{id} = cellfun(@num2str, num2cell(id * ones(size(g))), 'UniformOutput', false);
%             g = g+1;
%             x = x+1;
            
            if indiv_vis
                % visualize sample points
                beeswarm(g, x, 'corral_style', 'random', 'overlay_style', 'ci'); hold on;
            end
            
            % ==================================================================
            % ANOVA for var1 effect on var2 (within-subject)
            % ==================================================================
            [N_sig_sbj, storages] = anovan_indiv(id,  var1_name, var1_cond, var2_name, var2_cond, ... 
                g, x, N_sig_sbj, storages, ...
                'var1_unq', var1_unq);
            
            hold off;
            
    end
    
    if id==idsize
        
        switch vis_mode
            case 'time_series' % time series per sbj
                % ===== legend =====
                legend(stem_per_var1, var1_class_names)
                title_text2 = [];
                
            case 'distribution_per_condition' % distribution per sbj
                switch size(N_sig_sbj, 2)
                    case 1
                        fprintf('\n # of subject with significant var1 effect on var2: %d \n', ...
                            nansum(any(N_sig_sbj, 2)))
                        title_text2 = ['Sig. sbj N=' sprintf('%d', nansum(any(N_sig_sbj, 2)))];
                    case 3
                        fprintf('\n # of subject with significant var1 effect on var2: %d \n', ...
                            nansum(N_sig_sbj))
                        title_text2 = ['Sig. sbj N=' ... 
                            sprintf('%d, ', nansum(N_sig_sbj(:, 1))), ...
                            sprintf('%d, ', nansum(N_sig_sbj(:, 2))), ...
                            sprintf('%d', nansum(N_sig_sbj(:, 3)))];
                    case 6
                        fprintf('\n # of subject with significant var1 effect on var2: %d \n', ...
                            nansum(N_sig_sbj))
                        title_text2 = ['Sig. sbj N=' ... 
                            sprintf('%d, ', nansum(N_sig_sbj(:, 1))), ...
                            sprintf('%d, ', nansum(N_sig_sbj(:, 2))), ...
                            sprintf('%d, ', nansum(N_sig_sbj(:, 3))), ...
                            sprintf('%d, ', nansum(N_sig_sbj(:, 4))), ...
                            sprintf('%d, ', nansum(N_sig_sbj(:, 5))), ...
                            sprintf('%d', nansum(N_sig_sbj(:, 6)))];
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
        [~, p, ~, stats] = ttest(var2_means(:, 1), var2_means(:, 2));
        title_prefix = 'paired t-test p=';
        % figure
        % figure('Position', [0 1080/2-80 1920/8 1080/2]) % [left bottom width height], main figure 1 
        figure('Position', [0 1080/2-80 1920/6 1080/2]) % [left bottom width height] 
        % figure('Position', [0 1080/2-80 1920/8 1080/3]) % [left bottom width height] 
        
% % %         subplot(1,3,1:2)
% % % %         bar(var2_means); title({var2_name, [' per ' var1_name], title_prefix, p})
% % %         bar(var2_means); title([title_prefix, num2str(p)])
% % % %         legend(string(var1_unq))
% % %         legend(var1_class_names)
% % % %         ylim([0 .7])
% % %         subplot(1,3,3)
%         boxplot(var2_means, 'BoxStyle', 'filled');
%         plot(var2_means', 'LineWidth', 0.01, 'Marker', 'o', 'Color', [.9 .9 .9]); hold on; 
        %
        
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
                text(1.5 - .04 * length(sig_text), 1 + .07, sig_text, ...
                    'FontSize', 20)
            
                if contains(var2_name, 'ChoOptBinNaN')
                    % ylim([.3 1.1])
                    
                end
                
                if contains(var2_name, 'ChoConsist')
                    ylim([.5 1.1])
                end
                
                xlim([.7 2.3])
                % xlim([.6 2.4])
                % xlim([.5 2.5])
                % xticks(1:2); xticklabels({'Low', 'High'})
                xticks(1:2); xticklabels({'Early', 'Late'})
                
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
        
%         % 2023-04-30 original ~
%         % 2-way ANOVA: var1 effect & subject effect
%         var2mean_vec = reshape(var2_means, [], 1);
%         var1class_vec = reshape(repmat(1:N_var1_global, idsize, 1), [], 1);
%         id_vec = reshape(repmat(1:idsize, N_var1_global, 1)', [], 1); 
%         id_vec = cellfun(@num2str, num2cell(id_vec), 'UniformOutput', false);
%         [p, tbl, stats] = anovan(var2mean_vec, {var1class_vec, id_vec}, ...
%             'model', 'linear', ...
%             'varnames', {var1_name, 'Subject'}, ...
%             'display', 'off'); 
% %         [p, tbl, stats] = anovan(var2mean_vec, {var1class_vec, id_vec}, ...
% %             'model', 'interaction', ...
% %             'varnames', {var1_name, 'Subject'}, ...
% %             'display', 'off'); 
%         p = p(1); % main effect of var1
%         title_prefix = '2-way ANOVA p=';
%         % ~ 2023-04-30 original


        
        % #################################################################
        % ------------------------ Special cases ------------------------
        % #################################################################
        
        title_text = [];
        
        [title_text] = anovan_group(idsize, except_ids, N_var1_global, ...
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

if contains(var2_name, 'RT')
    ylim([.6 1.3])
end

% figure
% hist_per_var1 = [];
% for eli = 1:N_var1_global
%     h = histogram(var2_means(:, eli), nbins, ... 
%         'EdgeAlpha', 0.4, ... 
%         'FaceAlpha', FaceAlpha, ...
%         'FaceColor', color_map{eli});
%     box off; hold on
%     q = stem(nanmean(var2_means(:, eli)), -0.5, ...
%         'Marker', 'o', ...
%         'MarkerFaceColor', h.FaceColor, ...
%         'MarkerSize', 10, ...
%         'MarkerEdgeColor', 'none');
%     hist_per_var1 = [hist_per_var1, h];
% end
% hold off
% title([title_prefix num2str(p)])
% xlabel([var2_name var2_cond])
% legend(hist_per_var1, var1_class_names)

figure
multcompare(stats)

% ========================= Saving results =========================
switch N_var1_global
    case 2
        save(['ConditionedBehavior/' Exp '_' var1_name var1_cond '_' var2_name var2_cond '_'], ...
            'var2_means', 'N_sig_sbj', 'except_ids')

    otherwise
        save(['ConditionedBehavior/' Exp '_' var1_name var1_cond '_' var2_name var2_cond '_'], ...
            'var2_means', 'N_sig_sbj', 'except_ids', 'tbl', 'stats')
end

disp('end')

%% taskVar profile (categorical)

LoadDefaultSettings;
% LoadMyColor;
close all

disp('start')

% Setting =================================================================

% ===== Experiment =====
% Exp = 'Lee2014';  idsize = 22; except_id = [];
% Exp = 'Lee2014';  idsize = 22; except_id = [19 20];
% Exp = 'arb_virtual_episode100'; idsize = 22; except_id = [19 20];
% Exp = 'arb_virtual_episode1'; idsize = 22; except_id = [19 20];
% Exp = 'arb_virtual_episode1000'; idsize = 22; except_id = [19 20];
% Exp = 'Heo2021';  idsize = 28; except_id = [];
Exp = 'Kim2019';    idsize = 21; except_id = [];
% except_id = 19; % for state distribution
% except_id = [10, 15, 19];   % insufficient sessions


% Data load
% N_taskCond = nan * ones(idsize, 2, 2);
% N_taskCond = nan * ones(idsize, 2, 3);
% N_taskCond = nan * ones(idsize, 2, 4);
% N_taskCond = nan * ones(idsize, 3, 2);
% N_taskCond = nan * ones(idsize, 3, 3);
% N_taskCond = nan * ones(idsize, 3, 4);
% N_taskCond = nan * ones(idsize, 4, 2);
% N_taskCond = nan * ones(idsize, 4, 3);
% N_taskCond = nan * ones(idsize, 4, 4);
% N_taskCond = nan * ones(idsize, 5, 2);
% N_taskCond = nan * ones(idsize, 5, 4);
% N_taskCond = nan * ones(idsize, 6, 4);


% ===== var1 =====
% categorize1 = [];
categorize1 = 2;
% categorize1 = 3;
NameValue1 = {};

% --- name ---
% var1_name = 'Half'; Nclass_var1 = 2; var1_cond = [];
% var1_name = 'Session'; Nclass_var1 = 5; var1_cond = []; NameValue1 = {}; % for Lee 2014
% var1_name = 'Session'; Nclass_var1 = 5; var1_cond = 'OptW5G low'; NameValue1 = {}; % for Lee 2014
% var1_name = 'Session'; Nclass_var1 = 5; var1_cond = 'OptFW5G low'; NameValue1 = {}; % for Lee 2014
% var1_name = 'Session'; Nclass_var1 = 4; % for Kim 2019
% var1_name = 'blkCond'; Nclass_var1 = 4; var1_cond = [];
% var1_name = 'GoalCond'; Nclass_var1 = 2; var1_cond = [];
% var1_name = 'GoaltypePhase';
% var1_name = 'SpecGoalPhase2';
% var1_name = 'SpecGoalPhase3';
% var1_name = 'UncCond'; Nclass_var1 = 2; var1_cond = [];
% var1_name = 'UncCond'; Nclass_var1 = 2; var1_cond = 'G';
% var1_name = 'UncCond'; Nclass_var1 = 2; var1_cond = 'H';
% var1_name = 'UncPhase';
% var1_name = 'CmplxCond'; Nclass_var1 = 2;
% var1_name = 'Interaction';
% var1_name = 'Goal'; Nclass_var1 = 4; var1_cond = []; % Lee 2014
% var1_name = 'Goal'; Nclass_var1 = 4; var1_cond = []; NameValue1 = {'SimulationNumber', 1};
% var1_name = 'Goal'; Nclass_var1 = 3; var1_cond = []; % Kim 2019
% var1_name = 'Goal(6,7,8)'; Nclass_var1 = 3; var1_cond = [];
% var1_name = 'Goal(6,7,8)'; Nclass_var1 = 3; var1_cond = 'G';
% var1_name = 'Goal(6,7,8)'; Nclass_var1 = 3; var1_cond = 'G'; NameValue1 = {'SimulationNumber', 1};
% var1_name = 'Goal(6,7,8)'; Nclass_var1 = 3; var1_cond = 'uL';
% var1_name = 'Goal(6,7,8)'; Nclass_var1 = 3; var1_cond = 'uH';
% var1_name = 'Goal(6,7,8)xUC'; Nclass_var1 = 6; var1_cond = [];
% var1_name = 'Goal(6,7,8)xUCSwitch'; Nclass_var1 = 2;
% var1_name = 'Session x GS'; Nclass_var1 = 10; var1_cond = []; 
% var1_name = 'AchGoalS1RandA'; Nclass_var1 = 3;
% var1_name = 'AchGoalS2RandA'; Nclass_var1 = 3;
% var1_name = 'Block1Goal'; Nclass_var1 = 3;
% var1_name = 'Block2Goal'; Nclass_var1 = 3;
% var1_name = 'Block3Goal'; Nclass_var1 = 3;
% var1_name = 'Block4Goal'; Nclass_var1 = 3;
% var1_name = 'Block1AchGoalS1RandA'; Nclass_var1 = 3;
% var1_name = 'Goal1Block'; Nclass_var1 = 4;
% var1_name = 'Goal2Block'; Nclass_var1 = 4;
% var1_name = 'Goal3Block'; Nclass_var1 = 4;
% var1_name = 'S2'; Nclass_var1 = 4; var1_cond = [];
% var1_name = 'UC_H in S2 4'; Nclass_var1 = 2;
% var1_name = 'S3';
% var1_name = 'R';
% var1_name = 'relMB_G';
% var1_name = 'relMAX_G';
% var1_name = 'relDiff_G';
% var1_name = 'relComp2'; Nclass_var1 = 2; var1_cond = 'G';
var1_name = 'PMB'; Nclass_var1 = 2; var1_cond = [];
% var1_name = 'PMB_G'; Nclass_var1 = 2;
% var1_name = 'ChoOpt_old';
% var1_name = 'ChoOpt'; Nclass_var1 = 2;
% var1_name = 'binChoOpt2';
% var1_name = 'binChoOpt2 perC'; Nclass_var1 = 2;
% var1_name = 'binChoOpt perC'; Nclass_var1 = 2;
% var1_name = 'ChoOptBinNaN1n2'; Nclass_var1 = 3; var1_cond = 'G';
% var1_name = 'ChoOptBinNaN1n2'; Nclass_var1 = 2; var1_cond = 'G'; NameValue1 = {'MovingMean', 5, 'Mat2Row', 1, 'Categorize', 2};
% var1_name = 'binChoOptBinNaN1n2'; Nclass_var1 = 2; var1_cond = 'G';
% var1_name = 'ChoOpt_G';
% var1_name = 'ChoiceSwitch1'; Nclass_var1 = 2; var1_cond = 'G';
% var1_name = 'ChoiceSwitch2'; Nclass_var1 = 2; var1_cond = 'G';
% var1_name = 'ChoiceSwitch1n2(GoalSwitch)'; Nclass_var1 = 3; var1_cond = 'G';
% var1_name = 'ChoSwitch1n2(GoalSwitch)'; Nclass_var1 = 3; var1_cond = 'G';
% var1_name = 'ChoSwitch1n2(GoalSwitch)'; Nclass_var1 = 2; var1_cond = 'G'; NameValue1 = {'MovingMean', 5};
% var1_name = 'ChoConsist1n2'; Nclass_var1 = 3; var1_cond = 'G';
% var1_name = 'ChoConsist1n2'; Nclass_var1 = 2; var1_cond = 'G'; NameValue1 = {'MovingMean', 5};
% var1_name = 'binChoConsist1n2'; Nclass_var1 = 2; var1_cond = 'G';
% var1_name = 'binrelMAX';


% ===== var2 =====
categorize2 = [];
% categorize2 = 2;
% categorize2 = 3;

% --- bar colors ---
bar_colors = [];
% bar_colors = [0.8 0.8 0.8; 0.8 0.2 0; 0.2 0.6 1; 0.1 0.1 0.1]; % for var2 = Coin
% bar_colors = [0.8 0.8 0.8; 0.8 0.2 0; 0.2 0.6 1]; % for var2 = Goal (Kim 2019)
% bar_colors = {[0.8 0.8 0.8], [0.8 0.2 0], [0 0.8 0.2], [0.1 0.1 0.1]};
% bar_colors = [0.9 0.35 0.35; 0.35 0.35 0.9; 0.9 0.9 0.35; ...
%                 0.75*[0.9 0.35 0.35]; 0.75*[0.35 0.35 0.9]; 0.75*[0.9 0.9 0.35]]; % for var2 = Goal(6,7,8)xUC

% --- name ---
% var2_name = 'blkCond';      lgd_names = {'1', '2', '3', '4'}; Nclass_var2 = 4; var2_cond = [];
% var2_name = 'blkCond2';      lgd_names = {'Glow', 'Ghigh', 'Hlow', 'Hhigh'}; Nclass_var2 = 4; var2_cond = [];
% var2_name = 'GoalCond';     lgd_names = {'G', 'H'};
% var2_name = 'UncCond';      lgd_names = {'Low', 'High'}; Nclass_var2 = 2; var2_cond = [];
% var2_name = 'UncSwitch'; lgd_names = {'0', '1'}; Nclass_var2 = 2;
% var2_name = 'CmplxCond'; lgd_names = {'Low', 'High'}; Nclass_var2 = 2;
var2_name = 'Goal';     lgd_names = {'1', '2', '3'}; Nclass_var2 = 3; var2_cond = [];
% var2_name = 'Goal';     lgd_names = {'-1', '6', '7', '8'};
% var2_name = 'Goal(6,7,8)';     lgd_names = {'6', '7', '8'}; Nclass_var2 = 3; var2_cond = [];
% var2_name = 'Goal(6,7,8)xUC'; lgd_names = {'6L', '7L', '8L', '6H', '7H', '8H'}; Nclass_var2 = 6; var2_cond = [];
% var2_name = 'Goal(6,7,8)xUCd1'; lgd_names = {'6L', '7L', '8L', '6H', '7H', '8H'}; Nclass_var2 = 6; var2_cond = [];
% var2_name = 'Goal PMB1';     lgd_names = {'1', '2', '3'}; Nclass_var2 = 3;
% var2_name = 'GoalSwitch'; lgd_names = {'0', '1'}; Nclass_var2 = 2; var2_cond = [];
% var2_name = 'GoalSwitch'; lgd_names = {'0', '1'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'Goal(6,7,8)Switch'; lgd_names = {'0', '1'}; Nclass_var2 = 2;
% var2_name = 'Block1AchGoalS1RandA'; lgd_names = {'1', '2', '3'}; Nclass_var2 = 3;
% var2_name = 'Block2AchGoalS1RandA'; lgd_names = {'1', '2', '3'}; Nclass_var2 = 3;
% var2_name = 'Block3AchGoalS1RandA'; lgd_names = {'1', '2', '3'}; Nclass_var2 = 3;
% var2_name = 'Block4AchGoalS1RandA'; lgd_names = {'1', '2', '3'}; Nclass_var2 = 3;
% var2_name = 'AchGoalS1RandA'; lgd_names = {'1', '2', '3'}; Nclass_var2 = 3;
% var2_name = 'AchGoalS1OptA'; lgd_names = {'1', '2', '3'}; Nclass_var2 = 3;
% var2_name = 'AchGoalS2RandA'; lgd_names = {'1', '2', '3'}; Nclass_var2 = 3;
% var2_name = 'AchGoalS2OptA'; lgd_names = {'1', '2', '3'}; Nclass_var2 = 3;
% var2_name = 'S2';       lgd_names = {'2', '3'}; Nclass_var2 = 2;
% var2_name = 'S2';       lgd_names = {'2', '3', '4', '5'}; Nclass_var2 = 4; var2_cond = [];
% var2_name = 'S2';       lgd_names = {'2', '3', '4', '5'}; Nclass_var2 = 4; var2_cond = 'G';
% var2_name = 'A1'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = []; % for Lee 2014
% var2_name = 'A1'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = 'uL'; % for Lee 2014
% var2_name = 'A1'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = 'uH'; % for Lee 2014
% var2_name = 'A1'; lgd_names = {'R', 'L'}; Nclass_var2 = 2; % for Kim 2019
% var2_name = 'A2'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = [];
% var2_name = 'A2'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = 'uL';
% var2_name = 'A2'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = 'uH';
% var2_name = 'A2'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = 'state2';
% var2_name = 'A2'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = 'state3';
% var2_name = 'A2'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = 'state4';
% var2_name = 'A2'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = 'state5';
% var2_name = 'A2'; lgd_names = {'L1', 'R1', 'L2', 'R2'}; Nclass_var2 = 4; var2_cond = [];
% var2_name = 'A2'; lgd_names = {'L1', 'R1', 'L2', 'R2'}; Nclass_var2 = 4; var2_cond = 'uL';
% var2_name = 'A2'; lgd_names = {'L1', 'R1', 'L2', 'R2'}; Nclass_var2 = 4; var2_cond = 'uH';
% var2_name = 'A2 at state 4'; lgd_names = {'L', 'R'}; Nclass_var2 = 2; var2_cond = 'H';
% var2_name = 'S3';   lgd_names = {'6', '7', '8', '9'}; Nclass_var2 = 4; var2_cond = [];
% var2_name = 'S3';   lgd_names = {'6', '7', '8', '9'}; Nclass_var2 = 4; var2_cond = 'G';
% var2_name = 'Coin';   lgd_names = {'1', '2', '3', '4'}; Nclass_var2 = 4;
% var2_name = 'R';    lgd_names = {'0', '10', '20', '40'}; Nclass_var2 = 4; var2_cond = [];
% var2_name = 'Hit';    lgd_names = {'0', '1'};
% var2_name = 'episode';
% var2_name = 'binPMB'; lgd_names = {'low', 'high'}; Nclass_var2 = 2;
% var2_name = 'ChoOpt_old'; lgd_names = {'0', '1'}; Nclass_var2 = 2;
% var2_name = 'ChoOpt_G1';    % lgd_names = {'0', '1'};
% var2_name = 'ChoOpt_G2';
% var2_name = 'ChoOptBinNaN1n2';   lgd_names = {'0', '0.5', '1'}; Nclass_var2 = 3; var2_cond = 'G'; bar_colors = [0.2 0.2 0.2; 0.5 0.5 0.5; 0.8 0.8 0.8]; 
% var2_name = 'ChoOptBin1';   lgd_names = {'0', '0.5', '1'}; Nclass_var2 = 3; var2_cond = 'G'; bar_colors = [0.2 0.2 0.2; 0.5 0.5 0.5; 0.8 0.8 0.8]; 
% var2_name = 'ChoOptBin2';   lgd_names = {'0', '0.5', '1'}; Nclass_var2 = 3; var2_cond = 'G'; bar_colors = [0.2 0.2 0.2; 0.5 0.5 0.5; 0.8 0.8 0.8]; 
% var2_name = 'ChoiceSwitch1(GoalSwitch)'; lgd_names = {'0', '1'}; Nclass_var2 = 2; var2_cond = [];
% var2_name = 'ChoiceSwitch2(GoalSwitch)'; lgd_names = {'0', '1'}; Nclass_var2 = 2; var2_cond = [];
% var2_name = 'ChoCon1';
% var2_name = 'ChoCon2';
% var2_name = 'ChoCon1 per goal'; lgd_names = {'0', '1'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'ChoCon2 per goal'; lgd_names = {'0', '1'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'RPE2'; lgd_names = {'low', 'high'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'RPE3'; lgd_names = {'low', 'high'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'SPE2'; lgd_names = {'low', 'high'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'SPE3'; lgd_names = {'low', 'high'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'relDiff2'; lgd_names = {'low', 'high'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'relDiff3'; lgd_names = {'low', 'high'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'relComp2'; lgd_names = {'-', '+'}; Nclass_var2 = 2; var2_cond = 'G';
% var2_name = 'relComp3'; lgd_names = {'-', '+'}; Nclass_var2 = 2; var2_cond = 'G';

N_taskCond = nan * ones(idsize, Nclass_var1, Nclass_var2);

% --- class legend ---
% lgd_names = {};
% lgd_names = {'change', 'stay'};
% lgd_names = {'low', 'high'};
% lgd_names = {'1', '2', '3', '4'};
% lgd_names = {'-1', '6', '7', '8'};
% lgd_names = {'left', 'right'};
% lgd_names = {'non', 'opt'};
% lgd_names = {'1', '2', '3'};
% lgd_names = {'0', '10', '20', '40'};
% lgd_names = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16'};

% ===== N sample threshold =====
N_lab_threshold = 5;
% N_lab_threshold = 10;
% N_lab_threshold = 20;



freq_check = 1;

% figure('Position', [100 100 620 700]) 
figure % for time series plot

for id = 1:idsize
    if ismember(id, except_id); continue; end
    fprintf('%d ', id)
    
    switch var1_name
        
        case 'Half'
            var1 = arbMBMF_load_var(Exp, 'Session', id, var1_cond);
            var1(1:floor(length(var1)/2)) = 1;
            var1(ceil(length(var1)/2):end) = 2;
        case 'Block1Goal'
            var1 = arbMBMF_load_var(Exp, 'Goal', id, var1_cond);
            blk = arbMBMF_load_var(Exp, 'blkCond', id, var1_cond);
            var1(blk~=1) = nan;
        case 'Block2Goal'
            var1 = arbMBMF_load_var(Exp, 'Goal', id, var1_cond);
            blk = arbMBMF_load_var(Exp, 'blkCond', id, var1_cond);
            var1(blk~=2) = nan;
        case 'Block3Goal'
            var1 = arbMBMF_load_var(Exp, 'Goal', id, var1_cond);
            blk = arbMBMF_load_var(Exp, 'blkCond', id, var1_cond);
            var1(blk~=3) = nan;
        case 'Block4Goal'
            var1 = arbMBMF_load_var(Exp, 'Goal', id, var1_cond);
            blk = arbMBMF_load_var(Exp, 'blkCond', id, var1_cond);
            var1(blk~=4) = nan;
        case 'Goal1Block'
            var1 = arbMBMF_load_var(Exp, 'blkCond', id, var1_cond);
            goal = arbMBMF_load_var(Exp, 'Goal', id, var1_cond);
            var1(goal~=1) = nan;
        case 'Goal2Block'
            var1 = arbMBMF_load_var(Exp, 'blkCond', id, var1_cond);
            goal = arbMBMF_load_var(Exp, 'Goal', id, var1_cond);
            var1(goal~=2) = nan;
        case 'Goal3Block'
            var1 = arbMBMF_load_var(Exp, 'blkCond', id, var1_cond);
            goal = arbMBMF_load_var(Exp, 'Goal', id, var1_cond);
            var1(goal~=3) = nan;
        case 'binChoOpt2 perC'
            var2 = arbMBMF_load_var(Exp, var2_name, id, var1_cond); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            ChoOpt2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var1_cond);
            ChoOpt2 = ChoOpt2(2,:);
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                [binChoOpt_temp, cutoff] = regr2class(ChoOpt2(temp_v2_idx), 2);
%                 fprintf(cutoff.type)
%                 temp_nonopt_idx = temp_v2_idx(binChoOpt_temp==0);
%                 var1(temp_nonopt_idx) = 0;
                temp_opt_idx = temp_v2_idx(binChoOpt_temp==1);
                var1(temp_opt_idx) = 1;
            end
        case 'binChoOpt perC'
            var2 = arbMBMF_load_var(Exp, var2_name, id, var1_cond); var1 = zeros(size(var2));
            var2set = unq_elms(var2);
            ChoOpt = arbMBMF_load_var(Exp, 'ChoOpt', id, var1_cond);
            for i = 1:length(var2set)
                temp_v2_idx = find(var2 == var2set(i));
                [binChoOpt_temp, cutoff] = regr2class(ChoOpt(temp_v2_idx), 2);
                temp_opt_idx = temp_v2_idx(binChoOpt_temp==1);
                var1(temp_opt_idx) = 1;
            end
        case 'binrelMAX'
            var1 = arbMBMF_load_var(Exp, 'relMAX', id, var1_cond);
            var1 = regr2class(var1, 2);
        case 'binChoOpt2'
            var1 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var1_cond);
            var1 = regr2class(var1(2,:), 2);
        case 'binChoOpt'
            var1 = arbMBMF_load_var(Exp, 'ChoOpt', id, var1_cond);
            var1 = regr2class(var1, 2);
        case 'binChoOptBinNaN1n2'
            var1 = arbMBMF_load_var(Exp, 'ChoOptBinNaN1n2', id, var1_cond);
            var1 = mean(var1, 1, 'omitnan');
            var1(var1==.5) = 0;
        case 'binChoConsist1n2'
            var1 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var1_cond);
            var1 = mean(var1, 1, 'omitnan');
            var1(var1==.5) = 0;
        case 'ChoOpt_G'
            var1 = arbMBMF_load_var(Exp, 'ChoOpt', id, var1_cond);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var1_cond); % 1:G, 2:H
            var1(gt~=1) = nan;
        case 'ChoCon1'
            var1 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var1_cond);
            var1 = var1(1,:);
        case 'ChoCon2'
            var1 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var1_cond);
            var1 = var1(2,:);
        case 'UC_H in S2 4'
            var1 = arbMBMF_load_var(Exp, 'UncCond', id, var1_cond);
            S2 = arbMBMF_load_var(Exp, 'S2', id, var1_cond);
            GT = arbMBMF_load_var(Exp, 'GoalCond', id, var1_cond);
            var1(S2~=4) = nan; var1(GT~=2) = nan;
        otherwise
            if isempty(NameValue1)
                var1 = arbMBMF_load_var(Exp, var1_name, id, var1_cond);
            else
                var1 = arbMBMF_load_var(Exp, var1_name, id, var1_cond, ...
                    NameValue1{:});
            end
    end

    if size(var1, 1) == 2
        var1 = mean(var1, 1, 'omitnan');
    end

    if ~isempty(categorize1)   % continuous to categorical
        [var1, cutoff] = regr2class(var1, categorize1);
        [~, N_el] = unq_elms(var1);
        disp([cutoff.type ' ' num2str(N_el') ' ' num2str(cutoff.value)])
    end
    
    switch var2_name
        case 'Goal PMB1'
            var2 = arbMBMF_load_var(Exp, 'Goal', id, var2_cond);
            PMB = arbMBMF_load_var(Exp, 'PMB', id, var2_cond); binPMB = regr2class(PMB, 2);
            var2(binPMB~=1) = nan;
        case 'Block1AchGoalS1RandA'
            var2 = arbMBMF_load_var(Exp, 'AchGoalS1RandA', id, var2_cond);
            blk = arbMBMF_load_var(Exp, 'blkCond', id, var2_cond);
            var2(blk~=1) = nan;
        case 'Block2AchGoalS1RandA'
            var2 = arbMBMF_load_var(Exp, 'AchGoalS1RandA', id, var2_cond);
            blk = arbMBMF_load_var(Exp, 'blkCond', id, var2_cond);
            var2(blk~=2) = nan;
        case 'Block3AchGoalS1RandA'
            var2 = arbMBMF_load_var(Exp, 'AchGoalS1RandA', id, var2_cond);
            blk = arbMBMF_load_var(Exp, 'blkCond', id, var2_cond);
            var2(blk~=3) = nan;
        case 'Block4AchGoalS1RandA'
            var2 = arbMBMF_load_var(Exp, 'AchGoalS1RandA', id, var2_cond);
            blk = arbMBMF_load_var(Exp, 'blkCond', id, var2_cond);
            var2(blk~=4) = nan;
        case 'binPMB'
            PMB = arbMBMF_load_var(Exp, 'PMB', id, var2_cond);
            var2 = regr2class(PMB, 2);
        case 'ChoCon1'
            var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
            var2 = var2(1,:);
        case 'ChoCon2'
            var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
            var2 = var2(2,:);
        case 'ChoCon1 per goal'
            var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
            var2 = var2(1,:);
        case 'ChoCon2 per goal'
            var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
            var2 = var2(2,:);
        case 'ChoOpt_G1'
            var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(1,:);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2(gt~=1) = nan;
        case 'ChoOpt_G2'
            var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(2,:);
            gt = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
            var2(gt~=1) = nan;
        otherwise
            var2 = arbMBMF_load_var(Exp, var2_name, id, var2_cond);
    end

    if size(var2, 1) == 2
        var2 = mean(var2, 1, 'omitnan');
    end

    if ~isempty(categorize2)   % continuous to categorical
        [var2, cutoff] = regr2class(var2, categorize2);
        [~, N_el2] = unq_elms(var2);
        disp([cutoff.type ' ' num2str(N_el2') ' ' num2str(cutoff.value)])
    end    
    
% Time series visualization -----------------------------------------------
subplot(ceil(idsize/2), 2, id);
% ax = gca; ax.FontSize = 12;

% plot(1:length(var1), var1, 1:length(var1), var2)

stem(var1, 'LineStyle', 'none', 'MarkerFaceColor', 'auto'); hold on

% plot(1:length(var1), var2); hold off
stem(var2, 'Marker', '.'); hold off

legend([var1_name var1_cond], [var2_name var2_cond])
    
% ===============================================
    if freq_check
        uni1 = unq_elms(var1);
        uni2 = unq_elms(var2);
        
        if length(uni1) ~= Nclass_var1
            fprintf('var1 class lack: ')
            fprintf('(%d) ', uni1)
        end
        
        if length(uni2) ~= Nclass_var2
            fprintf('var2 class lack: ')
            fprintf('(%d) ', uni2)
        end
        
        for i1 = 1:length(uni1)
            for i2 = 1:length(uni2)
                N_taskCond(id, i1, i2) = sum(var1==uni1(i1) & var2==uni2(i2));
%                 N_taskCond(id, i1, i2) = sum(var1==uni1(i1) & var2==uni2(i2))... 
%                     /length(var1);
            end
        end
    end
    
end % id

lgd = legend([var1_name var1_cond], [var2_name var2_cond]); lgd.FontSize = 11;
disp(' ')

% Class frequency visualization -------------------------------------------

if freq_check
    % Visualization
    % ==============
    % ==============
    % ==============
    % horz_slot = 1;
    horz_slot = 2;
%     horz_slot = 3;
%     horz_slot = 5;
    % ==============
    % ==============
    % ==============
    
    figure('Position', [0 1080/2-80 1920 1080/2])
    % colororder(colors)
    % lineStyles = linspecer(10, 'qualitative');
    if str2double(erase(version('-release'),{'a', 'b'})) > 2018
        colororder(linspecer(Nclass_var2 + 3, 'qualitative'))
    end
%     figure
    %     for i1 = 1:length(uni1)
    for i1 = 1:Nclass_var1
%         if i1 <= length(uni1)
            subplot(ceil(length(uni1)/horz_slot), horz_slot, i1)
            barmat = squeeze(N_taskCond(:, i1, :));
            if isempty(bar_colors)
                B = bar([barmat; nanmean(barmat,1)], ...
                    'EdgeColor', 'none'); box off;
            else
                B = bar([barmat; nanmean(barmat,1)], ...
                    'EdgeColor', 'none'); box off;
                for i2 = 1:Nclass_var2
                    B(i2).FaceColor = bar_colors(i2, :);
                end
            end
            
            title(sprintf([var1_name var1_cond cell2mat(cellfun(@num2str, NameValue1, 'UniformOutput', false)) ... 
                ' %d ' var2_name var2_cond], uni1(i1)), ... 
                'FontSize', 14)
            %     legend(B(:), {'1', '2', '3', '4'})
            if ~isempty(lgd_names) && i1==Nclass_var1
                try
                    legend(B(:), lgd_names, 'FontSize', 12, 'Location', 'best');
                catch
                    legend(lgd_names{:}, 'FontSize', 12, 'Location', 'best');
                end
                hold off
            end
            ylabel('Frequency', 'FontSize', 15)
            xlabel('Subject ID', 'FontSize', 15)
            xlim([0 idsize+2])
            if length(diff(nanmean(barmat,1)))==1 % only for binary var2
                var2_diff_mean = abs(diff(nanmean(barmat,1)));
                hold on; line([0 idsize+3], [var2_diff_mean var2_diff_mean], ...
                    'LineStyle', '--')
                hold off
            end
            if N_lab_threshold
                hold on; line([0 idsize+3], [N_lab_threshold N_lab_threshold], ...
                    'LineStyle', '--')
                hold off
            end
%         end
    end
end

% contingency table & Cramer's V
% nRowVis = 2;
nRowVis = 3;
nColVis = ceil(idsize/nRowVis);
figure
for id = 1:idsize
    if ismember(id, except_id); continue; end
    subplot(nRowVis, nColVis, id)
    temp_table = squeeze(N_taskCond(id, :, :));
    imagesc(temp_table)
    stats = mestab(temp_table);
    ylabel([var1_name var1_cond])
    xlabel([var2_name var2_cond])
    if isfield(stats, 'cramerV')
        title(['cramerV ' num2str(stats.cramerV)])
    elseif isfield(stats, 'riskDifference')
        title(['RD ' num2str(stats.riskDifference)])
    end
end % id

disp('end');


%% Effect analysis (Categorical; one vs all SVM)
% 2021-11-17

LoadDefaultSettings;
close all

disp('start')

% Setting =================================================================

% ===== Experiment =====
% Exp = 'Lee2014';  idsize = 22;
% Exp = 'Heo2018';  idsize = 28;
Exp = 'Kim2019';    idsize = 21;

% ===== response =====
categorize_y = [];
% categorize_y = 2;
% categorize1 = 3;

% --- name ---
% y_name = 'Coin';  categ = 1;    Nclass_y_temp = 4;
% y_name = 'Coin_val';  categ = 1;  Nclass_y_temp = 3;
% y_name = 'Coin_val b1';  categ = 1;  Nclass_y_temp = 3;
% y_name = 'Coin_val b2';  categ = 1;  Nclass_y_temp = 3;
% y_name = 'Coin_val b3';  categ = 1;  Nclass_y_temp = 3;
y_name = 'Coin_val b4';  categ = 1;  Nclass_y_temp = 3;
% y_name = 'Goal';  categ = 1;
% y_name = 'ChoOpt';  categ = 0;
% y_name = 'ChoOpt1';  categ = 0;
% y_name = 'ChoOpt2';  categ = 0;

% --- N class ---
if ~isempty(categorize_y)
    Nclass_y = categorize_y;
    categ = 1;
else
    Nclass_y = Nclass_y_temp;
end

if Nclass_y == 2 || ~categ
    yloop = 1;
else
    yloop = Nclass_y;
end

% ===== feature =====
x_shuffle = 0;
% x_shuffle = 1;

% --- name ---
% X_names = {'CmplxCond', 'UncCond', 'Goal'}; Nclasses_X = [2 2 3]; Xticknames = {'CX', 'UC', 'G1', 'G2', 'G3'};
% X_names = {'CmplxCond', 'UncCond'}; Nclasses_X = [2 2]; Xticknames = {'CX', 'UC'};
X_names = {'Goal'}; Nclasses_X = 3; Xticknames = {'G1', 'G2', 'G3'};

% --- N class ---
N_coef = sum(Nclasses_X) - sum(Nclasses_X==2);

KFold = 5;
% X_coeffs = nan(idsize, Nclass_y, sum(Nclasses_X));
X_coeffs = nan(idsize, Nclass_y, N_coef);
CVaccs = nan(idsize, 1);
for id = 1:idsize
    fprintf('%d ', id)
    
    switch y_name
        case 'Coin_val'
            y = arbMBMF_load_var(Exp, 'Coin', id, []);
            y(y==4) = nan;
        case 'Coin_val b1'
            y = arbMBMF_load_var(Exp, 'Coin', id, []);
            b = arbMBMF_load_var(Exp, 'blkCond', id, []);
            y(y==4) = nan;
            y(b~=1) = nan;
        case 'Coin_val b2'
            y = arbMBMF_load_var(Exp, 'Coin', id, []);
            b = arbMBMF_load_var(Exp, 'blkCond', id, []);
            y(y==4) = nan;
            y(b~=2) = nan;
        case 'Coin_val b3'
            y = arbMBMF_load_var(Exp, 'Coin', id, []);
            b = arbMBMF_load_var(Exp, 'blkCond', id, []);
            y(y==4) = nan;
            y(b~=3) = nan;
        case 'Coin_val b4'
            y = arbMBMF_load_var(Exp, 'Coin', id, []);
            b = arbMBMF_load_var(Exp, 'blkCond', id, []);
            y(y==4) = nan;
            y(b~=4) = nan;
        case 'ChoOpt1'
            y = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, []);
            y = y(1,:);
        case 'ChoOpt2'
            y = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, []);
            y = y(2,:);
        otherwise
            y = arbMBMF_load_var(Exp, y_name, id, []);
    end
    % categorization of continous y
    if ~isempty(categorize_y)
        y = regr2class(y, categorize_y);
    end
        
    X = nan(N_coef, length(y));
%     X = nan(sum(Nclasses_X), length(y));
    onehot_idx = 0;
    % one hot X
    for xi = 1:length(X_names)
        
        x = arbMBMF_load_var(Exp, X_names{xi}, id, []);
        % shuffling for the null result
        if x_shuffle; x = x(randperm(length(x))); end
        xSet = unq_elms(x);
        % fprintf('(%d) ', xSet);
        
        % x class sanity check
        if length(xSet) ~= Nclasses_X(xi); fprintf('x class problem: '); fprintf('(%d) ', xSet); end
        
        % one hot encoding
        if Nclasses_X(xi) == 2
            x_onehot = (x == max(xSet));
            onehot_idx_new = onehot_idx + 1;
        else
            x_onehot = repmat(x, Nclasses_X(xi), 1);
            x_onehot = (x_onehot == xSet');
            onehot_idx_new = onehot_idx + Nclasses_X(xi);
        end
        X(onehot_idx+1:onehot_idx_new, :) = x_onehot;
        onehot_idx = onehot_idx_new;
    end
    
%     y_ = y(~isnan(y));
%     X_ = X(:, ~isnan(y));
    [y_, X_] = lab_pat_undersample(y, X);
    if categ
        CVMdl = fitcecoc(X_, y_, ...
            'Coding', 'onevsall', ...
            'Learners', templateLinear('Learner', 'svm', 'Regularization', 'ridge', 'Solver', 'dual'), ...
            'ObservationsIn', 'columns', ...
            'KFold', KFold);
    else
        CVMdl = fitrlinear(X_, y_, ...
            'Learner', 'svm', ...
            'ObservationsIn', 'columns', ...
            'Regularization', 'ridge', ...
            'Solver', 'dual', ...
            'KFold', KFold);
    end
    % generalized accuracy (CV acc)
    CVaccs(id) = 1 - kfoldLoss(CVMdl);    
    
    % y class sanity check
    if categ
        ySet = unq_elms(y);
        if length(ySet) ~= Nclass_y; fprintf('y class problem: '); fprintf('(%d) ', ySet); end
    end
    for yi = 1:yloop
        Beta_fold = nan(KFold, N_coef);
%         Beta_fold = nan(KFold, sum(Nclasses_X));
        for fold = 1:KFold
            if categ
                Beta_fold(fold, :) = CVMdl.Trained{fold}.BinaryLearners{yi}.Beta';
            else
                Beta_fold(fold, :) = CVMdl.Trained{fold}.Beta';
            end
        end
        X_coeffs(id, yi, :) = mean(Beta_fold, 1);
    end
    
end

disp('end')

% visualization
close all

figure
for yi = 1:yloop
    
%     stats = [CVaccs squeeze(X_coeffs(:, yi, :))];
%     px = reshape(repmat((1:1+N_coef), idsize, 1), [], 1);
    stats = squeeze(X_coeffs(:, yi, :));
    px = reshape(repmat((1:N_coef), idsize, 1), [], 1);
    py = reshape(stats, [], 1);
    
    % p vals
    hbox = nan(1, size(stats, 2));
    pbox = nan(1, size(stats, 2));
    for t = 1:size(stats, 2)
%         if t == 1
%             [h, p] = ttest(stats(:, 1), 1/Nclass_y);
%         else
            [h, p] = ttest(stats(:, t));
%         end
        hbox(t) = h;
        pbox(t) = p;
    end
    [~, p] = ttest(CVaccs, 1/Nclass_y);
    
    subplot(1, yloop, yi)
%     beeswarm(px, py, 'corral_style', 'gutter', 'overlay_style', 'ci');
    beeswarm(px, py, 'corral_style', 'random', 'overlay_style', 'ci');
    if categ
        title({[num2str(ySet(yi)) ' vs others'], num2str(hbox), ...
            ['Mean accuracy:' num2str(mean(CVaccs))], ... 
            ['SE:' num2str(std(CVaccs)/sqrt(idsize)) ', p=' num2str(p)]})
    else
        title({y_name, num2str(hbox), ... 
            ['Mean accuracy:' num2str(mean(CVaccs))], ... 
            ['SE:' num2str(std(CVaccs)/sqrt(idsize)) ', p=' num2str(p)]})
    end
    xticks(1:N_coef)
    xticklabels(Xticknames)
end

disp('end');


%% Effect analysis (continous; fitlm)
% 2021-11-17

% LoadDefaultSettings;
close all

disp('start')

% Setting =================================================================

% ===== Experiment =====
% Exp = 'Lee2014';  idsize = 22;
% Exp = 'Heo2018';  idsize = 28;
Exp = 'Kim2019';    idsize = 21;
except_id = [];
% except_id = 19; % for state distribution
% except_id = [10, 15, 19];   % insufficient sessions

% ===== response =====
categorize_y = [];
% categorize_y = 2;
% categorize1 = 3;

% --- name ---
% y_name = 'Coin';
% y_name = 'Coin_val';
% y_name = 'Goal';
% y_name = 'ChoOpt';  categ = 0;
% y_name = 'ChoOpt1';  categ = 0;
y_name = 'ChoOpt2';  categ = 0;

% --- N class ---
if ~isempty(categorize_y)
    Nclass_y = categorize_y;
    categ = 1;
else
    % Nclass_y = 4;
    Nclass_y = 3;
end
if Nclass_y == 2 || ~categ
    yloop = 1;
else
    yloop = Nclass_y;
end

% ===== feature =====
% --- name ---
X_names = {'CmplxCond', 'UncCond', 'Goal'}; Nclasses_X = [2 2 3];
% X_names = {'CmplxCond', 'UncCond'}; Nclasses_X = [2 2];

% --- N class ---
N_coef = sum(Nclasses_X) - sum(Nclasses_X==2);

KFold = 5;
% X_coeffs = nan(idsize, Nclass_y, sum(Nclasses_X));
X_coeffs = nan(idsize, Nclass_y, N_coef+1);
for id = 1:idsize
    fprintf('%d ', id)
    
    switch y_name
        case 'Coin_val'
            y = arbMBMF_load_var(Exp, 'Coin', id, []);
            y(y==4) = nan;
        case 'ChoOpt1'
            y = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, []);
            y = y(1,:);
        case 'ChoOpt2'
            y = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, []);
            y = y(2,:);
        otherwise
            y = arbMBMF_load_var(Exp, y_name, id, []);
    end
    % categorization of continous y
    if ~isempty(categorize_y)
        y = regr2class(y, categorize_y);
    end
        
    X = nan(N_coef, length(y));
%     X = nan(sum(Nclasses_X), length(y));
    onehot_idx = 0;
    % one hot X
    for xi = 1:length(X_names)
        
        x = arbMBMF_load_var(Exp, X_names{xi}, id, []);
        xSet = unq_elms(x);
        % fprintf('(%d) ', xSet);
        
        % x class sanity check
        if length(xSet) ~= Nclasses_X(xi); fprintf('x class problem: '); fprintf('(%d) ', xSet); end
        
        % one hot encoding
        if Nclasses_X(xi) == 2
            x_onehot = (x == max(xSet));
            onehot_idx_new = onehot_idx + 1;
        else
            x_onehot = repmat(x, Nclasses_X(xi), 1);
            x_onehot = (x_onehot == xSet');
            onehot_idx_new = onehot_idx + Nclasses_X(xi);
        end
        X(onehot_idx+1:onehot_idx_new, :) = x_onehot;
        onehot_idx = onehot_idx_new;
    end
    
    y_ = y(~isnan(y));
    X_ = X(:, ~isnan(y));
    if categ
        CVMdl = fitcecoc(X_, y_, ...
            'Coding', 'onevsall', ...
            'Learners', templateLinear('Learner', 'svm', 'Regularization', 'ridge', 'Solver', 'dual'), ...
            'ObservationsIn', 'columns', ...
            'KFold', KFold);
    else
        mdl = fitlm(X_', y_');
    end
    % generalized accuracy (CV acc)
    mdl
    
    for yi = 1:yloop
        X_coeffs(id, yi, :) = mdl.Coefficients.Estimate;
    end
    
end

disp('end')

%% visualization
close all

figure
for yi = 1:yloop
    
    stats = squeeze(X_coeffs(:, yi, :));
    px = reshape(repmat((1:1+N_coef), idsize, 1), [], 1);
    py = reshape(stats, [], 1);
    
    % p vals
    hbox = nan(1, size(stats, 2));
    pbox = nan(1, size(stats, 2));
    for t = 1:size(stats, 2)
        if t == 1
            [h, p] = ttest(stats(:, 1), 1/Nclass_y);
        else
            [h, p] = ttest(stats(:, t));
        end
        hbox(t) = h;
        pbox(t) = p;
    end
    
    subplot(1, yloop, yi)
%     beeswarm(px, py, 'corral_style', 'gutter', 'overlay_style', 'ci');
    beeswarm(px, py, 'corral_style', 'random', 'overlay_style', 'ci');
    if categ
        title({[num2str(ySet(yi)) ' vs others'], num2str(hbox)})
    else
        title({y_name, num2str(hbox)})
    end
%     xticklabels({''})
end

disp('end');


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

function [title_text] = anovan_group(idsize, except_ids, N_var1_global, ...
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

else
    title_text = [];

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
    case 'ChoCon1 G'
        var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
        GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
        var2 = var2(1,:);
        var2(GT~=1) = nan;
    case 'ChoCon2 G'
        var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
        GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
        var2 = var2(2,:);
        var2(GT~=1) = nan;
    case 'ChoCon ignoring context G'
        var2 = arbMBMF_load_var(Exp, 'ChoCon ignoring context', id, var2_cond);
        GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
        var2(:, GT~=1) = nan;
    case 'ChoCon1 per goal'
        var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
        var2 = var2(1,:);
    case 'ChoCon2 per goal'
        var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
        var2 = var2(2,:);
    case 'ChoCon1 per goal G'
        var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
        GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
        var2 = var2(1,:);
        var2(GT~=1) = nan;
    case 'ChoCon2 per goal G'
        var2 = arbMBMF_load_var(Exp, 'ChoConsist1n2', id, var2_cond);
        GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
        var2 = var2(2,:);
        var2(GT~=1) = nan;
    case 'StrictChoCon G'
        var2 = arbMBMF_load_var(Exp, 'StrictChoCon', id, var2_cond);
        GT = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
        var2(GT~=1) = nan;
    case 'ChoOpt1'

        switch variant
            case 'default'

                % standard approach (default parameters)
                var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(1,:);

            case 'OptimalA2SoftmaxA1blkCond'

                % special case for comparison with T_est inspection - OptimalA2, context-dependent softmax T
                load(['C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace\GoalSettingSoftmaxT\Kim2019 24\' variant '.mat'])
                blkCond = arbMBMF_load_var(Exp, 'blkCond', id, []);
                var2 = [blkCond; blkCond];
                for blk = 1:4
                    selector = ismember(blkCond, blk);
                    temp_var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, selector, ...
                        'Policy', 'optimal', ...
                        'SoftmaxT', paramset(id, blk), ...
                        'DoubleInComplex', false);
                    var2(:, selector) = temp_var2(:, selector);
                end
                var2 = var2(1,:);
        end

    case 'ChoOptBin1'
        var2 = arbMBMF_load_var(Exp, 'ChoOptBin1n2', id, var2_cond); var2 = var2(1,:);

    case 'ChoOptBinNew1'
        var2 = arbMBMF_load_var(Exp, 'ChoOptBinNew1n2', id, var2_cond); var2 = var2(1,:);

    case 'ChoOptTransEst1'

        switch variant
            case 'default'

                % standard approach (default parameters)
                var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2TransEst', id, var2_cond); var2 = var2(1,:);

            case 'OptimalA2SoftmaxA1TransEstblkCond'

                % special case for T_est inspection - OptimalA2, context-dependent softmax T
                load(['C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace\GoalSettingSoftmaxT\Kim2019 24\' variant '.mat'])
                blkCond = arbMBMF_load_var(Exp, 'blkCond', id, []);
                var2 = [blkCond; blkCond];
                for blk = 1:4
                    selector = ismember(blkCond, blk);
                    temp_var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2TransEst', id, selector, ...
                        'Policy', 'optimal', ...
                        'SoftmaxT', paramset(id, blk), ...
                        'DoubleInComplex', false);
                    var2(:, selector) = temp_var2(:, selector);
                end
                var2 = var2(1,:);
        end

    case 'ChoOpt2'

        switch variant
            case 'default'

                % standard approach (default parameters)
                var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(2,:);

            case 'OptimalA2SoftmaxA2blkCond'

                % special case for comparison with T_est inspection - OptimalA2, context-dependent softmax T
                load(['C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace\GoalSettingSoftmaxT\Kim2019 24\' variant '.mat'])
                blkCond = arbMBMF_load_var(Exp, 'blkCond', id, []);
                var2 = [blkCond; blkCond];
                for blk = 1:4
                    selector = ismember(blkCond, blk);
                    temp_var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, selector, ...
                        'Policy', 'optimal', ...
                        'SoftmaxT', paramset(id, blk), ...
                        'DoubleInComplex', false);
                    var2(:, selector) = temp_var2(:, selector);
                end
                var2 = var2(2,:);
        end

    case 'ChoOptBin2'
        var2 = arbMBMF_load_var(Exp, 'ChoOptBin1n2', id, var2_cond); var2 = var2(2,:);

    case 'ChoOptBinNew2'
        var2 = arbMBMF_load_var(Exp, 'ChoOptBinNew1n2', id, var2_cond); var2 = var2(2,:);
    case 'ChoOptTransEst2'
        switch variant
            case 'default'

                % standard approach (default parameters)
                var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2TransEst', id, var2_cond); var2 = var2(2,:);

            case 'OptimalA2SoftmaxA2TransEstblkCond'

                % special case for T_est inspection - OptimalA2, context-dependent softmax T
                load(['C:\Users\User\Desktop\Research\Dimensionality_control_for_prefrontal_meta_RL\workspace\GoalSettingSoftmaxT\Kim2019 24\' variant '.mat'])
                blkCond = arbMBMF_load_var(Exp, 'blkCond', id, []);
                var2 = [blkCond; blkCond];
                for blk = 1:4
                    selector = ismember(blkCond, blk);
                    temp_var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2TransEst', id, selector, ...
                        'Policy', 'optimal', ...
                        'SoftmaxT', paramset(id, blk), ...
                        'DoubleInComplex', false);
                    var2(:, selector) = temp_var2(:, selector);
                end
                var2 = var2(2,:);
        end

    case 'ChoOpt_G1'
        var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(1,:);
        gt = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
        var2(gt~=1) = nan;
    case 'ChoOpt_G2'
        var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond); var2 = var2(2,:);
        gt = arbMBMF_load_var(Exp, 'GoalCond', id, var2_cond);
        var2(gt~=1) = nan;
    case 'ChoOpt_cL1'
        var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond, 'DoubleInComplex', false);
        var2 = var2(1,:);
        cx = arbMBMF_load_var(Exp, 'CmplxCond', id, var2_cond);
        var2(cx~=1) = nan;
    case 'ChoOpt_cL2'
        var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond, 'DoubleInComplex', false);
        var2 = var2(2,:);
        cx = arbMBMF_load_var(Exp, 'CmplxCond', id, var2_cond);
        var2(cx~=1) = nan;
    case 'ChoOpt_cH1'
        var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond, 'DoubleInComplex', false);
        var2 = var2(1,:);
        cx = arbMBMF_load_var(Exp, 'CmplxCond', id, var2_cond);
        var2(cx~=2) = nan;
    case 'ChoOpt_cH2'
        var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2', id, var2_cond, 'DoubleInComplex', false);
        var2 = var2(2,:);
        cx = arbMBMF_load_var(Exp, 'CmplxCond', id, var2_cond);
        var2(cx~=2) = nan;
    case 'ChoOpt1AchGoalOpt'
        var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2AchGoalOpt', id, var2_cond); var2 = var2(1,:);
    case 'ChoOpt1AchGoalRan'
        var2 = arbMBMF_load_var(Exp, 'ChoOpt1n2AchGoalRan', id, var2_cond); var2 = var2(1,:);
    case 'absRPE'
        var2 = arbMBMF_load_var(Exp, 'RPE', id, var2_cond);
        var2 = abs(var2);
    case 'absRPE3'
        var2 = arbMBMF_load_var(Exp, 'RPE3', id, var2_cond);
        var2 = abs(var2);
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
