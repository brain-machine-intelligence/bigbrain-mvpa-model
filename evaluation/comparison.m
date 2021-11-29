clear;

LIST_SBJ={'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
    'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};
% rejecting abnormal(different session length)
list_sbj={LIST_SBJ{2:15} LIST_SBJ{17:24}};
list_sbj={list_sbj{1:14} list_sbj{16:end}};
BS = cell(1, length(list_sbj));
for subind = 1 : length(list_sbj)
    matname = [pwd '\result_save\' 'SBJ_' list_sbj{subind}  '*DPON*.mat'];
    matlist = dir(matname);
    setmodel = cell(length(matlist),3);
    for i = 1 : length(matlist)
        load([pwd '\result_save\' matlist(i).name]);
        resb = length(SBJtot);
        RESULT = cell(resb,1);
        RESULT_SB = zeros(resb,1);
        for ii = 1 : resb
            SBSBJ = SBJtot{ii,1};
            RESULT_SB(ii) = SBSBJ{1,1}.model_BayesArb.val;
            RESULT{ii,1} = SBSBJ{1,1}.model_BayesArb.param;
        end
        MINMIN = min(RESULT_SB);
        findindmin = find(RESULT_SB == MINMIN);
        OPTIMPARAM = RESULT{findindmin(1),1};
        setmodel{i,1} = SBJtot{findindmin(1),1}{1, 1}.model_BayesArb.val;
        setmodel{i,2} = SBJtot{findindmin(1),1}{1, 1}.model_BayesArb.param;
        setmodel{i,3} = SBJtot{findindmin(1),1};
    end
    BS{1,subind} = setmodel;
end

resbs = length(BS);
valset = cell(resbs,1);
paramset = cell(resbs,1);
minRef = 999;
SBJset = cell(resbs,1);
for i = 1 : resbs
    ind = size(BS{1,i},1);
    valset{i,1} = minRef;
    for in = 1 : ind
        mintemp = BS{1,i}{in,1};
        if mintemp < valset{i,1}
            valset{i,1} = mintemp;
            indexs = in;
        end
    end
    paramset{i,1} = BS{1,i}{indexs,2};
    SBJset = BS{1,i}{indexs,3};
end
paramset = cell2mat(paramset);

%% Model comparison
% load for both comparison
load('SBJ_structure_each_exp_BESTcollectionJan29.mat')
SBJ_ori = {SBJ{1:14} SBJ{16:22}};
valset2 = cell(length(SBJ_ori),1);
numK = zeros(length(SBJ_ori),1);
for subind = 1 : length(list_sbj)
    valset2{subind} = SBJ_ori{1, subind}.model_BayesArb.val;
    numK(subind) = SBJ_ori{1, subind}.num_data;
end

% value to loglikelihood sum
DPval = -cell2mat(valset);
Orival = -cell2mat(valset2);

% 1. loglikelihood sum spm_BMS test
[alpha,exp_r,xp,pxp,bor] = spm_BMS([Orival DPval]);
[h_ll, p_ll] = ttest(Orival, DPval, 'Alpha', 0.01, 'Tail', 'left'); % This means that we believe DPval is greater than Orival.

% 2. BIC based spm_BMS test
[h_bic,p_bic] = ttest(-Orival* 2 + 6 * log(numK) , -DPval* 2 + 4 * log(numK), 'Alpha', 0.01, 'Tail', 'right');

