% Description : Bayesian-IGMM artibtration model that use one-by-one
% approach.
%% valid subject
function [] = BATCH_ORI_FULL(sis,MAXD,maxcore)

    t = clock;
    LIST_SBJ={'Oliver', 'Hao', 'Breanna', 'Derek', 'Timothy', 'Teagan', 'Jeffrey', 'Seung', 'Carole', 'Tony', 'Surendra', 'Lark',...
        'Joaquin', 'DavidB', 'Christopher', 'Gjergji', 'Charles', 'Erin', 'Connor', 'Domenick', 'Thao', 'Arin', 'Pauline', 'Tho'};
    % rejecting abnormal(different session length)
list_sbj={LIST_SBJ{2:15} LIST_SBJ{17:24}};
    list_sbj= { list_sbj{sis:sis} };
%     parpool('local', maxcore);


    %% param_in
    % original parameter description
    % Description : param1 = zero prediction error threshold of SPE
    %               param2 = learning rate for the estimate of absolute reward
    %                           prediction error
    %               param3 = the amplitude of a transition rate function(mb-mf
    %               param4 = the amplitude fo a transition rate function(mf-mb
    %               param5 = inverse softmax temparature
    %               param6 = learning rate of the model based and the model
    %                      free (both of them use same value of learning rate)

    % DPNMM parameter description
    % We don't need any of DPNMM parameter because we just want the optimized
    % cluster model. we don't need to be fast.
    % Description : param1 = param 3 of original model
    %                           the amplitude of a transition rate function(mb-mf
    %               param2 = param 4 of original model
    %                           the amplitude of a transition rate function(mf-mb
    %               param3 = param 5 of original model
    %                           inverse softmax temparature
    %               param4 = param 6 of original model
    %                           learning rate of the model based and the model

    mode.param_BoundL = [0.3 0.01  0.1 0.1 0.01 0.01];
    mode.param_BoundU = [0.8 0.35 10 10 0.5 0.2];
    param_init= zeros(1,6);
    for i = 1 : 6
        rand(floor(mod(sum(clock*10),10000)));
        param_init(i) = rand  * (mode.param_BoundU(i) - mode.param_BoundL(i))  + mode.param_BoundL(i);
    end



    %% mode
    mode.param_length = 6 ;%size(param_init,2);
    mode.total_simul = 10; %30
    mode.experience_sbj_events=[1 1]; % ones(1,2); % experience_sbj_events(2) = 1 일때 subject의 state transition을 이용.
    mode.USE_FWDSARSA_ONLY=0; % 이게 1이면 fwd only고 이게 2이면 sarsa only. 둘다 아니게 하려고 이런 정보를 넣었지만 아예 그 if문을 둘다 빼버리는 방법도 있음.
    mode.USE_BWDupdate_of_FWDmodel=1; % BWD update 사용하는 방식으로 학습하기로.
    mode.DEBUG_Q_VALUE_CHG=0;
    mode.simul_process_display=0;
    mode.out=1;
    mode.opt_ArbModel = 1; % % 1: naive model(m1_wgt) . 2: posterior model(posterior)
    mode.boundary_12 = 0.1; 
    mode.boundary_21 = 0.01;
    mode.max_iter=200;
    maxd=MAXD;


    %% data_in
    % Load from list_sbj.
    maxi=size(list_sbj,2);% because it is for 1 subject. maxi=size(list_sbj,2);
    data_in_tag{1,1}='data in의 열들은 subject 넘버. 그리고 그 내부의 HIST behavior info내부는 각각 세션을 의미';
    % outputval={};

    % storing
    PreBehav = {};
    PreBlck = {};
    MainBehav = {};
    MainBlck = {};
    SBJ = cell(1, maxi);

    % in the dirrectory?


    for i = 1 : maxi
        TEMP_PRE=load([pwd '/result_save/' list_sbj{i} '_pre_1.mat']);
        PreBehav{i}=TEMP_PRE.HIST_behavior_info{1,1};
        PreBlck{i}=TEMP_PRE.HIST_block_condition{1,1};


        tt = dir([pwd '/result_save']);
        tt = {tt.name};
        maxsess = sum(cell2mat(strfind(tt,[list_sbj{i} '_fmri_']))) - 1;

        try
            for ii = 1 : maxsess
                TEMP_MAIN=load([pwd '/result_save/' list_sbj{i} '_fmri_' num2str(ii) '.mat']);
                MainBehav{i,ii}=TEMP_MAIN.HIST_behavior_info0;
                MainBlck{i,ii}=TEMP_MAIN.HIST_block_condition{1,ii};
            end
        catch
            warning('MAX session error');
        end

    end
    save([pwd '/result_save/F_DATA.mat'], 'PreBehav', 'PreBlck','MainBehav', 'MainBlck');
    disp('DATA STORING DONE!');

    SBJtot = cell(maxd,1);
    for d = 1 : maxd
        %% PE, DPNMM

        [SBJ]=ArbBat_Ori(maxi,mode,PreBehav,PreBlck,MainBehav,MainBlck,list_sbj,param_init)
        SBJtot{d,1} = SBJ;
    end
    list_month={'Jan','Feb','Mar','Apr', 'May', 'Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
    save([pwd '/result_save/SBJ_' list_sbj{1} '_DPON_' num2str(t(1)) list_month{t(2)} num2str(t(3)) '_' num2str(t(4)) num2str(t(5)) num2str(floor(t(6))) '.mat'], 'SBJtot');
end

