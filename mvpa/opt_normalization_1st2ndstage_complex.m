function [outputArg1,outputArg2] = opt_normalization_1st2ndstage_complex(second_state, action, token_value, ind, condi, varargin)
% modified by ydSung (22-01-11)

% disp(varargin)
% disp(length(varargin))

% default input parameters (from the 2019 original settings)
options = struct('Policy', 'random', ...
                'SoftmaxT', [], ...
                'DoubleInComplex', true);
% varargin validity
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(opt_normalization_1st2ndstage_complex) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(opt_normalization_1st2ndstage_complex) %s is not a recognized parameter name', pair{1})
    end
end

% for action 1 optimality computation
switch options.Policy
    case 'random'
        policy = @mean; % random policy for action 2, 2019 original
    case 'optimal'
        policy = @max; % optimal policy for action 2
end

% softmax temperature
if ~isempty(options.SoftmaxT)
    if length(options.SoftmaxT) < 2
        T1 = options.SoftmaxT;
        T2 = options.SoftmaxT;
    else
        T1 = options.SoftmaxT(1);
        T2 = options.SoftmaxT(2);
    end
end

expect_value = []; % expected value for all actions, regarding token value. L1, R1, L2, R2 respectively
if second_state == 1
    % ************ action 1 : Right // 2 : Left ************
    switch ind
        case 1
            state2 = policy([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2)]);
            state3 = policy([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3)]);
            expect_value = [0.9*state2 + 0.1*state3, 0.9*state3 + 0.1*state2];
        case 2
            state2 = policy([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2), 0.9*token_value(3) + 0.1*0, 0.9*token_value(2) + 0.1*token_value(1)]);
            state3 = policy([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3), 0.9*0 + 0.1*token_value(1), 0.9*token_value(3) + 0.1*token_value(2)]);
            expect_value = [0.9*state2 + 0.1*state3, 0.9*state3 + 0.1*state2];
        case 3
            state2 = policy([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2)]);
            state3 = policy([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3)]);
            expect_value = [0.5*state2 + 0.5*state3, 0.5*state3 + 0.5*state2];
        case 4
            state2 = policy([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2), 0.5*token_value(3) + 0.5*0, 0.5*token_value(2) + 0.5*token_value(1)]);
            state3 = policy([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3), 0.5*0 + 0.5*token_value(1), 0.5*token_value(3) + 0.5*token_value(2)]);
            expect_value = [0.5*state2 + 0.5*state3, 0.5*state3 + 0.5*state2];    
    end
    
    % choice likelihood: 1) normalized value, 2) softmax(value)
    if isempty(options.SoftmaxT)
        likelihood = expect_value ./ sum(expect_value); % normalized value (2019 original)
    else
        likelihood = exp(expect_value/T1) ./ sum(exp(expect_value/T1));
    end
    
    if (sum(expect_value) == 0)
%         error('(opt_normalization_1st2ndstage_complex) zero expect_value')
        likelihood = ones(size(expect_value));
        likelihood = likelihood/sum(likelihood);
    end
    
    if strcmp (condi, 'likelihood')
        outputArg1 = likelihood(action);
        outputArg2 = find(likelihood == max(likelihood));
        
    elseif strcmp (condi, 'binary')
        if sum(likelihood == max(likelihood)) == length(likelihood)
            outputArg1 = 1/length(likelihood);
        else
            outputArg1 = ismember(action, find(likelihood == max(likelihood))) / sum(likelihood == max(likelihood));
        end
        
    elseif strcmp(condi, 'leftbias')
        outputArg1 = likelihood(2) - likelihood(1);
        
    elseif strcmp(condi, 'optchoice')
        M = containers.Map([1 2], [-4 4]); % left bias colormap range
        max_vec = find(likelihood == max(likelihood));
        outputArg1 = nan(size(max_vec));
        N_max = sum(likelihood == max(likelihood));
        for maxi = 1:N_max
            outputArg1(maxi) = M(max_vec(maxi));
        end
        outputArg1 = mean(outputArg1);    
        
    end
    
    
else
    switch ind
        case 1 % LL
            %                 opt(1,2) = opt(1,2)+1;
            if second_state == 2 % a1 should be L1. L1 : NOTHING, L2 : BLUE, R2 : RED, R1 : SILVER
                expect_value = [0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2)];
            else % second_state == 3 L1 : SILVER, L2 : NOTHING, R2 : BLUE, R1 : RED
                expect_value = [0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3)];
            end
            
        case 2 %HL
            %                 opt(2,2) = opt(2,2)+1;
            if second_state == 2
                expect_value = [0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2), 0.9*token_value(3) + 0.1*0, 0.9*token_value(2) + 0.1*token_value(1)];
            else
                expect_value = [0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3), 0.9*0 + 0.1*token_value(1), 0.9*token_value(3) + 0.1*token_value(2)];
            end
            
        case 3 % LH
            %                 opt(3,2) = opt(3,2)+1;
            if second_state == 2 % a1 should be L1. L1 : NOTHING, L2 : BLUE, R2 : RED, R1 : SILVER
                expect_value = [0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2)];
            else % second_state == 3 L1 : SILVER, L2 : NOTHING, R2 : BLUE, R1 : RED
                expect_value = [0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3)];
            end
            
        case 4 %HH
            %                 opt(4,2) = opt(4,2)+1;
            if second_state == 2
                expect_value = [0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2), 0.5*token_value(3) + 0.5*0, 0.5*token_value(2) + 0.5*token_value(1)];
                
                %             I = find(expect_value == max(expect_value));
                %             if ~isempty(find(sbj_data(4) == I))
                %                 cc_HH(sub,2) = cc_HH(sub,2) + 1;
                %                 if isequal(I,[2,4])
                %                     cc_HH(sub,1) = cc_HH(sub,1) + 1;
                %                 end
                %             end
            else
                
                expect_value = [0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3), 0.5*0 + 0.5*token_value(1), 0.5*token_value(3) + 0.5*token_value(2)];
            end
            
    end
    
    % choice likelihood: 1) normalized value, 2) softmax(value)
    if isempty(options.SoftmaxT)
        likelihood = expect_value ./ sum(expect_value); % normalized value (2019 original)
    else
        likelihood = exp(expect_value/T2) ./ sum(exp(expect_value/T2));
    end
    
    if (sum(expect_value) == 0)
%         error('(opt_normalization_1st2ndstage_complex) zero expect_value')
        likelihood = ones(size(expect_value));
        likelihood = likelihood/sum(likelihood);
    end
        
    if strcmp (condi, 'likelihood')
        if (ind == 2) || (ind == 4)
            if options.DoubleInComplex
                % doubling the likelihood in the complexity high condition
                likelihood = likelihood*2;
            end
        end
        outputArg1 = likelihood(action);
        outputArg2 = find(likelihood == max(likelihood));
    elseif strcmp (condi, 'uncertainty')
        expect_expect = [];
        
        switch ind
            case 1
                expect_expect = [0.9*mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2)]) + 0.1*mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3)]),...
                    0.1*mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2)]) + 0.9*mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3)])];
            case 2
                expect_expect = [0.9*mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2), 0.9*token_value(3) + 0.1*0, 0.9*token_value(2) + 0.1*token_value(1)]) + 0.1*mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3), 0.9*0 + 0.1*token_value(1), 0.9*token_value(3) + 0.1*token_value(2)]),...
                    0.1*mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2), 0.9*token_value(3) + 0.1*0, 0.9*token_value(2) + 0.1*token_value(1)]) + 0.9*mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3), 0.9*0 + 0.1*token_value(1), 0.9*token_value(3) + 0.1*token_value(2)])];
            case 3
                expect_expect = [0.5*mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2)]) + 0.5*mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3)]),...
                    0.5*mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2)]) + 0.5*mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3)])];
            case 4
                expect_expect = [0.5*mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2), 0.5*token_value(3) + 0.5*0, 0.5*token_value(2) + 0.5*token_value(1)]) + 0.5*mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3), 0.5*0 + 0.5*token_value(1), 0.5*token_value(3) + 0.5*token_value(2)]),...
                    0.5*mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2), 0.5*token_value(3) + 0.5*0, 0.5*token_value(2) + 0.5*token_value(1)]) + 0.5*mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3), 0.5*0 + 0.5*token_value(1), 0.5*token_value(3) + 0.5*token_value(2)])];
        end
        likelihood = expect_expect ./ sum(expect_expect);
        outputArg1 = likelihood(action);
        
    elseif strcmp (condi, 'binary')
        if sum(likelihood == max(likelihood)) == length(likelihood)
            outputArg1 = 1/length(likelihood);
        else
            outputArg1 = ismember(action, find(likelihood == max(likelihood))) / sum(likelihood == max(likelihood));
        end        
        
    elseif strcmp(condi, 'leftbias')
        if length(likelihood) == 2
            outputArg1 = likelihood(1) - likelihood(2);
        elseif length(likelihood) == 4
            outputArg1 = sum(likelihood([1 3])) - sum(likelihood([2 4]));
        end
        
    elseif strcmp(condi, 'optchoice')
        M = containers.Map([1 2 3 4], [4 -4 2 -2]); % left bias
        max_vec = find(likelihood == max(likelihood));
        outputArg1 = nan(size(max_vec));
        N_max = sum(likelihood == max(likelihood));
        for maxi = 1:N_max
            outputArg1(maxi) = M(max_vec(maxi));
        end
        outputArg1 = mean(outputArg1);        
    end
end

end

