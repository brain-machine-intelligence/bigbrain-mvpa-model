function [dim_msr, pat_sizes, N_nonzero_eigs] = ana_dim(boldpat, varargin)
% unsupervised dimensionality of multi-voxel pattern 

% input: boldpat [N_vox x N_observation]

% output: dimensionality measure (participation ratio / number of top-n-PCs)
% if 'Labels' is empty -> a scalar number
% if 'Labels' = n-class vector -> [n x 1] vector

% =========================================================================

% default input parameters
options = struct('Labels', [], ...
                'DimMethod', 'both', ... 'PR' / 'nPC' / 'both' / 'entropy' / 'ER'
                'ExpVarPercent', 90, ...
                'Normalization', 'none', ... 'none' / 'basic' / 'voxels' / 'trials' / 'undersample'
                'UndersampleNumber', 1, ...
                'RandomSeed', 2022, ...
                'ShuffleRepeat', [1 0]); % [iterator shuffle_number], if shuffle_number = 0 then do nothing
% read input parameters
option_names = fieldnames(options);
if mod(length(varargin),2) == 1
    error('(arbMBMF_boldpat_dimension) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, option_names))
        options.(pair{1}) = pair{2};
    else
        error('(arbMBMF_boldpat_dimension) %s is not a recognized parameter name', pair{1})
    end
end

% Random seed
if options.ShuffleRepeat(2)
    RandSeed = options.RandomSeed + 10000 * mod(options.ShuffleRepeat(1), options.ShuffleRepeat(2));
%     fprintf(' t1') % test
else
    RandSeed = options.RandomSeed;
%     fprintf(' t2') % test
end
rng(RandSeed);

% parameter assignment
label = options.Labels;
method = options.DimMethod;
expvar = options.ExpVarPercent;
normalization = options.Normalization;

%% Dimensionality measure

if isempty(label)
    label = ones(1, size(boldpat, 2));
end

if strcmp(normalization, 'undersampleWithoutNaN')
    labelboldpat = [label; boldpat];
    labelboldpat = labelboldpat(:, ~any(isnan(boldpat), 1));
    label = labelboldpat(1, :);
    boldpat = labelboldpat(2:end, :);
end

[lab_class, class_lengths] = unq_elms(label); N_min_class = min(class_lengths);
N_class = length(lab_class);
N_vox = size(boldpat, 1);

N_nonzero_eigs = nan(N_class, 1);
pat_sizes = nan(N_class, 1);

if contains(method, 'both')
    dim_msr = zeros(N_class, 2);
elseif contains(method, 'all')
    dim_msr = zeros(N_class, 3);
else
    dim_msr = zeros(N_class, 1);
end

if isempty(options.UndersampleNumber)
    options.UndersampleNumber = 1;
end

for ui = 1:options.UndersampleNumber
    
    dim_msr_temp = nan(size(dim_msr));
    
    for c = 1:N_class
        
        temp_pat = boldpat(:, label==lab_class(c));
        N_trial_c = size(temp_pat, 2);
        
        if strcmp(normalization, 'undersample') || strcmp(normalization, 'undersampleWithoutNaN')
            rng(RandSeed + ui);
%             fprintf(' %d', RandSeed + ui) % test
            temp_pat = temp_pat(:, randperm(N_trial_c, N_min_class));
        end
        rng(RandSeed);
        
        pat_sizes(c) = size(temp_pat, 2);

        if contains(method, 'disp')

            % Mean pattern of class c
            mu_c = mean(temp_pat, 2);

            % Distances of each sample from mu_c
            distVec = sqrt( sum( (temp_pat - mu_c).^2, 1 ) );
            % distVec is [1 x #samples_in_class]

            % Average distance
            avgDist = mean(distVec);

            dim_msr_temp(c, 1) = avgDist;

        else

            %     [~, ~, latent] = pca(temp_pat');
            [~, ~, latent, ~, explained] = pca(temp_pat'); explained = explained / 100;
            %     [~, ~, ~, ~, explained] = pca(temp_pat', 'Economy', false);

            %     fprintf('disp(size(temp_pat)) ')
            %     disp(size(temp_pat))
            %     fprintf('disp(size(explained)) ')
            %     disp(size(explained))

            %     fprintf('disp(sum(explained~=0)): %d', sum(explained~=0))
            %         fprintf('%d ', sum(explained~=0))
            N_nonzero_eigs(c) = sum(explained~=0);

            if contains(method, 'PR')
                % participation ratio
                dim_msr_temp(c, 1) = sum(explained)^2/sum(explained.^2);
            elseif contains(method, 'nPC')
                % number of top PCs
                dim_msr_temp(c, 1) = sum( cumsum(explained)/sum(explained) <= expvar/100 );
            elseif contains(method, 'entropy')
                dim_msr_temp(c, 1) = explained' * -log(explained);
            elseif contains(method, 'ER')
                dim_msr_temp(c, 1) = exp(explained' * -log(explained));
            elseif contains(method, 'MV') % mean variance
                dim_msr_temp(c, 1) = mean(latent);
            elseif contains(method, 'both')
                dim_msr_temp(c, 1) = sum(explained)^2/sum(explained.^2);
                dim_msr_temp(c, 2) = sum( cumsum(explained)/sum(explained) <= expvar/100 );
            elseif contains(method, 'all')
                dim_msr_temp(c, 1) = sum(explained)^2/sum(explained.^2);
                %         dim_msr(c, 2) = - explained' * log(explained);
                dim_msr_temp(c, 2) = exp(- explained' * log(explained));
                dim_msr_temp(c, 3) = sum( cumsum(explained)/sum(explained) <= expvar/100 );
            end

        end
        
        if strcmp(normalization, 'none')
            % do nothing
        elseif strcmp(normalization, 'undersample') || strcmp(normalization, 'undersampleWithoutNaN')
            % do nothing
        elseif strcmp(normalization, 'basic')
            %         N_pc_nonzero = sum(explained~=0);
            %         normalizer = min([N_vox, N_trial_c, N_pc_nonzero]);
            normalizer = min([N_vox, size(temp_pat, 2)]);
            dim_msr_temp(c, :) = dim_msr_temp(c, :) / normalizer;
        elseif strcmp(normalization, 'voxels')
            dim_msr_temp(c, :) = dim_msr_temp(c, :) / N_vox;
        elseif strcmp(normalization, 'trials')
            dim_msr_temp(c, :) = dim_msr_temp(c, :) / N_trial_c;
        else
            error('(arbMBMF_boldpat_dimension) wrong Name-Value pair for Normalization')
        end
        
    end
    
%     disp(reshape(dim_msr_temp, 1, [])) % temp test
    dim_msr = (dim_msr * (ui-1) + dim_msr_temp) / ui;
end


