function [subj] = load_spm_mask_2(subj,new_maskname,filename,varargin)

% Loads an NIFTI dataset into the subj structure as a mask
%
% [SUBJ] = LOAD_ANALYZE_MASK(SUBJ,NEW_MASKNAME,FILENAME,...)
%
% Adds the following objects:
% - mask object called NEW_MASKNAME
%
% License:
%=====================================================================
%
% This is part of the Princeton MVPA toolbox, released under
% the GPL. See http://www.csbmb.princeton.edu/mvpa for more
% information.
% 
% The Princeton MVPA toolbox is available free and
% unsupported to those who might find it useful. We do not
% take any responsibility whatsoever for any problems that
% you have related to the use of the MVPA toolbox.
%
% ======================================================================

% default
args = struct('mask_num', 0, ...
              'hemisphere', []);
args_names = fieldnames(args);
if mod(length(varargin), 2) == 1
    error('(load_spm_mask_2) Unmatched name-value pair')
end
for pair = reshape(varargin, 2, [])
    if any(strcmp(pair{1}, args_names))
        args.(pair{1}) = pair{2};
    else
        error('(load_spm_mask_2) %s is not a recognized parameter name', pair{1})
    end
end

% Initialize the new mask
subj = init_object(subj,'mask',new_maskname);

% Create a volume
vol = spm_vol(filename); % struct with fields: fname, dim, dt, ...

V = spm_read_vols(vol); % [79 x 95 x 79 double]

% Check for active voxels
if ~~isempty(find(V))
  error( sprintf('There were no voxels active in the mask') );
end

V(find(isnan(V))) = 0;

% Does this consist of solely ones and zeros?W
if length(find(V)) ~= (length(find(V==0))+length(find(V==1)))

    if args.mask_num == 0 
        % No specific restriction of mask number (ID)
        V(find(V ~= 0)) = 1; % All nonzero values make a valid mask
    else
        % using the specific mask number value(s) for a valid mask
        V(~ismember(V, args.mask_num)) = 0;
        V(ismember(V, args.mask_num)) = 1;
    end
    
    if ~isempty(args.hemisphere)
        [xdim, ydim, zdim] = size(V);
        if vol.mat(1,1) < 0
            changed_LR = setdiff({'left', 'right'}, args.hemisphere);
            args.hemisphere = [changed_LR{:}];
        end
        switch args.hemisphere
            case 'left'
                % nullify right idx
                V(ceil(xdim/2):xdim,:,:) = 0;
            case 'right'
                % nullify left idx
                V(1:floor(xdim/2),:,:) = 0;
        end
    end
    
end

% Store the data in the new mask structure
subj = set_mat(subj,'mask',new_maskname,V);

% Add the AFNI header to the patterns
hist_str = sprintf('Mask ''%s'' created by load_analyze_pattern',new_maskname);
subj = add_history(subj,'mask',new_maskname,hist_str,true);

% Add information to the new mask's header, for future reference
subj = set_objsubfield(subj,'mask',new_maskname,'header', ...
			 'vol',vol,'ignore_absence',true);

% Record how this mask was created
created.function = 'load_analyze_mask';
subj = add_created(subj,'mask',new_maskname,created);