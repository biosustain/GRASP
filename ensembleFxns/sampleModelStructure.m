function y = sampleModelStructure(idx,ksamples,weights)
%--------------------------------------------------------------------------
% Weighted sample of a model structure with replacement
%
% Inputs:       ensemble (structure) , workerIdx (double)
%
% Outputs:      initial sampled ensemble structure
%--------------------- Pedro Saa 2016 -------------------------------------
dim   = find(size(idx)~=1,1);
sumw  = sum(weights);
p     = weights(:)'/sumw;
edges = min([0,cumsum(p)],1);            % protect against accumulated round-off
edges(end) = 1;                          % get the upper edge exact
[~, i] = histc(rand(1,ksamples),edges);

% Use the index vector to sample from the data.
if ismatrix(idx)
    if dim == 1
        y = idx(i,:);
    elseif dim == 2
        y = idx(:,i);
    else
        reps = [ones(1,dim-1) ksamples];
        y = repmat(idx,reps);
    end
else % N-D
    subs = repmat({':'},1,max(ndims(idx),dim));
    subs{dim} = i;
    y = idx(subs{:});
end