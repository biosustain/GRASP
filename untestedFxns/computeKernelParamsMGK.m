function ensemble = computeKernelParamsMGK(ensemble,popIdx,prevParams)
%--------------------------------------------------------------------------
% Compute kernel parameters 
%
% Inputs:       ensemble (structure) , popIdx (double)
%
% Outputs:      ensemble kernel parameters updated (structure)
%--------------------- Pedro Saa 2016 -------------------------------------
% Loop through all the structures and rxns in the model
for strucIdx = 1:ensemble.numStruct
	
	% Check whether the structure has any surving particles
	if isempty(prevParams)
		continue;
	
	% If so, then compute sigma as twice the weighted cov matrix of the previous population
	else
	    
		% Extract alive particles from the current model structure
		idxParticleStructure = find((ensemble.populations(popIdx-1).strucIdx==strucIdx)&(ensemble.populations(popIdx-1).weights~=0));
		weights              = robustWeightNorm(ensemble.populations(popIdx-1).weights(idxParticleStructure));											% Robust normalize weigths
		
		% Compute kernel parameters
		paramArray = prevParams(strucIdx).paramArray';
		muX        = weights'*paramArray;
		diffX      = paramArray-muX(ones(size(paramArray,1),1),:);
		diffW      = weights(:,ones(1,size(paramArray,2))).*diffX;
		covX       = diffW'*diffX/(1-sum(weights.^2));
		
		% Save kernel parameters
		ensemble.populations(popIdx).probParams(strucIdx).globalKernelParam = 2*covX; 	% twice the weighted sample covariance
	end
end
disp('Kernel parameters updated.');