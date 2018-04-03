function log_p = logp_logitBeta(z,alpha)
log_p = alpha(:,1).*z(:) - sum(alpha,2).*log(1 + exp(z(:))) - (sum(gammaln(alpha),2) - gammaln(sum(alpha,2)));