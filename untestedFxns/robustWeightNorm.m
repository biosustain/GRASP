function w = robustWeightNorm(w)
w = log(w)-max(log(w));
w = exp(w)/sum(exp(w));