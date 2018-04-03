function f = logsumexp(x)
the_max = max(x);
x       = x - the_max;
f       = the_max + log(sum(exp(x)));