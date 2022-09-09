function l=likelihood_function(Y,X,h,sigma)

% this algorithm returns the likelihood function L(y|x)

a=Y-h_function(h,X);
exponent=-a.*a/2./(sigma.*sigma)';
p=1/sqrt(2*pi)./sigma'.*exp(exponent);  % the Likelihood for each channel.
l=prod(p);