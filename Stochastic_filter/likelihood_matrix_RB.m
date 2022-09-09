function l=likelihood_matrix_RB(Y,X,h,sigma)

a=Y-h_function_vector_form(h,X);
exponent=-a.*a/2./(sigma.*sigma)';
l=1/sqrt(2*pi)./sigma'.*exp(exponent);  % the Likelihood for each channel.