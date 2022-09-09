function Likelihood=likelihood_matrix(Y,h,sigma,truncated_size)

Likelihood=zeros(truncated_size');
Likelihood=reshape(Likelihood,[],1); % transform it into a column vector;

size_of_the_matrix=prod(truncated_size);

parfor LI=1:size_of_the_matrix
    MI=Linear_index_2_matrix_index(truncated_size,LI);
    if size(MI,1)==1
       MI=MI';    % it should be column vector 
    end
    X=double(MI-1);
    Likelihood(LI)=likelihood_function(Y,X,h,sigma);
end

if size(truncated_size,1)>1
   Likelihood=reshape(Likelihood,truncated_size');
end