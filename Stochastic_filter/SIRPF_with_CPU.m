function pi_PF=SIRPF_with_CPU(TY,Y,k,phi,reactant_vector,product_vector,...
                                IniDistribution_Mean,...
                                IniDistribution_Variance,...
                                h,sigma,sample_size)
                            
% This function returns the solution of the particle filter for SRNs

% Preparation
X=Initial_samples(IniDistribution_Mean,IniDistribution_Variance,...
                  sample_size);
W=ones(1,sample_size)/sample_size;
pi_PF=cell(2,max(size(TY))); 


% Iteration 
for i=1:max(size(TY))
    % determine dt
    if i==1
        dt=TY(1);
    else
        dt=TY(i)-TY(i-1);
    end
    
    parfor j=1:sample_size
        % simulate particles
        [t_temp,X_temp]=next_reaction_method_single_trajectory...
               (X(:,j),k,phi,dt,reactant_vector,product_vector);
        X(:,j)=X_temp(:,end);
        %update weights
        W(j)=W(j)*likelihood_function(Y(:,i),X(:,j),h,sigma);
    end
    
    % save the result 
    pi_PF{1,i}=X;
    pi_PF{2,i}=W/sum(W);
    
    %resample
    [X,W]=resampling_step(X,W); 
end