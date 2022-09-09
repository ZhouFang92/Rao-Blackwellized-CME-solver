function [TX,X]=next_reaction_method_multi_trajectories_no_GPU...
                (k,phi,FT,reactant_vector,product_vector,...
                 IniDistribution_Mean,IniDistribution_Variance,...
                 sample_size)
             
             
% This function uses the next reaction method to generate
% multi_trajectories for the network. 
% GPUs are not used


% Initialization
X=Initial_samples(IniDistribution_Mean,IniDistribution_Variance,...
                  sample_size);

%simulation
parfor i=1:sample_size
    [t_temp,X_temp]=next_reaction_method_single_trajectory...
               (X(:,i),k,phi,FT,reactant_vector,product_vector);
    TX(i)=t_temp(end);
    X(:,i)=X_temp(:,end);
end