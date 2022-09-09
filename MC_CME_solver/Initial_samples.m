function [X]=Initial_samples(IniDistribution_Mean,IniDistribution_Variance,...
                             sample_size)

 num_species=size(IniDistribution_Mean,1);
 
 for i=1:num_species
     X(i,:)= IniDistribution_Mean(i)-IniDistribution_Variance(i)...
             +poissrnd(IniDistribution_Variance(i),1,sample_size);
 end