function d=Initial_distribution(IniDistribution_Mean,...
                                IniDistribution_Variance,...
                                truncated_size)
                            
% This function generates initial distribution for the FSP method. 

% generate d
if size(truncated_size,1)==1 %only one species
    d=zeros(truncated_size,1); 
else
    d=zeros(truncated_size');
end

% prepare the marginal distribution
for i=1:size(truncated_size,1)
    Marginal{i}=zeros(truncated_size(i),1);
    for j=1:truncated_size(i)
        x=double(j-1); % the value of the state
        if x>=IniDistribution_Mean(i)-IniDistribution_Variance(i)
           x=x-IniDistribution_Mean(i)+IniDistribution_Variance(i);
           x=double(x);
           Marginal{i}(j)= poisspdf(x,IniDistribution_Variance(i));
        end
    end
end



% compute the initial distribution
num_species=size(truncated_size,1);
for LI=1:prod(truncated_size)  %LI is the linear index of the matrix
    MI=Linear_index_2_matrix_index(truncated_size,LI);
    value=1;
    for i=1:num_species
        value=value*Marginal{i}(MI(i));
    end
    d(LI)=value;
end