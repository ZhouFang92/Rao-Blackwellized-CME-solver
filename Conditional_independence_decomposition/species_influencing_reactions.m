function v=species_influencing_reactions(num_species,propensity)

% This function return a cell array the j-th of which shows those species
% influencing the j-th reaction.
v{1}=[];
for j=1:size(propensity,2)
    v{j}=[];
    for i=1:num_species
        substring=strcat("X(",string(i),")");
        if contains(propensity(j),substring)==1 % check if the i-th species 
                                             % influence the reaction
           v{j}=[v{j},i];
        end
    end
end