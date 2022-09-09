function v=species_measured_by_observation_channels(num_species,h)

% This function return a cell array the j-th of which shows the species
% measured by the j-th observation channel. 
v{1}=[];
for j=1:size(h,2)
    v{j}=[];
    for i=1:num_species
        substring=strcat("X(",string(i),")");
        if contains(h(j),substring)==1 % check if the i-th species 
                                             % influence the reaction
           v{j}=[v{j},i];
        end
    end
end