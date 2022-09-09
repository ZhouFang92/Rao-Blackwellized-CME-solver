function PS=Mass_action_propensity_generator(reaction_index,Rvector)

PS=strcat("k(",string(reaction_index),")");

for i=1:size(Rvector,1)
    for j=1:Rvector(i,1)
        if j==1
           PS=strcat(PS,"*X(",string(i),")");
        else
           PS=strcat(PS, "*(X(" , string(i), ")-", string(j-1), ")");
        end
    end
end