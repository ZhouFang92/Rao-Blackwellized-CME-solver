function T=Prepare_network_information_for_display(num_reactions,reactant_vector,product_vector,reaction_type,k,propensity)

       for i=1:num_reactions
                   Reactants{i,1}=reactant_product_vector_to_string(reactant_vector(:,i));
                   Products{i,1}=reactant_product_vector_to_string(product_vector(:,i));
                   if reaction_type(i)==1
                       Reaction_type{i,1}="mass action";
                       Reaction_constants{i,1}=k(i);
                   else
                       Reaction_type{i,1}="non-mass action";
                       Reaction_constants{i,1}=[];
                   end
                   Propensity_function{i,1}=propensity(i);
       end
       index=[1:num_reactions]';
       index=strcat("R",string(index));
       T=table(index, Reactants,Products,Reaction_type,Reaction_constants,Propensity_function);