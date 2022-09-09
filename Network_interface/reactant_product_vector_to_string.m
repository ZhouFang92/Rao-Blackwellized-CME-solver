function RPstring=reactant_product_vector_to_string(RPvector)

% RPvector is a column vector. This function convert it into a string
% telling the reactant or the prodcut complex.

RPstring='';
for i=1:size(RPvector,1)
    if RPvector(i)~=0 % check if the i-th species participate in the reaction
        num_molecules=string(RPvector(i));
        if size(RPstring,2)>0 % check if it is the first species appearing in the complex
            RPstring=strcat(RPstring,'+'); % add a plus sign
        end
        if  RPvector(i)~=1
            RPstring=strcat(RPstring,num_molecules); % add the number of molecules
        end
        RPstring=strcat(RPstring,'S',int2str(i));
    end
end

if size(RPstring,2)==0
   RPstring='0';
end

%RPstring=convertCharsToStrings(RPstring);