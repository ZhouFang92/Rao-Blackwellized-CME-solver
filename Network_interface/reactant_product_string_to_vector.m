function [RPvector,check]=reactant_product_string_to_vector(num_reactions,RPstring)

% convert a string showing a reactant or product complex into a vector
% check==1 succeed, ortherwise 0

RPvector=zeros(num_reactions,1);
strV=strsplit(RPstring,'+'); % obtain 4Si

if strcmp(RPstring,"0")==1 || strcmp(RPstring,"")==1
   check=1;
   return;
end

for i=1:size(strV,2)
    v=strsplit(strV(i),'S'); % obtain 4 i;
    if size(v,2)~=2
       check =0;
       return;
    end
    % obtain the index of the species and the number of molecules
    species_index=str2num(v(2)); 
    if strcmp(v(1),"")==1
       num_molecules=1;
    else
       num_molecules=str2num(v(1));
    end
    RPvector(species_index,1)=num_molecules;
end

check=1;