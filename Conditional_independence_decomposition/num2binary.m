function v=num2binary(num_species,n)

% convert a number n to a binary array v

v=zeros(1,num_species);

for i=1:num_species
    v(i)= mod(n,2);
    n=(n-v(i))/2;
end