function [t,X]=deterministic_process...
               (X0,k,phi,FT,reactant_vector,product_vector)
         
% Initialization
dt=0.01;
reaction_vector=product_vector-reactant_vector;
num_reactions=size(reaction_vector,2);
X=zeros(size(X0,1),FT/dt+1);
X(:,1)=X0;
t=[0:dt:FT];

for i=2:size(t,2)
    f=0;
    for j=1:num_reactions
        f=f+propensity_function(j,X(:,i-1),k,phi)*reaction_vector(:,j);
    end
    dX=f*dt;
    X(:,i)=X(:,i-1)+dX;
end
