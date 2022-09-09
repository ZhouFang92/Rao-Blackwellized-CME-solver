function d=Distribution_PF(W,X,truncated_size)

% This algorithm returns the distribution by the particles from a particle
% filter. 

W=W/sum(W);
sample_size=size(X,2); % number of particles
d=zeros(truncated_size'); % initialize the distribution
if size(truncated_size,1)==1
   d=zeros(truncated_size,1); 
end
X=gather(X); % store the particles on CPU memory

for j=1:sample_size
      i=X(:,j)+1; % get the index for the j-th particle
      % update the distribution
      if all(i<=truncated_size)==1
        i=Matrix_index_2_linear_index(truncated_size,i);
        d(i)=d(i)+W(j);
      end
end

