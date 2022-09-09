function d=plot_RB_leader_plus_follower(i,j,X,P,truncated_size,follower_subsystems,...
                        species_decomposition_LF)
% i is the follower, j is the leader

sample_size=size(X,2); 

d=zeros(truncated_size(i),truncated_size(j));
                    
% indentify which subsystem do they belongs to 
for index=1:size(follower_subsystems,2)
    if ismember(j,follower_subsystems{index}) ==1
       follower_index2=index;
       break;
    end
end

for sample_index=1:sample_size
    if X(i,sample_index)+1>truncated_size(i)
       continue; % it is outside the truncated space;
    end
    for LI2=1:size(P{follower_index2,1},1)  % tranverse all the probability vector
           truncated_size_follower_subsystem=truncated_size(follower_subsystems{follower_index2});
           MI_j=Linear_index_2_matrix_index(truncated_size_follower_subsystem,LI2);
           a=find(follower_subsystems{follower_index2}==j);
           index_j=MI_j(a);
           value2=P{follower_index2,sample_index}(LI2);
           d(X(i,sample_index)+1,index_j)=...
               d(X(i,sample_index)+1,index_j)+value2/sample_size;
    end
 end