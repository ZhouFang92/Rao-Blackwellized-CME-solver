%% Description
% 1. type in the name of the network
% 2. set the data size as you wish
% 3. change function summarizing the distribution for RB method.
% 4. Go to the directory RB-SRN App

%%
addpath(genpath('./')); % add path

%input file name
name='Linear_metabolite_pathway_2'; %'Genetic_toggle_switch'; %'Linear_metabolite_pathway_3'  
file_name=strcat('./Networks/',name,'.mat');
load(file_name); % load the network selected;
data_size=[10^4, 2*10^4, 4*10^4, 10^5];

% generate Matlab script
propensity_function_generator_for_multi_trajectories(propensity,name);
propensity_function_generator(propensity,name);
h_function_generator(h,name);
propensity_generator_for_RB(propensity,name);

[species_decomposition_LF,follower_subsystems,reaction_decomposition_LF,...
           reaction_decomposition_F,Size_of_truncated_space]...
          =conditional_independence_decomposition(num_species,reactant_vector,...
           product_vector,propensity,h,truncated_size,threshold_follower_subsystem);
       
%% FSP
       
%  d_FSP=fd(1:truncated_size(1),1:truncated_size(2),1:truncated_size(3),1:truncated_size(4));
%  d_FSP=gather(d_FSP);
%  t_FSP=consumed_time;
tic
    d_FSP=FSP_Power_Series(k,phi,reactant_vector,product_vector,...
                                 FT,...
                                 IniDistribution_Mean,...
                                 IniDistribution_Variance,...
                                 truncated_size);
t_FSP=toc
%% Monte-Carlo method and the RB-CME solvers

parfor i=1:1 i; end;  % connected to workers

for i=1:1%size(data_size,2)
    sample_size=data_size(i);
    
tic
        [TX,X]=next_reaction_method_multi_trajectories_no_GPU...
                (k,phi,FT,reactant_vector,product_vector,...
                 IniDistribution_Mean,IniDistribution_Variance,...
                 10*sample_size);
t_MC(i)=toc
        d_MC=Distribution_Monte_Carlo(X,truncated_size);

 


tic;
   [X,P]=RB_CME_solver...
                 (k,phi,FT,num_species,reactant_vector,product_vector,...
                  follower_subsystems,...
                  reaction_decomposition_LF,reaction_decomposition_F,...
                  species_decomposition_LF,...
                  IniDistribution_Mean,IniDistribution_Variance,...
                  sample_size,truncated_size,propensity);
t_RB(i)=toc

%d_RB=distribution_from_RB_AIC(X,P,truncated_size);
%d_RB=distribution_from_RB_Genetic_toggle_switch(X,P,truncated_size);
%d_RB=Distribution_RB(X,P,truncated_size,follower_subsystems,species_decomposition_LF);
                            
error1=d_FSP-d_MC;

error2=d_FSP-d_RB;

error_MC(i)=sum(abs(error1),'all')
error_RB(i)=sum(abs(error2),'all')
end

%%
filename=strcat('./Saved_distributions/',name,'_Result');
save(filename,'name','num_species','num_reactions','reactant_vector','product_vector',...
                 'k','size_phi','phi','reaction_type','propensity','h','sigma','truncated_size',...
                 'threshold_follower_subsystem',"FT","observation_period","X0","IniDistribution_Mean",...
                 "IniDistribution_Variance",'sample_size','d_FSP',"d_MC","d_RB",...
                 't_MC','t_RB','error_MC','error_RB','data_size','t_FSP');
             
             
%% plot

 plot_error_computational_cost(t_FSP,10*data_size,t_MC,error_MC,...
                                        data_size,t_RB,error_RB);
