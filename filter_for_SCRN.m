%% Description
% 1. type in the name of the network
% 2. set the data size as you wish
% 3. change function summarizing the distribution for RB method.
% 4. Go to the directory RB-SRN App

%%
addpath(genpath('./')); % add path

%input file name
name='Genetic_toggle_switch';
file_name=strcat(name,'_filtering','.mat');
load(file_name); % load the network selected;
data_size=[10^4, 2*10^4, 4*10^4, 10^5];
%TY=[2:2:10];
%Y=[-0.5 29 13 28 39];
sigma=1;

%% generate Matlab script
propensity_function_generator_for_multi_trajectories(propensity,name);
propensity_function_generator(propensity,name);
h_function_generator(h,name);
h_function_vector_form_generator(h,name)
propensity_generator_for_RB(propensity,name);

[species_decomposition_LF,follower_subsystems,reaction_decomposition_LF,...
           reaction_decomposition_F,Size_of_truncated_space]...
          =conditional_independence_decomposition(num_species,reactant_vector,...
           product_vector,propensity,h,truncated_size,threshold_follower_subsystem);
           
       
%% generate observation
% 
% [TX,X_trajectory]=next_reaction_method_single_trajectory...
%                (X0,k,phi,FT,reactant_vector,product_vector);
% [TY,Y]=observations(TX,X_trajectory,h,sigma,observation_period);
% 
% plot(TX,X_trajectory);
% hold on
% plot(TY,Y);

%% different methods.
tic
pi_FSP=FSP_filter(TY,Y,k,phi,reactant_vector,product_vector,...
                                IniDistribution_Mean,...
                                IniDistribution_Variance,...
                                truncated_size,....
                                h,sigma);
toc

for index=1:size(data_size,2)

sample_size=data_size(index);
    
tic
pi_PF=SIRPF_with_CPU(TY,Y,k,phi,reactant_vector,product_vector,...
                                IniDistribution_Mean,...
                                IniDistribution_Variance,...
                                h,sigma,10*sample_size);
t_PF(index)=toc

tic;
pi_RBPF=RBPF_with_CPU...                 
                 (TY,Y,k,phi,FT,num_species,reactant_vector,product_vector,...
                  follower_subsystems,...
                  reaction_decomposition_LF,reaction_decomposition_F,...
                  species_decomposition_LF,...
                  IniDistribution_Mean,IniDistribution_Variance,...
                  sample_size,truncated_size,propensity,h,sigma);
t_RBPF(index)=toc

% error analysis
error1=0;
error2=0;
error3=0;
error4=0;
for i=1:size(Y,2)
  X=pi_PF{1,i};
  W=pi_PF{2,i};
  d_PF=Distribution_PF(W,X,truncated_size);
  X=pi_RBPF{1,i};
  P=pi_RBPF{2,i};
  W=pi_RBPF{3,i};
  d_RBPF=distribution_from_RBPF_Genetic_toggle_switch(W,X,P,truncated_size);
  error1=error1+sum(abs(d_PF-pi_FSP{i}),'all');
  error2=error2+sum(abs(d_RBPF-pi_FSP{i}),'all');
  error3=error3+error_of_marginal_distribution_GTS(pi_FSP{i},d_PF,truncated_size);
  error4=error4+error_of_marginal_distribution_GTS(pi_FSP{i},d_RBPF,truncated_size);
end
error_PF(index)=error1/size(Y,2);
error_RBPF(index)=error2/size(Y,2);
marginal_error_PF(:,index)=error3/size(Y,2);
marginal_error_RBPF(:,index)=error4/size(Y,2);



end 

%% save the result

filename=strcat('./Saved_distributions/',name,'_filter_Result_single_observation');
save(filename,'name','num_species','num_reactions','reactant_vector','product_vector',...
                 'k','size_phi','phi','reaction_type','propensity','h','sigma','truncated_size',...
                 'threshold_follower_subsystem',"FT","observation_period","X0","IniDistribution_Mean",...
                 "IniDistribution_Variance",'sample_size',...
                 "TX","X_trajectory","TY","Y",...
                 "data_size",...
                 "t_PF","t_RBPF","pi_RBPF","pi_FSP","pi_PF",...
                 "error_PF","error_RBPF","marginal_error_PF","marginal_error_RBPF");

%% plot

 plot_filter_error_computational_cost(0,10*data_size,t_PF,error_PF,...
                                        data_size,t_RBPF,error_RBPF);

