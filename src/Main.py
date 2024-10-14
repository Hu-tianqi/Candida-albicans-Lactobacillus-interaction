import cobra
import glob
import pickle
import copy
import cplex
import math
from cobra import Model, Reaction, Metabolite
import pandas as pd
from functions import *


#import C. albicans model
candida_albicans = cobra.io.read_sbml_model('D:\Candida-albicans-Lactobacillus-interaction\models\candida_albicans.xml')


#read metabolite mapping file
id_mapping = pd.read_csv('D:\Candida-albicans-Lactobacillus-interaction\dat\id_mapping.csv',index_col=0)


#make C. albicans model with lumen reactions
candida_albicans = join_lumen_rxns(candida_albicans,id_mapping,Metabolite,Reaction)


#import diets

#based on C. albicans lumen exchange ids. Also works for joint model lumen exchange reactions
diet_candida = list()
diet_candida.append(pickle.load(open("D:\Candida-albicans-Lactobacillus-interaction/dat/calb_dmem3_media.pckl","rb")))



#based on agora lumen exchange ids
diet_agora = list()
diet_agora.append(pickle.load(open("D:\Candida-albicans-Lactobacillus-interaction/dat/agora_dmem3_media.pckl","rb")))





#import simulated models
simulated_models = pd.read_csv("dat/bac_info.csv")
simulated_models = list(simulated_models["Model_ID"])

#Agora models 
considered_models = pd.read_csv("dat/bac_info.csv")
considered_models = list(considered_models["Model_ID"])

#take the model directories
model_dir_agora = glob.glob('models/agora/*.xml')
model_dir = model_dir_agora

#create a list of model source
model_source = ['agora']*len(model_dir_agora)

##implement whole pairwise simulations
#index 0 -> DMEM
gr_calb_indv_all = list()
gr_calb_pairs_all = list()
gr_bacs_pairs_all = list()
gr_bacs_indv_all = list()
solution_status_all = list()
coefs_calb_all = list()
coefs_bacs_all = list()
int_types_all = list()
int_types_all_considered_models = list()
roles_all = list()
individual_flux_distribution = pd.DataFrame()
pairwise_flux_distribution = pd.DataFrame()
ca_indiv_flux_distribution = pd.DataFrame()

#set diet on C. albicans' exchange constraints
candida_albicans = constrain_model(candida_albicans,diet_candida[0],"EX_u")

candida_albicans.solver = 'cplex'

#do FBA for C. albicans alone
candida_albicans.objective = 'scerbiomasspseudoreaction'
fba_solution = candida_albicans.optimize()
gr_calb_indv = fba_solution.objective_value
#do pFBA for C. albicans alone
sol_ca = cobra.flux_analysis.pfba(candida_albicans)
ca_indiv_flux_distribution = pd.DataFrame(index=sol_ca.fluxes.index,data = {candida_albicans.id:sol_ca.fluxes})

#iterate for pairwise simulation with all the models
    
gr_calb_pairs = list()
gr_bacs_pairs = list()
gr_bacs_indv = list()
bacs_model_ids = list()
solutions = list()
solution_status = list()
coefs_calb = list()
coefs_bacs = list()
int_types = list()
int_types_considered_models = list()
roles = list()
    
for ii in range(len(model_dir)):
       

        #read sbml models
        bac = cobra.io.read_sbml_model(model_dir[ii])

        #skip the model if it is not in desired list of models
        #if bac.id not in simulated_models:
            #continue
        
        bac.solver = 'cplex'

        
        #set diet on model's exchange constraints
        #if model_source =='agora':
        bac = constrain_model(bac,diet_agora[0],"EX_")
      
        # else:
           # bac = constrain_model(bac,diet_carveme[0],"EX_")

        #do FBA for the model alone
        sol = bac.optimize()
        gr_bac_indv = sol.objective_value
        print(gr_bac_indv)
        
        #check if model is able to grow on diet. Go for the next one if it is unable
        if gr_bac_indv<0.0001:
            gr_calb_pairs.append('')
            gr_bacs_pairs.append('')
            gr_bacs_indv.append(0)
            bacs_model_ids.append(bac.id)
            solutions.append('')
            solution_status.append('')
            coefs_calb.append('')
            coefs_bacs.append('')
            roles.append('')
            int_types.append('')
            continue

        #do pFBA for feasible models
        sol = cobra.flux_analysis.pfba(bac)

        #merge dataframes to expand the table of pFBA-derived flux distribution for bacterial models
        df2merge = pd.DataFrame(index=sol.fluxes.index,data = {bac.id:sol.fluxes})
        individual_flux_distribution = pd.merge(individual_flux_distribution,df2merge,how='outer',left_index=True,right_index=True)

        bacs_model_ids.append(bac.id)
        gr_bacs_indv.append(gr_bac_indv)

        
        #change IDs
        bac = change_ids(bac)

        [bac,bac_lumen_ext_rxns,bac_lumen_ext_mets] = tune_bac_lumen_rxns(bac,candida_albicans,Metabolite,Reaction,model_source[ii])

        joint_model = make_joint_model(candida_albicans,bac,bac_lumen_ext_rxns,bac_lumen_ext_mets,diet_candida[0],copy,model_source[ii])


        #constrain the joint model by using the diet
        joint_model = constrain_model(joint_model,diet_candida[0],"EX_u")
        
        with open('results/joint_model_problem.lp', 'w') as out:
            out.write(str(joint_model.solver))

    
        cpx = cplex.Cplex("results/joint_model_problem.lp")


        cpx,loc_bio_1,loc_bio_2,first_forward_vars,first_backward_vars,second_forward_vars,second_backward_vars = apply_coupling_constraints(cpx,cplex,model_source[ii])


        cpx = pFBA_joint_model(cpx,loc_bio_1,loc_bio_2,cplex,first_forward_vars,first_backward_vars,second_forward_vars,second_backward_vars)


        
        solution_status.append(cpx.solution.get_status())
        
        gr_calb_pair = cpx.solution.get_values(loc_bio_1)
        gr_bac_pair = cpx.solution.get_values(loc_bio_2)
        
        gr_calb_pairs.append(gr_calb_pair)
        gr_bacs_pairs.append(gr_bac_pair)

        coef_calb = math.log2(gr_calb_pair/gr_calb_indv)
        coef_bac = math.log2(gr_bac_pair/gr_bac_indv)
        
        coefs_calb.append(coef_calb)
        coefs_bacs.append(coef_bac)

        #unify the cplex variables and merge dataframes to expand the table of pFBA-derived flux distribution for joined models
        cplex_vars = cpx.variables.get_names()
        reverse_loc = list(map(lambda x:x.find('reverse'),cplex_vars))
        revised_cplex_vars = list(map(lambda x,y:x[:y+7] if y!=-1 else x,cplex_vars,reverse_loc))
        hashtag_loc = list(map(lambda x:x.find('#'),revised_cplex_vars))
        revised_cplex_vars = list(map(lambda x,y:x[:y] if y!=-1 else x,revised_cplex_vars,hashtag_loc))
        
        df2merge = pd.DataFrame(index=revised_cplex_vars,data = {'C_albicans-'+bac.id:cpx.solution.get_values()})
        pairwise_flux_distribution = pd.merge(pairwise_flux_distribution,df2merge,how='outer',left_index=True,right_index=True)

        #interaction type
        int_type,role = interaction_type(coef_calb,coef_bac,math)

        int_types.append(int_type)
        roles.append(role)
        
        if bac.id in considered_models:
            int_types_considered_models.append(int_type)
    
gr_calb_indv_all.append([gr_calb_indv]*len(bacs_model_ids))
gr_calb_pairs_all.append(gr_calb_pairs)
gr_bacs_pairs_all.append(gr_bacs_pairs)
gr_bacs_indv_all.append(gr_bacs_indv)
coefs_calb_all.append(coefs_calb)
coefs_bacs_all.append(coefs_bacs)
solution_status_all.append(solution_status)
int_types_all.append(int_types)
int_types_all_considered_models.append(int_types_considered_models)
roles_all.append(roles)
individual_flux_distribution = individual_flux_distribution.transpose()
pairwise_flux_distribution = pairwise_flux_distribution.transpose()
    
#End of pairwise simulations
#export pFBA-derived flux distribution for the individual models
ca_indiv_flux_distribution.to_csv('results/pfba_fluxes_indivdual_ca_dmem3.csv')
#export pFBA-derived flux distribution for the individual models
individual_flux_distribution.to_csv('results/pfba_fluxes_indivdual_bacs_dmem3.csv')


#export pFBA-derived flux distribution for the joined models
pairwise_flux_distribution.to_csv('results/pfba_fluxes_pairwise_dmem3.csv')


#export pair wise simulation results
results_dmem3 = pd.DataFrame({'c_albicans_single_growth':gr_calb_indv_all[0],'c_albicans_paired_growth':gr_calb_pairs_all[0],
                                'bacteria_single_growth':gr_bacs_indv_all[0],'bacteria_paired_growth':gr_bacs_pairs_all[0],
                                'solution_status':solution_status_all[0],'c_albicans_coefficient':coefs_calb_all[0],
                                'bacteria_coefficient':coefs_bacs_all[0],'interaction_type':int_types_all[0],
                                'bacteria_interaction_role':roles_all[0]},index = bacs_model_ids)


                                 


results_dmem3.to_csv('results/results_dmem3.csv')



#interaction distribution
dist = list()
dist.append(int_types_all[0].count('Parasitism [-/+]'))
dist.append(int_types_all[0].count('Commensalism [o/+]'))
dist.append(int_types_all[0].count('Amensalism [o/-]'))
dist.append(int_types_all[0].count('Neutralism [o/o]'))
dist.append(int_types_all[0].count('Parasitism [+/-]'))
dist.append(int_types_all[0].count('Commensalism [+/o]'))
dist.append(int_types_all[0].count('Mutualism [+/+]'))
dist.append(int_types_all[0].count('Amensalism [-/o]'))
dist.append(int_types_all[0].count('Competition [-/-]'))

category = ['Parasitism [-/+]',
            'Commensalism [o/+]',
            'Amensalism [o/-]',
            'Neutralism [o/o]',
            'Parasitism [+/-]',
            'Commensalism [+/o]',
            'Mutualism [+/+]',
            'Amensalism [-/o]',
            'Competition [-/-]'] 

Media =['dmem3 diet']*9
interaction_distribution = pd.DataFrame({'count':dist,'category':category,'Media':Media})

interaction_distribution.to_csv('results/interaction_distribution.csv')




