import subprocess
import fasta_plus_mcc_to_mcc
import main_code_generation_experiment



main_code_generation_experiment.main(num_seeds=100)

chain_length, log_every, echo = 10000000, 10000, 1000000
probas_file = open('true_probabilities.txt', 'r')
probas_lines = probas_file.readlines()
probas = [int(float(proba_line[:-1])*100) for proba_line in probas_lines]
exps_file = open('experiments.txt', 'r')
exps_lines = exps_file.readlines()
exps = [int(exp_line[:-1]) for exp_line in exps_lines]
proba_index = 0

for exp in exps:
    # Generate the nexus file from the fasta file (this is dumb) for the chosen true probability in experiment exp
    subprocess.run(f'~/Documents/BEASTGen_v0.3pre_thorney/bin/beastgen  fasta_to_nexus.template  experiment_beast/code_{probas[proba_index]}_exp_{exp}.fasta experiment_beast/code_{probas[proba_index]}_exp_{exp}.nex', shell = True)
    # Generate BEAST file from nexus
    subprocess.run(f'~/Documents/BEASTGen_v0.3pre_thorney/bin/beastgen -D "chain_length={chain_length},log_every={log_every},echo={echo}" code_proba_exp_exp.template experiment_beast/code_{probas[proba_index]}_exp_{exp}.nex experiment_beast/code_{probas[proba_index]}_exp_{exp}.xml', shell = True)
    # Run BEAST
    subprocess.run(f"~/Documents/BEASTv1.10.5pre_thorney_0.1.2/bin/beast -overwrite -working experiment_beast/code_{probas[proba_index]}_exp_{exp}.xml", shell=True)
    # Compute annotated tree
    subprocess.run(f"~/Documents/BEASTv1.10.5pre_thorney_0.1.2/bin/treeannotator -heights ca experiment_beast/code_{probas[proba_index]}_exp_{exp}.trees experiment_beast/beast_tree_{probas[proba_index]}_exp_{exp}.nex", shell=True)
    
    for proba in probas[:proba_index] + probas[proba_index + 1:]:
        # Generate nexus again for all other true probabilities
        subprocess.run(f'~/Documents/BEASTGen_v0.3pre_thorney/bin/beastgen  fasta_to_nexus.template  experiment_beast/code_{proba}_exp_{exp}.fasta experiment_beast/code_{proba}_exp_{exp}.nex', shell = True)
    # 
    fasta_plus_mcc_to_mcc.generate_mcc_multiples(probas, exp, proba_index=proba_index)



# exp = 1
