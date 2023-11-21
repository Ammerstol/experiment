def write_mcc_to_nex(proba_x_100, lines_mcc):
    codes_file = open(f'experiment_beast/code_{proba_x_100}.fasta', 'r')
    lines_codes = codes_file.readlines()
    new_names = {name[1:-3]: name[1:-1] for name in lines_codes[::2]}
    new_mcc_file = open(f'experiment_beast/beast_tree_{proba_x_100}.nex', 'w')
    for line in lines_mcc:
        key_label = line[3:-4]
        key_translate = line[5:-5]
        if key_label in new_names:
            new_mcc_file.write(f"\t\t'{new_names[key_label]}'\n")
        elif key_translate in new_names:
            new_mcc_file.write(line[:5] + f"{new_names[key_translate]}',\n")
        elif key_translate[1:] in new_names:
            new_mcc_file.write(line[:6] + f"{new_names[key_translate[1:]]}',\n")
        elif key_translate[2:] in new_names:
            new_mcc_file.write(line[:7] + f"{new_names[key_translate[2:]]}',\n")
        elif key_translate[3:] in new_names:
            new_mcc_file.write(line[:8] + f"{new_names[key_translate[3:]]}',\n")
        else:
            new_mcc_file.write(line)
    new_mcc_file.close()


probas_file = open('true_probabilities.txt', 'r')
lines_mcc = probas_file.readlines()
probas = [int(float(proba_line[:-1])*100) for proba_line in lines_mcc]


# Using readlines()
mcc_5_file = open(f'experiment_beast/beast_tree_{probas[0]}.nex', 'r')
lines_mcc = mcc_5_file.readlines()


for proba in probas[1:]:
    write_mcc_to_nex(proba,lines_mcc)