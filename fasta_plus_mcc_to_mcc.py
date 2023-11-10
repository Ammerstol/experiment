def write_mcc_to_nex(proba_x_100, mcc_lines, exp, n_tax_digits):
    codes_file = open(f'experiment_beast/code_{proba_x_100}_exp_{exp}.fasta', 'r')
    lines_codes = codes_file.readlines()
    new_names = {name[1:-3]: name[1:-1] for name in lines_codes[::2]}
    new_mcc_file = open(f'experiment_beast/beast_tree_{proba_x_100}_exp_{exp}.nex', 'w')
    for line in mcc_lines:
        key_label = line[3:-4]
        key_translate = line[5:-5]
        if key_label in new_names:
            new_mcc_file.write(f"\t\t'{new_names[key_label]}'\n")
        else:
            for i in range(n_tax_digits):
                if key_translate[i:] in new_names:
                    new_mcc_file.write(line[:5 + i] + f"{new_names[key_translate[i:]]}',\n")
                    break
            else:
                new_mcc_file.write(line)
    new_mcc_file.close()

def generate_mcc_multiples(probas, exp, proba_index = 0):
    mcc_first_file = open(f'experiment_beast/beast_tree_{probas[proba_index]}_exp_{exp}.nex', 'r')
    mcc_lines = mcc_first_file.readlines()

    for i, symbol in enumerate(mcc_lines[3]):
        if symbol == "=":
            n_tax_digits = len(str(mcc_lines[3][i+1:-2]))
            break

    for proba in probas[:proba_index] + probas[proba_index + 1:]:
        write_mcc_to_nex(proba, mcc_lines, exp, n_tax_digits)