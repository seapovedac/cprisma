import random

                                    ############## Matrix Generator ##############

def message_error1_in_matrix_generator_pair(met_ref1,lenseq,met_ref):
    print(f"Attention!!! Method selected by the user is '{met_ref1}' but the number of sequences (={lenseq}) is not a even number. The method return to '{met_ref}'.\n\n")
    ali_log.write(f"Attention!!! Method selected by the user is '{met_ref1}' but the number of sequences (={lenseq}) is not a even number. The method return to '{met_ref}'.\n\n")

def message_error2_in_matrix_generator_multiple(met_ref1,dic_ref,met_ref):
    ali_log.write(f"Attention!!! Method selected by the user is '{met_ref1}' but there are not any array in the dictionary {dic_ref}. The method return to '{met_ref}'.\n\n")
    print((f"Attention!!! Method selected by the user is '{met_ref1}' but there are not any array in the dictionary {dic_ref}. The method return to '{met_ref}'."))

def message_error3_in_matrix_generator_multiple(met_ref1,met_ref,no_exist,protein):
    ali_log.write(f"Attention!!! Method selected by the user is '{met_ref1}' but the number sequence {no_exist} no exist in {protein}: {len(protein)} != {no_exist}. The method return to '{met_ref}'.\n\n")
    print(f"Attention!!! Method selected by the user is '{met_ref1}' but the number sequence {no_exist} no exist in {protein}: {len(protein)} != {no_exist}. The method return to '{met_ref}'.")

def message_error4_in_matrix_generator_multiple(met_ref1,met_ref,key):
    ali_log.write(f"Attention!!! Method selected by the user is '{met_ref1}' but there is missing a reference in {key}. The method return to '{met_ref}'.\n\n")
    print(f"Attention!!! Method selected by the user is '{met_ref1}' but there is missing a reference in {key}. The method return to '{met_ref}'.")

def message_error5_in_matrix_generator_multiple(met_ref1,met_ref,dic_ref):
    ali_log.write(f"Attention!!! Method selected by the user is '{met_ref1}' but there is a main key name missing in the array {dic_ref}. This is not recommended! The method return to '{met_ref}'.\n\n")
    print(f"Attention!!! Method selected by the user is '{met_ref1}' but there is a main key name missing in the array {dic_ref}. This is not recommended! The method return to '{met_ref}'.")


                                    ################   Features   ################

br_op='{'
br_co='}'

def message_error_in_main_key(ali_log,suffix_feature,key,final_sentece):
    ali_log.write(f"\nAttention!!! The main key '{key}' is missing in dict_{suffix_feature}. The arrays dict_{suffix_feature} != dict_ref or the name of the main key is not the same. {final_sentece}.\n\n")
    print(f"Attention!!! The main key '{key}' is missing in dict_{suffix_feature}. The arrays dict_{suffix_feature} != dict_ref or the name of the main key is not the same. {final_sentece}.")

def message_error_in_comparison(ali_log,suffix_feature,key,final_sentece,len_dr,len_do):
    ali_log.write(f"\nAttention!!! The main key '{key}' does not have the same number of comparisons {len_dr} (dict_ref) != {len_do} (dict_{suffix_feature}). {final_sentece}.\n\n")
    print(f"Attention!!! The main key '{key}' does not have the same number of comparisons {len_dr} (dict_ref) != {len_do} (dict_{suffix_feature}). {final_sentece}.")

### Operation ###

def message_error1_in_operation(ali_log,suffix_feature,key,final_sentece,f):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{f}' and this is invalid or missing! It is only possible ['d','da','m','n'] inside brackets '{br_op}{br_co}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{f}' and this is invalid or missing! It is only possible ['d','da','m','n'] inside brackets '{br_op}{br_co}'. {final_sentece}.")

def message_error2_in_operation(ali_log,suffix_feature,key,final_sentece,f):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{list(f)[0]}' and this is invalid or missing! Maybe your descriptor need aditional brackets '{br_op}{br_co}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{list(f)[0]}' and this is invalid or missing! Maybe your descriptor need aditional brackets '{br_op}{br_co}'. {final_sentece}.")

def message_error3_in_operation(ali_log,suffix_feature,key,final_sentece,f):
    ali_log.write(f"\nAttention!!! In the main key '{key}' of dict_{suffix_feature}, the operation type '{list(f)[0] }' does not exist. It is only possible ['d','da','m','n']. {final_sentece}.\n\n")
    print(f"Attention!!! In the main key '{key}' of dict_{suffix_feature}, the operation type '{list(f)[0] }' does not exist. It is only possible ['d','da','m','n']. {final_sentece}.")

### Maximum ###

def message_error1_in_maximum(ali_log,suffix_feature,key,final_sentece):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' the brackets '{br_op}{br_co}' are empty. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' the brackets '{br_op}{br_co}' are empty. {final_sentece}.")

def message_error2_in_maximum(ali_log,suffix_feature,key,final_sentece,f):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{list(f)[0]}' and this is invalid or missing! It is only possible ['nr','rt','ra','rm','ram','rs','rsa','rsm','rsam','rspi','rspim','rpi'] or the same descriptors with aditional letter 'Y' (e.g. nrY) inside brackets '{br_op}{br_co}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{list(f)[0]}' and this is invalid or missing! It is only possible ['nr','rt','ra','rm','ram','rs','rsa','rsm','rsam','rspi','rspim','rpi'] or the same descriptors with aditional letter 'Y' (e.g. nrY) inside brackets '{br_op}{br_co}'. {final_sentece}.")

def message_error3_in_maximum(ali_log,suffix_feature,key,final_sentece,f):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' has missing parameters for the descriptor '{list(f)[0]}'! {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' has missing additional parameters for the descriptor '{list(f)[0]}'! {final_sentece}.")

def message_error4_in_maximum(ali_log,suffix_feature,key,final_sentece,sub_k):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor {br_op}'{sub_k}'{br_co} and this is invalid or missing! It is only possible ['nr','rt','ra','rm','ram','rs','rsa','rsm','rsam','rspi','rspim','rpi'] or the same descriptors with aditional letter 'Y' (e.g. nrY). {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor {br_op}'{sub_k}'{br_co} and this is invalid or missing! It is only possible ['nr','rt','ra','rm','ram','rs','rsa','rsm','rsam','rspi','rspim','rpi'] or the same descriptors with aditional letter 'Y' (e.g. nrY). {final_sentece}.")

def message_error5_in_maximum(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has unnecessary parameters: {sub_v}! {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has unnecessary parameters: {sub_v}! {final_sentece}.")

def message_error6_in_maximum_rt(ali_log,suffix_feature,key,final_sentece,sub_k):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have only two parameters one numeric (e.g. 1.0) and another as symbol of order relation ['>=', '<=', '>', '<', '==', '!=']. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have only two parameters one numeric (e.g. 1.0) and another as symbol of order relation ['>=', '<=', '>', '<', '==', '!=']. {final_sentece}.")

def message_error7_in_maximum_rt(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the first parameter of the descriptor '{sub_k}' of the main key '{key}' should have an integer or float number (e.g. 1.0) and not '{sub_v[0]}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the first parameter of the descriptor '{sub_k}' of the main key '{key}' should have an integer or float number (e.g. 1.0) and not '{sub_v[0]}'. {final_sentece}.")

def message_error8_in_maximum_rt(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the second parameter ({sub_v[1]}) of the descriptor '{sub_k}' of the main key '{key}' should be a symbol of order relation ['>=', '<=', '>', '<', '==', '!=']. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the second parameter ({sub_v[1]}) of the descriptor '{sub_k}' of the main key '{key}' should be a symbol of order relation ['>=', '<=', '>', '<', '==', '!=']. {final_sentece}.")

def message_error9_in_maximum_rpi(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,index):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of numbers (e.g. [{random.choice(index)}]) and not '{sub_v}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of numbers (e.g. [{random.choice(index)}]) and not '{sub_v}'. {final_sentece}.")

def message_error10_in_maximum_rpi(ali_log,suffix_feature,key,final_sentece,sub_k,index):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of numbers belonging to the processed data and not ' '. Numbers should be in a range of {min(list(index))}-{max(list(index))}. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of numbers belonging to the processed data and not ' '. Numbers should be in a range of {min(list(index))}-{max(list(index))}. {final_sentece}.")

def message_error11_in_maximum_rpi(ali_log,suffix_feature,key,final_sentece,sub_k,index,lpi):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of numbers belonging to the processed data and not '{lpi}'. Numbers should be in a range of {min(list(index))}-{max(list(index))}. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of numbers belonging to the processed data and not '{lpi}'. Numbers should be in a range of {min(list(index))}-{max(list(index))}. {final_sentece}.")

def message_error12_in_maximum_ra(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of amino acids (e.g. ['D', 'R', 'T']) and not '{sub_v}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of amino acids (e.g. ['D', 'R', 'T']) and not '{sub_v}'. {final_sentece}.")

def message_error13_in_maximum_ra(ali_log,suffix_feature,key,final_sentece,sub_k,tr_r,lam):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have as elements amino acids with one letter format (e.g. [{random.choice(list(tr_r)[0])}]) and not ['{lam}']. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have as elements amino acids with one letter format (e.g. [{random.choice(list(tr_r)[0])}]) and not ['{lam}']. {final_sentece}.")

def message_error14_in_maximum_ra(ali_log,suffix_feature,key,final_sentece,sub_k,tr_r,lam):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have as elements amino acids with one letter format (e.g. [{random.choice(list(tr_r)[0])}]) and related to target-residues of the tuple {tr_r}, and not ['{lam}']. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have as elements amino acids with one letter format (e.g. [{random.choice(list(tr_r)[0])}]) and related to target-residues of the tuple {tr_r}, and not ['{lam}']. {final_sentece}.")

def message_error15_in_maximum_rs(ali_log,suffix_feature,key,final_sentece,sub_k,list_gen):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a list with at least one number sequence similar and *not different* like those reported in dict_ref {list_gen} (e.g. '{sub_k}': [at least one number sequence of dict_ref]). {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a list with at least one number sequence similar and *not different* like those reported in dict_ref {list_gen} (e.g. '{sub_k}': [at least one number sequence of dict_ref]). {final_sentece}.")

def message_error16_in_maximum_rsa(ali_log,suffix_feature,key,final_sentece,sub_k,list_gen,tr_r):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a dictionary with at least one number sequence similar and *not different* like those reported in dict_ref {list_gen} in the 'third key' followed by a 'third value' as a list of any amino acids related to target-residues.\n")
    ali_log.write(f"            1- e.g. '{sub_k}': [{br_op}'third key' with at least one number sequence of dict_ref: ['third value' as a list of any amino acids related to target-residues] {br_co}].\n")
    ali_log.write(f"            2- e.g. '{sub_k}': [{br_op} '{random.choice(list_gen)}': {random.choice(list(tr_r))} {br_co}].\n")
    ali_log.write(f"            {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a dictionary with at least one number sequence similar and *not different* like those reported in dict_ref {list_gen} in the 'third key' followed by a 'third value' as a list of any amino acids related to target-residues.")
    print(f"            1- e.g. '{sub_k}': [{br_op}'third key' with at least one number sequence of dict_ref: ['third value' as a list of any amino acids related to target-residues] {br_co}].")
    print(f"            2- e.g. '{sub_k}': [{br_op} '{random.choice(list_gen)}': {random.choice(list(tr_r))} {br_co}].")
    print(f"            {final_sentece}.")

def message_error17_in_maximum_rspi(ali_log,suffix_feature,key,final_sentece,sub_k,list_gen,index):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a dictionary with at least one number sequence similar and *not different* like those reported in dict_ref {list_gen} in the 'third key' followed by a 'third value' as a list of position index numbers.\n")
    ali_log.write(f"            1- e.g. '{sub_k}': [{br_op}'third key' with at least one number sequence of dict_ref: ['third value' as a list of position index numbers] {br_co}].\n")
    ali_log.write(f"            2- e.g. '{sub_k}': [{br_op} '{random.choice(list_gen)}': [{random.choice(index)}] {br_co}].\n")
    ali_log.write(f"            {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a dictionary with at least one number sequence similar and *not different* like those reported in dict_ref {list_gen} in the 'third key' followed by a 'third value' as a list of position index numbers.")
    print(f"            1- e.g. '{sub_k}': [{br_op}'third key' with at least one number sequence of dict_ref: ['third value' as a list of position index numbers] {br_co}].")
    print(f"            2- e.g. '{sub_k}': [{br_op} '{random.choice(list_gen)}': [{random.choice(index)}] {br_co}].")
    print(f"            {final_sentece}.")

def message_error18_in_maximum_rspi(ali_log,suffix_feature,key,final_sentece,sub_k,index,sv):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of numbers (e.g. [{random.choice(index)}]) in the 'third key' and not '{sv}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should be a list of numbers (e.g. [{random.choice(index)}]) in the 'third key' and not '{sv}'. {final_sentece}.")

def message_error_in_special_case(ali_log,suffix_feature,descriptor,met_ref):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{descriptor}', it is not supported for the combination between this descriptor and the reference method '{met_ref}'. All the descriptors will be equal to 'nr'.\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{descriptor}', it is not supported for the combination between this descriptor and the reference method '{met_ref}'. All the descriptors will be equal to 'nr'.")

### Color ###

def message_error19_in_color(ali_log,suffix_feature,key,final_sentece):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' the brackets '{br_op}{br_co}' are empty. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' the brackets '{br_op}{br_co}' are empty. {final_sentece}.")

def message_error20_in_color(ali_log,suffix_feature,key,final_sentece,f):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{list(f)[0]}' and this is invalid or missing! It is only possible ['nc','ssc','fsc','fac','fmac,'pic','pimc','tc','tmc','pkac'] inside brackets '{br_op}{br_co}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{list(f)[0]}' and this is invalid or missing! It is only possible ['nc','ssc','fsc','fac','fmac,'pic','pimc','tc','tmc','pkac'] inside brackets '{br_op}{br_co}'. {final_sentece}.")

def message_error21_in_color(ali_log,suffix_feature,key,final_sentece,f):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' has missing parameters for the descriptor '{list(f)[0]}'! {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' has missing additional parameters for the descriptor '{list(f)[0]}'! {final_sentece}.")

def message_error22_in_color(ali_log,suffix_feature,key,final_sentece,sub_k):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{sub_k}' and this is invalid or missing! It is only possible ['nc','ssc','fsc','fac','fmac,'pic','pimc','tc','tmc','pkac']. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' has the descriptor '{sub_k}' and this is invalid or missing! It is only possible ['nc','ssc','fsc','fac','fmac,'pic','pimc','tc','tmc','pkac']. {final_sentece}.")

def message_error23_in_color(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has unnecessary parameters: {sub_v}! {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has unnecessary parameters: {sub_v}! {final_sentece}.")

def message_error24_in_color_free(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,name_color_file):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a list of index numbers related to specific colors of {name_color_file} or 'random' option, and not '{sub_v}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a list of index numbers related to specific colors of {name_color_file} or 'random' option, and not '{sub_v}'. {final_sentece}.")

def message_error25_in_color_free(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,ld):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' have a list in dict_{suffix_feature} with different quantity of elements compared to the dict_ref list (target-sequences, reference sequence is not taking into account): {len(sub_v)} {sub_v}) != {len(ld)} {ld}. The 'random' option is possible to do not specify a list. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' have a list in dict_{suffix_feature} with different quantity of elements compared to the dict_ref list (target-sequences, reference sequence is not taking into account): {len(sub_v)} {sub_v} != {len(ld)} {ld}. The 'random' option is possible to do not specify a list. {final_sentece}.")

def message_error26_in_color_free(ali_log,suffix_feature,key,final_sentece,sub_k,sv,name_color_file,list_colors):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' have a list with a number index color ('{sv}') out of range colors reported in {name_color_file} ({min(list_colors)}-{max(list_colors)}) or is not a variable type 'int'. The 'random' option is possible to do not specify a list. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' have a list with a number index color ('{sv}') out of range colors reported in {name_color_file} ({min(list_colors)}-{max(list_colors)}) or is not a variable type 'int'. The 'random' option is possible to do not specify a list. {final_sentece}.")

def message_error27_in_color_defined(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should has a boolean variable as a value and not '{sub_v}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should has a boolean variable as a value and not '{sub_v}'. {final_sentece}.")

def message_error28_in_color_defined(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,tr_r):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' have a list in dict_{suffix_feature} with different quantity of elements compared to the target-residues tuple: {len(sub_v)} {sub_v} != {len(tr_r)} {tr_r}. The 'random' option is possible to do not specify a list. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' have a list in dict_{suffix_feature} with different quantity of elements compared to the target-residues tuple: {len(sub_v)} {sub_v} != {len(tr_r)} {tr_r}. The 'random' option is possible to do not specify a list. {final_sentece}.")

def message_error29_in_color_defined(ali_log,suffix_feature,key,final_sentece,sub_k,ch_r):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} to apply the descriptor '{sub_k}' of the main key '{key}' your tuple of target-residues should have at least one ionizable residue like {ch_r}. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} to apply the descriptor '{sub_k}' of the main key '{key}' your tuple of target-residues should have at least one ionizable residue like {ch_r}. {final_sentece}.")

def message_error30_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should has a list of dictionaries with color position index conditions instead of '{sub_v}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should has a list of dictionaries with color position index conditions instead of '{sub_v}'. {final_sentece}.")

def message_error31_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has an element in the list of dictionaries not valid: '{supi}'. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has an element in the list of dictionaries not valid: '{supi}'. {final_sentece}.")

def message_error32_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,ksupi,key_seq):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has an invalid sequence ('{ksupi}') in '{supi}'. Sequence should be equal to '{key_seq}' in the subkey. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has an invalid sequence ('{ksupi}') in '{supi}'. Sequence should be equal to '{key_seq}' in the subkey. {final_sentece}.")

def message_error33_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,ksupi,tr_seq):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has an invalid sequence in the subkey ('{ksupi}'), check: '{supi}'. Number sequence should belong to '{tr_seq}'. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has an invalid sequence in the subkey ('{ksupi}'), check: '{supi}'. Number sequence should belong to '{tr_seq}'. {final_sentece}.")

def message_error34_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,max_col):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has problem in the first two parameters of the list '{vsupi}'. First element should be a color number (int) <= {max_col} and second element should be a format parameter like ['n', 'b', 'i', 'd', 'u']. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has problem in the first two parameters of the list '{vsupi}'. First element should be a color number (int) <= {max_col} and second element should be a format parameter like ['n', 'b', 'i', 'd', 'u']. {final_sentece}.")

def message_error35_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,max_ind):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has at least one no appropriate interval in the list '{vsupi}'. The first value in the interval must be < than the second one (e.g. 17-26). Besides both of them must be <= {max_ind} and >= 0. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has at least one no appropriate interval in the list '{vsupi}'. The first value in the interval must be < than the second one (e.g. 17-26). Besides both of them must be <= {max_ind} and >= 0. {final_sentece}.")

def message_error36_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,nsupi,max_ind):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has a no appropriate value ('{nsupi}') in the list '{vsupi}'. The value must be <= {max_ind} and >= 0. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' has a no appropriate value ('{nsupi}') in the list '{vsupi}'. The value must be <= {max_ind} and >= 0. {final_sentece}.")

def message_error37_in_color_threshold(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a list with three elements: one color number, one float/int number, and one symbol of order relation ['>=', '<=', '>', '<', '==', '!=']. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}' should have a lits with three: one color number, one float/int number, and one symbol of order relation ['>=', '<=', '>', '<', '==', '!=']. {final_sentece}.")

def message_error38_in_color_threshold(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,color,first):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}', the list {sub_v} should have at the first position a color number between 0 - {color} and not '{first}'. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}', the list {sub_v} should have at the first position a color number between 0 - {color} and not '{first}'. {final_sentece}.")

def message_error39_in_color_threshold(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,color,second):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}', the list {sub_v} should have at the second position a float/int number and not '{second}'. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}', the list {sub_v} should have at the second position a float/int number and not '{second}'. {final_sentece}.")

def message_error40_in_color_threshold(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,color,third):
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}', the list {sub_v} should have at the third position a symbol of order relation ['>=', '<=', '>', '<', '==', '!='] and not '{third}'. {final_sentece}.\n\n")
    print(f"\nAttention!!! In dict_{suffix_feature} the descriptor '{sub_k}' of the main key '{key}', the list {sub_v} should have at the third position a symbol of order relation ['>=', '<=', '>', '<', '==', '!='] and not '{third}'. {final_sentece}.")

### visualization ###

def message_error30_in_visualization(ali_log,suffix_feature,key,final_sentece):
    list_Re=['ReY','ReN']
    list_De=['DeY','DeN']
    list_Le=['LeY','LeN']
    ali_log.write(f"\nAttention!!! In dict_{suffix_feature} the main key '{key}' has an invalid descriptor or missing! It is only possible three different parameters simultaneously (e.g. '{random.choice(list_Re)}','{random.choice(list_De)}','{random.choice(list_Le)}') inside brackets '{br_op}{br_co}'. {final_sentece}.\n\n")
    print(f"Attention!!! In dict_{suffix_feature} the main key '{key}' has an invalid descriptor or missing! It is only possible three different parameters simultaneously (e.g. '{random.choice(list_Re)}','{random.choice(list_De)}','{random.choice(list_Le)}') inside brackets '{br_op}{br_co}'. {final_sentece}.")

                                    ################   Data Array Relationship   ################

def message1_error_in_data_array_relationship(ali_log,flag_des,rf_n):
    ali_log.write(f"Attention!!! For the descriptor '{flag_des}' will no be considered maximum restriction, because you are using the same sequence number as a reference [{rf_n}].\n")
    print(f"Attention!!! For the descriptor '{flag_des}' will no be considered maximum restriction, because you are using the same sequence number as a reference [{rf_n}].")

def message2_error_in_data_array_relationship(ali_log,descriptor,mutation_col):
    ali_log.write(f"Attention!!! For the descriptor '{descriptor}' the color for mutation {mutation_col+2} it is not accepted, only values <= 142 could be considered. The color for mutation now is gray.\n")
    print(f"Attention!!! For the descriptor '{descriptor}' the color for mutation {mutation_col+2} it is not accepted, only values <= 142 could be considered. The color for mutation now is gray.")

def message3_error_in_data_array_relationship(ali_log,descriptor,sequence_col):
    ali_log.write(f"Attention!!! For the descriptor '{descriptor}' the color for sequence {sequence_col+2} it is not accepted, only values <= 142 could be considered. The color for sequence now is red.\n")
    print(f"Attention!!! For the descriptor '{descriptor}' the color for sequence {sequence_col+2} it is not accepted, only values <= 142 could be considered. The color for sequence now is red.")
