from cprisma.visualization import final_z

    ######## Treatment to the file alignment input ########

def treatment_alignment(turn_com,name,ali_log,ch_r,hy_r,ck_ns,name_sequence,counter):

    if counter == 0:
        print(f"Processing {name}.")
        ali_log.write(f"{' '*35}--- Data from {name} ---\n\n")          # Writing log file

    with open (name) as ali:                                    # Treatment sequence alignment

        try:
            data=[d for d in ali.read().split('\n') if d != '']     # Split alignment file
            len_d=len(data[0])
            dic_ali={'conserv':''}              # Name or key of the conservation line
            for d in data:                      # Reading the lines
                if any(c in d[:len_d] for c in ['*','.',':']):
                    dic_ali['conserv']+=d[16:]
                else:
                    aux=d.split()
                    if aux[0] not in dic_ali.keys():
                        dic_ali[aux[0]]=aux[1]          # Name or key of each sequence
                    else:
                        dic_ali[aux[0]]+=aux[1]         # Rest of the sequence
        except IndexError:
            ali_log.write(f"The program stop! Improper alignment format. Check ClustalW format at https://www.ebi.ac.uk/Tools/msa/muscle.\n")
            print(f"The program stop! Improper alignment format. Check ClustalW format at https://www.ebi.ac.uk/Tools/msa/muscle.")
            quit()


    #### Giving names to each sequence ####

    l = list(dic_ali.keys())[1:]                    # Defining the name of each sequence
    cc=0
    if ck_ns:
        if len(name_sequence) == len(l):
            for n in name_sequence:
                old=l[cc]
                dic_ali[n] = dic_ali.pop(old)
                cc+=1
            if counter == 0:
                ali_log.write(f"Variable check_name_seq = {ck_ns}, so the sequences' names have been changed.\n\n")
                print(f"Variable check_name_seq = {ck_ns}, so the sequences' names have been changed.")
        else:
            if counter == 0:
                ali_log.write(f"Variable check_name_seq = {ck_ns} but the sequences' number does not match with tuple input {name_sequence}.\n\n")
                print(f"Variable check_name_seq = {ck_ns} but the sequences' number does not match with tuple input: {name_sequence}.")
    else:
        if counter == 0:
            ali_log.write(f"Variable check_name_seq = {ck_ns}, so the sequences' names of {name} have been conserved.\n\n")
            print(f"Variable check_name_seq = {ck_ns}, so the sequences' names of {name} have been conserved.")
    l = list(dic_ali.keys())[1:]

    #### Putting alignment input to log file ####

    list_big,max_z=final_z(l) # Function to do not have problem with spaces

    cc=0
    for key in l:
        z=list_big[cc]
        if counter == 0 and turn_com == True:
            ali_log.write(f"{key}{' '*z}\t{dic_ali[key]}\n")
        cc+=1
    if counter == 0 and turn_com == True:
        ali_log.write(f"{' '*(int(sum([len(x) for x in l])/len(l)))}\t{dic_ali['conserv']}\n\n")

    #### Counting residues ####

    space=len(key)
    if counter == 0:
        ali_log.write(f"{' '*30} --- Number of residues in alignment ---\n\n")
        ali_log.write(f"{' '*space}{' '*z}\t\t\tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW\tY\tV\tZ\t\n")

    dic_composition={}
    cc=0
    for key in l:

        c=dic_ali[key]
        z=list_big[cc]
        if counter == 0:
            ali_log.write(f"{key}{' '*z}\t\t\t{c.count('A')}\t{c.count('R')}\t{c.count('N')}\t{c.count('D')}\t{c.count('C')}\t{c.count('Q')}\t{c.count('E')}\t{c.count('G')}\t{c.count('H')}\t{c.count('I')}\t{c.count('L')}\t{c.count('K')}\t{c.count('M')}\t{c.count('F')}\t{c.count('P')}\t{c.count('S')}\t{c.count('T')}\t{c.count('W')}\t{c.count('Y')}\t{c.count('V')}\t{c.count('Z')}\t\n")
        cc+=1

        ch=0
        hy=0
        dic_composition[key]=''

        for b in list(c):
            if any(d in ch_r for d in b): # Just titratable residues
                ch+=1
            if any(d in hy_r for d in b): # Just hydrophobic residues
                hy+=1
        dic_composition[key]+='charged='
        dic_composition[key]+=str(ch)
        dic_composition[key]+=' '
        dic_composition[key]+=str(ch_r)
        dic_composition[key]+=', hydrophobic='
        dic_composition[key]+=str(hy)
        dic_composition[key]+=' '
        dic_composition[key]+=str(hy_r)

    #### Observing composition of residues ####

    if counter == 0:
        ali_log.write(f"\n")
        ali_log.write(f"{' '*30} --- Kind of residues in alignment ---\n\n")
    cc=0
    for x,y in dic_composition.items():
        z=list_big[cc]

        if counter == 0:
            ali_log.write(f"{x}{' '*z} {y}\n")

        cc+=1
    if counter == 0:
        ali_log.write(f"\n")

    return dic_ali, l
