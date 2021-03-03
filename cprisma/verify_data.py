from cprisma.interpret_residue import interpret_residue
import pandas as pd

    ########## First verification ##########

def verify_data(ali_log,df_dat,protein,dic_ali,name_dat,bol,counter):

    is_NaN=(str(df_dat.isnull().iloc[0])).find('True') # To detect NaN in first row
    is_col=(len(df_dat.columns))%2 # To verify if number of columns are ok
    match1=(len(df_dat.columns))/2-len(protein) # Data and sequences should match

    if is_NaN < 0 and is_col == 0 and match1 == 0:
        if counter == 0:
            ali_log.write(f"Number of columns in {name_dat} is ok!\n\n")
            print(f"Number of columns in {name_dat} is ok!")
        bol=True
    else:
        ali_log.write(f"The program stop! Number of columns are not complete in {name_dat}.\n\n")
        print(f"The program stop! Number of columns are not complete in {name_dat}.")
        quit()

    return bol

    ##########  Second verification ##########

def compare_data(dic_aa,dic_da,protein,dic_ali,tr_r,ck_tr,name_dat,name_ali,bol,ali_log,counter):

    bol=False

    #### Converting all values in dic_aa and dic_da to a one letter aminoacid name ####

    cc=0
    for p in protein:
        l_aa=(list(dic_aa.values())[cc])
        for a in l_aa:
            am=interpret_residue(a)
            l_aa=[t.replace(str(a), am) for t in l_aa]
        dic_aa.update({list(dic_aa.keys())[cc]:l_aa})
        cc+=1

    #### Counting target residues in dic_ali ####

    if ck_tr:           # If exist a set of target residues specific
        if (len(tr_r)) == 0:
            ali_log.write(f"The program stop! The variable to check the target-residues is {ck_tr}, but the tuple is equal to 0.")
            print(f"The program stop! The variable to check the target-residues is {ck_tr}, but the tuple is equal to 0.")
            quit()
    else:
        tr_r=('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')

    dic_tar={}
    for t in tr_r: dic_tar[t]=''    # It's created a new dictionary for target residues
    cc=0
    for t in dic_tar:
        list_t=[]
        for p in protein:
            c=dic_ali[p]
            ct=c.count(t)
            list_t.append(ct)
        dic_tar.update({list(dic_tar.keys())[cc]:list_t})
        cc+=1
    #print(dic_tar)

    #### Comparing if target residues match with input data ####

    for x,y in dic_aa.items():
        for z in y:
            if z != 'X':
                if not any(d in tr_r for d in z):
                    ali_log.write(f"The program stop! The residue {z} [{x}] is appearing in {name_dat} but is missing in the target-residues: {tr_r}.")
                    print(f"The program stop! The residue {z} [{x}] is appearing in {name_dat} but is missing in the target-residues: {tr_r}.")
                    quit()

    #### Comparing data of alignment and input data ####

    cc=0
    for t in dic_tar:
        dd=0
        for value in dic_aa:
            cdat=list(dic_aa.values())[dd]
            cdat=cdat.count(t)
            cali=list(dic_tar.values())[cc][dd]
            if cali != cdat:
                ali_log.write(f"The program stop! Number of target-residues [{value}] is not the same in the input files. Residue {t}: {cdat} ({name_dat}) != {cali} ({name_ali}).")
                print(f"The program stop! Number of target-residues [{value}] is not the same in the input files. Residue {t}: {cdat} ({name_dat}) != {cali} ({name_ali}).")
                quit()
            dd+=1
        cc+=1

    bol=True

    #### Data frame to show correctly target residues count in log file ####

    dfp=pd.DataFrame((protein),columns =[""])
    dft=pd.DataFrame.from_dict(dic_tar)
    frames = [dfp, dft]
    frames=(pd.concat(frames, axis=1, sort=False))
    frames=frames.to_string(index=False)
    if counter == 0:
        ali_log.write(f"Target-residues {tr_r} match in {name_ali} and {name_dat}!\n\n")
        ali_log.write(f"{frames}\n\n")
        print(f"Target-residues {tr_r} match in {name_ali} and {name_dat}!")

    return bol, tr_r
