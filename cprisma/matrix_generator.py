from cprisma.arrays_default import array_default
from cprisma.print_array import *
from cprisma.alert_array import *
import pandas as pd
import numpy as np

    ########## Creating new dictionaries to relate reference and target sequences ##########

def matrix_generator(turn_com,ck_ref,dic_ref,met_ref,dic_da,df_da,protein,ali_log,counter):

    if counter == 0:
        ali_log.write(f"{' '*20}--- Generating arrays to compare ---\n\n")

    if ck_ref:                  # If you want to compare sequences wich method you prefer default, pair or multiple

        #### Adjusting array for default method ####

        if met_ref == 'default':
            dic_ref=array_default(protein)

        #### Verifying and creating array for pair method ####

        elif met_ref == 'pair':
            lenseq=len(protein)
            div_p=lenseq%2

            if div_p == 0:
                dic_ref={}
                dic_sub={}
                for i in range(lenseq):
                    cc=i%2
                    list_ref=[]
                    list_seq=[]

                    if cc == 0:
                        dic_sub[str(i)]=[]
                        join='main_'+str(i)
                        dic_ref[join]=[dic_sub]
                    else:
                        list_seq.append(i)
                        dic_sub[str(i-1)]=list_seq  # Adjusting array to pair mode
                        dic_sub={}
            else:
                met_ref1=met_ref
                met_ref='default'
                message_error1_in_matrix_generator_pair(met_ref1,lenseq,met_ref)
                dic_ref=array_default(protein)       # Adjusting array to default mode

        #### Verifying array for multiple method ####

        elif met_ref == 'multiple':

            if bool(dic_ref): # Checking if there is array
                detect=True

                ## Reading main key of dictionary ##
                for key,list_dic in dic_ref.items():

                    lk=len(key)
                    if lk == 0: # Checking if main key have string greater than 0
                        met_ref1=met_ref
                        met_ref='default'
                        message_error5_in_matrix_generator_multiple(met_ref1,met_ref,dic_ref)
                        dic_ref=array_default(protein)       # Adjusting array to default mode
                        break

                    ## Reading the list of dictionaries from main key ##
                    for d in list_dic:

                        ## Reading dictionaries with comparison between reference and target sequences ##
                        for x,y in d.items():

                            lenseq=len(protein)-1

                            try:                    # Checking if there are missing references in array
                                if int(x) > lenseq: # Checking if the reference sequence is greater than number of sequences analyzed
                                    detect=False
                                    no_exist=int(x)
                                    break
                                for z in y:
                                    if z > lenseq:  # Detecting if the target sequences are greater than number of sequences analyzed
                                        detect=False
                                        no_exist=z
                                        break
                            except ValueError:
                                met_ref1=met_ref
                                met_ref='default'
                                message_error4_in_matrix_generator_multiple(met_ref1,met_ref,key)
                                dic_ref=array_default(protein)       # Adjusting array to default mode
                                break

                if detect:
                    pass
                else:
                    met_ref1=met_ref
                    met_ref='default'
                    message_error3_in_matrix_generator_multiple(met_ref1,met_ref,no_exist,protein)
                    dic_ref=array_default(protein)       # Adjusting array to default mode
            else:
                met_ref1=met_ref
                met_ref='default'
                message_error2_in_matrix_generator_multiple(met_ref1,dic_ref,met_ref)
                dic_ref=array_default(protein)       # Adjusting array to default mode

        ################ Generating matrices for comparisons ################

        title='Method used '+"'"+met_ref+"'"+' for the array: '
        if counter == 0:
            print_array(title,dic_ref,ali_log)              # Function just to print the array separating the main keys

        df_bas_seq=pd.DataFrame()                           # Creating a new dataframe to add target sequences of interest

        for key,list_dic in dic_ref.items():
            df_auxs=pd.DataFrame()

            for d in list_dic:

                for rf_n,sq_n in d.items():
                    dic_trs1={}                             # We define a new dictionary for target sequences

                    for s in sq_n:                          # Generating dataframe for sequences to compare
                        # Generating new dictionary based on an element in the list of the dictionary
                        dic_trs1[str(s)]=[]
                        # Select column with data to compare in dataframe df_da
                        sl=df_da.iloc[0:, int(s)]
                        # Converting new dataframe to list
                        sl=sl.values.tolist()
                        # Adding the list with values in the dictionary
                        dic_trs1[str(s)]=sl
                        # Converting the dictionary to dataframe
                        dfs=pd.DataFrame.from_dict(dic_trs1)
                        # Creating array between df_auxs and dataframe created
                        df_auxs=[df_auxs, dfs]
                        # Concatenating the data frames
                        df_auxs=(pd.concat(df_auxs, axis=1, sort=False))
                        # Removing target sequences repeated
                        df_auxs=df_auxs.iloc[0:, -len(sq_n): ]
                    # Adding to the main or base dataframe df_base_seq
                    if len(sq_n) != 0:
                        df_bas_seq=pd.concat([df_bas_seq,df_auxs],sort=False)
                        #print(df_bas_seq)
                    else:    # If the reference sequence does not have any sequence to compare
                        index=df_da.index
                        dfnull = pd.DataFrame(np.nan, index=index,columns=['null'])
                        dfnull['null'] = dfnull['null'].replace(np.nan, 0)
                        df_bas_seq=pd.concat([df_bas_seq,dfnull],sort=False)

    ################ If it is not considered data comparison ################

    else:
        if counter == 0:
            ali_log.write(f"No comparison among the data.\n\n")
            print(f"No comparison among the data.")

        lenseq=len(protein)
        dic_ref={}
        list_seq=[]
        for i in range(lenseq):
            dic_sub={}
            dic_sub[str(i)]=[]
            join='main_'+str(i)
            dic_ref[join]=[dic_sub]
        df_bas_seq=pd.DataFrame()           # It is created an 'artificial' dataframe of sequences to compare just to avoid problems later

    ########## Defining reference sequences ##########

    dic_ref1={}
    df_bas_ref=pd.DataFrame()               # Creating a new dataframe to add reference sequences of interest

    for key,list_dic in dic_ref.items():

        for d in list_dic:

            for rf_n,sq_n in d.items():       # Generating dataframe for reference sequences similar process like before
                                                    # Note: This part is common for all methods
                if not ck_ref:                              # If there are not comparisons
                    if len(sq_n) == 0:                      # If the reference sequence does not have any sequence to compare
                        index=df_da.index
                        dfnull = pd.DataFrame(np.nan, index=index,columns=['null'])
                        dfnull['null'] = dfnull['null'].replace(np.nan, 0)
                        df_bas_seq=pd.concat([df_bas_seq,dfnull],sort=False)

                dic_ref1[str(rf_n)]=[]
                rl=df_da.iloc[0:, int(rf_n)]
                rl=rl.values.tolist()
                dic_ref1[str(rf_n)]=rl
                dfr=pd.DataFrame.from_dict(dic_ref1)
                dfr=dfr[str(rf_n)]
                dfr=pd.DataFrame(dfr)
                df_bas_ref=pd.concat([df_bas_ref,dfr],sort=False)

    #### Reading the new matrices of comparisons ####

    index=df_da.index
    numrws=len(index)                   # Number of rows
    cc=0
    dd=1
    ee=0

    if counter == 0:
        print('Comparisons: ')

    for key,list_dic in dic_ref.items():

        if counter == 0:
            ali_log.write(f"\t-- Set of comparisons: {key} --\n\n")

        for d in list_dic:

            for rf_n,sq_n in d.items():

                v=len(sq_n)

                ref=df_bas_ref.iloc[ ee:numrws*dd , : ]         # Selecting references
                ref=ref.dropna(axis=1, how='all')

                if v != 0:
                    seq=df_bas_seq.iloc[ ee:numrws*dd , : ]      # Selecting targets
                    seq=seq.dropna(axis=1, how='all')
                    aux=str(protein[int(rf_n)])+' X '
                else:
                    seq=df_bas_seq.iloc[ ee:numrws*dd , : ]
                    seq=seq['null']
                    aux=str(protein[int(rf_n)])+' X no comparison '

                for s in sq_n: aux+=str(protein[s])+' '

                if counter == 0:
                    print(f"\t{aux}")
                    ali_log.write(f"\t{aux}\n\n")

                aux+='\n'
                ssqq=[str(sq) for sq in sq_n]
                seq1=seq[ssqq]
                dfcomp=[ref, seq1]
                dfcomp=(pd.concat(dfcomp, axis=1, sort=False))        # Concatenating
                dfcompvi=dfcomp.to_string()

                cc+=1
                dd+=1
                ee=numrws+ee

    if counter == 0:
        ali_log.write(f"Matrices with data input were done!\n\n")
        print(f"Matrices with data input were done!")

    return df_bas_ref, df_bas_seq, dic_ref
