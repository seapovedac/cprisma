import pandas as pd
import numpy as np
import copy

def interpret_array(turn_com, tr_r, dic_ref, df_aa, df_da, df_bas_ref, df_bas_seq, df_array, feature, protein, df_bas_feature, accuracy, ali_log):

    index=df_da.index
    numrws=len(index)                   # Number of rows
    cc=0
    cc_col=cc
    dd=1
    ee=0
    hola=1

    df_auxs=pd.DataFrame()
    if feature == 'operation':
        # DataFrame for operations
        df_bas_ope=pd.DataFrame()
    elif feature == 'maximum':
        df_bas_ope=df_bas_feature
        # DataFrame for maximum restriction
        df_bas_mxr=pd.DataFrame(columns=['main','minimum (abs)','maximum (abs)','lenght'])
    elif feature == 'color':
        # It is passed a list of additional parameters
        list_op_col=df_bas_feature
        # DataFrame of operation
        df_bas_ope=list_op_col[2]
        df_bas_col=pd.DataFrame()

    list_t=[]
    type_list=type(list_t)
    str_t='hola'
    type_str=type(str_t)
    bool_t=True
    type_bool=type(bool_t)
    int_t=17
    type_int=type(int_t)
    float_t=17.0
    type_float=type(float_t)

    for key,list_dic in dic_ref.items():

        if turn_com:
            ali_log.write(f"\t-- Set of comparisons: {key} --\n\n")

        for d in list_dic:

            list_protein=[]    # List

            for rf_n,sq_n in d.items():

                v=len(sq_n)

                # Reference numeric data
                ref=df_bas_ref.iloc[ ee:numrws*dd , : ]
                ref=ref.dropna(axis=1, how='all')

                # Numeric data to compare
                if v != 0:                                      # If we have more than one sequence
                    seq=df_bas_seq.iloc[ ee:numrws*dd , : ]
                    seq=seq.dropna(axis=1, how='all')
                    aux=str(protein[int(rf_n)])+' X '
                    # List of protein involved as a reference with comparison
                    list_protein.append(protein[int(rf_n)])
                else:
                    seq=df_bas_seq.iloc[ ee:numrws*dd , : ]
                    seq=seq['null']
                    aux=str(protein[int(rf_n)])+' X no comparison '
                    # List of protein involved as a reference without comparison
                    list_protein.append(protein[int(rf_n)])

                for s in sq_n:
                    # List of protein involved as a target
                    list_protein.append(protein[s])
                    # To print later the target sequences
                    aux+=str(protein[s])+' '

                aux+='\n'

                # Data frame with residues of target sequences
                df_aux_protein=df_aa[list_protein]

                # Features

                if feature == 'operation':

                    from cprisma.data_array_relationship import operations

                    # According to a descriptor
                    df_bas_ope,df_auxs,descriptor=operations(ref,seq,sq_n,numrws,v,df_array,cc,df_auxs,df_bas_ope,ali_log,aux,turn_com)
                     # Creating DataFrame to separate columns with and without operations
                    df_separator=pd.DataFrame('->', index=range(numrws), columns=[descriptor])
                    df_bas_feature=copy.deepcopy(df_bas_ope)

                elif feature == 'maximum':

                    from cprisma.data_array_relationship import maximum_restriction_check

                    # Selecting data from operation DataFrame
                    seq_ope=df_bas_ope.iloc[ ee:numrws*dd , : ]
                    seq_ope=seq_ope.dropna(axis=1, how='all')

                    if turn_com:
                        ali_log.write(f"\t{aux}\n")

                    df_bas_mxr=maximum_restriction_check(key,df_bas_mxr,seq_ope,df_array,tr_r,type_bool,type_list,type_str,accuracy,rf_n,sq_n,numrws,cc,v,seq,protein,df_aux_protein,list_protein,turn_com,ali_log)
                    df_bas_feature=copy.deepcopy(df_bas_mxr)

                if feature == 'color':

                    from cprisma.data_array_relationship import color

                    if turn_com:
                        ali_log.write(f"\t{aux}\n")

                    df_bas_col,cc_col,hola=color(protein,hola,df_bas_col,list_op_col,df_array,df_aux_protein,df_bas_ope,tr_r,cc_col,dd,ee,numrws,sq_n,type_bool,type_int,type_float,type_list,accuracy,turn_com,ali_log)
                    df_bas_feature=copy.deepcopy(df_bas_col)

                # Visualization of comparison
                if turn_com:

                    if feature == 'operation':
                        # Removing column names
                        df_auxs.columns = [''] * len(df_auxs.columns)
                        if v != 0:
                            ssqq=[str(sq) for sq in sq_n]
                            seq1=seq[ssqq]
                        else:
                            seq1=copy.deepcopy(seq)
                        dfcomp=[ref, seq1, df_separator, df_auxs]
                        # Concatenating all df
                        dfcomp=(pd.concat(dfcomp, axis=1, sort=False))
                        dfbas=dfcomp.to_string()
                         # Saving in log file
                        ali_log.write(f"\t{aux}\n")
                        ali_log.write(f"{dfbas}\n\n")

                cc+=1
                #cc_col+=1
                dd+=1
                ee=numrws+ee

    if feature == 'maximum':
        df_bas_mxr=df_bas_mxr.to_string(index=False)
        ali_log.write(f"{df_bas_mxr}\n\n")

    ali_log.write(f"Matrices for {feature} were done!\n\n")
    print(f"Matrices for {feature} were done!")

    return df_bas_feature
