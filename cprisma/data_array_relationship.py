import pandas as pd
import numpy as np

    ##################################################################################################
    #                                    Dataframe for Operations                                    #
    ##################################################################################################

def operations(ref,seq,sq_n,numrws,v,df_array,cc,df_auxs,df_bas_ope,ali_log,aux,turn_com):

    a=np.array(ref)
    a=a.reshape(numrws,1)

    if v != 0:
        ssqq=[str(sq) for sq in sq_n]
        seq=seq[ssqq]
        b=np.array(seq)
    else:
        b=np.array(seq)
        b=b.reshape(numrws,1)

    descriptor=df_array['descriptor'].iloc[cc]
                            # Broadcasting
    # None
    if descriptor == 'n':
        if len(sq_n) != 0:
            z=b
        else:
            z=a
    # Delta
    elif descriptor == 'd':
        if len(sq_n) != 0:
            z=b-a
        else:
            ali_log.write(f"Attention!!! Because next comparison does not have any target data, the operation 'd' will not be applied.\n\n")
            if turn_com == False:
                ali_log.write(f"\t{aux}\n")
            z=a
    # Absolute delta
    elif descriptor == 'da':
        z=abs(b-a)
    # Multiply
    elif descriptor == 'm':
        if len(sq_n) != 0:
            z=b*a
        else:
            ali_log.write(f"Attention!!! Because next comparison does not have any target data, it was assumed that the reference data should be multiply by 1 for the descriptor 'm'.\n\n")
            if turn_com == False:
                ali_log.write(f"\t{aux}\n")
            z=a


    df_auxs=pd.DataFrame(z)
    #print(seq,df_auxs)
    df_bas_ope=pd.concat([df_bas_ope,df_auxs],sort=False)

    return df_bas_ope, df_auxs, descriptor

    ##################################################################################################
    #                                    Dataframe for Maximum 1                                     #
    ##################################################################################################

def maximum_restriction(cc,v,seq,df_aux_protein,protein,list_protein,seq_ope,df_array,rf_n,sq_n,type_bool,type_list,type_str,accuracy,ali_log):

    from cprisma.alert_array import message1_error_in_data_array_relationship
    import copy
    import math

    # To follow the same header order of the ref and seq df we should do some fittings
    # Fitting the order of the columns through their numbers for target sequence
    try:
        # Converting headers of DataFrame to a list
        list_column_num=seq.columns.tolist()
    except AttributeError:
        list_column_num=['null']

    # Creating list of name from list of numbers
    list_column_name=[]

    # Adding the header of the reference first to the list_column_name
    list_column_name.append(list(df_aux_protein.columns.tolist())[0])

    for lcn in list_column_num:
        if lcn != 'null':
            list_column_name.append(protein[int(lcn)])
        else:
            list_column_name=['null']

    # Changing the order of the columns to the correct one
    if list_column_name[0] != 'null':
        df_aux_protein=df_aux_protein[list_column_name]

    dfcomp=[df_aux_protein, seq_ope]
    # Concatenating values of operations and residues
    dfcomp=(pd.concat(dfcomp, axis=1, sort=False))

    # Selecting the comparison (it can be more than one subcomparisons)
    sub_array=df_array.loc[df_array['comparison'] == cc+1] # cc+1 to select correctly the comparison

    #print(sub_array)

    # Index of df_array for one specific comparison
    indx_des=sub_array.index
    indx_des=len(indx_des)
    list_sub_arr=[]
    list_max=[]

    # Loop to see the subcomparisons and calculate the maximum restrictions
    for ii in range (indx_des):

        descriptor=sub_array.iloc[ii]
        flag_des=descriptor['descriptor']

        # nr
        if flag_des == 'nr' or flag_des == 'nrY':
            dfcomp=dfcomp.iloc[ : , len(list_protein): ]
            dfcomp=np.transpose(dfcomp)
            dfcomp=dfcomp.values.tolist()
            list_max=[sj for si in dfcomp for sj in si]

        # rsa, rsam, rspi, rspim
        if flag_des == 'rsa' or flag_des == 'rsam' or flag_des == 'rspi' or flag_des == 'rspim' or flag_des == 'rsaY' or flag_des == 'rsamY' or flag_des == 'rspiY' or flag_des == 'rspimY':

            sequence=descriptor['sequence']
            sequence=int(sequence)
            col_pos=dfcomp.columns.get_loc(protein[sequence])
            col_pos=v+col_pos
            # Getting the residues of the reference
            sub_res_ref=dfcomp[protein[int(rf_n)]]
            # Getting the residues of the target squence
            sub_res_seq=dfcomp[protein[sequence]]
            # Position of the target sequence reported in the descriptor of the dict_array
            sub_ope_seq=dfcomp.iloc[ : , col_pos ]
            # Concatenating everything
            sub_dfcomp=[sub_res_ref, sub_res_seq, sub_ope_seq]
            sub_dfcomp=(pd.concat(sub_dfcomp, axis=1, sort=False))

            # Removing data when mutation happen
            mutation=descriptor['mutation']
            if type(mutation) == type_bool:
                # rsam and rspim
                sub_dfcomp=sub_dfcomp[sub_dfcomp[protein[int(rf_n)]].eq(sub_dfcomp[protein[sequence]])]

            # Removing specific kind of residues
            if flag_des == 'rsa' or flag_des == 'rsam' or flag_des == 'rsaY' or flag_des == 'rsamY':
                residues=descriptor['residue']
                # Removing residues based in the list of the descriptor
                sub_dfcomp=sub_dfcomp[~sub_dfcomp[:].isin(residues)]
                sub_dfcomp=sub_dfcomp.dropna()
                # Selecting just numerical values
                sub_dfcomp=sub_dfcomp.iloc[ : , -1 : ]

            # Removing positions
            if flag_des == 'rspi' or flag_des == 'rspim' or flag_des == 'rspiY' or flag_des == 'rspimY':
                position=descriptor['position']
                # Selecting just the numerical values
                sub_dfcomp=sub_dfcomp.iloc[ : , -1 : ]
                # If rspim is probable that my index change so we should evaluate and update the position in my descriptor
                index_sub_dfcomp=sub_dfcomp.index
                position=set(position).intersection(set(index_sub_dfcomp))
                # Removing the position based on index
                sub_dfcomp=sub_dfcomp.drop(position)

            # Transposing the new data set
            sub_dfcomp=np.transpose(sub_dfcomp)
            # Converting to a list
            list_sub_dfcomp=sub_dfcomp.values.tolist()
            # More adjustements
            list_sub_dfcomp=[j for i in list_sub_dfcomp for j in i]
            # Appending to list_sub_arr
            list_sub_arr.append(list_sub_dfcomp)

            # We should add the sequences that are not part of df_array
            if ii == indx_des-1:
                # Select sequences not declared in df_array for rsa, rsam, rspi, rspim... etc
                list_fast=[]
                # Creating list with the sequences of df_array
                for iii in range (indx_des):
                    descriptor=sub_array.iloc[iii]
                    sequence=descriptor['sequence']
                    list_fast.append(int(sequence))
                # Getting just the sequences that are not in df_array
                no_rs=list(set(sq_n) - set(list_fast))

                # Reading the new list
                for nrs in no_rs:
                    col_pos=dfcomp.columns.get_loc(protein[nrs])
                    col_pos=v+col_pos
                    sub_ope_seq_aux=dfcomp.iloc[ : , col_pos ]
                    sub_dfcomp=[sub_ope_seq_aux]
                    sub_dfcomp=(pd.concat(sub_dfcomp, axis=1, sort=False))
                    sub_dfcomp=sub_dfcomp.iloc[ : , -1 : ]
                    # Transposing the new data set
                    sub_dfcomp=np.transpose(sub_dfcomp)
                    # Converting to a list
                    list_sub_dfcomp=sub_dfcomp.values.tolist()
                    # More adjustements
                    list_sub_dfcomp=[j for i in list_sub_dfcomp for j in i]
                    # Appending to list_sub_arr
                    list_sub_arr.append(list_sub_dfcomp)

                ####

        # rm, ra, ram, rs, rsm
        if flag_des == 'rm' or flag_des == 'ra' or flag_des == 'ram' or flag_des == 'rsm' or flag_des == 'rs' or flag_des == 'rmY' or flag_des == 'raY' or flag_des == 'ramY' or flag_des == 'rsmY' or flag_des == 'rsY':

            # Removing data when mutation happen
            mutation=descriptor['mutation']

            if type(mutation) == type_bool:
                ref_name=list_column_name[0]
                for seq_name in (list_column_name):
                    # rm, ram, rsm
                    if seq_name != 'null':
                        dfcomp=dfcomp[dfcomp[ref_name].eq(dfcomp[seq_name])]

            # Removing specific kind of residues
            if flag_des == 'ra' or flag_des == 'ram' or flag_des == 'raY' or flag_des == 'ramY':

                residues=descriptor['residue']
                dfcomp=dfcomp[~dfcomp[:].isin(residues)]
                dfcomp=dfcomp.dropna()

            # Removing specific sequence
            if flag_des == 'rs' or flag_des == 'rsm' or flag_des == 'rsY' or flag_des == 'rsmY':

                dfcomp_aux=copy.deepcopy(dfcomp)
                sequence=descriptor['sequence']

                list_aux_seq=[]
                # If reference is in the list of descriptor it is removed this element
                for csq in sequence:
                    if csq != int(rf_n):
                        list_aux_seq.append(csq)

                sequence=list_aux_seq
                list_0=[]
                # Just target squences will be considered!
                if sequence != list_0:

                    # Sorting number list of sequences in the descriptor respect to list of the array of comparison.
                    sequence=sorted(sequence, key=lambda x: sq_n.index(x))
                    scc=0

                    for dsq in sequence:

                        try:                    # If the user put the sequence by duplicate
                            # Checking exact position of the protein to remove
                            col_pos1=dfcomp_aux.columns.get_loc(protein[dsq])
                            # Using as a reference col_pos1 to select the numeric data of that protein
                            col_pos2=v+col_pos1
                            # Excluding that column
                            dfcomp1=dfcomp_aux.iloc[  : , 0 : col_pos1  ]
                            dfcomp2=dfcomp_aux.iloc[  : , col_pos1+1 : col_pos2-scc ]
                            dfcomp3=dfcomp_aux.iloc[  : , col_pos2+1-scc :  ]
                            dfcomp_aux=[dfcomp1, dfcomp2, dfcomp3]
                            dfcomp_aux=(pd.concat(dfcomp_aux, axis=1, sort=False))
                        except KeyError:
                            pass
                        scc+=1
                    # Fitting data to have the numerical data of interest
                    num=((list(dfcomp_aux.shape)[1]-1)/2)+1
                    dfcomp=dfcomp_aux
                    dfcomp=dfcomp.iloc[ : , int(num) : ]
                    dfcomp=np.transpose(dfcomp)
                    dfcomp=dfcomp.values.tolist()
                    list_max=[sj for si in dfcomp for sj in si]

                else:   # if the the list sequence of the descriptor does not have any element
                    message1_error_in_data_array_relationship(ali_log,flag_des,rf_n)
                    list_max=[]
                    pass

            if flag_des == 'rm' or flag_des == 'ra' or flag_des == 'ram' or flag_des == 'rmY' or flag_des == 'raY' or flag_des == 'ramY':
                dfcomp=dfcomp.iloc[ : , len(list_protein): ]
                dfcomp=np.transpose(dfcomp)
                dfcomp=dfcomp.values.tolist()
                list_max=[sj for si in dfcomp for sj in si]

        # rpi
        if flag_des == 'rpi' or flag_des == 'rpiY':

            position=descriptor['position']

            if type(position) == type_list:
                # position with list
                rpi=float(list(position)[0])
            else:
                # position with NaN
                rpi=position
            if math.isnan(rpi): position=[]

            dfcomp=dfcomp.iloc[ : , len(list_protein): ]
            dfcomp=dfcomp.drop(position)
            dfcomp=np.transpose(dfcomp)
            dfcomp=dfcomp.values.tolist()
            list_max=[sj for si in dfcomp for sj in si]

        # rt
        if flag_des == 'rt' or flag_des == 'rtY':

            def threshold_zeros(dfcomp,dfcomp_bool):

                # Converting based on threshold to zero
                for key, value in dfcomp_bool.iteritems():
                    idfb=0
                    for vvv in value:
                        if vvv == False:
                            value_tr=(dfcomp[key].iloc[ idfb ])
                            dfcomp=(dfcomp.replace({ key : float(value_tr)}, 0.0))
                        idfb+=1

                return dfcomp

            # selecting just the values numerical values
            dfcomp=dfcomp.iloc[ : , len(list_protein): ]
            threshold=descriptor['threshold']

            try:
                rt1=threshold[1]
                rt2=threshold[0]
                rt2=round(rt2,accuracy)
            except TypeError:
                rt1='no_rt'
            #print(dfcomp)

            if rt1 == '==':
                # Creating a df of Booleans
                dfcomp_bool=dfcomp.select_dtypes(include=['number']) == float(rt2)
                dfcomp=threshold_zeros(dfcomp,dfcomp_bool)

            elif rt1 == '<=':
                # Creating a df of Booleans
                dfcomp_bool=dfcomp.select_dtypes(include=['number']) <= float(rt2)
                dfcomp=threshold_zeros(dfcomp,dfcomp_bool)

            elif rt1 == '>=':
                # Creating a df of Booleans
                dfcomp_bool=dfcomp.select_dtypes(include=['number']) >= float(rt2)
                dfcomp=threshold_zeros(dfcomp,dfcomp_bool)

            elif rt1 == '<':
                # Creating a df of Booleans
                dfcomp_bool=dfcomp.select_dtypes(include=['number']) < float(rt2)
                dfcomp=threshold_zeros(dfcomp,dfcomp_bool)

            elif rt1 == '>':
                # Creating a df of Booleans
                dfcomp_bool=dfcomp.select_dtypes(include=['number']) > float(rt2)
                dfcomp=threshold_zeros(dfcomp,dfcomp_bool)

            elif rt1 == '!=':
                # Creating a df of Booleans
                dfcomp_bool=dfcomp.select_dtypes(include=['number']) != float(rt2)
                dfcomp=threshold_zeros(dfcomp,dfcomp_bool)

            dfcomp=np.transpose(dfcomp)
            dfcomp=dfcomp.values.tolist()
            list_max=[sj for si in dfcomp for sj in si]

    list_0=[]
    if list_sub_arr != list_0: # If we have sub_arrays
        list_max=[sj for si in list_sub_arr for sj in si]

    list_max_no_abs=copy.deepcopy(list_max)
    list_max_no_abs=[round(item,accuracy) for item in list_max_no_abs] # Applying abs function

    # Finding the maximum
    if len(list_max) != 0:
        # Because there is not "negative" color we should be equivalent with negative numbers
        list_max=[abs(item) for item in list_max] # Applying abs function
        max_restriction=max(list_max)
        if max_restriction == 0:
            max_restriction=round(1.0,accuracy)
        minimum=min(list_max)
    else:
        max_restriction=round(1.0,accuracy)
        minimum=round(0.0,accuracy)
    list_max=[round(float(item),accuracy) for item in list_max]

    return list_max, list_max_no_abs, max_restriction, minimum

    ##################################################################################################
    #                                    Dataframe for Maximum 2                                     #
    ##################################################################################################

def maximum_restriction_check(key,df_bas_mxr,seq_ope,df_array,tr_r,type_bool,type_list,type_str,accuracy,rf_n,sq_n,numrws,cc,v,seq,protein,df_aux_protein,list_protein,turn_com,ali_log):

    import copy

    list_max=[]
    max_restriction=[]
    minimum=[]

    sub_array=df_array.loc[df_array['comparison'] == cc+1] # cc+1 to select correctly the comparison
    # Index of df_array for one specific comparison
    indx_des=sub_array.index
    indx_des=len(indx_des)

    for ii in range (indx_des):
        separate=sub_array['separate'].iloc[ii]

    # If it's separated the maximums of each target residue
    if type(separate) == type_bool:

        # For case with comparisons
        if len(sq_n) != 0:

            # loop of target_residues
            for ttrr in tr_r:
                seq_ope_res=copy.deepcopy(seq_ope)
                # loop for columns
                for ccisq in range(len(sq_n)):
                    # loop for rows
                    for isq in range(numrws):
                        # getting residue to evaluate
                        res_tar=df_aux_protein[protein[sq_n[ccisq]]].iloc[isq]
                        if res_tar == ttrr:
                            pass
                        else:
                            seq_ope_res.loc[  isq , ccisq]=0.0
                list_max_tr,list_max_no_abs,max_restriction_tr,minimum_tr=maximum_restriction(cc,v,seq,df_aux_protein,protein,list_protein,seq_ope_res,df_array,rf_n,sq_n,type_bool,type_list,type_str,accuracy,ali_log)
                max_restriction.append(round(max_restriction_tr,accuracy))
                minimum.append(round(minimum_tr,accuracy))
                if turn_com:
                    ali_log.write(f'list_max[{ttrr}] = {list_max_no_abs}\n\n')
            df_bas_mxr=df_bas_mxr.append({'main' : key , 'minimum (abs)': minimum, 'maximum (abs)' : max_restriction}, ignore_index=True)

        # For case without comparisons
        else:
            # loop of target_residues
            for ttrr in tr_r:
                seq_ope_res=copy.deepcopy(seq_ope)
                for isq in range(numrws):
                    # getting residue to evaluate
                    res_tar=df_aux_protein[protein[int(rf_n)]].iloc[isq]
                    if res_tar == ttrr:
                        pass
                    else:
                        seq_ope_res.loc[ isq ]=0.0
                list_max_tr,list_max_no_abs,max_restriction_tr,minimum_tr=maximum_restriction(cc,v,seq,df_aux_protein,protein,list_protein,seq_ope_res,df_array,rf_n,sq_n,type_bool,type_list,type_str,accuracy,ali_log)
                max_restriction.append(round(max_restriction_tr,accuracy))
                minimum.append(round(minimum_tr,accuracy))
                if turn_com:
                    ali_log.write(f'list_max[{ttrr}] = {list_max_no_abs}\n\n')
            df_bas_mxr=df_bas_mxr.append({'main' : key , 'minimum (abs)': minimum, 'maximum (abs)' : max_restriction}, ignore_index=True)

    # If it's not separated the maximums of each target residue
    else:
        list_max,list_max_no_abs,max_restriction,minimum=maximum_restriction(cc,v,seq,df_aux_protein,protein,list_protein,seq_ope,df_array,rf_n,sq_n,type_bool,type_list,type_str,accuracy,ali_log)
        df_bas_mxr=df_bas_mxr.append({'main' : key , 'minimum (abs)': minimum, 'maximum (abs)' : round(max_restriction,accuracy), 'lenght': len(list_max)}, ignore_index=True)
        if turn_com:
            ali_log.write(f'list_max = {list_max_no_abs}\n\n')

    return df_bas_mxr

    ##################################################################################################
    #                                      Dataframe for Color                                       #
    ##################################################################################################

def color(protein,hola,df_bas_col,list_op_col,df_array,df_aux_protein,df_bas_ope,tr_r,cc,dd,ee,numrws,sq_n,type_bool,type_int,type_float,type_list,accuracy,turn_com,ali_log):

    import random
    import copy
    import os

    def lenght_str_col(color):

        len_color=len(color)
        if len_color < 17:
            len_color=17-len_color
            color+=' '*len_color

        return color

    def random_color(target_col,list_ran_col,list_colors,elec):

        while target_col in list_ran_col or target_col in [116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142]:

            if elec != 0:
                target_col=random.choice(list_colors)

        list_ran_col.append(target_col)

        return target_col, list_ran_col

    def color_residues(residue_col,colors,df_auxs,name_col):

        color=colors['RGB'].iloc[residue_col]
        color=lenght_str_col(color)
        df_auxs=df_auxs.append({name_col:color} , ignore_index=True)

        return df_auxs

    def color_pka(column_col_ope,rcope,colors,df_auxs,name_col):

        if column_col_ope[rcope] > 0:
            pka_col=94
            df_auxs=color_residues(pka_col,colors,df_auxs,name_col)
        # Positive pKa
        elif column_col_ope[rcope] < 0:
            pka_col=6
            df_auxs=color_residues(pka_col,colors,df_auxs,name_col)
        # Neutral pKa
        elif column_col_ope[rcope] == 0:
            pka_col=116
            df_auxs=color_residues(pka_col,colors,df_auxs,name_col)

        return df_auxs

    def call_var_pic(df_array,row_pi_res,ccpi):

        sequence=df_array['sequence'].iloc[row_pi_res]
        # Choosing list of residues
        pic_res=df_array['residue'].iloc[row_pi_res]
        pic_res=pic_res[1]
        # Color
        sequence_col=df_array['residue'].iloc[row_pi_res][0]
        color=colors['RGB'].iloc[sequence_col]
        # Checking the lenght of that list
        len_pic_res=len(pic_res)
        # Choosing the elements of that list based on ccpi counter
        in_pic_res=pic_res[ccpi]
        # Splitting the interval
        spic=in_pic_res.split("-")

        return sequence, sequence_col, color, spic, len_pic_res


    def add_list_pic(df_array,row_pi_res,tar_col_num):

        seg_col_pic=df_array['residue'].iloc[row_pi_res]
        seg_col_pic=seg_col_pic[0]
        tar_col_num.append(seg_col_pic)

        return tar_col_num

    def final_pic(col_c, num_col_seq,rcope,numrws,descriptor,tar_col_num,mutation_col,turn_com,ali_log):

        if col_c == num_col_seq-1 and rcope == numrws-1:
            if descriptor == 'pimc':
                tar_col_num.append(mutation_col)
            if turn_com:
                ali_log.write('')

    # Common variables

    tar_col_num=[]

    # Colors based on https://htmlcolorcodes.com/color-names/
    dir_run=os.path.dirname(os.path.abspath(__file__))
    color_dir=dir_run+"/colors.csv"
    name_color_file=color_dir
    colors=pd.read_csv(name_color_file)
    index_colors=colors.index
    list_colors=list(index_colors)
    color_neutral=colors['RGB'].iloc[116]

    mutation_col=list_op_col[0]
    sequence_col=list_op_col[1]

    # Selecting the descriptor
    row_pi_res=cc # Variable for pic and pimc
    descriptor=df_array['descriptor'].iloc[row_pi_res]
    # Number of columns in df_aux_protein
    num_col_seq=df_aux_protein.shape[1]
    df_auxs_bas=pd.DataFrame()

    # Verifyng if there are mutations
    mutated_res=df_array['mutation'].iloc[row_pi_res]

    if mutated_res == True:
        if mutation_col > (len(index_colors)-1):
            from cprisma.alert_array import message2_error_in_data_array_relationship
            message2_error_in_data_array_relationship(ali_log,descriptor,mutation_col)
            mutation_col=137

    if descriptor == 'ssc':
        if sequence_col > (len(index_colors)-1):
            from cprisma.alert_array import message3_error_in_data_array_relationship
            message3_error_in_data_array_relationship(ali_log,descriptor,sequence_col)
            sequence_col=6

    ######################## Selecting color for residue methods ########################

    residues=df_array['residue'].iloc[row_pi_res]

    # nc
    if descriptor == 'nc':

        signal=0
        comparison=df_array['comparison'].iloc[row_pi_res]

        color=color_neutral
        if turn_com:
            ali_log.write(f'Descriptor: {descriptor} - Color (sequence): White - all\n')
            ali_log.write('\n')
        tar_col_num=[116]

    if descriptor == 'pkac':

        signal=0
        comparison=df_array['comparison'].iloc[row_pi_res]

        # Selecting data from operation DataFrame to later use as a base to put color
        seq_ope=df_bas_ope.iloc[ ee:numrws*dd , : ]
        seq_ope=seq_ope.dropna(axis=1, how='all')

        if turn_com:
            if mutated_res == False:
                ali_log.write(f'Descriptor: {descriptor} - Color (residue): Red - acid\n')
                ali_log.write(f'Descriptor: {descriptor} - Color (residue): Blue - basic\n')
            else:
                ali_log.write(f'Descriptor: {descriptor} - Color (residue): Red - acid\n')
                ali_log.write(f'Descriptor: {descriptor} - Color (residue): Blue - basic\n')
                ali_log.write(f'Descriptor: {descriptor} - Color (residue): Yellow - acid mutation\n')
                ali_log.write(f'Descriptor: {descriptor} - Color (residue): Green - basic mutation\n')
            ali_log.write('\n')

        tar_col_num=[94,6,54,22]

    # Defining the method for color residue: 'fac', 'fmac'
    if descriptor == 'fac' or descriptor == 'fmac':

        signal=0
        comparison=df_array['comparison'].iloc[row_pi_res]

        # fac, fmac
        if residues[0] == 'random': # If it's selected to build a list with random colors
            list_ran_residues=[]
            elec=0 # Variable to select one time the first random color in residue_rcol
            for ttrr in range(len(tr_r)):
                if elec == 0:
                    residue_rcol=random.choice(list_colors)
                    elec=1
                residue_rcol,list_ran_residues=random_color(residue_rcol,list_ran_residues,list_colors,elec)
            residues=list_ran_residues
            goal_res=tr_r

            tar_col_num=copy.deepcopy(residues)
            if descriptor == 'fmac':
                tar_col_num.append(mutation_col)

        elif type(residues[0])== type_int:  # else it's defined a list of colors
            residues=[xc for xc in residues]
            goal_res=tr_r
            tar_col_num=copy.deepcopy(residues)
            if descriptor == 'fmac':
                tar_col_num.append(mutation_col)

        if turn_com:
            for col_name in range(len(residues)):
                rcn=residues[col_name]
                rcn_tr=tr_r[col_name]
                rcn=colors['Name'].iloc[rcn]
                ali_log.write(f'Descriptor: {descriptor} - Color (residue): {rcn} - {rcn_tr}\n')
            if mutated_res == True:
                rcn=colors['Name'].iloc[mutation_col]
                ali_log.write(f'Descriptor: {descriptor} - Color (residue): {rcn} - mutation\n')
            ali_log.write('\n')

    ###########################################################################
    #                         Starting reading...                             #
    ###########################################################################

    col_c=0
    list_ran_col=[]

    #### Reading column ####

    while col_c < num_col_seq:

    ######################## Selecting color for sequence methods ########################

        # Defining the method for color sequence: ssc, fsc, pic, pimc

        try:
            sequence=df_array['sequence'].iloc[row_pi_res]
        except IndexError:
            pass

        elec=0 # Variable to select one time the first random color in sequence_col
        shift=0

        if type(sequence) != type_float:

            # ssc
            if sequence[0] == 'same' and descriptor == 'ssc':

                signal=0
                comparison=df_array['comparison'].iloc[row_pi_res]
                sequence=df_array['sequence'].iloc[row_pi_res]
                color=colors['RGB'].iloc[sequence_col]
                # Function to have simmetry in the spaces
                color=lenght_str_col(color)
                tar_col_num=[sequence_col]

            # fsc['random']
            elif sequence[0] == 'random' and descriptor == 'fsc':

                signal=0
                comparison=df_array['comparison'].iloc[row_pi_res]

                if elec == 0:
                    # Giving the first random color to a given sequence
                    sequence_col=random.choice(list_colors)
                    elec=1
                # Function to do not repeat colors or select white color
                sequence_col,list_ran_col=random_color(sequence_col,list_ran_col,list_colors,elec)
                # Selecting the color based on the file colors.csv
                color=colors['RGB'].iloc[sequence_col]
                # Function to have simmetry in the spaces
                color=lenght_str_col(color)
                tar_col_num=list_ran_col

            # fsc with list and pic/pimc
            elif type(sequence[0])== type_int:

                if descriptor == 'fsc':

                    signal=0
                    comparison=df_array['comparison'].iloc[row_pi_res]

                    if col_c != 0:
                        sequence_col=sequence[col_c-1]
                        color=colors['RGB'].iloc[sequence_col]
                        color=lenght_str_col(color)
                    else: # Auxiliar
                        #sequence_col=sequence[col_c-1]
                        color=colors['RGB'].iloc[116]
                        color=lenght_str_col(color)

                    color_neutral=colors['RGB'].iloc[116]
                    sss=[xc for xc in sequence]#
                    sss = [116] + sss # It is added a auxiliar color representing the color of the ref
                    tar_col_num=sss
                    ccpi=0

                elif  descriptor == 'pic' or  descriptor == 'pimc':       # pic or pimc

                    signal=0

                    aux_tar_col_num=[]
                    try:
                        sequence=df_array['sequence'].iloc[row_pi_res]
                        descriptor=df_array['descriptor'].iloc[row_pi_res]
                        comparison=df_array['comparison'].iloc[row_pi_res]
                    except IndexError:
                        pass
                    color_neutral=colors['RGB'].iloc[116]
                    ccpi=0

                elif descriptor == 'tc' or descriptor == 'tmc':           # tc or tmc

                    signal=0
                    # Color parameters
                    comparison=df_array['comparison'].iloc[row_pi_res]
                    sequence=df_array['sequence'].iloc[row_pi_res]
                    sequence_col=sequence[0]
                    tar_col_num=[sequence_col]

                    if descriptor == 'tmc':
                        tar_col_num.append(mutation_col)

                    # Threshold parameters
                    threshold_number=float(df_array['threshold'].iloc[row_pi_res][0])
                    threshold_symbol=df_array['threshold'].iloc[row_pi_res][1]

                    # Operation parameters
                    seq_ope=df_bas_ope.iloc[ ee:numrws*dd , : ]
                    seq_ope=seq_ope.dropna(axis=1, how='all')

        #######################################################################################

        df_auxs=pd.DataFrame()

        # Reading column
        column_col=df_aux_protein.iloc[ : , col_c ]
        column_col_ref=df_aux_protein.iloc[ : , 0 ]

        # Getting the name of the column
        name_col=list(pd.DataFrame(column_col).to_dict().keys())[0]

        # Sequence color in log file
        if turn_com:

            if type(sequence) == type_list and descriptor == 'ssc' or descriptor == 'fsc' or descriptor == 'tc' or descriptor == 'tmc' and comparison == hola:

                if col_c != 0:
                    scn=colors['Name'].iloc[sequence_col]
                    ali_log.write(f'Descriptor: {descriptor} - Color (sequence): {scn} - {name_col}\n')

                    if descriptor == 'tmc' and col_c == num_col_seq-1: # Mutation
                        rcn=colors['Name'].iloc[mutation_col]
                        ali_log.write(f'Descriptor: {descriptor} - Color: {rcn} - mutation\n')

                    if col_c == num_col_seq-1:
                        ali_log.write('\n')

                else: # Case of no comparison but that need color

                    if col_c != 0 or col_c == 0 and len(sq_n) == 0:
                        scn=colors['Name'].iloc[sequence_col]
                        ali_log.write(f'Descriptor: {descriptor} - Color (sequence): {scn} - {name_col}\n')

                        if descriptor == 'tmc' and col_c == num_col_seq-1: # Mutation
                            rcn=colors['Name'].iloc[mutation_col]
                            ali_log.write(f'Descriptor: {descriptor} - Color: {rcn} - mutation\n\n')
                        else:
                            ali_log.write('\n')


            elif type(sequence) == type_list and descriptor == 'ssc' or descriptor == 'fsc' or descriptor == 'tc' or descriptor == 'tmc' and comparison != hola:
                if col_c != 0 and col_c == num_col_seq-1:
                    ali_log.write('\n')

        rcref=0 # counter for reference residues
        rcope=0 # counter for values of operation df (only for pkac descriptor)

        #### Reading rows ####

        for row_col in column_col:

            # Putting color to target sequence or if reference does not have target sequence
            if col_c != 0 or col_c == 0 and len(sq_n) == 0:

                # Putting colors to target residues!
                if row_col != 'X': # The residue should not be a gap

                    # nc
                    if descriptor == 'nc':
                        shift=1
                        df_auxs=df_auxs.append({name_col:color} , ignore_index=True)

                    # 'fac', 'fmac'
                    elif descriptor == 'fac' or descriptor == 'fmac':
                        shift=1
                        # fac
                        if type(mutated_res) != type_bool:
                            residue_col=residues[goal_res.index(row_col)]
                            df_auxs=color_residues(residue_col,colors,df_auxs,name_col)
                        # fmac
                        else:
                            res_ref=column_col_ref[rcref] # Residue in reference
                            if row_col == res_ref: # If it is the same like in reference
                                residue_col=residues[goal_res.index(row_col)]
                                df_auxs=color_residues(residue_col,colors,df_auxs,name_col)
                            else: # If it is a mutation
                                residue_col=mutation_col
                                df_auxs=color_residues(residue_col,colors,df_auxs,name_col)

                    # 'pkac'
                    elif descriptor == 'pkac':
                        shift=1
                        # Getting values of operation df
                        column_col_ope=seq_ope.iloc[: , col_c-1]

                        if mutated_res == False:
                            df_auxs=color_pka(column_col_ope,rcope,colors,df_auxs,name_col)
                        else:
                            res_ref=column_col_ref[rcref] # Residue in reference
                            if row_col == res_ref: # If it is the same like in reference
                                df_auxs=color_pka(column_col_ope,rcope,colors,df_auxs,name_col)
                            else: # Mutations acids or basics
                                if row_col == 'D' or row_col == 'E':
                                    mutation_col=22 # Acid one
                                elif row_col == 'R' or row_col == 'H' or row_col == 'K' or row_col == 'Y' or row_col == 'C':
                                    mutation_col=54 # Basic one
                                elif row_col == 'A' or row_col == 'N' or row_col == 'Q' or row_col == 'G' or row_col == 'I' or row_col == 'L' or row_col == 'M' or row_col == 'F' or row_col == 'P' or row_col == 'S' or row_col == 'T' or row_col == 'W' or row_col == 'V' or row_col == 'Z' or row_col == 'X' or row_col == '0':
                                    mutation_col=116 # Basic one
                                pka_col=mutation_col
                                df_auxs=color_residues(pka_col,colors,df_auxs,name_col)

                    # tc, tmc
                    if descriptor == 'tc' or descriptor == 'tmc' and comparison == hola:

                        shift=1

                        if mutated_res != True:
                            color=colors['RGB'].iloc[sequence_col]
                            color=lenght_str_col(color)

                        else:

                            res_ref=column_col_ref[rcref] # Residue in reference

                            if row_col == res_ref: # If it is the same like in reference
                                color=colors['RGB'].iloc[sequence_col]
                                color=lenght_str_col(color)
                            else:                  # Mutations
                                color=colors['RGB'].iloc[mutation_col]
                                color=lenght_str_col(color)

                        # Getting values of operation df
                        column_col_ope=seq_ope.iloc[: , col_c-1]
                        tope=round(column_col_ope[rcope], accuracy)

                        if threshold_symbol == '==':
                            if tope == threshold_number:
                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                            else:
                                df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                        if threshold_symbol == '!=':
                            if tope != threshold_number:
                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                            else:
                                df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                        if threshold_symbol == '>=':
                            if tope >= threshold_number:
                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                            else:
                                df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                        if threshold_symbol == '<=':
                            if tope <= threshold_number:
                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                            else:
                                df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                        if threshold_symbol == '>':
                            if tope > threshold_number:
                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                            else:
                                df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                        if threshold_symbol == '<':
                            if tope < threshold_number:
                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                            else:
                                df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                    elif descriptor == 'tc' or descriptor == 'tmc' and comparison != hola:
                        shift=0
                        df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                    # ssc, fsc
                    if descriptor == 'ssc' or descriptor == 'fsc' and comparison == hola:
                        shift=1
                        df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                    elif descriptor == 'ssc' or descriptor == 'fsc' and comparison != hola:
                        shift=0
                        df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                    # pic, pimc
                    if descriptor == 'pic' or descriptor == 'pimc' and comparison == hola:

                        try:

                            # Calling variables
                            sequence,sequence_col,color,spic,len_pic_res=call_var_pic(df_array,row_pi_res,ccpi)
                            signal=0
                            # Auxiliar variable
                            aux_col=0

                            if len(spic) == 2:    # If it's an interval num1-num2 123
                                # Creating list based on intervals
                                list_inpi = [ pia for pia in range(int(spic[0]), int(spic[1])+1) ]

                                if any (v_pi in [rcope] for v_pi in list_inpi):

                                    aux_col=1
                                    signal=1

                                    if protein[sequence[0]] == name_col and comparison == hola:

                                        # pic
                                        if type(mutated_res) != type_bool:
                                            style=df_array['style'].iloc[row_pi_res]
                                            color=style+'_'+color
                                            df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                                        # pimc
                                        else:
                                             res_ref=column_col_ref[rcref] # Residue in reference
                                             if row_col == res_ref: # If it is the same like in reference
                                                style=df_array['style'].iloc[row_pi_res]
                                                color=style+'_'+color
                                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                                             else:
                                                color=colors['RGB'].iloc[mutation_col]
                                                style=df_array['style'].iloc[row_pi_res]
                                                color=style+'_'+color
                                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)

                                    else:
                                        df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                                    if rcope == list_inpi[-1]:
                                        ccpi+=1

                                        if ccpi == len_pic_res and protein[sequence[0]] == name_col:

                                            tar_col_num=add_list_pic(df_array,row_pi_res,tar_col_num)
                                            row_pi_res+=1
                                            ccpi=0  # we should restart the counter

                                            if turn_com and protein[sequence[0]] == name_col and comparison == hola:
                                                scn=colors['Name'].iloc[sequence_col]
                                                ali_log.write(f'Descriptor: {descriptor} - Color (sequence {sequence[0]}): {scn} - {name_col}\n')

                            elif len(spic) == 1:   # Else it's just one number 123

                                if rcope == int(spic[0]):
                                    aux_col=1
                                    signal=1
                                    ccpi+=1

                                    if protein[sequence[0]] == name_col:

                                        # pic
                                        if type(mutated_res) != type_bool:
                                            style=df_array['style'].iloc[row_pi_res]
                                            color=style+'_'+color
                                            df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                                        # pimc
                                        else:
                                             res_ref=column_col_ref[rcref] # Residue in reference
                                             if row_col == res_ref: # If it is the same like in reference
                                                style=df_array['style'].iloc[row_pi_res]
                                                color=style+'_'+color
                                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)
                                             else:
                                                color=colors['RGB'].iloc[mutation_col]
                                                style=df_array['style'].iloc[row_pi_res]
                                                color=style+'_'+color
                                                df_auxs=df_auxs.append({name_col:color} , ignore_index=True)

                                    else:
                                        df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                                    if ccpi == len_pic_res and protein[sequence[0]] == name_col:

                                        tar_col_num=add_list_pic(df_array,row_pi_res,tar_col_num)
                                        row_pi_res+=1
                                        ccpi=0  # we should restart the counter

                                    if turn_com and protein[sequence[0]] == name_col and comparison == hola:
                                        scn=colors['Name'].iloc[sequence_col]
                                        ali_log.write(f'Descriptor: {descriptor} - Color (sequence {sequence[0]}): {scn} - {name_col}\n')

                            if aux_col == 0:
                                df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                            final_pic(col_c, num_col_seq,rcope,numrws,descriptor,tar_col_num,mutation_col,turn_com,ali_log)

                        except TypeError:
                            df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)
                            final_pic(col_c, num_col_seq,rcope,numrws,descriptor,tar_col_num,mutation_col,turn_com,ali_log)

                        except IndexError:
                            df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)
                            final_pic(col_c, num_col_seq,rcope,numrws,descriptor,tar_col_num,mutation_col,turn_com,ali_log)

                    # For sequences without assign color but with descriptor pic or pimc
                    elif descriptor == 'pic' or descriptor == 'pimc':
                        df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)
                        if rcope == numrws-1 and signal == 1:
                            row_pi_res+=1
                            descriptor='pic'

                ######## Gaps will be colored like white... these gaps are just related to the input data gaps. #####

                else:

                    if descriptor == 'ssc' or descriptor == 'fsc' or descriptor == 'fac' or descriptor == 'fmac' or  descriptor == 'pkac' or descriptor == 'nc' or descriptor == 'tc' or descriptor == 'tmc' and comparison == hola:
                        shift=1
                        df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)
                    elif descriptor == 'ssc' or descriptor == 'fsc' or descriptor == 'fac' or descriptor == 'fmac' or  descriptor == 'pkac' or descriptor == 'nc' or descriptor == 'tc' or descriptor == 'tmc' and comparison != hola:
                        shift=0
                        df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                    if descriptor == 'pic' or descriptor == 'pimc' and comparison == hola:

                        color_neutral=colors['RGB'].iloc[116]

                        try:

                            # Calling variables
                            sequence,sequence_col,color,spic,len_pic_res=call_var_pic(df_array,row_pi_res,ccpi)
                            signal=0
                            # Auxiliar variable
                            aux_col=0

                            if len(spic) == 2:  # If it's an interval num1-num2 123

                                list_inpi = [ pia for pia in range(int(spic[0]), int(spic[1])+1) ]
                                # Creating list based on intervals
                                if any (v_pi in [rcope] for v_pi in list_inpi):
                                    signal=1
                                    aux_col=1
                                    df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                                    if rcope == list_inpi[-1]:
                                        ccpi+=1

                                        if ccpi == len_pic_res and protein[sequence[0]] == name_col:

                                            tar_col_num=add_list_pic(df_array,row_pi_res,tar_col_num)
                                            row_pi_res+=1
                                            ccpi=0  # we should restart the counter

                                            if turn_com and protein[sequence[0]] == name_col and comparison == hola:
                                                scn=colors['Name'].iloc[sequence_col]
                                                ali_log.write(f'Descriptor: {descriptor} - Color (sequence {sequence[0]}): {scn} - {name_col}\n')

                            elif len(spic) == 1:        # Else it's just one number 123

                                if rcope == int(spic[0]):
                                    signal=1
                                    aux_col=1
                                    ccpi+=1
                                    df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                                    if ccpi == len_pic_res and protein[sequence[0]] == name_col:

                                        tar_col_num=add_list_pic(df_array,row_pi_res,tar_col_num)
                                        row_pi_res+=1
                                        ccpi=0          # We should restart the counter

                                        if turn_com and protein[sequence[0]] == name_col and comparison == hola:
                                            scn=colors['Name'].iloc[sequence_col]
                                            ali_log.write(f'Descriptor: {descriptor} - Color (sequence {sequence[0]}): {scn} - {name_col}\n')

                            if aux_col == 0:
                                df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

                            final_pic(col_c, num_col_seq,rcope,numrws,descriptor,tar_col_num,mutation_col,turn_com,ali_log)

                        except TypeError:
                            df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)
                            final_pic(col_c, num_col_seq,rcope,numrws,descriptor,tar_col_num,mutation_col,turn_com,ali_log)

                        except IndexError:
                            df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)
                            final_pic(col_c, num_col_seq,rcope,numrws,descriptor,tar_col_num,mutation_col,turn_com,ali_log)

                    # For sequences without assign color but with descriptor pic or pimc
                    elif descriptor == 'pic' or descriptor == 'pimc':
                        df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)
                        if rcope == numrws-1 and signal == 1:
                            row_pi_res+=1

            else:   # Reference will be colored as white

                df_auxs=df_auxs.append({name_col:color_neutral} , ignore_index=True)

            rcref+=1
            rcope+=1


        df_auxs=pd.DataFrame(df_auxs)
        df_auxs_bas=[df_auxs_bas, df_auxs]
        df_auxs_bas=(pd.concat(df_auxs_bas, axis=1, sort=False))

        col_c+=1

    dfbas=df_auxs_bas.to_string(index=False)

    if descriptor == 'fsc':

        if len(sq_n) != 0 and shift == 1:
            tar_col_num=tar_col_num[1:]

    #print(tar_col_num) # List of residues to use for scale
    df_auxs_bas=df_auxs_bas.append({name_col:tar_col_num}, ignore_index=True)
    #print(df_auxs_bas)
    df_auxs_bas=[df_bas_col, df_auxs_bas]
    df_bas_col=(pd.concat(df_auxs_bas, axis=1, sort=False))

    hola+=1 # Counter of comparisons

    if descriptor == 'pic' or descriptor == 'pimc':
        cc=row_pi_res
        if turn_com and mutated_res == True:
            rcn=colors['Name'].iloc[mutation_col]
            ali_log.write(f'Descriptor: {descriptor} - Color {rcn} - mutation\n\n')
    else:
        if shift == 1:
            cc=row_pi_res+1
        else:
            cc=row_pi_res

    # If you want to verify how is going the color DataFrame
    #ali_log.write(f"\t{dfbas}\n\n")
    return df_bas_col, cc, hola
