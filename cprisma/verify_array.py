from cprisma.arrays_default import *
from cprisma.print_array import *
from cprisma.alert_array import *
import pandas as pd
import random

def verify_array_feature(default_parameter,dic_array,ck_dic_arr,feature,ck_ref,dic_ref,met_ref,ch_r,tr_r,ck_tr,df_da,ali_log,counter):

    ali_log.write('------------------------------------------------------------------------------------------------------------------------------------------\n\n')
    print('------------------------------------------------------------------------------------------------------------------------------------------')

    ali_log.write(f"{' '*40}--- Verifying array based on feature '{feature}' ---\n\n")

                ######## Common variables for all methods ########

    if feature == 'operation':
        suffix_feature='ope'
        default_parameter_aux={'n'}
        column_names=["main", "comparison", "descriptor"]

    if feature == 'maximum':
        suffix_feature='mxr'
        default_parameter_aux={'nr'}
        column_names=["main", "comparison", "descriptor", "sequence", "residue", "position", "mutation", "threshold", "separate"]

    if feature == 'color':
        suffix_feature='col'
        default_parameter_aux={'nc'}
        column_names=["main", "comparison", "descriptor", "sequence", "residue", "mutation", "style", "threshold"]

    if  feature == 'visualization':
        suffix_feature='vis'
        default_parameter_aux={ 'ReY','DeN','LeY' }
        column_names=["main", "comparison", "reference", "degraded", "letters"]

    df_array=pd.DataFrame(columns = column_names)
    dict_r={}
    type_dict=type(dict_r)
    set_r={''}
    type_set=type(set_r)
    num_r=17
    type_int=type(num_r)
    num_r=17.0
    type_float=type(num_r)
    lis_r=[]
    type_list=type(lis_r)
    bool_r=True
    type_bool=type(bool_r)

    index=df_da.index
    array_error=False
    check_pic=False

    final_sentece="All the descriptors for "+"'"+feature+"'"+" feature were changed to "+"'"+str(default_parameter_aux)+"'"+' in dict_'+suffix_feature

    ali_log.write(f"Comparing dict_ref with dict_{suffix_feature}:\n\n")
    print(f"Comparing dict_ref with dict_{suffix_feature}:")

                    #### Verifying array for default method ####

    if met_ref == 'default':

        ck_dic_arr=True

        if not ck_ref:           # If it is not compared nothing so ck_dic_arr should be False too
            ck_dic_arr=False
        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter,ali_log,suffix_feature)

                #### Verifying array for pair method ####

    elif met_ref == 'pair':

        ck_dic_arr=True

        if not ck_ref:           # If it is not compared nothing so ck_dic_arr should be False too
            ck_dic_arr=False

        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter,ali_log,suffix_feature)

                #### Verifying array for multiple method ####

    elif met_ref == 'multiple':

        if not ck_ref:           # If it is not compared nothing so ck_dic_arr should be False too
            ck_dic_arr=False

        if not ck_dic_arr:       # If it is not, all comparison will have the default_parameter selected

            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter,ali_log,suffix_feature)

    ######################################## GENERAL CHECKING ########################################

    pos_comparison=1 # To count the total number of comparison

    for key,list_dic in dic_ref.items():

                 #### Verification of main key. Note: common to all features

                if key in dic_array:

                    len_dr=len(dic_ref[key])
                    len_do=len(dic_array[key])

                     #### Verification of number of comparisons. Note: common to all features

                    if len_dr != len_do:

                        array_error=True
                        message_error_in_comparison(ali_log,suffix_feature,key,final_sentece,len_dr,len_do)
                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                        break

                    ############################################################################################

                    ############################ Verification for operation feature ##############################

                    ############################################################################################

                    if feature == 'operation':

                        for ld in range(len(list_dic)):

                            f=dic_array[key][ld]
                            drv=list_dic[ld]

                            try:  # If try f is ok will be iterated, else the brackets are empty

                                ali_log.write(f"\t\t\t\t{drv}........ {list(f)[0]}\n")
                                print('\t\t\t\t',drv,'........', list(f)[0])

                                if type(f) == type_set:            # If it is <class 'set'>, like {'da'}

                                    if list(f)[0] in ['d','da','m','n']:
                                        if array_error == False:
                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : list(f)[0]} , ignore_index=True)

                                    else:
                                        message_error3_in_operation(ali_log,suffix_feature,key,final_sentece,f)
                                        array_error=True
                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                        break

                                    pos_comparison+=1  # Counting position of the comparisons

                                else:
                                    message_error2_in_operation(ali_log,suffix_feature,key,final_sentece,f)
                                    array_error=True
                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                    break

                            except TypeError:
                                message_error1_in_operation(ali_log,suffix_feature,key,final_sentece,f)
                                default_parameter='n'
                                array_error=True
                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                break

                            except IndexError:
                                message_error1_in_operation(ali_log,suffix_feature,key,final_sentece,f)
                                default_parameter='n'
                                array_error=True
                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                break

                    ############################################################################################

                    ############################ Verification for maximum feature ##############################

                    ############################################################################################

                    if feature == 'maximum':

                        for ld in range(len(list_dic)):

                            f=dic_array[key][ld]
                            drv=list_dic[ld]

                            if bool(f):       # If the subarray have parameters inside pass

                                ali_log.write(f"\t\t\t\t{drv}........ {dic_array[key][ld]}\n")
                                print('\t\t\t\t',drv,'........', dic_array[key][ld])

                                if type(f) == type_dict:            # If it is <class 'dict'>, like {'rt': [0.01,'>=']}

                                    for sub_k,sub_v in f.items():

                                        if sub_k in ['rt','ra','ram','rs','rsa','rsm','rsam','rspi','rspim','rpi',
                                                     'rtY','raY','ramY','rsY','rsaY','rsmY','rsamY','rspiY','rspimY','rpiY']: #### 1      # Descriptors with addtional conditions

                                            if sub_k == 'rt' or sub_k == 'rtY':                                  # Restriction by threshold

                                                if len(sub_v) == 2:                 # We need at least two parameters a number and a symbol of order relation

                                                    if type(sub_v[0])==type_int or type(sub_v[0])==type_float:
                                                        pass
                                                    else:
                                                        message_error7_in_maximum_rt(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v)
                                                        array_error=True
                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                        break

                                                    if sub_v[1] in ['>=', '<=', '>', '<', '==', '!=']:
                                                        if array_error == False:
                                                            if sub_k == 'rt':
                                                                df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'threshold': sub_v} , ignore_index=True)
                                                            else:
                                                                df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'threshold': sub_v, 'separate': True} , ignore_index=True)
                                                    else:
                                                        message_error8_in_maximum_rt(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v)
                                                        array_error=True
                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                        break

                                                else:
                                                    message_error6_in_maximum_rt(ali_log,suffix_feature,key,final_sentece,sub_k)
                                                    array_error=True
                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                    break


                                            elif sub_k == 'rpi' or sub_k == 'rpiY':                               # Restriction by position index

                                                if type(sub_v) == type_list:      # Detecting if it is a list

                                                    list_pi=sub_v
                                                    array_error,dic_array,df_array=check_position_df(array_error,df_array,ck_dic_arr,met_ref,ck_ref,dic_array,list_pi,index,suffix_feature,sub_k,key,final_sentece,dic_ref,default_parameter,default_parameter_aux,ali_log)

                                                    if array_error == False:

                                                        if sub_k == 'rpi':
                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'position': sub_v} , ignore_index=True)
                                                        else:
                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'position': sub_v, 'separate': True} , ignore_index=True)

                                                else:
                                                    message_error9_in_maximum_rpi(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,index)
                                                    array_error=True
                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                    break

                                            elif sub_k == 'ra' or sub_k == 'ram' or sub_k == 'raY' or sub_k == 'ramY':              # Restriction by amino acid or by amino acid and mutation

                                                if type(sub_v) == type_list:

                                                    list_am=sub_v           # Verifying list of amino acids of ra and ram
                                                    array_error,dic_array,df_array=check_aminoacid_list(array_error,df_array,ck_dic_arr,met_ref,ck_ref,dic_array,list_am,ck_tr,default_parameter,default_parameter_aux,suffix_feature,sub_k,key,tr_r,dic_ref,final_sentece,ali_log)

                                                    if array_error == False:
                                                        if sub_k == 'ra':
                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'residue': sub_v} , ignore_index=True)
                                                        elif sub_k == 'ram':
                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'residue': sub_v, 'mutation': True} , ignore_index=True)
                                                        elif sub_k == 'raY':
                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'residue': sub_v, 'separate': True} , ignore_index=True)
                                                        elif sub_k == 'ramY':
                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'residue': sub_v, 'mutation': True, 'separate': True} , ignore_index=True)

                                                else:
                                                    message_error12_in_maximum_ra(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v)
                                                    array_error=True

                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                    break

                                            elif sub_k == 'rs' or sub_k == 'rsm' or sub_k == 'rsY' or sub_k == 'rsmY':              # Restriction by sequence or by sequence and mutation

                                                import copy
                                                begin_dic_ref=copy.deepcopy(dic_ref) # Just to avoid problems with dic_ref later

                                                for sref, sseq in drv.items():          # Creating list with all sequences of dict_ref

                                                    list_gen=sseq
                                                    list_gen.append(int(sref))                   # Concatenation of the reference and target sequences of the dict_ref to compare if the list sub_v has any sequence different
                                                    ck_rs=all(item in list_gen for item in sub_v)

                                                    dic_ref=begin_dic_ref
                                                    list_0=[]

                                                    if ck_rs == True and sub_v != list_0:   # If there is at least one sequence similar and not different in the comparison of dict_mxr and dict_ref it will be accepted

                                                        if array_error == False:
                                                            if sub_k == 'rs':
                                                                df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sub_v} , ignore_index=True)
                                                            elif sub_k == 'rsm':
                                                                df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sub_v, 'mutation': True} , ignore_index=True)
                                                            elif sub_k == 'rsY':
                                                                df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sub_v, 'separate': True} , ignore_index=True)
                                                            elif sub_k == 'rsmY':
                                                                df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sub_v, 'mutation': True, 'separate': True} , ignore_index=True)

                                                    else:
                                                        message_error15_in_maximum_rs(ali_log,suffix_feature,key,final_sentece,sub_k,list_gen)
                                                        array_error=True
                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                        break

                                            elif sub_k == 'rsa' or sub_k == 'rsam' or sub_k == 'rsaY' or sub_k == 'rsamY':            # Restriction by sequence and amino acid and restriction by sequence, amino acid and mutation

                                                list_0=[]
                                                dic_0={}

                                                import copy
                                                begin_dic_ref=copy.deepcopy(dic_ref) # Just to avoid problems with dic_ref later
                                                ck_d,list_rsa,list_gen=check_match_sequence(sub_v,drv)
                                                dic_ref=begin_dic_ref
                                                ck_rd=all(it in list_gen for it in list_rsa)

                                                if ck_d == True and list_rsa != list_0:  # If there is at least one sequence similar and not different in the comparison of dict_mxr and dict_ref it will be accepted

                                                    for sv in sub_v:

                                                        sv_key=list(sv.keys())[0]

                                                        for list_am in sv.values():         # Verifying list of amino acids of rsa and rsam

                                                                array_error,dic_array,df_array=check_aminoacid_list(array_error,df_array,ck_dic_arr,met_ref,ck_ref,dic_array,list_am,ck_tr,default_parameter,default_parameter_aux,suffix_feature,sub_k,key,tr_r,dic_ref,final_sentece,ali_log)

                                                                if array_error == False:
                                                                    if sub_k == 'rsa':
                                                                        df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sv_key, 'residue': list_am} , ignore_index=True)
                                                                    elif sub_k == 'rsam':
                                                                        df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sv_key, 'residue': list_am, 'mutation': True} , ignore_index=True)
                                                                    elif sub_k == 'rsaY':
                                                                        df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sv_key, 'residue': list_am, 'separate': True} , ignore_index=True)
                                                                    elif sub_k == 'rsamY':
                                                                        df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sv_key, 'residue': list_am, 'mutation': True, 'separate': True} , ignore_index=True)

                                                else:
                                                    message_error16_in_maximum_rsa(ali_log,suffix_feature,key,final_sentece,sub_k,list_gen,tr_r)
                                                    array_error=True
                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                    break

                                            elif sub_k == 'rspi' or sub_k == 'rspim' or sub_k == 'rspiY' or sub_k == 'rspimY':            # Restriction by position index and restriction by sequence and position index

                                                try:

                                                    list_0=[]
                                                    import copy
                                                    begin_dic_ref=copy.deepcopy(dic_ref) # Just to avoid problems with dic_ref later
                                                    ck_d,list_rspi,list_gen=check_match_sequence(sub_v,drv)
                                                    dic_ref=begin_dic_ref
                                                    ck_d=all(item in list_gen for item in list_rspi)

                                                    if ck_d == True and list_rspi != list_0:  # If there is at least one sequence similar and not different in the comparison of dict_mxr and dict_ref it will be accepted

                                                        for sv in sub_v:

                                                            sv_key=list(sv.keys())[0]

                                                            for list_am in sv.values():         # Verifying list of amino acids of rsa and rsam

                                                                if type(list_am) == type_list:       # Detecting if it is a list

                                                                    list_pi=list_am
                                                                    array_error,dic_array,df_array=check_position_df(array_error,df_array,ck_dic_arr,met_ref,ck_ref,dic_array,list_pi,index,suffix_feature,sub_k,key,final_sentece,dic_ref,default_parameter,default_parameter_aux,ali_log)

                                                                    if array_error == False:
                                                                        if sub_k == 'rspi':
                                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sv_key, 'position': list_am} , ignore_index=True)
                                                                        if sub_k == 'rspim':
                                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sv_key, 'position': list_am, 'mutation': True} , ignore_index=True)
                                                                        if sub_k == 'rspiY':
                                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sv_key, 'position': list_am, 'separate': True} , ignore_index=True)
                                                                        if sub_k == 'rspimY':
                                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence': sv_key, 'position': list_am, 'mutation': True, 'separate': True} , ignore_index=True)

                                                                else:
                                                                    message_error18_in_maximum_rspi(ali_log,suffix_feature,key,final_sentece,sub_k,index,sv)
                                                                    array_error=True
                                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                    break
                                                    else:
                                                        message_error17_in_maximum_rspi(ali_log,suffix_feature,key,final_sentece,sub_k,list_gen,index)
                                                        array_error=True
                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                        break

                                                except UnboundLocalError:           # Error when is used 'rsa' or 'rsam' format instead of 'rspi' or 'rspim' format
                                                    message_error17_in_maximum_rspi(ali_log,suffix_feature,key,final_sentece,sub_k,list_gen)
                                                    array_error=True
                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                    break


                                        else:  #### We put two more filter to detect if the descriptors 'rm' and 'rt' have unnecessary parameters or if it's missing the descriptor

                                            if sub_k in ['nr','rm', 'nrY','rmY']:
                                                message_error5_in_maximum(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v)
                                                array_error=True
                                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                break
                                            else:
                                                message_error4_in_maximum(ali_log,suffix_feature,key,final_sentece,sub_k)
                                                array_error=True
                                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                break

                                else:                               # Else it is <class 'set'>, like {'nr'}

                                    if list(f)[0] in ['nr','rm','nrY','rmY']:                       # No restriction and restriction by mutation

                                        if list(f)[0] == 'rm' or list(f)[0] == 'rmY':
                                            if array_error == False:
                                                if list(f)[0] == 'rm':
                                                    df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : 'rm', 'mutation': True} , ignore_index=True)
                                                else:
                                                    df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : 'rmY', 'mutation': True, 'separate': True} , ignore_index=True)

                                        elif list(f)[0] == 'nr' or list(f)[0] == 'nrY':
                                            if array_error == False:
                                                if list(f)[0] == 'nr':
                                                    df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : 'nr'} , ignore_index=True)
                                                else:
                                                    df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : 'nrY', 'separate': True} , ignore_index=True)

                                    else:                           # If the descriptors 'rt','ra','rm','ram','rs','rsa','rsm','rsam','rspi','rspim', and 'rpi' need aditional parameters or it is missing a descriptor
                                        if list(f)[0] in ['rt','ra','ram','rs','rsa','rsm','rsam','rspi','rspim','rpi','rtY','raY','ramY','rsY','rsaY','rsmY','rsamY','rspiY','rspimY','rpiY']:
                                            message_error3_in_maximum(ali_log,suffix_feature,key,final_sentece,f)
                                            array_error=True
                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                            break
                                        else:
                                            message_error2_in_maximum(ali_log,suffix_feature,key,final_sentece,f)
                                            array_error=True
                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                            break

                                pos_comparison+=1  # Counting position of the comparisons

                            else: # Else does not have nothing inside
                                message_error1_in_maximum(ali_log,suffix_feature,key,final_sentece)
                                default_parameter='nr'
                                array_error=True
                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                break


                    ############################################################################################

                    ############################# Verification for color feature ###############################

                    ############################################################################################

                    if feature == 'color':

                        import os

                        dir_run=os.path.dirname(os.path.abspath(__file__))
                        color_dir=dir_run+"/colors.csv"
                        name_color_file=color_dir
                        colors=pd.read_csv(name_color_file)
                        index_colors=colors.index
                        list_colors=list(index_colors)

                        for ld in range(len(list_dic)):

                            f=dic_array[key][ld]
                            drv=list_dic[ld]

                            if bool(f):       # If the subarray have parameters inside pass

                                ali_log.write(f"\t\t\t\t{drv}........ {dic_array[key][ld]}\n")
                                print('\t\t\t\t',drv,'........', dic_array[key][ld])

                                if type(f) == type_dict:            # If it is <class 'dict'>, like {'aic': True}

                                    for sub_k,sub_v in f.items():

                                        #print(sub_k)

                                        if sub_k in ['fsc','fac','fmac', 'pkac', 'pic', 'pimc', 'tc', 'tmc']:

                                            list_0=[]

                                            if sub_k == 'fsc'or sub_k == 'fac' or sub_k == 'fmac': # Free sequence color, free amino acid color and free and mutation amino acid color

                                                try:        # Because sometimes it is introduced a boolean in this variables we should try

                                                    if sub_v != list_0: # Checking values of the descriptor have full list

                                                        if sub_k == 'fsc': ### 1 Free sequence color

                                                            ld=list(drv.values())[0]

                                                            if len(sub_v) == len(ld) or sub_v[0] == 'random':

                                                                for sv in sub_v:

                                                                    if sv in list_colors or sv == 'random':
                                                                        pass
                                                                    else:
                                                                        message_error26_in_color_free(ali_log,suffix_feature,key,final_sentece,sub_k,sv,name_color_file,list_colors)
                                                                        array_error=True
                                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                        break

                                                                if array_error == False:    # Dataframe fsc
                                                                    df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence' : sub_v} , ignore_index=True)

                                                            else:
                                                                message_error25_in_color_free(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,ld)
                                                                array_error=True
                                                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                break

                                                        else:           #### 1 Free amino acid color or free amino acid and mutation color ### fac fmac

                                                            if len(sub_v) == len (tr_r) or sub_v[0] == 'random':

                                                                for sv in sub_v:

                                                                    if sv in list_colors or sv == 'random':
                                                                        pass
                                                                    else:
                                                                        array_error=True
                                                                        message_error26_in_color_free(ali_log,suffix_feature,key,final_sentece,sub_k,sv,name_color_file,list_colors)
                                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                        break

                                                                if array_error == False: # Dataframe fac, fmac
                                                                    if sub_k == 'fac':
                                                                        df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'residue' : sub_v} , ignore_index=True)
                                                                    else: # fmac
                                                                        df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'residue' : sub_v, 'mutation': True} , ignore_index=True)

                                                            else:
                                                                array_error=True
                                                                message_error28_in_color_defined(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,tr_r)
                                                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                break

                                                    else:
                                                        message_error24_in_color_free(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,name_color_file)
                                                        array_error=True
                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                        break

                                                except TypeError:
                                                    message_error24_in_color_free(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,name_color_file)
                                                    array_error=True
                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                    break

                                            elif sub_k == 'pkac': # pKa condition

                                                ck_res_ion=all(item in ch_r for item in list(tr_r)) # checking if tr_r have ionizable residues

                                                if ck_res_ion: # pkac just work with at least one ionizable residue in your list of target residues

                                                    if type(sub_v) == type_bool:
                                                        pass
                                                    else:
                                                        message_error27_in_color_defined(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v)
                                                        array_error=True
                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                        break

                                                    if array_error == False:    # Dataframe 'pkac
                                                        df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'residue' : ['ionizable'], 'mutation': sub_v} , ignore_index=True)

                                                else:
                                                    message_error29_in_color_defined(ali_log,suffix_feature,key,final_sentece,sub_k,ch_r)
                                                    array_error=True
                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

                                            elif sub_k == 'pic' or sub_k == 'pimc': # Position index color and position index and mutation color

                                                if type(sub_v) == type_list:

                                                    for supi in sub_v:

                                                        if type(supi) == type_dict:

                                                            for ksupi,vsupi in dict.items(supi):

                                                                if list(drv.values())[0] ==  list_0: # if we don't comparison 123

                                                                    if ksupi == list(drv.keys())[0]: # we confirm that it's a valid sequences
                                                                        col_ts=vsupi[0]
                                                                        stl_ts=vsupi[1]

                                                                        try: # Because you can put an int in the list

                                                                            if col_ts <= (max(list(index_colors))) and col_ts >= 0 and any(format_seq in stl_ts for format_seq in ['n', 'b', 'i', 'd', 'u']):

                                                                                for lvsupi in vsupi[2:]:

                                                                                    nsupi=lvsupi.split("-")

                                                                                    list_two=[1,2]

                                                                                    if len(nsupi) == len(list_two): # if we have an interval

                                                                                        if int(nsupi[1]) > int(nsupi[0]) and int(nsupi[0]) >= 0 and int(nsupi[0]) <= max(list(index)) and int(nsupi[1]) >= 0 and int(nsupi[1]) <= max(list(index)):

                                                                                            check_pic=True

                                                                                        else:
                                                                                            message_error35_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,max(list(index)))
                                                                                            array_error=True
                                                                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                                            break

                                                                                    else:
                                                                                                                # else we have just one number
                                                                                        if int(nsupi[0]) >= 0 and int(nsupi[0]) <= max(list(index)):

                                                                                            check_pic=True

                                                                                        else:
                                                                                            message_error36_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,nsupi[0],max(list(index)))
                                                                                            array_error=True
                                                                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                                            break

                                                                                if check_pic:

                                                                                    if array_error == False:

                                                                                        if  sub_k == 'pic':
                                                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence' : [int(ksupi)] ,'residue' : [col_ts, vsupi[2:]], 'style': stl_ts }, ignore_index=True)
                                                                                        else:
                                                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence' : [int(ksupi)] ,'residue' : [col_ts, vsupi[2:]], 'mutation': True, 'style': stl_ts }, ignore_index=True)

                                                                            else:
                                                                                message_error34_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,max(list(index_colors)))
                                                                                array_error=True
                                                                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                                break

                                                                        except TypeError:
                                                                            message_error34_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,max(list(index_colors)))
                                                                            array_error=True
                                                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                            break

                                                                    else:
                                                                        message_error32_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,ksupi,list(drv.keys())[0])
                                                                        array_error=True
                                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                        break

                                                                else:                               # else we have comparison 123
                                                                    tr_seq=list(drv.values())[0]
                                                                    int_ksupi=int(ksupi)

                                                                    if any(trs in [int_ksupi] for trs in tr_seq):
                                                                        col_ts=vsupi[0]
                                                                        stl_ts=vsupi[1]

                                                                        try: # Because you can put an int in the list

                                                                            if col_ts <= (max(list(index_colors))) and col_ts >= 0 and any(format_seq in stl_ts for format_seq in ['n', 'b', 'i', 'd', 'u']):

                                                                                for lvsupi in vsupi[2:]:

                                                                                    nsupi=lvsupi.split("-")

                                                                                    list_two=[1,2]

                                                                                    if len(nsupi) == len(list_two): # if we have an interval

                                                                                        if int(nsupi[1]) > int(nsupi[0]) and int(nsupi[0]) >= 0 and int(nsupi[0]) <= max(list(index)) and int(nsupi[1]) >= 0 and int(nsupi[1]) <= max(list(index)):

                                                                                            check_pic=True

                                                                                        else:
                                                                                            message_error35_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,max(list(index)))
                                                                                            array_error=True
                                                                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                                            break

                                                                                    else:
                                                                                                                # else we have just one number
                                                                                        if int(nsupi[0]) >= 0 and int(nsupi[0]) <= max(list(index)):

                                                                                            check_pic=True

                                                                                        else:
                                                                                            message_error36_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,nsupi[0],max(list(index)))
                                                                                            array_error=True
                                                                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                                            break

                                                                                if check_pic: # Putting parameters

                                                                                    if array_error == False:

                                                                                        if  sub_k == 'pic':
                                                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence' : [int(ksupi)] ,'residue' : [col_ts, vsupi[2:]], 'style': stl_ts }, ignore_index=True)
                                                                                        else:
                                                                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence' : [int(ksupi)] ,'residue' : [col_ts, vsupi[2:]], 'mutation': True, 'style': stl_ts }, ignore_index=True)

                                                                            else:
                                                                                message_error34_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,max(list(index_colors)))
                                                                                array_error=True
                                                                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                                break

                                                                        except TypeError:
                                                                            message_error34_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,vsupi,max(list(index_colors)))
                                                                            array_error=True
                                                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                            break

                                                                    else:
                                                                        message_error33_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi,ksupi,tr_seq)
                                                                        array_error=True
                                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                                        break

                                                        else:
                                                            message_error31_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,supi)
                                                            array_error=True
                                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                            break

                                                else:
                                                    message_error30_in_color_position(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v)
                                                    array_error=True
                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                    break

                                            elif sub_k == 'tc' or sub_k == 'tmc': # Threshold color and threshold and mutation color

                                                if len(sub_v) == 3:

                                                    col_tc=max(list(index_colors))
                                                    first_tc=sub_v[0]

                                                    try: # Verifying first position

                                                        if int(first_tc) <= col_tc: # First element should be a color number

                                                            second_tc=sub_v[1]

                                                            if type(second_tc) == type_int or type(second_tc) == type_float: # Verifying second position

                                                                third_tc=sub_v[2]

                                                                if third_tc in ['>=', '<=', '>', '<', '==', '!=']:  # Verifying third position

                                                                    if sub_k == 'tc':
                                                                        df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence' : [first_tc], 'threshold': [second_tc, third_tc]} , ignore_index=True)
                                                                    else:
                                                                        df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : sub_k, 'sequence' : [first_tc], 'threshold': [second_tc, third_tc], 'mutation': True} , ignore_index=True)

                                                                else:
                                                                    message_error40_in_color_threshold(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,col_tc,third_tc)
                                                                    array_error=True
                                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

                                                            else:
                                                                message_error39_in_color_threshold(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,col_tc,second_tc)
                                                                array_error=True
                                                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

                                                        else:
                                                            message_error38_in_color_threshold(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,col_tc,first_tc)
                                                            array_error=True
                                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

                                                    except ValueError:
                                                        message_error38_in_color_threshold(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,col_tc,first_tc)
                                                        array_error=True
                                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

                                                else:
                                                    message_error37_in_color_threshold(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v)
                                                    array_error=True
                                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

                                        else:

                                            if sub_k in ['nc','ssc']:                       # They should not have additional parameters
                                                message_error23_in_color(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v)
                                                array_error=True
                                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                break
                                            else:
                                                message_error22_in_color(ali_log,suffix_feature,key,final_sentece,sub_k)
                                                array_error=True
                                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                                break

                                else:                               # Else it is <class 'set'>, like {'ssc'}

                                    if list(f)[0] in ['nc','ssc']:                       # No color and same sequence color

                                        if array_error == False:    # Dataframe 'nc','ssc'
                                            if list(f)[0] == 'nc':
                                                df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : list(f)[0]} , ignore_index=True)
                                            elif list(f)[0] == 'ssc':
                                                df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : list(f)[0], 'sequence': ['same']} , ignore_index=True)

                                    else:                           # If the descriptors 'fsc','fmac','fac','pkac','pic', 'tc' need aditional parameters or it is missing a descriptor
                                        if list(f)[0] in ['fsc','fmac','fac','pkac','pic','pimc', 'tc', 'tmc']:
                                            message_error21_in_color(ali_log,suffix_feature,key,final_sentece,f)
                                            array_error=True
                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                            break

                                        else:
                                            message_error20_in_color(ali_log,suffix_feature,key,final_sentece,f)
                                            array_error=True
                                            array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                            break

                                pos_comparison+=1  # Counting position of the comparisons

                            else: # Else does not have nothing inside

                                message_error19_in_color(ali_log,suffix_feature,key,final_sentece)
                                default_parameter='nc'
                                array_error=True
                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                break

                    ############################################################################################

                    ######################### Verification for visualization feature ###########################

                    ############################################################################################

                    if feature == 'visualization':

                        for ld in range(len(list_dic)):

                            f=dic_array[key][ld]
                            drv=list_dic[ld]

                            try:  # If try f is ok will be iterated, else the brackets are empty

                                ali_log.write(f"\t\t\t\t{drv}........ {list(f)[0:]}\n")
                                print('\t\t\t\t',drv,'........', list(f)[0:])

                                if type(f) == type_set:            # If it is <class 'set'>

                                    list_gen=['ReY', 'DeY', 'LeY', 'ReN', 'DeN', 'LeN']

                                    ck_vis=all(item in list_gen for item in list(f)[0:])

                                    if ck_vis:

                                        if array_error == False:

                                            if list(f)[0] == 'ReY' or list(f)[1] == 'ReY' or list(f)[2] == 'ReY': # Reference visible or not
                                                ref_visual=True
                                            else:
                                                ref_visual=False

                                            if list(f)[0] == 'DeY' or list(f)[1] == 'DeY' or list(f)[2] == 'DeY': # Degraded visible or not
                                                deg_visual=True
                                            else:
                                                deg_visual=False

                                            if list(f)[0] == 'LeY' or list(f)[1] == 'LeY' or list(f)[2] == 'LeY': # Letters visible or not
                                                let_visual=True
                                            else:
                                                let_visual=False

                                            ## Dataframe for visualization ##

                                            df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'reference' : ref_visual, 'degraded' : deg_visual,'letters' : let_visual} , ignore_index=True)

                                    else:
                                        message_error30_in_visualization(ali_log,suffix_feature,key,final_sentece)
                                        array_error=True
                                        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                        break

                                    pos_comparison+=1  # Counting position of the comparisons

                                else:
                                    message_error30_in_visualization(ali_log,suffix_feature,key,final_sentece)
                                    array_error=True
                                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                                    break

                            except TypeError:
                                message_error30_in_visualization(ali_log,suffix_feature,key,final_sentece)
                                default_parameter={ 'ReY','DeN','LeY' }
                                array_error=True
                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

                            except IndexError:
                                message_error30_in_visualization(ali_log,suffix_feature,key,final_sentece)
                                default_parameter={ 'ReY','DeN','LeY' }
                                array_error=True
                                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

                    ############################################################################################



                    ############################################################################################


                #### Verifying if there are same name or main key. Note: common to all features

                else:

                    array_error=True
                    message_error_in_main_key(ali_log,suffix_feature,key,final_sentece)
                    array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                    break

    ##################################################################################################

                                    ### Printing array ###

    ali_log.write(f'\n')
    title='Array generated for the method '+"'"+met_ref+"'"+' and feature '+"'"+feature+"'"+': '
    #print_array(title,dic_array,ali_log)                  # Function just to print the array separating the main keys
    ali_log.write(f"The array of the feature '{feature}' dict_{suffix_feature} is compatible with the array of comparison sequences dict_ref!\n\n")
    print(f"The array of the feature '{feature}' dict_{suffix_feature} is compatible with the array of comparison sequences dict_ref!")

                                     ### Printing df ###



    df_array_v=df_array.to_string(index=False)
    ali_log.write(f"{df_array_v}\n\n")
    ali_log.write('------------------------------------------------------------------------------------------------------------------------------------------\n\n')
    print(f"{df_array_v}")
    print('------------------------------------------------------------------------------------------------------------------------------------------')

    return dic_array, df_array, dic_ref
