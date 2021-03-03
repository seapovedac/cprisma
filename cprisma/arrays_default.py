from cprisma.alert_array import *
import random

##########################################################################################################################

def array_default(protein):

    lenseq=len(protein)
    dic_ref={}
    list_ref=[]
    list_seq=[]
    for i in range(lenseq):
        if i == 0:
            dic_sub={}
            dic_sub[str(i)]=[]
            join='main_'+str(i)
            dic_ref[join]=[dic_sub]
        else:
            list_seq.append(i)
    dic_sub['0']=list_seq

    return dic_ref

##########################################################################################################################

def special_case_maximum (met_ref,descriptor,list_d,default_parameter,ali_log,suffix_feature,d,dict_aux,list_aux,pos_comparison,array_error):

    if descriptor in ['rs', 'rsm']:   # Maximum decriptors
      default_parameter={'nr'}
      # it cannot support restricting the values for the same sequence so restart all the array
      array_error=True
      list_d.append(default_parameter)
      message_error_in_special_case(ali_log,suffix_feature,descriptor,met_ref)

    elif descriptor in ['rsa', 'rsam', 'rspi', 'rspim']:  # Maximum decriptors
        # Creating sub-sub-sub_values
      list_sub_sub_values=(list((list(list(default_parameter.values())[0])[0]).values())[0])
        # Creating sub-sub-sub_key
      sub_sub_keys=str(list(d.keys())[0])
        # Creating sub-sub-sub_dictionary
      dict_aux[sub_sub_keys]=list_sub_sub_values
        # Creating a list for each sequence without comparison
      list_aux.append(dict_aux)
        # Restarting dict_aux
      dict_aux={}
        # Adding the new sub-sub-sub condition
      dict_aux[descriptor]=(list_aux)
        # Appending in the list of the descriptor

      if array_error == False:      # if the array is ok pass
        list_d.append(dict_aux)
      else:                         # if the array is bad, should restart all the comparison
        df_array=df_array_null(df_array,key,pos_comparison)
        list_d.append(dict_aux)

    else: # For ra, ram, rt for dict_mxr
      list_d.append(default_parameter)

##########################################################################################################################

def array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter,ali_log,suffix_feature):

    import pandas as pd

    dic_array={}
    dict_r={}
    type_dict=type(dict_r)

    ###### Functions if Error in array #####

    def header_df_array():
        # Headers for each feature
        df_array=pd.DataFrame()

        if suffix_feature == 'ope':
            column_names=["main", "comparison", "descriptor"]
        elif suffix_feature == 'mxr':
            column_names=["main", "comparison", "descriptor", "sequence", "residue", "position", "mutation", "threshold", "separate"]
        elif suffix_feature == 'col':
            column_names=["main", "comparison", "descriptor", "sequence", "residue", "mutation", "style", "threshold"]
        elif suffix_feature == 'vis':
            column_names=["main", "comparison", "reference", "degraded", "letters"]

        df_array=pd.DataFrame(columns = column_names)
        pos_comparison=1 # To count the total number of comparison

        return df_array, pos_comparison

    def df_array_null(df_array,key,pos_comparison):
        # Putting default parameters for all descriptors
        if suffix_feature == 'mxr':
          df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : 'nr'} , ignore_index=True)
        elif suffix_feature == 'ope':
          df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : 'n'} , ignore_index=True)
        elif suffix_feature == 'col':
          df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'descriptor' : 'nc'} , ignore_index=True)
        elif suffix_feature == 'vis':
          df_array=df_array.append({'main' : key , 'comparison': pos_comparison, 'reference' : True, 'degraded' : False,'letters' : True} , ignore_index=True)

        return df_array

    no_more=0 # Variable of especial case!
    for key,list_dic in dic_ref.items():

      dic_array[key]=''
      list_d=[]

      for d in list_dic:

          if ck_dic_arr == True:

              if met_ref == 'default':  ##### 1 For 'default'

                  df_array,pos_comparison=header_df_array()
                  list_d.append(default_parameter)

              elif met_ref == 'pair':  ##### 1 For 'pair'

                  descriptor=list(default_parameter)[0]
                  # for pair reference method is only possible for cases different to 'rs','rsa', 'rsam', 'rspi', 'rspim'
                  if descriptor in ['rs','rsa', 'rsam', 'rspi', 'rspim']:  # Maximum decriptors where sequence is involved will return error
                      array_error=True
                      list_d.append(default_parameter)
                      message_error_in_special_case(ali_log,suffix_feature,descriptor,met_ref)
                  else:
                      df_array,pos_comparison=header_df_array()
                      list_d.append(default_parameter)

              elif met_ref == 'multiple': ##### 1 For 'multiple' reference method

                  list_d.append(default_parameter)

          else: # Special case!!! For the condition: 'no reference' and 'no dict array' just descriptor_mxr will be consider to built an special array

              if no_more == 0:
                  df_array,pos_comparison=header_df_array()
                  no_more+=1

              if suffix_feature == 'mxr':

                  dict_aux={}
                  list_aux=[]

                  if type(default_parameter) == type_dict:      # For all decriptor with dictionary structure for dict_mxr

                      descriptor=list(default_parameter.keys())[0]
                      special_case_maximum(met_ref,descriptor,list_d,default_parameter,ali_log,suffix_feature,d,dict_aux,list_aux,pos_comparison,array_error)

                  else:  # For all decriptor with set structure for dict_mxr
                      list_d.append(default_parameter)

              elif suffix_feature == 'ope' or suffix_feature == 'col' or suffix_feature == 'vis':   # For features like color, operation and visualization
                  list_d.append(default_parameter)

      dic_array[key]=list_d # Generating the new array

     # If the array has problem should restart its data frame
    if array_error:
        df_array,pos_comparison=header_df_array()
        for erk, erv in dic_array.items():
            for cc in erv:
                df_array=df_array_null(df_array,erk,pos_comparison)
                pos_comparison+=1

    return array_error, dic_array, df_array

##########################################################################################################################

def check_aminoacid_list(array_error,df_array,ck_dic_arr,met_ref,ck_ref,dic_array,list_am,ck_tr,default_parameter,default_parameter_aux,suffix_feature,sub_k,key,tr_r,dic_ref,final_sentece,ali_log):

    global lam

    list_0=[]

    if list_am != list_0:

        for lam in list_am:

          if ck_tr: # If we have target residues ck_tr=True and not empty list
              if lam in tr_r:
                  pass
              else:
                  message_error14_in_maximum_ra(ali_log,suffix_feature,key,final_sentece,sub_k,tr_r,lam)
                  array_error=True
                  array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                  break

          else:    # Else it is compared with all amino acids type ck_tr=False
              if lam in ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','Z']:
                  pass
              else:
                  message_error13_in_maximum_ra(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,tr_r,lam)
                  array_error=True
                  array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                  break

    else:     # If the list is empty
        lam=''
        message_error13_in_maximum_ra(ali_log,suffix_feature,key,final_sentece,sub_k,sub_v,tr_r,lam)
        array_error=True
        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

    return array_error,dic_array,df_array

##########################################################################################################################

def check_match_sequence(sub_v,drv):

    list_sd=[]

    for sv in sub_v:                                # Creating list of sequences of the descriptor

        try:
            seq_d=int(list(sv.keys())[0])
        except ValueError:
            seq_d=''
        except AttributeError:
            seq_d=''
        except IndexError:
            seq_d=''
        except UnboundLocalError:
            seq_d=''

        list_sd.append(seq_d)

    for sref, sseq in drv.items():          # Creating list with all sequences of dict_ref
        list_gen=sseq
        list_gen.append(int(sref))                   # Concatenation of the reference and target sequences of the dict_ref to compare if the list sub_v has any sequence different

    ck_d=all(it in list_gen for it in list_sd)

    return ck_d, list_sd, list_gen

##########################################################################################################################

def check_position_df(array_error,df_array,ck_dic_arr,met_ref,ck_ref,dic_array,list_pi,index,suffix_feature,sub_k,key,final_sentece,dic_ref,default_parameter,default_parameter_aux,ali_log):

    list_0=[]
    if list_pi != list_0:

        for lpi in list_pi:

            if lpi in list(index):
                pass
            else:
                message_error11_in_maximum_rpi(ali_log,suffix_feature,key,final_sentece,sub_k,index,lpi)
                array_error=True
                array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)
                break
    else:
        message_error10_in_maximum_rpi(ali_log,suffix_feature,key,final_sentece,sub_k,index)
        array_error=True
        array_error,dic_array,df_array=array_default_feature(array_error,ck_dic_arr,ck_ref,dic_ref,met_ref,default_parameter_aux,ali_log,suffix_feature)

    return array_error, dic_array, df_array
