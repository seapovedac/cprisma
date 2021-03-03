from cprisma.visualization import *
import pandas as pd
#import os
import shutil

    ##########################################################################################################################

    ##########################################################################################################################

def compile(met_ref, dict_ali, tr_r, df_aa, protein, dict_ref, df_ope, df_mxr, df_col, df_col_array, df_vis, factor_color, first_num, join, conservation, ali_method, tridimensional_col, check_dict_vis, res_per_col, ali_log):

    ali_log.write('------------------------------------------------------------------------------------------------------------------------------------------\n\n')
    print('------------------------------------------------------------------------------------------------------------------------------------------')

    ali_log.write(f"\t--- Visualization features ---\n\n")

    if ali_method ==1:
        ali_log.write(f" Alignment method = {ali_method}\n")
    else:
        ali_log.write(f" Alignment method = {ali_method}\n Join = {join}\n")
    ali_log.write(f" First number = {first_num}\n Residues per column = {res_per_col}\n Factor color = {factor_color}\n\n")

    index=df_aa.index
    numrws=len(index)                   # Number of rows
    cc=0
    dd=1
    ee=0
    numcols1=0
    numcols2=0

    # Cleaning directory if pymol script is generated
    try:
        shutil.rmtree('3D_representation')
    except FileNotFoundError:
        pass

    len_ali=(len(dict_ali[protein[0]]))
    dict_main={}
    t_list=[]
    type_list=type(t_list)

    for key,list_dic in dict_ref.items():

        dict_main[key]={}
        dict_comparison={}

        for d in list_dic:

            i=0

            for rf_n,sq_n in d.items():

                # Maximum restriction
                max=df_mxr['maximum (abs)'].iloc[ cc ]
                # Visibility of degraded
                vis_deg=df_vis['degraded'].iloc[ cc ]
                # Visibility of letters
                vis_let=df_vis['letters'].iloc[ cc ]
                # comparison
                comparison=df_vis['comparison'].iloc[ cc ]

                dict_comparison[comparison]={}

                # Numbers of columns per comparison
                comp=(len([rf_n])+len(sq_n))-1

                # Column counter
                numcols2=numcols1+comp+1
                # Reading columns in df_col
                color_data=df_col.iloc[ : , numcols1 : numcols2 ]

                # Reading columns in df_ope
                seq_ope=df_ope.iloc[ ee:numrws*dd , : ]
                seq_ope=seq_ope.dropna(axis=1, how='all')

                # Extracting all names of comparison
                name_col=list(pd.DataFrame(color_data).to_dict().keys())[0:]

                hh=0 # Row counter for df_col and df_ope
                cre=0
                dict_protein={}
                # Reading the aligment input
                for i in range(len_ali):

                    list_r=[]
                    # Creating a list with aminoacids for each set of proteins in a specific comparison
                    for nncc in name_col:
                        list_r.append(list(dict_ali[nncc])[i])

                    gg=0

                    # Loop for each protein
                    for nncc in name_col:

                        # Creating a new key with name sequence
                        if cre < len(name_col):
                            dict_protein[nncc]=''
                            cre+=1

                        # Filtering just the target residues based on the list created
                        if any(a in tr_r for a in (list_r)):

                            # Calling this function to correctly read the rows of df_col and df_ope
                            hh=read_rows(hh,gg,name_col,df_aa)
                            # Getting the rows in df_col
                            row_color=color_data.iloc[hh-1]
                            # Getting the rows in df_ope
                            row_operation=seq_ope.iloc[hh-1]

                            ck_tr=True
                            ff=0

                            # Loop for the aminoacids in a specific row
                            for a in list_r:

                                if gg == ff:

                                    # Target sequence
                                    if ff != 0:
                                        val_ope=list(row_operation)[ff-1]

                                        # If there is not maximum per residue
                                        if type(max) != type_list:
                                            color_target_3d=True
                                            HTML=color_residue(key,comparison,nncc,tridimensional_col,color_target_3d,ck_tr,row_color,val_ope,ff,vis_deg,vis_let,tr_r,a,factor_color,max,i,first_num)

                                        # Else there is maximum per residue
                                        else:
                                            if a in tr_r:
                                                pos_max=list(tr_r).index(a)
                                                max_r=max[pos_max]
                                            else:
                                                max_r=1.0
                                            color_target_3d=True
                                            HTML=color_residue(key,comparison,nncc,tridimensional_col,color_target_3d,ck_tr,row_color,val_ope,ff,vis_deg,vis_let,tr_r,a,factor_color,max_r,i,first_num)

                                        dict_protein[nncc]+=HTML

                                    # Reference sequence
                                    else:
                                        # If there is not maximum per residue
                                        if type(max) != type_list:
                                            val_ope=list(row_operation)[ff-1]
                                            # If it is a reference with target sequences
                                            if len(name_col) != 1:
                                                val_ope=max
                                            color_target_3d=True
                                            HTML=color_residue(key,comparison,nncc,tridimensional_col,color_target_3d,ck_tr,row_color,val_ope,ff,vis_deg,vis_let,tr_r,a,factor_color,max,i,first_num)

                                        # Else there is maximum per residue
                                        else:
                                            val_ope=list(row_operation)[ff-1]
                                            # If it is a reference with target sequences
                                            if len(name_col) != 1:
                                                max_r=1.0
                                                val_ope=max_r
                                            # If the reference does not have comparisons
                                            else:
                                                if a in tr_r:
                                                    pos_max=list(tr_r).index(a)
                                                    max_r=max[pos_max]
                                                else:
                                                    max_r=1.0
                                            color_target_3d=True
                                            HTML=color_residue(key,comparison,nncc,tridimensional_col,color_target_3d,ck_tr,row_color,val_ope,ff,vis_deg,vis_let,tr_r,a,factor_color,max_r,i,first_num)

                                        dict_protein[nncc]+=HTML
                                ff+=1

                        # Amino acids not involved as a target residue
                        else:
                            ff=0
                            # Loop for the aminoacids in a specific row
                            for a in list_r:
                                if gg == ff:
                                    ck_tr=False
                                    row_color='rgb(255, 255, 255'
                                    val_ope=1.0
                                    color_target_3d=False
                                    HTML=color_residue(key,comparison,nncc,tridimensional_col,color_target_3d,ck_tr,row_color,val_ope,ff,vis_deg,vis_let,tr_r,a,factor_color,max,i,first_num)
                                    dict_protein[nncc]+=HTML
                                ff+=1
                        gg+=1

                # Adding dict_protein to dict_comparison
                dict_comparison[comparison]=dict_protein
                # Column counter df_col
                numcols1=numcols1+comp+1
                # Column counter df_ope
                dict_main[key]=dict_comparison
                cc+=1
                dd+=1
                ee=numrws+ee

    ##########################################################################################################################

    ##########################################################################################################################

    # Converting descriptor pic or pimc to fsc in df_array_col

    column_names=["main", "comparison", "descriptor", "sequence", "residue", "mutation", "style", "threshold"]
    df_auxs_col=pd.DataFrame(columns = column_names)

    # Getting number of last comparison
    last_com=df_col_array['comparison'].iloc[ -1 ]

    for cpa in range(1,last_com+1):

        for des in range (2):

            if des == 0:
                val_des='pic'
            else:
                val_des='pimc'

            new1=df_col_array[df_col_array['descriptor'].str.contains(val_des)]
            new1=new1.loc[new1['comparison'] == cpa]
            ind_ne1=len(new1.index)

            list_pic=[]

            for cpic in range(ind_ne1):
                col_pic=new1['residue'].iloc[cpic]
                col_pic=col_pic[0]
                list_pic.append(col_pic)

            if len(list_pic) != 0:
                key=new1['main'].iloc[cpic]
                if val_des == 'pic':
                    df_auxs_col=df_auxs_col.append({'main' : key , 'comparison': cpa, 'descriptor' : 'fsc', 'sequence' : list_pic}, ignore_index=True)
                else:
                    df_auxs_col=df_auxs_col.append({'main' : key , 'comparison': cpa, 'descriptor' : 'fsc', 'sequence' : list_pic, 'mutation': True}, ignore_index=True)

    # Droping all rows with pic or pimc
    df = df_col_array[~df_col_array['descriptor'].isin(['pic','pimc'])]

    df_auxs_col=df_auxs_col.append(df, ignore_index=True)
    df_auxs_col=df_auxs_col.sort_values(by=['comparison'])
    #print(df_auxs_col)

    ##########################################################################################################################

    ##########################################################################################################################

    html=open('alignment.html', 'w+')
    html.write('<!DOCTYPE html>\n<html>\n<body>\n<pre>\n<font face="Courier">\n<small>\n')
    m=len('<del><font style="background-color:rgb(128, 128, 128, 1.0);">K</font></del>')
    cc=0
    dict_aux={}
    con1=list(dict_ali.keys())[0]
    con2=list(dict_ali.values())[0]
    dict_aux[con1]=con2
    list_big,max_z=final_z(protein)
    dict_num=number_sequence(dict_ali,protein,first_num)
    t_bool=True
    type_bool=type(t_bool)
    l=res_per_col
    n=l*m

    if ali_method == 1:

        method_visualization1(html,dict_main,met_ref,dict_aux,dict_num,df_auxs_col,df_col,df_mxr,df_vis,cc,list_big,max_z,protein,n,l,type_bool,type_list,tr_r,conservation,ali_method)

    elif ali_method == 2:

        method_visualization2(html,dict_main,met_ref,dict_aux,dict_num,protein,n,l,list_big,max_z,df_auxs_col,df_col,df_mxr,df_vis,join,type_list,type_bool,tr_r,conservation,check_dict_vis,ali_method)

    html.write('</small>\n</font>\n</pre>\n</body>\n</html>')
