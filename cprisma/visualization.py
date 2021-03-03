    ##  Function to create a pymol script ##

def pymol_3d(color_target_3d,key,comparison,name_protein,RGB,alpha,i):

    import os

    try:
        os.mkdir('3D_representation')
    except FileExistsError:
        pass

    if color_target_3d:

        try:
            name_main='3D_representation/'+key
            os.mkdir(name_main)
        except FileExistsError:
            pass

        name_script_pml=name_main+'/'+key+'_comparison'+str(comparison)+'_'+name_protein+'.pml'
        pymol_3d=open(name_script_pml, 'a+')

        if i == 0:

            pymol_3d.write(f"# {name_protein} - {key} comparison {comparison}\n")
            pymol_3d.write(f"protein=''\ncmd.color('white' ,protein)\ncmd.show('cartoon' , protein)\nhide lines, all\n\n")

        spi_pymolRGB=RGB.split("rgb(")

        red_spi=int(spi_pymolRGB[1].split(",")[0])
        green_spi=int(spi_pymolRGB[1].split(",")[1])
        blue_spi=int(spi_pymolRGB[1].split(",")[2])

        # Converting RGBA to RGB
        red_RGB=int(255-(255*float(alpha) - red_spi*float(alpha)))
        green_RGB=int(255-(255*float(alpha) - green_spi*float(alpha)))
        blue_RGB=int(255-(255*float(alpha) - blue_spi*float(alpha)))

        pymol_3d.write(f"show sticks, res {i+1}\n")
        pymol_3d.write(f"set_color color_{i+1}, [{red_RGB},{green_RGB},{blue_RGB}]\n")
        pymol_3d.write(f"color color_{i+1}, res {i+1}\n")

    ## Function to put HTML code color per residues ##

def color_residue(key,comparison,name_protein,tridimensional_col,color_target_3d,ck_tr,row_color,val_ope,ff,vis_deg,vis_let,tr_r,a,factor_color,max,i,first_num):

    if ck_tr == True:
        RGB=list(row_color)[ff]
    else:
        RGB=row_color

    # if it is to visualize with degrade
    if vis_deg == True:

        if a in tr_r:
            val_ope=abs(val_ope)
            try:
                alpha=(val_ope/max)*factor_color
            except ZeroDivisionError:
                alpha=0.0
            alpha=repr(round(alpha,1))
            if float(alpha) > 1:
                alpha=1.0
        else:
            alpha=1.0
    else:
        alpha=1.0

    # if letters display is disabled
    if vis_let == False:
        a=' '

    spiRGB=RGB.split("_")

    if len(spiRGB) == 1:

        HTML='<font style="background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'
        space=75-len(HTML)
        HTML='<font style='+' '*space+'"background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'

        if tridimensional_col:
            if any (py_res in [a] for py_res in tr_r):
                pymol_3d(color_target_3d,key,comparison,name_protein,RGB,alpha,i)

    else: # For pic and pimc descriptors

        style=spiRGB[0]
        RGB=spiRGB[1]

        if style == 'b':
            let_beg='<b>'
            let_end='</b>'
            HTML=let_beg+'<font style="background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'+let_end
            space=75-len(HTML)
            HTML=let_beg+'<font style='+' '*space+'"background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'+let_end
        elif style == 'i':
            let_beg='<i>'
            let_end='</i>'
            HTML=let_beg+'<font style="background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'+let_end
            space=75-len(HTML)
            HTML=let_beg+'<font style='+' '*space+'"background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'+let_end
        elif style == 'd':
            let_beg='<del>'
            let_end='</del>'
            HTML=let_beg+'<font style="background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'+let_end
            space=75-len(HTML)
            HTML=let_beg+'<font style='+' '*space+'"background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'+let_end
        elif style == 'u':
            let_beg='<ins>'
            let_end='</ins>'
            HTML=let_beg+'<font style="background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'+let_end
            space=75-len(HTML)
            HTML=let_beg+'<font style='+' '*space+'"background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'+let_end
        elif style == 'n':
            HTML='<font style="background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'
            space=75-len(HTML)
            HTML='<font style='+' '*space+'"background-color:'+RGB+', '+str(alpha)+');">'+a+'</font>'

        if tridimensional_col:
            if any (py_res in [a] for py_res in tr_r):
                pymol_3d(color_target_3d,key,comparison,name_protein,RGB,alpha,i)

    return HTML

    ## Function to correctly read the rows in df_col and df_ope ##

def read_rows(hh,gg,name_col,df_aa):

    if gg == 0:
        list_vnc=[]
        ccvnc=0
        for vnc1 in name_col:
            list_vnc.append(df_aa[vnc1].iloc[hh])
            ccvnc+=1
            # Evaluation to increase one row in df_col and df_ope when df_aa have 'X'
            if ccvnc == len(name_col):
                while list_vnc.count('X') == len(name_col):
                    hh+=1
                    list_vnc=[]
                    for vnc2 in name_col:
                        list_vnc.append(df_aa[vnc2].iloc[hh])
        hh+=1

    return hh

    ## Function to read features ##

def read_feature(df_col_array,df_col,df_mxr,df_vis,cc):

    # Descriptor of color

    descriptor_col=df_col_array['descriptor'].iloc[ cc ]

    # List of colors for each comparison

    if descriptor_col == 'nc':
        color_list=df_col.iloc[-1].dropna()
    else:
        color_list=df_col.iloc[-1].dropna()[cc]

    # Mutation
    mutation=df_col_array['mutation'].iloc[ cc ]
    # Maximum restriction
    max=df_mxr['maximum (abs)'].iloc[ cc ]
    # Minimum restriction
    min=df_mxr['minimum (abs)'].iloc[ cc ]
    # Visibility of degraded
    vis_deg=df_vis['degraded'].iloc[ cc ]
    # Visibility of reference
    vis_ref=df_vis['reference'].iloc[ cc ]

    return descriptor_col, color_list, mutation, max, min, vis_deg, vis_ref

    ## Function to put number in the aligment ##

def number_sequence(dict_ali,protein,first_num):

    def space_num(len_n1,u,len_n2):
        global space
        if len_n1 == 1:
            space=""
        if (u-len_n2) > 0 and (u != 1):
            space=" "*(u-len_n2)
        elif (u-len_n2) > 0 and (u == 1):
            space=""
        if (u-len_n2) < 0:
            space=" "*(21-len_n2)
        space=str(space)
        return space

                        # Variables to define position of number above the alignment
    n=1                 # Amino acid counter
    char_per_line=80    # Number of characters per line
    o=char_per_line/4   # "Distance" between each value that will appear above the alignment
    p=0                 # Variable to select the exact position above the alignment
    q=0                 # Variable to oscillate between zero and char_per_line every o variable
    r=0                 # Variable to exclude the last position per set of alignment division
    s=1                 # Number to appear above the alignment!
    t=0                 # Variable to include the first position for set of alignment division
    u=0                 # Space counter and oscillator to correctly fit the number above the alignment
    len_n2=0
    dict_num={}
    dict_num['number']=''

    for a in dict_ali[protein[0]]:

        u=u+1
        if n == o or p == 0 or n == p:
            if n !=r:
                s=first_num+n-1
                space=space_num(len(str(s)),u,len_n2)
                dict_num['number']+=space
                dict_num['number']+=str(s)
            len_n2=len(str(s))
            u=0
            p=o+p
            q=o+q
        if n == (char_per_line*t+1):
            if n!=1:
                s=first_num+n-1
                space=space_num(len(str(s)),u,len_n2)
                dict_num['number']+=space
                dict_num['number']+=str(s)
                len_n2=len(str(s))
                u=0
        if q == char_per_line:
            q=0
            r=r+char_per_line
            t=t+1
        n=n+1

    return dict_num

    ## Function to avoid problems of space between name sequence and data residues in output files ##

def final_z(l):
    list_big=[]
    for key in l: list_big.append(len(key))
    z=max(list_big)
    list_big=[]
    for key in l: list_big.append(z-len(key))
    return list_big, z

    ## Auxiliar function ##

def max_sc(dict_scale):

    list_s=[]
    for sc in dict_scale.values():
        list_s.append(len(sc))
    max_sc=max(list_s)

    return max_sc

    ## Fucntion for scale ##

def scale(protein,vis_deg,descriptor_col,type_list,min,max,max_z,color_list,html,tr_r,mutation,type_bool,ali_method):

    import pandas as pd
    import copy
    import os

    list_sc=[] # list for second method

    if vis_deg == True:

        dir_run=os.path.dirname(os.path.abspath(__file__))
        color_dir=dir_run+"/colors.csv"
        name_color_file=color_dir
        colors=pd.read_csv(name_color_file)
        aux_scale=''

        if descriptor_col == 'pkac':

            if type(max) != type_list:

                scale=-1.0
                for sc in range (-11,10):
                    if sc < 0:
                     target_color=color_list[1]
                    else:
                     target_color=color_list[0]
                    RGB=colors['RGB'].iloc[target_color]
                    alpha=round(abs(scale),1)
                    HTML='<font style="background-color:'+RGB+', '+str(alpha)+');"> </font>'
                    aux_scale+=HTML
                    scale+=0.1

                if ali_method == 1:

                    html.write(f"{' '*max_z}\t-{max} {aux_scale} {max}\n")

                elif ali_method == 2:

                    pklen=len(str(max))
                    if pklen > 3:
                        ax=pklen-3
                        pklen0=pklen-ax
                    else:
                        pklen0=pklen

                    space=' '*int(pklen0/2)
                    sc_vis='\t\t<small>'+space+'-'+str(max)+' '+aux_scale+' '+str(max)+'</small>'
                    list_sc.append(sc_vis)

            else:

                # To have good spaces
                sptr=[]
                for ttrr in tr_r:
                    pos_max=list(tr_r).index(ttrr)
                    max_r=max[pos_max]
                    sptr.append(len(str(max_r)))

                max_sptr=0
                for mr in sptr:
                    if mr > max_sptr:
                        max_sptr=mr # Number with maximum characters

                for ttrr in tr_r:

                    scale=-1.0
                    aux_scale=''
                    for sc in range (-11,10):
                        if sc < 0:
                            target_color=color_list[1]
                        else:
                            target_color=color_list[0]
                        pos_max=list(tr_r).index(ttrr)
                        max_r=max[pos_max]
                        aux_space=max_sptr-len(str(max_r))
                        RGB=colors['RGB'].iloc[target_color]
                        alpha=round(abs(scale),1)
                        HTML='<font style="background-color:'+RGB+', '+str(alpha)+');"> </font>'
                        aux_scale+=HTML
                        scale+=0.1

                    if ali_method == 1:
                        html.write(f"{' '*max_z}\t{ttrr} -{max_r}{' '*aux_space} {aux_scale} {max_r}{' '*aux_space}\n")
                    elif ali_method == 2:
                        aux_space=' '*aux_space
                        sc_vis='\t\t'+ttrr+' -'+str(max_r)+aux_space+' '+aux_scale+' '+str(max_r)+aux_space
                        list_sc.append(sc_vis)

            if type (mutation) == type_bool:
                target_color=color_list[3]
                RGB=colors['RGB'].iloc[target_color]
                HTML1='<font style="background-color:'+RGB+', '+'1.0);"> </font>'
                target_color=color_list[2]
                RGB=colors['RGB'].iloc[target_color]
                HTML2='<font style="background-color:'+RGB+', '+'1.0);"> </font>'

                if ali_method == 1:

                    html.write(f"\n")
                    html.write(f"{' '*max_z}\t  acid mutation {HTML1} basic mutation {HTML2}\n\n")

                elif ali_method == 2:

                    if type(max) != type_list:

                        if pklen > 3:
                            ax=pklen-3
                            pklen1=pklen+ax
                        else:
                            pklen1=pklen

                        space=' '*int(pklen1/2)
                        sc_vis='\t\t<small>'+space+'         acid mutation '+HTML1+'</small>'+'\t'
                        list_sc.append(sc_vis)

                        space=' '*int(pklen1/2)
                        sc_vis='\t\t<small>'+space+'        basic mutation '+HTML2+'</small>'+'\t'
                        list_sc.append(sc_vis)

                    else:

                        space='  '
                        sc_vis='\t<small>  '+space+'                     acid mutation '+HTML1+'</small>'+'\t'
                        list_sc.append(sc_vis)
                        sc_vis='\t<small>  '+space+'                    basic mutation '+HTML2+'</small>'+'\t'
                        list_sc.append(sc_vis)

        if descriptor_col == 'ssc' or descriptor_col == 'tc' or descriptor_col == 'tmc':

            target_color=color_list[0]
            RGB=colors['RGB'].iloc[target_color]
            scale=0

            for sc in range (11):
                alpha=round(abs(scale),1)
                HTML='<font style="background-color:'+RGB+', '+str(alpha)+');"> </font>'
                aux_scale+=HTML
                scale+=0.1

            if type(max) != type_list:

                if ali_method == 1:
                    html.write(f"{' '*max_z}\t{min} {aux_scale} &#177;{max}\n\n")
                elif ali_method == 2:
                    sc_vis='\t\t'+str(min)+' '+aux_scale+' &#177;'+str(max)+'\t'
                    list_sc.append(sc_vis)

            else:
                for tt in range(len(tr_r)):

                    if ali_method == 1:
                        html.write(f"{' '*max_z}\t{tr_r[tt]} {min[tt]} {aux_scale} &#177;{max[tt]}\n")
                    elif ali_method == 2:
                        sc_vis='\t\t'+tr_r[tt]+' '+str(min[tt])+' '+aux_scale+' &#177;'+str(max[tt])+'\t'
                        list_sc.append(sc_vis)

            # for tmc
            if type (mutation) == type_bool:

                target_color=color_list[1]
                RGB=colors['RGB'].iloc[target_color]
                HTML='<font style="background-color:'+RGB+', '+'1.0);"> </font>'

                if ali_method == 1:
                    html.write(f"\n")

                    if type(min) == type_list: # with target-residues-separated
                        tar=tr_r[tt]+' '
                        space=' '*len(str(tar))+' '*len(str(min[tt-1]))+' '*1
                    else:
                        tar=''
                        space=' '*len(str(tar))+' '*len(str(min))+' '*1
                    html.write(f"{' '*max_z}\t{space} mutation {HTML}\n\n")

                else: # method 2

                    if type(min) == type_list: # with target-residues-separated
                        space=' '*11
                        sc_vis='\t<small>'+space+'        mutation  '+HTML+'</small>'+'\t'
                        list_sc.append(sc_vis)
                    else:
                        space=' '*11
                        sc_vis='\t<small>'+space+'      mutation  '+HTML+'</small>'+'\t\t'
                        list_sc.append(sc_vis)

        if descriptor_col == 'fsc' or descriptor_col == 'fac' or descriptor_col == 'fmac':

            tt=0

            if descriptor_col == 'fmac':
                aux_mut=1 # To avoid go out of range
            else:
                aux_mut=0

            if descriptor_col == 'fsc' and type (mutation) == type_bool:

                color_aux_fsc=copy.deepcopy(color_list)
                color_list=color_list[:-1]

            for cl in range(len(color_list)-aux_mut):

                target_color=color_list[cl]

                if descriptor_col != 'fsc':
                    tar=tr_r[tt]+' '
                else:
                    tar=''

                RGB=colors['RGB'].iloc[target_color]
                scale=0
                aux_scale=''

                for sc in range (11):
                    alpha=round(abs(scale),1)
                    HTML='<font style="background-color:'+RGB+', '+str(alpha)+');"> </font>'
                    aux_scale+=HTML
                    scale+=0.1

                # When it is not showing residues
                if type(max) != type_list:

                    if ali_method == 1:
                        html.write(f"{' '*max_z}\t{tar}{min} {aux_scale} &#177;{max}\n")

                    elif ali_method == 2:
                        sc_vis='\t\t'+tar+str(min)+' '+aux_scale+' &#177;'+str(max)+'\t'
                        list_sc.append(sc_vis)

                else:

                    # 'fac' and 'fmac'
                    if descriptor_col != 'fsc':

                        if ali_method == 1:
                            html.write(f"{' '*max_z}\t{tar}{min[tt]} {aux_scale} &#177;{max[tt]}\n")
                        elif ali_method == 2:
                            sc_vis='\t\t'+tar+str(min[tt])+' '+aux_scale+' &#177;'+str(max[tt])+'\t'
                            list_sc.append(sc_vis)

                    # 'fsc'
                    else:

                        for ttt in range(len(tr_r)):
                            tar=tr_r[ttt]+' '
                            if ali_method == 1:
                                html.write(f"{' '*max_z}\t{tar}{min[tt]} {aux_scale} &#177;{max[ttt]}\n")
                            if ali_method == 2:
                                sc_vis='\t\t'+tar+str(min[tt])+' '+aux_scale+' &#177;'+str(max[ttt])+'\t'
                                list_sc.append(sc_vis)

                tt+=1
            # When mutation
            if type (mutation) == type_bool:

                # Remember: pimc was converted to fsc!
                if descriptor_col != 'fsc':
                    target_color=color_list[-1]
                else:
                    target_color=color_aux_fsc[-1]

                RGB=colors['RGB'].iloc[target_color]
                HTML='<font style="background-color:'+RGB+', '+'1.0);"> </font>'

                if ali_method == 1:
                    html.write(f"\n")

                    if type(min) == type_list:
                        space=' '*len(str(tar))+' '*len(str(min[tt-1]))+' '*1
                    else:
                        space=' '*len(str(tar))+' '*len(str(min))+' '*1
                    html.write(f"{' '*max_z}\t{space} mutation {HTML}\n\n")

                elif ali_method == 2:

                    if type(max) != type_list:

                        if descriptor_col != 'fsc':
                            space=' '*11
                            sc_vis='\t<small>'+space+'        mutation  '+HTML+'</small>'+'\t\t'
                            list_sc.append(sc_vis)
                        else:
                            space=' '*11
                            sc_vis='\t<small>'+space+'      mutation  '+HTML+'</small>'+'\t\t'
                            list_sc.append(sc_vis)
                    else:
                        space=' '*11
                        sc_vis='\t<small>'+space+'        mutation  '+HTML+'</small>'+'\t\t'
                        list_sc.append(sc_vis)

    else:
        # To avoid big white spaces
        if ali_method == 2:
            white='<font style="background-color:rgb(255, 255, 255);"></font>'
            for i in range(len(tr_r)*len(protein)+5):
                list_sc.append(f'{white}')

    return list_sc

    ## Function for method of visualization 1 ##

def method_visualization1(html,dict_main,met_ref,dict_aux,dict_num,df_col_array,df_col,df_mxr,df_vis,cc,list_big,max_z,protein,n,l,type_bool,type_list,tr_r,conservation,ali_method):

    import copy

    for main, dic_com in dict_main.items():

        for comparison, aligment in dic_com.items():

            if met_ref == 'multiple' or met_ref == 'pair':
                html.write(f'<b><p style="font-size:17px;">Comparison: {comparison}</p></b>\n')

            ref=(list(aligment.keys())[0])

            # Dictionaries for conservation and numbers
            dict_conser=copy.deepcopy(dict_aux)
            dict_num1=copy.deepcopy(dict_num)

            while(len(aligment[ref])>0):

                descriptor_col,color_list,mutation,max,min,vis_deg,vis_ref=read_feature(df_col_array,df_col,df_mxr,df_vis,cc)

                dd=0
                for key in aligment:

                    z=list_big[protein.index(key)]

                    # Reference
                    if dd == 0:

                        aux=dict_num1['number'][:l]
                        dict_num1['number']=dict_num1['number'][l:]
                        html.write(f"{' '*max_z}\t{aux}\n")

                        if vis_ref == True:
                            aux=aligment[key][:n]
                            aligment[key]=aligment[key][n:]
                            html.write(f"{key}{' '*z}\t{aux}\n")
                        else:
                            aux=aligment[key][:n]
                            aligment[key]=aligment[key][n:]
                    # Target
                    else:
                        aux=aligment[key][:n]
                        aligment[key]=aligment[key][n:]
                        html.write(f"{key}{' '*z}\t{aux}\n")

                    dd+=1

                if conservation == False:
                    aux=dict_conser['conserv'][:l]
                    dict_conser['conserv']=dict_conser['conserv'][l:]
                    html.write(f"{' '*max_z}\t{aux}\n\n")
                else:
                    html.write(f"\n")
            cc+=1

            scale(protein,vis_deg,descriptor_col,type_list,min,max,max_z,color_list,html,tr_r,mutation,type_bool,ali_method)
            html.write(f"\n\n")

    ## Function for method of visualization 2 ##

def method_visualization2(html,dict_main,met_ref,dict_aux,dict_num,protein,n,l,list_big,max_z,df_col_array,df_col,df_mxr,df_vis,join,type_list,type_bool,tr_r,conservation,check_dict_vis,ali_method):

    import copy

    aligment_base=list(list(dict_main.values())[0].values())[0]
    key_let=list(list(list(dict_main.values())[0].values())[0].keys())[0]

    # Dictionaries for conservation and numbers
    dict_conser=copy.deepcopy(dict_aux)
    dict_num1=copy.deepcopy(dict_num)

    while(len(aligment_base[key_let])>0):

        cc=0

        aux=dict_num1['number'][:l]
        dict_num1['number']=dict_num1['number'][l:]
        html.write(f"{' '*max_z}\t{aux}\n")

        for main, dic_com in dict_main.items():

            for comparison, aligment in dic_com.items():

                descriptor_col,color_list,mutation,max,min,vis_deg,vis_ref=read_feature(df_col_array,df_col,df_mxr,df_vis,cc)

                dd=0

                for key in aligment:

                    z=list_big[protein.index(key)]

                    # Reference
                    if dd == 0:
                        if vis_ref == True:
                            aux=aligment[key][:n]
                            aligment[key]=aligment[key][n:]
                            html.write(f"{key}{' '*z}\t{aux}  <b>comparison {comparison}</b>\n")
                        else:
                            aux=aligment[key][:n]
                            aligment[key]=aligment[key][n:]
                    # Target
                    else:
                        if vis_ref == True:
                            aux=aligment[key][:n]
                            aligment[key]=aligment[key][n:]
                            html.write(f"{key}{' '*z}\t{aux}\n")
                        else:
                            if dd != 1:
                                aux=aligment[key][:n]
                                aligment[key]=aligment[key][n:]
                                html.write(f"{key}{' '*z}\t{aux}\n")
                            else:
                                aux=aligment[key][:n]
                                aligment[key]=aligment[key][n:]
                                html.write(f"{key}{' '*z}\t{aux}  <b>comparison {comparison}</b>\n")
                    dd+=1
                cc+=1

                if join == False:
                    html.write(f"\n")

        if conservation == False:
            aux=dict_conser['conserv'][:l]
            dict_conser['conserv']=dict_conser['conserv'][l:]
            html.write(f"{' '*max_z}\t{aux}\n\n")
        else:
            if join == False:
                pass
            else:
                html.write(f"\n")

    html.write(f"\n")

    ## Scale part ##

    cc=0
    dict_scale={}
    for main, dic_com in dict_main.items():

        for comparison, aligment in dic_com.items():

            ref=(list(aligment.keys())[0])
            dict_num1=copy.deepcopy(dict_num)
            descriptor_col,color_list,mutation,max,min,vis_deg,vis_ref=read_feature(df_col_array,df_col,df_mxr,df_vis,cc)
            dict_scale[cc]=''
            scale_aux=scale(protein,vis_deg,descriptor_col,type_list,min,max,max_z,color_list,html,tr_r,mutation,type_bool,ali_method)
            dict_scale[cc]=scale_aux
            cc+=1

    # To distributed each scale

    # Reading the máximum number of rows in the scale
    rsc=max_sc(dict_scale)

    if check_dict_vis:
        factor_sp=32
    else:
        factor_sp=32

    for s in range(rsc):
        # Reading each scale of each comparison
        for sc in dict_scale.values():
            try:
                html.write(f"{sc[s]}")
            # Because some scale does not have the same number of rows it is added a space
            except IndexError:
                html.write(f"{' '*factor_sp}")
                #html.write(f"\t\t\t\t")
        html.write(f"\n")
