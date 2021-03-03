from cprisma.interpret_residue import interpret_residue3
import numpy as np
import pandas as pd
import copy

    ########## Creating dictionary of input data ##########

def treatment_data(name_dat,df_dat,df_stat,protein,ali_log,counter):

    dic_aa={}
    dic_da={}
    list_he=[]

    #### Creating dictionary of data_input and change name of header ####

    cc=0
    for p in protein:
        pdat=p+'_dat'
        unm="Unnamed: "
        unm+=str(cc)
        df_dat[unm] = df_dat[unm].replace(np.nan, 'XXX') # Convert nan to XXX
        df_dat.rename(columns={unm:p}, inplace=True)
        cc+=1
        unm="Unnamed: "
        unm+=str(cc)
        df_dat[unm] = df_dat[unm].replace(np.nan, 0) # Convert nan to 0
        df_dat.rename(columns={unm:pdat}, inplace=True)
        cc+=1
        list_he.append(p)
        list_he.append(pdat)
        dic_aa[p]=''
        dic_da[pdat]=''

    cc=1
    for p in protein:
        pdat=p
        unm="Unnamed: "
        unm+=str(cc)
        df_stat.rename(columns={unm:p}, inplace=True)
        cc+=2

    head_all_df=list(df_dat.head(0))   # Selecting head data frame

    if counter == 0:
        ali_log.write(f"{' '*15}--- Information in {name_dat} ---\n\n")

    df_all=df_dat.to_string(index=False)

    if counter == 0:
        ali_log.write(f" {df_all} \n\n")

    #### Adding values to dictionary ####

    cc=0
    dd=0
    for l in list_he:
        for d in dic_aa:
            if l==d:
                df_ls=df_dat[l].tolist()
                dic_aa.update({list(dic_aa.keys())[cc]:df_ls})
                cc+=1
        for d in dic_da:
            if l==d:
                df_ls=df_dat[l].tolist()
                dic_da.update({list(dic_da.keys())[dd]:df_ls})
                dd+=1

    return dic_aa, dic_da, head_all_df, df_stat

    ##########  Aligning input data dictionary using alignment and target residues as basis ##########

def adjusting_data(bol,dic_ali,tr_r,dic_aa,dic_da,acc,protein,head_all_df,df_stat,ali_log,counter):

    if bol:

        i=0
        j=0
        k=0
        list_yy=[]

        for i in range (0,(len(list(dic_ali.values())[1]))):
            list_y=[]
            for x,y in dic_ali.items():
                list_y.append((y)[i])
            if any(a in tr_r for a in (list_y)):            # List of residues from alignment
                i1=i+1
                j=j+1
                list_vv=[]
                list_bb=[]
                list_yy=list_y[1:]
                for key in dic_aa:
                    list_v=list(dic_aa[key])
                    list_vv.append(list_v[k])
                for v in range(0,len(list_vv)):             # Fitting data values
                    am=list_vv[v]
                    if am != list_yy[v]:
                        list_z = (list(dic_aa.values())[v][:k] + list('X') + list(dic_aa.values())[v][k:])
                        dic_aa.update({list(dic_aa.keys())[v]:list_z})
                        list_z = (list(dic_da.values())[v][:k] + list('0') + list(dic_da.values())[v][k:])
                        list_z=[round(float(item),acc) for item in list_z]
                        dic_da.update({list(dic_da.keys())[v]:list_z})
                k=k+1

        #### We find the minimum lenght in dic_aa ####

        list_min=[]
        for d in dic_aa.values():
            list_min.append(len(d))
        mind=min(list_min)

        if counter == 0:
            ali_log.write(f"{' '*25}--- Data processed ---\n\n")

        #### With the new minimum lenght it is fitted each dictionary to have the same lenght for dataframe later ####

        v=0
        aux="Gaps generated for "
        for d in dic_aa.values():
            list_n=d[:mind]
            cx=list_n.count('X') # Function to count the number of gaps
            protein1=protein[v]
            if v != len(list(dic_aa.keys()))-1:
                aux+=protein1+"="+str(cx)+", "
            dic_aa.update({list(dic_aa.keys())[v]:list_n})
            v+=1
        aux+=protein1+"="+str(cx)+".\n\n"

        v=0
        for d in dic_da.values():
            list_n=d[:mind]
            dic_da.update({list(dic_da.keys())[v]:list_n})
            v+=1

        #### Creating new dictionary using format of three letters aminoacid ####

        dic_aa3={}
        cc=0
        for p in protein:
            dic_aa3[p]=[]
            l_aa=copy.copy(list(dic_aa.values())[cc])       # Function 'copy' to avoid simultaneously  change in dic_aa
            dd=0
            for a in l_aa:
                am=interpret_residue3(a)
                l_aa[dd]=am
                dd+=1
            dic_aa3.update({list(dic_aa3.keys())[cc]:l_aa})
            cc+=1

        #### Creating data frame for all dictionaries ####

        df_da=pd.DataFrame.from_dict(dic_da)
        df_aa=pd.DataFrame.from_dict(dic_aa)
        df_aa3=pd.DataFrame.from_dict(dic_aa3)

        frames = [df_aa3, df_da]
        frames=(pd.concat(frames, axis=1, sort=False))
        frames=frames[head_all_df]      # To put data in the same order of input data file
        frames=frames.to_string(index=True)

        if counter == 0:
            ali_log.write(f"{frames}\n\n")
            ali_log.write(f"{aux}")

            ali_log.write(f"{' '*18}--- Statistics ---\n\n")
            df_stat1=df_stat.to_string(index=True)
            ali_log.write(f"{df_stat1}\n\n")
            print("Data input processed!")

        return dic_aa, dic_aa3, dic_da, df_aa, df_aa3, df_da, frames
