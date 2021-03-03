from cprisma.treatment_alignment import *
from cprisma.treatment_data import *
from cprisma.verify_data import *
from cprisma.matrix_generator import *
from cprisma.verify_array import verify_array_feature
from cprisma.interpret_array import interpret_array
import pandas as pd
import copy

                #### MAIN CLASS ####

class Information:
    pass
    def __init__(self, turn_com, ali_log, check_tar_res, target_residues, accuracy, check_name_seq, name_sequence):

                                            # Private
        # Name input data
        self.__name_dat='./data_input.csv'
        # Reading input data as a DataFrame
        self.__df_dat=pd.read_csv(self.__name_dat)
        # General statistics description of input data
        self.__df_stat=self.__df_dat.describe()
        # Alignment from https://www.ebi.ac.uk/Tools/msa/muscle/ in output forma ClustalW
        self.__name_ali='./alignment.dat'
        # To verify later the quality of all information
        self.__state_dat=False
        # Checking if is introduced a change in the name
        self.__ck_ns=check_name_seq
                                                # Public
        # Output data
        self.ali_log=ali_log
        # Turn visualization of comparison in log file
        self.turn_com=turn_com
        # Target residues to evaluate
        self.tr_r=target_residues
        # Checking presence of target residues
        self.ck_tr=check_tar_res
        # Accuracy of my data
        self.acc=accuracy
        # Name of each protein
        self.name_sequence=name_sequence
        # Ionizable residues
        self.ch_r=['D','E','R','H','K','Y','C']
        # Hydrophobic residues
        self.hy_r=['A','I','L','M','F','P','W','V']

    def treatment(self,counter):
        # Reading alignment
        self.dict_ali,self.protein=treatment_alignment(self.turn_com,
                                                       self.__name_ali,
                                                       self.ali_log,
                                                       self.ch_r,
                                                       self.hy_r,
                                                       self.__ck_ns,
                                                       self.name_sequence,
                                                       counter)

        if counter == 0:
            print('Sequences: ',self.protein,'.')
        # Verification of data input
        self.__state_dat=verify_data(self.ali_log,
                                     self.__df_dat,
                                     self.protein,
                                     self.dict_ali,
                                     self.__name_dat,
                                     self.__state_dat,
                                     counter)

        # If data input is ok the program continue
        if self.__state_dat:

            # Creating dictionaries base on the data input
            self.dict_aa,self.dict_da,self.head_all_df,self.__df_stat=treatment_data(self.__name_dat,
                                                                                     self.__df_dat,
                                                                                     self.__df_stat,
                                                                                     self.protein,
                                                                                     self.ali_log,
                                                                                     counter)

            # Comparing data between aligment and data input using as a bases target sequence
            self.__state_dat,self.tr_r=compare_data(self.dict_aa,
                                                    self.dict_da,
                                                    self.protein,
                                                    self.dict_ali,
                                                    self.tr_r,
                                                    self.ck_tr,
                                                    self.__name_dat,
                                                    self.__name_ali,
                                                    self.__state_dat,
                                                    self.ali_log,
                                                    counter)

            # Aligning input data dictionary using alignment and target residues as basis
            (self.dict_aa,self.dict_aa3,self.dict_da,
             self.df_aa,self.df_aa3,self.df_da,self.frames)=adjusting_data(self.__state_dat,
                                                                           self.dict_ali,
                                                                           self.tr_r,
                                                                           self.dict_aa,
                                                                           self.dict_da,
                                                                           self.acc,
                                                                           self.protein,
                                                                           self.head_all_df,
                                                                           self.__df_stat,
                                                                           self.ali_log,
                                                                           counter)

    def reference(self,check_reference,dict_ref,method_reference,counter):

        self.ck_ref=check_reference
        self.dict_ref=dict_ref
        self.met_ref=method_reference
        self.df_bas_ref,self.df_bas_seq,self.dict_ref=matrix_generator(self.turn_com,
                                                                       self.ck_ref,
                                                                       self.dict_ref,
                                                                       self.met_ref,
                                                                       self.dict_da,
                                                                       self.df_da,
                                                                       self.protein,
                                                                       self.ali_log,
                                                                       counter)
        return self.df_bas_ref, self.df_bas_seq, self.dict_ref

                #### VERIFY ARRAY CLASS ####

class VerifyArray(Information):

    def __init__(self,turn_com,ali_log,check_tar_res,target_residues,accuracy,check_name_seq,name_sequence,
                 feature, default_parameter, dict_array, ck_dict_arr):
        super().__init__(turn_com,ali_log,check_tar_res,target_residues,accuracy,check_name_seq,name_sequence)
        self.feature=feature
        self.default_parameter=default_parameter
        self.dict_array=dict_array
        self.ck_dict_arr=ck_dict_arr

    def reference(self,check_reference,dict_ref,method_reference,counter):

        super().reference(check_reference,dict_ref,method_reference,counter)
        self.dict_ref_alter=copy.deepcopy(self.dict_ref)
        self.dict_array,self.df_array,self.dict_ref=verify_array_feature(self.default_parameter,
                                                                         self.dict_array,
                                                                         self.ck_dict_arr,
                                                                         self.feature,
                                                                         self.ck_ref,
                                                                         self.dict_ref_alter,
                                                                         self.met_ref,
                                                                         self.ch_r,
                                                                         self.tr_r,
                                                                         self.ck_tr,
                                                                         self.df_da,
                                                                         self.ali_log,
                                                                         counter)

        return self.dict_array, self.df_array, self.dict_ref

                #### OPERATION CLASS ####

class Operation(VerifyArray):
    pass

    def matrix_operation(self,df_bas_ope):

        self.df_bas_ope=df_bas_ope
        if self.feature == 'operation':
            self.df_bas_ope=interpret_array(self.turn_com,
                                            self.tr_r,
                                            self.dict_ref,
                                            self.df_aa,
                                            self.df_da,
                                            self.df_bas_ref,
                                            self.df_bas_seq,
                                            self.df_array,
                                            self.feature,
                                            self.protein,
                                            self.df_bas_ope,
                                            self.acc,
                                            self.ali_log)
        return self.df_bas_ope

                #### MAXIMUM CLASS ####

class Maximum(Operation):
    pass

    def matrix_operation(self,df_bas_ope):

        super().matrix_operation(df_bas_ope)
        self.df_bas_ope=df_bas_ope

    def matrix_maximum(self):

        self.df_bas_mxr=interpret_array(self.turn_com,
                                        self.tr_r,
                                        self.dict_ref,
                                        self.df_aa,
                                        self.df_da,
                                        self.df_bas_ref,
                                        self.df_bas_seq,
                                        self.df_array,
                                        self.feature,
                                        self.protein,
                                        self.df_bas_ope,
                                        self.acc,
                                        self.ali_log)
        return self.df_bas_mxr

        #### COLOR CLASS ####

class Color(Operation):
    pass
    def matrix_operation(self,df_bas_ope):

        super().matrix_operation(df_bas_ope)
        self.df_bas_ope=df_bas_ope

    def matrix_color(self, mutation_col, sequence_col):

        self.list_op_col=[mutation_col, sequence_col, self.df_bas_ope]
        self.df_bas_col=interpret_array(self.turn_com,
                                        self.tr_r,
                                        self.dict_ref,
                                        self.df_aa,
                                        self.df_da,
                                        self.df_bas_ref,
                                        self.df_bas_seq,
                                        self.df_array,
                                        self.feature,
                                        self.protein,
                                        self.list_op_col,
                                        self.acc,
                                        self.ali_log)

        return self.df_bas_col, self.df_array

        #### VISUALIZATION CLASS ####

class Visualization(VerifyArray):
    pass
    def matrix_visualization(self):

        self.ali_log.write(f"Matrices for {self.feature} were done!\n\n")
        print(f"Matrices for {self.feature} were done!")

        return self.dict_ref, self.df_array

        #### ALIGNMENT CLASS ####

class Alignment(Information):

    def __init__(self, turn_com, ali_log, check_tar_res, target_residues, accuracy, check_name_seq, name_sequence,
                 method_reference, dict_ref, df_ope, df_mxr, df_col, df_col_array, df_vis):

        super().__init__(turn_com,ali_log,check_tar_res,target_residues,accuracy,check_name_seq,name_sequence)

        self.met_ref=method_reference
        self.dict_ref=dict_ref
        self.df_ope=df_ope
        self.df_mxr=df_mxr
        self.df_col=df_col
        self.df_col_array=df_col_array
        self.df_vis=df_vis

    def add_color(self, factor_color, first_num, join, conservation, res_per_col, tridimensional_col, check_dict_vis, ali_method):

        from cprisma.compile import compile

        self.factor_color=factor_color
        self.first_num=first_num
        self.join=join
        self.conservation=conservation
        self.res_per_col=res_per_col
        self.tridimensional_col=tridimensional_col
        self.check_dict_vis=check_dict_vis
        self.ali_method=ali_method
        self.aligment_col=compile(self.met_ref,
                                  self.dict_ali,
                                  self.tr_r,
                                  self.df_aa,
                                  self.protein,
                                  self.dict_ref,
                                  self.df_ope,
                                  self.df_mxr,
                                  self.df_col,
                                  self.df_col_array,
                                  self.df_vis,
                                  self.factor_color,
                                  self.first_num,
                                  self.join,
                                  self.conservation,
                                  self.ali_method,
                                  self.tridimensional_col,
                                  self.check_dict_vis,
                                  self.res_per_col,
                                  self.ali_log)

        print(f'Process completed!')

        from time import gmtime, strftime
        var_time=strftime("%a, %d %b %Y %H:%M:%S", gmtime())
        self.ali_log.write(f'Process completed! [{var_time}]\n\n')
