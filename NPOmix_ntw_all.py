#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import re
import numpy as np
import glob
import os
from collections import defaultdict
import networkx
from networkx.algorithms.components.connected import connected_components
from packages import npomix
import time
from datetime import datetime
from sklearn.metrics import jaccard_score
import csv #new
import matplotlib.pyplot as plt #new



start = time.time()


#GNPS files
    #networkedges_selfloop
edges_path = "/home/andres/Documents/Lab_Actino/UCSD_work/1000_Genomes/NPOmix/NPOmix_data/GNPS_321/networkedges_selfloop/ca909324641045abbfa26b6d21c74c1d..selfloop"
    #clusterinfosummarygroup file
nodes_path = "/home/andres/Documents/Lab_Actino/UCSD_work/1000_Genomes/NPOmix/NPOmix_data/GNPS_321/clusterinfosummarygroup_attributes_withIDs_withcomponentID/d659f4f36de74a7ebdefe626d58024d3.clustersummary"
    #result_specnets_DB
hits_path = pd.read_csv("/home/andres/Documents/Lab_Actino/UCSD_work/1000_Genomes/NPOmix/NPOmix_data/GNPS_321/result_specnets_DB/b527f46471ad41cca3a503816f82d3ef.tsv",sep='\t')

#BGCs files
    #BiG-SCAPE network
input_bigscape_net = "/home/andres/Documents/Lab_Actino/UCSD_work/1000_Genomes/NPOmix/NPOmix_data/Bigscape/output_bigscape_mibig31_filteredclusters321/bigscape_all_c030_321strains.txt"
    #BiG-SCAPE annotation table
input_bigscape_ann = "/home/andres/Documents/Lab_Actino/UCSD_work/1000_Genomes/NPOmix/NPOmix_data/Bigscape/output_bigscape_mibig31_filteredclusters321/Network_Annotations_Full.tsv" #new
    #antiSMASH results folder
antismash_folder = "/home/andres/Documents/Lab_Actino/UCSD_work/1000_Genomes/NPOmix/NPOmix_data/antismash_only_gbk/"

#NPOmix files
ena_df_file = "/home/andres/Documents/Lab_Actino/UCSD_work/1000_Genomes/NPOmix/NPOmix_python/ena_dict-210315.csv"
    #matched MIBiG BGCs with spec lib
mibig_df = pd.read_csv("/home/andres/Documents/Lab_Actino/UCSD_work/1000_Genomes/NPOmix/NPOmix_python/matched_mibig_gnps_update.tsv",sep='\t')
    #MIBiG 3.1 main product table
mibig_prod = pd.read_csv("MIBiG_products.csv",sep=',') #new

#Specify results folder
results_folder = '/home/andres/Documents/Lab_Actino/UCSD_work/1000_Genomes/NPOmix/NPOmix_data/20230127_notderep_321_edited'

current_date = datetime.today().strftime('%Y%m%d')

if not os.path.isdir(results_folder):
    os.mkdir(results_folder)


nodes_df = pd.read_csv(nodes_path,sep='\t')

nodes_df[:5]


clusterindex_list = []

for i,r in nodes_df.iterrows():
# remove this line to use everything
#    if type(r['LibraryID']) == float:
        clusterindex_list.append(r['cluster index'])
        
clusterindex_list


edges_df = pd.read_csv(edges_path,sep='\t')

def get_neighbors(target,dataframe,column1,column2):
    subset1 = dataframe[(dataframe[column1]==target)]
    subcat = subset1.append(dataframe[(dataframe[column2]==target)])
    temp_list = []
    for index,row in subcat.iterrows():
        temp_list.append(subcat[column1][index])
        temp_list.append(subcat[column2][index])
    temp_list = list(np.unique(temp_list))
    return temp_list

def to_edges(l):
    it = iter(l)
    last = next(it)
    for current in it:
        yield last, current
        last = current

def to_graph(l):
    G = networkx.Graph()
    for part in l:
        G.add_nodes_from(part)
        G.add_edges_from(to_edges(part))
    return G

def get_family_dict(components_list,dataframe,dictionary,column1,column2,column3):
    count = 0
    for family in list(components_list):
        count += 1
        for fam_member in family:
            dictionary['MF%s'%count].append(fam_member)
    return dictionary

def main_get_families(gnps_df):
    targets_list = np.unique([gnps_df.CLUSTERID1,gnps_df.CLUSTERID2])
    neighbors_list = []
    for target in targets_list:
        neighbors_list.append(get_neighbors(target,gnps_df,'CLUSTERID1','CLUSTERID2'))
    G = to_graph(neighbors_list)
    C = connected_components(G)
    mf_dict = defaultdict(list)
    mf_dict = get_family_dict(C,gnps_df,mf_dict,'CLUSTERID1','CLUSTERID2','Cosine')
    return mf_dict

mf_dict = main_get_families(edges_df)


for ion in clusterindex_list:
    print(ion)


#Get a table with MFs and cluster indexes
mf_df = pd.DataFrame.from_dict(mf_dict, orient='index') #convert from dict list to df
mf_dfwide = pd.melt(mf_df, ignore_index=False, id_vars=None) #melt to long format
#mf_dfwide = mf_dfwide.drop(['data','index','index_names']) #remove rows with data, index and index_names
mf_dfwide = mf_dfwide.drop(columns=['variable']) #remove columns
mf_dfwide = mf_dfwide.dropna() #remove rows with NaN
mf_dfwide = mf_dfwide.sort_values(by=['value']) #sort by cluster index
mf_dfwide.reset_index(inplace=True) #set index as column
mf_dfwide = mf_dfwide.rename(columns={'index':'MFs','value':'cluster index'}) #rename columns


nodes_df2= nodes_df[["cluster index","DefaultGroups","UniqueFileSources","parent mass","sum(precursor intensity)","LibraryID"]]

#Concatenate Node tables with MFs
nodes_mf = pd.merge(nodes_df2, mf_dfwide, right_index = True, left_index = True)
nodes_mf  = nodes_mf.drop(columns=['cluster index_y']) #remove columns
nodes_mf  = nodes_mf.rename(columns={'cluster index_x':'metabolite_ID','UniqueFileSources':'mzML_file_sources','sum(precursor intensity)':'sum_intensity','parent mass':'parent_mass'}) #rename columns
nodes_mf = nodes_mf[["metabolite_ID","MFs","DefaultGroups","mzML_file_sources","parent_mass","LibraryID"]]
#include GNPS compound classes
hits_path2 = hits_path[["#Scan#","npclassifier_superclass","npclassifier_class","npclassifier_pathway"]] #select only scan# and npclassifier columns
hits_path2  = hits_path2.rename(columns={'#Scan#':'metabolite_ID'}) #rename columns
nodes_mf = pd.merge(nodes_mf, hits_path2, how = 'left' ,on='metabolite_ID')
nodes_mf.to_csv(os.path.join(results_folder,'NodesMF_dic-NPOmix1.0-%s.txt'%(current_date)))



lcms_file_list = []

for item in nodes_df['UniqueFileSources']:
    for lcms_file in item.split('|'):
        if lcms_file not in lcms_file_list:
            print(lcms_file)
            lcms_file_list.append(lcms_file)


len(lcms_file_list)


def get_presence_files(nodes_df):
    for unique_list in list(nodes_df[nodes_df['cluster index'] == ion]['UniqueFileSources']):
        return list(unique_list.split('|'))

all_rows_list,testing_indexes_list = [],[]
for ion in clusterindex_list:
    subset_edges = edges_df[(edges_df.CLUSTERID1 == ion) | (edges_df.CLUSTERID2 == ion)]
    cosine = round(max(subset_edges['Cosine']),2)
    presence_list = get_presence_files(nodes_df)
    single_row_list = []
    for lcms_file in lcms_file_list:
        if lcms_file in presence_list:
            single_row_list.append(cosine)
        else:
            single_row_list.append(0)
    all_rows_list.append(single_row_list)
    testing_indexes_list.append(ion)
    
all_rows_list

pre_testing_df = pd.DataFrame(all_rows_list,index=testing_indexes_list,columns=lcms_file_list)

pre_testing_df


def get_final_df(training_df,testing_df,neighbors_array,results_folder):
    final_df = pd.DataFrame(columns=('metabolite_ID','predicted_GCFs','max_jaccard','jaccard_scores','parent_mass','sum_intensity')) #added jaccard_scores
    for i,ccms_id in enumerate(testing_df.index):
        jaccard_scores = []
        for j in range(0,len(neighbors_array[i])):
            query_bgc = neighbors_array[i][j]
            bgc_fp = training_df[training_df['label'] == query_bgc].iloc[0]
            bgc_fp = bgc_fp.drop("label")
            ms_fp = testing_df.loc[ccms_id]
            bgc_binary = npomix.get_binary(bgc_fp)
            ms_binary = npomix.get_binary(ms_fp)
            jaccard_scores.append(jaccard_score(bgc_binary,ms_binary))
        max_jaccard = round(max(jaccard_scores),2)
        parentmass = nodes_df[nodes_df['cluster index'] == ccms_id]['precursor mass'].item()
        intensity = nodes_df[nodes_df['cluster index'] == ccms_id]['sum(precursor intensity)'].item()
        final_df.loc[i] = ccms_id,neighbors_array[i],max_jaccard,jaccard_scores,parentmass,intensity #added jaccard_scores
    #final_df.to_csv("%sfinal_df-noderep-NPOmix1.0-%s.txt"%(results_folder,current_date),sep="\t",index_label=False)
    return final_df


k_value = 3

merged_ispec_mat = npomix.get_merged_ispec_mat(pre_testing_df)
merged_ispec_mat = npomix.renaming_merged_ispec_mat(ena_df_file,merged_ispec_mat)
print('Obtaining BiG-SCAPE dataframe and BiG-SCAPE dictionary')
bigscape_df,bigscape_dict = npomix.get_bigscape_df(ena_df_file,input_bigscape_net)
bigscape_df,bigscape_dict2 = npomix.rename_bigscape_df(antismash_folder,bigscape_df,bigscape_dict)
print('BiG-SCAPE create with %s GCFs'%len(bigscape_dict))
strain_list,bgcs_list = npomix.get_strain_list(bigscape_df)
print('Creating training dataframe')
affinity_df,affinity_bgcs = npomix.get_pre_training_df(bigscape_df,bigscape_dict2,strain_list,bgcs_list)
affinity_df = npomix.renaming_affinity_df(affinity_df)
networked_cols = npomix.get_networked_cols(merged_ispec_mat,affinity_df)
training_df,training_bgcs = npomix.get_training_df(affinity_df,networked_cols,results_folder,affinity_bgcs)
bgcs_df = pd.DataFrame(training_bgcs, columns=['bgcs'])
testing_df = npomix.get_testing_df(merged_ispec_mat,networked_cols,results_folder)
print('Running KNN using k equals to %s'%k_value)
neighbors_array = npomix.running_knn(training_df,testing_df,k_value)
print('Creating final dataframe')
final_df = get_final_df(training_df,testing_df,neighbors_array,results_folder)


final_df = final_df[final_df['max_jaccard'] > 0]
final_df = final_df.sort_values(by='max_jaccard',ascending=False)
final_df = final_df.reset_index(drop=True)

testing_df.shape


def getListFiles(dirName):
    listFile = os.listdir(dirName)
    allFiles = list()
    for entry in listFile:
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles

BGC_list= getListFiles(antismash_folder)

#convert from BGC list to df
BGC_df = pd.DataFrame(BGC_list,columns=['BGCs'])
BGC_df = BGC_df.BGCs.str.split("/", expand = True) #separate all the folders in columns
BGC_df = BGC_df.iloc[:,-2:] #select only the last two columns, strain folder and the BGCs present in that folder
BGC_df = BGC_df.rename(columns={10:'strain', 11:'BGC'}) #rename columns Needs to be improved so it changes for any column name
BGC_df = BGC_df[["BGC","strain"]]

#Get a table with the folder containing BGCs
strain_df = BGC_df.strain.str.split(".", expand = True) #separate folder name by a point
strain_df = strain_df.rename(columns={0:'strain'}) #rename columns
strain_df = strain_df[['strain']]

#Merge table with BGCs and respective strains
BGC_df = BGC_df[["BGC"]]
BGC_df = pd.merge(BGC_df, strain_df, right_index = True, left_index = True)
BGC_df = BGC_df.replace('.gbk','',regex=True)

#Get a table with GCFss and BGCs
GCF_df = pd.DataFrame.from_dict(bigscape_dict, orient='index') #convert from dict list to df
GCF_dfwide = pd.melt(GCF_df, ignore_index=False, id_vars=None) #melt to long format
GCF_dfwide = GCF_dfwide.drop(columns=['variable']) #remove columns
GCF_dfwide = GCF_dfwide.dropna() #remove rows with NaN
GCF_dfwide = GCF_dfwide.sort_values(by=['value']) #sort by cluster index
GCF_dfwide.reset_index(inplace=True) #set index as column
GCF_dfwide = GCF_dfwide.rename(columns={'index':'GCFs','value':'BGC'}) #rename columns

#concatenate GCFs with BiG-SCAPE annotation table
ann_df = pd.read_table(input_bigscape_ann, delimiter='\t')
ann_df = ann_df.rename(columns={'BiG-SCAPE class':'BiG-SCAPE_class'})
ann_df = ann_df[["BGC","BiG-SCAPE_class"]]
ann_df = pd.merge(GCF_dfwide, ann_df, how = 'right' ,on='BGC')
ann_df = pd.merge(ann_df, mibig_prod, how = 'left' ,on='BGC')
ann_df = ann_df[["BGC","GCFs","MIBiG_main_product","BiG-SCAPE_class"]]

#Merge with table with BGCs and strain info
GCF_BGC = pd.merge(ann_df, BGC_df, how = 'outer' ,on='BGC')
GCF_BGC.GCFs.fillna('singleton', inplace=True)
GCF_BGC.to_csv(os.path.join(results_folder,'Bigscape_dic-NPOmix1.0-%s.txt'%(current_date)))

mibig_name_dict = dict(zip(mibig_df['mibig_id'],mibig_df['mibig_name']))

ccmsid_mibig_dict = dict(zip(mibig_df['# mgf_spectrum_id'],mibig_df['mibig_id']))

ccmsid_mibig_dict

subset_final_df = final_df[final_df['max_jaccard'] >= 0.7]

#Concatenate Node tables with predicted_GCFs
mf_gcf = pd.merge(nodes_mf, subset_final_df, left_on='metabolite_ID', right_on='metabolite_ID')
mf_gcf  = mf_gcf.drop(columns=['parent_mass_y']) #remove columns
mf_gcf  = mf_gcf.rename(columns={'parent_mass_x':'parent_mass'}) #rename columns
mf_gcf.to_csv(os.path.join(results_folder,'mf_gcf_all_dic-NPOmix1.0-%s.txt'%(current_date)))

#Define a table without the training dataset (G1)
mf_gcf_g2 = mf_gcf
mf_gcf_g2 = mf_gcf_g2[mf_gcf_g2.DefaultGroups != "G1"]
mf_gcf_g2.to_csv(os.path.join(results_folder,'mf_gcf_g2_dic-NPOmix1.0-%s.txt'%(current_date)))

#Create a table with MFs counts
MFs_count = nodes_mf['MFs'].value_counts()
MFs_count = pd.DataFrame.from_dict(MFs_count) #convert from dict list to df
MFs_count.reset_index(inplace=True) #set index as column
MFs_count = MFs_count.rename(columns={'MFs':'count','index':'MFs'}) #rename columns

#Define a list of MFs
MFs_list = []
for i,r in nodes_mf.iterrows():
    if type(r['LibraryID']) != float:
        MFs_list.append(r['MFs'])
#Define a list of Library ID
Lib_list = []        
for i,r in nodes_mf.iterrows():
    if type(r['LibraryID']) != float:        
        Lib_list.append(r['LibraryID'])

#Convert to dataframe and merge MFs and LibraryIDs        
MFs_list = pd.DataFrame.from_dict(MFs_list) #convert from dict list to df 
MFs_list = MFs_list.rename(columns={0:'MFs'})
Lib_list = pd.DataFrame.from_dict(Lib_list) #convert from dict list to df 
Lib_list = Lib_list.rename(columns={0:'LibraryIDs'})
MFs_list = pd.concat([MFs_list, Lib_list], axis=1).drop_duplicates() #concatenate and drop duplicates

#Group MFs with more than one annotation using a ;
MFs_list2 = MFs_list.groupby(['MFs'])['LibraryIDs'].apply('; '.join).reset_index()

#Merge with MFs counts
MFs_list2 = pd.merge(MFs_count, MFs_list2, how='left', on='MFs')

#Create a table with information about MFs and default GNPS groups. G1 training dataset, G2 other dataset
Groups = nodes_mf[["MFs","DefaultGroups"]].drop_duplicates()
Groups = Groups.groupby(['MFs'])['DefaultGroups'].apply('; '.join).reset_index()

#Merge with MFs counts
MFs_list2 = pd.merge(MFs_list2, Groups, how='left', on='MFs')
#Replace everything other than G1 or G2 with G1,G2
MFs_list2['DefaultGroups'] = np.where(MFs_list2['DefaultGroups'].isin(['G1','G2']), MFs_list2['DefaultGroups'], 'G1,G2')
#Add a column that indicates that these are MF
MFs_list2.insert(4,'Type','MF')
MFs_list2 = MFs_list2.rename(columns={'LibraryIDs':'Annotation'})

MFs_list2.to_csv(os.path.join(results_folder,'Network_metadata_MFs_list_all-NPOmix1.0-%s.txt'%(current_date)))

#Create a table with GCFs counts
GCFs_count = GCF_BGC['GCFs'].value_counts()
GCFs_count = pd.DataFrame.from_dict(GCFs_count) #convert from dict list to df
GCFs_count.reset_index(inplace=True) #set index as column
GCFs_count = GCFs_count.rename(columns={'GCFs':'count','index':'GCFs'}) #rename columns
GCFs_count = GCFs_count[GCFs_count['GCFs'] != 'singleton']
#Define a list of GCFs with MIBiG hits
GCFs_list = []
for i,r in GCF_BGC.iterrows():
    if type(r['MIBiG_main_product']) != float:
        GCFs_list.append(r['GCFs'])
#Define a list of MIBiG hits
MIBiG_list = []        
for i,r in GCF_BGC.iterrows():
    if type(r['MIBiG_main_product']) != float:        
        MIBiG_list.append(r['MIBiG_main_product'])
        
#Convert to dataframe and merge GCFs and MIBiG main products        
GCFs_list = pd.DataFrame.from_dict(GCFs_list) #convert from dict list to df 
GCFs_list = GCFs_list.rename(columns={0:'GCFs'})
MIBiG_list = pd.DataFrame.from_dict(MIBiG_list) #convert from dict list to df 
MIBiG_list = MIBiG_list.rename(columns={0:'Annotation'})
GCFs_list = pd.concat([GCFs_list, MIBiG_list], axis=1).drop_duplicates() #concatenate and drop duplicates       
GCFs_list = GCFs_list[GCFs_list['GCFs'] != 'singleton']

#Group GCFs with more than one annotation using a ;
GCFs_list = GCFs_list.groupby(['GCFs'])['Annotation'].apply('; '.join).reset_index()
#Merge with MFs counts
GCFs_list = pd.merge(GCFs_count, GCFs_list, how='left', on='GCFs')
GCFs_list.insert(3,'Type','GCF')


#Counts by BGCs classes
#NRPS
NRPS_BGC = GCF_BGC[GCF_BGC['BiG-SCAPE_class'] == 'NRPS'] #filter only NRPSs from table
NRPS_count = NRPS_BGC['GCFs'].value_counts().to_frame() #count
NRPS_count.reset_index(inplace=True) #set index as column
NRPS_count = NRPS_count.rename(columns={'GCFs':'NRPS','index':'GCFs'}) #rename columns
NRPS_count = NRPS_count[NRPS_count['GCFs'] != 'singleton'] #remove singletons
GCFs_list = pd.merge(GCFs_list,NRPS_count, how='left', on='GCFs')

#others
Others_BGC = GCF_BGC[GCF_BGC['BiG-SCAPE_class'] == 'Others'] #filter only NRPSs from table
Others_count = Others_BGC['GCFs'].value_counts().to_frame() #count
Others_count.reset_index(inplace=True) #set index as column
Others_count = Others_count.rename(columns={'GCFs':'Others','index':'GCFs'}) #rename columns
Others_count = Others_count[Others_count['GCFs'] != 'singleton'] #remove singletons
GCFs_list = pd.merge(GCFs_list,Others_count, how='left', on='GCFs')

#PKS-NRP_Hybrids
PKSNRP_Hybrids_BGC = GCF_BGC[GCF_BGC['BiG-SCAPE_class'] == 'PKS-NRP_Hybrids'] #filter only NRPSs from table
PKSNRP_Hybrids_count = PKSNRP_Hybrids_BGC['GCFs'].value_counts().to_frame() #count
PKSNRP_Hybrids_count.reset_index(inplace=True) #set index as column
PKSNRP_Hybrids_count = PKSNRP_Hybrids_count.rename(columns={'GCFs':'PKS_NRP_Hybrids','index':'GCFs'}) #rename columns
PKSNRP_Hybrids_count = PKSNRP_Hybrids_count[PKSNRP_Hybrids_count['GCFs'] != 'singleton'] #remove singletons
GCFs_list = pd.merge(GCFs_list,PKSNRP_Hybrids_count, how='left', on='GCFs')

#PKSI
PKSI_BGC = GCF_BGC[GCF_BGC['BiG-SCAPE_class'] == 'PKSI'] #filter only NRPSs from table
PKSI_count = PKSI_BGC['GCFs'].value_counts().to_frame() #count
PKSI_count.reset_index(inplace=True) #set index as column
PKSI_count = PKSI_count.rename(columns={'GCFs':'PKSI','index':'GCFs'}) #rename columns
PKSI_count = PKSI_count[PKSI_count['GCFs'] != 'singleton'] #remove singletons
GCFs_list = pd.merge(GCFs_list,PKSI_count, how='left', on='GCFs')

#PKSother
PKSother_BGC = GCF_BGC[GCF_BGC['BiG-SCAPE_class'] == 'PKSother'] #filter only NRPSs from table
PKSother_count = PKSother_BGC['GCFs'].value_counts().to_frame() #count
PKSother_count.reset_index(inplace=True) #set index as column
PKSother_count = PKSother_count.rename(columns={'GCFs':'PKSother','index':'GCFs'}) #rename columns
PKSother_count = PKSother_count[PKSother_count['GCFs'] != 'singleton'] #remove singletons
GCFs_list = pd.merge(GCFs_list,PKSother_count, how='left', on='GCFs')

#RiPPs
RiPPs_BGC = GCF_BGC[GCF_BGC['BiG-SCAPE_class'] == 'RiPPs']
RiPPs_count = RiPPs_BGC['GCFs'].value_counts().to_frame()
RiPPs_count.reset_index(inplace=True) #set index as column
RiPPs_count = RiPPs_count.rename(columns={'GCFs':'RiPPs','index':'GCFs'}) #rename columns
RiPPs_count = RiPPs_count[RiPPs_count['GCFs'] != 'singleton']
GCFs_list = pd.merge(GCFs_list,RiPPs_count, how='left', on='GCFs')

#Terpene
Terpene_BGC = GCF_BGC[GCF_BGC['BiG-SCAPE_class'] == 'Terpene']
Terpene_count = Terpene_BGC['GCFs'].value_counts().to_frame()
Terpene_count.reset_index(inplace=True) #set index as column
Terpene_count = Terpene_count.rename(columns={'GCFs':'Terpene','index':'GCFs'}) #rename columns
Terpene_count = Terpene_count[Terpene_count['GCFs'] != 'singleton']
GCFs_list = pd.merge(GCFs_list,Terpene_count, how='left', on='GCFs')
GCFs_list = GCFs_list.rename(columns={'MIBiG_main_product':'Annotation'})

GCFs_list.to_csv(os.path.join(results_folder,'Network_metadata_GCFs_list-NPOmix1.0-%s.txt'%(current_date)))

#Create a network that connects MFs with the three GCFs
ntwks = mf_gcf[["MFs","predicted_GCFs","jaccard_scores"]]

#Loop to make a df with each GCF according to k_value
for i in range(1,k_value+1):
    exec(f'ntwks{i}=mf_gcf[["predicted_GCFs"]]')
    exec(f'ntwks{i}["predicted_GCFs"] = ntwks{i}["predicted_GCFs"].str.get({i-1})')

nt_GCF=pd.concat([ntwks1, ntwks2, ntwks3], axis=1)
nt_GCF.columns=["GCF"+str(i) for i in range(1,k_value+1)]

#Loop to make a df with each jaccard according to k_value
for i in range(1,k_value+1):
    exec(f'ntwks{i}=mf_gcf[["jaccard_scores"]]')
    exec(f'ntwks{i}["jaccard_scores"] = ntwks{i}["jaccard_scores"].str.get({i-1})')

#Concatenate jaccard and GCFs. Needs improvement cause only works with k=3
nt_jac=pd.concat([ntwks1, ntwks2, ntwks3], axis=1)
nt_jac.columns=["jaccard_scores"+str(i) for i in range(1,k_value+1)]

nt_GCFs = pd.concat([nt_GCF, nt_jac], axis=1)
ntwks = pd.concat([ntwks.MFs, nt_GCFs], axis=1)
ntwks = ntwks.set_index('MFs')

#melt to long format each GCF with its jaccard score
for i in range(1,k_value+1):
    exec(f'var{i}=ntwks[["GCF{i}","jaccard_scores{i}"]]')
    exec(f'var_wide{i}=pd.melt(var{i}, ignore_index=False, id_vars=["GCF{i}"], value_vars=["jaccard_scores{i}"])')
    exec(f'var_wide{i}.columns=["GCF","variable","jaccard_score"]')
    
ntwks_wide=pd.concat([var_wide1, var_wide2, var_wide3]) #concatenate again. Needs improvement too
ntwks_wide.reset_index(inplace=True)
ntwks_wide = ntwks_wide[["MFs","GCF","jaccard_score"]]
ntwks_wide = ntwks_wide[ntwks_wide['jaccard_score'] >= 0.7] #filter jaccard 
ntwks_wide = ntwks_wide.drop_duplicates() #remove duplicates
ntwks_wide.to_csv(os.path.join(results_folder,'MF_GCF_network_all-NPOmix1.0-%s.txt'%(current_date)))

#Filter MFs from training dataset 
ntwks_wideG2 = pd.merge(ntwks_wide, MFs_list2, how='left', on='MFs')
ntwks_wideG2 = ntwks_wideG2[ntwks_wideG2.DefaultGroups != "G1"]
ntwks_wideG2 = ntwks_wideG2.drop(columns=['count','Annotation','DefaultGroups','Type']) #remove columns
ntwks_wideG2.to_csv(os.path.join(results_folder,'MF_GCF_network_G2-NPOmix1.0-%s.txt'%(current_date)))
ntwks_wideG2

ax = subset_final_df['max_jaccard'].plot(x='metabolite_ID',figsize=(10,6))

ax.set_ylim(0.65, 1.01)
ax.set_xlabel('%s metabolites'%len(subset_final_df['max_jaccard']),size=20)
ax.set_ylabel('Maximum jaccard similarity for k = %s'%k_value,size=20)

ax

#save figure
plt.savefig(os.path.join(results_folder,'Jaccard_plot-NPOmix1.0-%s.png'%(current_date)))

boxplot = pd.DataFrame(subset_final_df['parent_mass']).boxplot(figsize=[16,8],rot=90,fontsize='xx-large')
#save figure
plt.savefig(os.path.join(results_folder,'parentmass_boxplot-NPOmix1.0-%s.png'%(current_date))) #save plot

boxplot = pd.DataFrame(subset_final_df['sum_intensity']).boxplot(figsize=[16,8],rot=90,fontsize='xx-large')
#save figure
plt.savefig(os.path.join(results_folder,'sum_intensity_boxplot-NPOmix1.0-%s.png'%(current_date))) #save plot

k_list = [1,3,5,10,25,50,100]
precision_scores = [100,80,90.1,85.7,60,52.17,60]
annotation_count = [3,10,11,14,20,23,24]
recall_list = []
for item in annotation_count:
    recall_list.append(round(item/max(annotation_count),2))
    
labels = ['k_{0}'.format(i) for i in k_list]

plt.figure(figsize=(10, 10), dpi=120)
plt.scatter(recall_list, precision_scores, s=k_list*np.repeat(25,len(k_list)))
for label, x, y in zip(labels, recall_list, precision_scores):
    plt.annotate(
        label,
        xytext=(-5,5),
        textcoords='offset points', ha='right', va='bottom',
        xy=(x, y))
#plt.show()

#save figure
plt.savefig(os.path.join(results_folder,'k_list-NPOmix1.0-%s.png'%(current_date))) #save plot


end = time.time()
hours, rem = divmod(end-start, 3600)
minutes, seconds = divmod(rem, 60)
run_time = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds)
print(run_time)
