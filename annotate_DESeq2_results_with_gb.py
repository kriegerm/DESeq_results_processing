#!/usr/bin/env python3

#annotate_DESeq2_results.py
#UDPATED: 12/30/22
#It's very annoying how DESEq2 usually gives a mixture of locus tags and gene names, and it's even more annoying when you're trying to use the OLD locus tags for something (like kEGG analysis!)
#So this code takes in a csv of DEGs with the "Gene	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj" header. If you have the counts or anyhting extra on there, it'll be lost.
#It'll take in a NewGenbank file with the fancy new identifiers that they have for species like S. mutans, and then give you a csv
#With all the important info from the DEG results as well as  all the informatio you could ever want from the genbank file!

#NOTES ABOUT THIS CODE OUTPUT:
#####-Column A "Gene and column B "Gene_name_check" should be identical. This is a sanity check.
#####-If there is a "*" next to the gene name in column B, that means that there were TWO or more entries for that gene name in the genbank file! If you used a gb to map the reads,
######which you did, then sometimes the mapper uses the gene names (like "gyrA") instead of the unique identifiers (like "SMU_001"). There can be multiple copies of "gyrA" in a genome,
######so if you are looking at enrichmnet in operons this might mess things up because it's unclear which "gyrA" the reads are mapping to. So just be careful with these!
######It's also good to remember that the protein ID and product of these genes may be incorrect, because it's not clear if the reads mapped to this particular gene copy or not.

import sys
from Bio import SeqIO
import csv
import pandas as pd
import numpy as np

def annotate_DESeq2_results(DEGS_csv, NewGenbank):

    old_locus_tags_list = []
    log2FC_list = []
    pval_list = []
    padj_list = []

    with open(DEGS_csv, 'r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)  # skip the headers
        for line in reader:
            name = line[0]
            log2FC = line[2]
            pval = line[5]
            padj = line [6]

            old_locus_tags_list.append(name.strip())
            log2FC_list.append(log2FC.strip())
            pval_list.append(pval.strip())
            padj_list.append(padj.strip())

#Now we have a list of locgus tags, Log2FC, pvals, and adjusted p-vals. Lets just make sure they're all the same length and things can run smoothly.
    if len(old_locus_tags_list) != len(log2FC_list) != len(pval_list) != len(padj_list):
        print("WARNING: You need to examine your input csv to make sure all the columns have the same number of entries!")
        exit()

#Build a dictionary from the gb gb csvfile
    Genbank_dict = {}

    with open(NewGenbank) as input_handle:
        for gb_record in SeqIO.parse(input_handle, "genbank"):

            #Look through each record in the gb
            for f in gb_record.features:
                if f.type == "CDS": #Want to start with the CDS records because those have the most information

                    locus_tag = str(f.qualifiers.get('locus_tag')) #This is the NEW locus tag
                    locus_tag = locus_tag.lstrip("['").rstrip("']'")

                    gene_info = dict(f.qualifiers)
                    Genbank_dict[locus_tag] = gene_info #Genbank_dict is a dictionary of dictionaries!

                if f.type != "CDS": #Want to put things in there that don't have CDS entries
                    locus_tag = str(f.qualifiers.get('locus_tag')) #This is the NEW locus tag
                    locus_tag = locus_tag.lstrip("['").rstrip("']'")

                    if locus_tag not in Genbank_dict.keys():
                            gene_info = dict(f.qualifiers)
                            Genbank_dict[locus_tag] = gene_info


#I noticed that some of the genes might have duplicate entries for their gene names, and I just want to make note of that.
#For example, if there are two "galE" genes, then my code will only annotate the DESeq2 results with the first one.
#If they have the same gene name they are probably the same gene, but if I go back and add coordinates or something they might be incorrect!
#Ultimately this will depend on the mapper - if there are two galE genes in the gb file, who knows which one the reads were mapping to.
    gene_list = []
    for entry in Genbank_dict:
        entry_dict = Genbank_dict[entry]
        if str("gene") in entry_dict.keys():
            gene = entry_dict["gene"]
            gene_list.append(gene)

    uniqueList = []
    duplicateList = []

    for i in gene_list:
        if i not in uniqueList:
            uniqueList.append(i)
        elif i not in duplicateList:
            duplicateList.append(i)

    duplicateList2 = [] #This just makes the duplicate list more manageable, because right now it's a list of lists and that's not helpful
    if len(duplicateList) != 0:
        for thing in duplicateList:
            duplicateList2.append(thing[0])

#Optional warning message if you really care about looking through your genbank file
    #    print("WARNING: Your genbank file contains multiple annotations for the following genes. Make sure you examine your output closely to ensure that you are actually looking at the gene you think you are. Next time, it's better to use unique gene names to map your reads!")
    #    for item in duplicateList:
    #        print(item[0])


#Now we are ready to look through the list of old locus tags and extract more info about them from the Genbank dictionary.
    results_list = []
    locus_tag_list = []
    original_locus_tag_list = []
    protein_id_list = []
    product_list = []
    GO_function_list = []
    GO_pathway_list = []

    for item in old_locus_tags_list: #For each item in the list of genes DESeq2 is reporting
        count = 0

        if item in Genbank_dict.keys(): #If there's an entry for that item in the genbank dictionary, which there should be
            count = count + 1
            entry_dict = Genbank_dict[item] #We are going inception here to access the dictionary within the dictionary!

            results_list.append(item)

            if str("locus_tag") in entry_dict.keys():
                locus_tag_fromdict = []
                locus_tag_fromdict = list(entry_dict["locus_tag"])
                locus_tag_list.append(locus_tag_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
            if str("locus_tag") not in entry_dict.keys():
                locus_tag_list.append("No new locus tag found for " + item)

            if str("old_locus_tag") in entry_dict.keys():
                old_locus_tag_fromdict = []
                old_locus_tag_fromdict = list(entry_dict["old_locus_tag"])
                if len(old_locus_tag_fromdict) == 2:
                    original_locus_tag_list.append(old_locus_tag_fromdict[1]) #The SMU_ is always second to SMU. in the list, and this one is better. So We are going ot use that if it's an option!
                if len(old_locus_tag_fromdict) == 1:
                    original_locus_tag_list.append(old_locus_tag_fromdict[0])
            if str("old_locus_tag") not in entry_dict.keys():
                original_locus_tag_list.append("No old locus tag found for " + item)

            if str("product") in entry_dict.keys():
                product_fromdict = []
                product_fromdict = list(entry_dict["product"])
                product_list.append(product_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
            if str("product") not in entry_dict.keys():
                product_list.append("No product found")

            if str("protein_id") in entry_dict.keys():
                protein_id_fromdict = []
                protein_id_fromdict = list(entry_dict["protein_id"])
                protein_id_list.append(protein_id_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
            if str("protein_id") not in entry_dict.keys():
                protein_id_list.append("No protein ID found")

            if str("GO_function") in entry_dict.keys():
                GO_function_fromdict = []
                GO_function_fromdict = list(entry_dict["GO_function"])
                GO_function_list.append(GO_function_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
            if str("GO_function") not in entry_dict.keys():
                GO_function_list.append("No GO function found")

            if str("GO_pathway") in entry_dict.keys():
                GO_pathway_fromdict = []
                GO_pathway_fromdict = list(entry_dict["GO_pathway"])
                GO_pathway_list.append(GO_pathway_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
            if str("GO_pathway") not in entry_dict.keys():
                GO_pathway_list.append("No GO pathway found")

#Things get more complicated if you have a gene name as input, like "gyrA". So we have to go into the dictionary of dictionary and find the locus_tag for the gene namne.
        if item not in Genbank_dict.keys():

            for entry in Genbank_dict:
                entry_dict = Genbank_dict[entry] #Looking through each dictionary entry within the dictionary

                if str("gene") in entry_dict.keys():
                    gene = entry_dict["gene"]
                    if gene[0] == item and count == 0: #if the gene name matches up with what's in the dictionary AND we haven't already looked at an entry with that gene name
                        count = count + 1

                        if item in duplicateList2:
                            item2 = str(item+"*")
                            results_list.append(item2)
                        if item not in duplicateList2:
                            results_list.append(item)

                        if str("locus_tag") in entry_dict.keys():
                            locus_tag_fromdict = []
                            locus_tag_fromdict = list(entry_dict["locus_tag"])
                            locus_tag_list.append(locus_tag_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
                        if str("locus_tag") not in entry_dict.keys():
                            locus_tag_list.append("No new locus tag found for " + item)

                        if str("old_locus_tag") in entry_dict.keys():
                            old_locus_tag_fromdict = []
                            old_locus_tag_fromdict = list(entry_dict["old_locus_tag"])
                            if len(old_locus_tag_fromdict) == 2:
                                original_locus_tag_list.append(old_locus_tag_fromdict[1]) #The SMU_ is always second to SMU. in the list, and this one is better. So We are going ot use that if it's an option!
                            if len(old_locus_tag_fromdict) == 1:
                                original_locus_tag_list.append(old_locus_tag_fromdict[0])
                        if str("old_locus_tag") not in entry_dict.keys():
                            original_locus_tag_list.append("No KEGG locus tag found for " + item)

                        if str("product") in entry_dict.keys():
                            product_fromdict = []
                            product_fromdict = list(entry_dict["product"])
                            product_list.append(product_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
                        if str("product") not in entry_dict.keys():
                            product_list.append("No product found")

                        if str("protein_id") in entry_dict.keys():
                            protein_id_fromdict = []
                            protein_id_fromdict = list(entry_dict["protein_id"])
                            protein_id_list.append(protein_id_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
                        if str("protein_id") not in entry_dict.keys():
                            protein_id_list.append("No protein ID found")

                        if str("GO_function") in entry_dict.keys():
                            GO_function_fromdict = []
                            GO_function_fromdict = list(entry_dict["GO_function"])
                            GO_function_list.append(GO_function_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
                        if str("GO_function") not in entry_dict.keys():
                            GO_function_list.append("No GO function found")

                        if str("GO_pathway") in entry_dict.keys():
                            GO_pathway_fromdict = []
                            GO_pathway_fromdict = list(entry_dict["GO_pathway"])
                            GO_pathway_list.append(GO_pathway_fromdict[0]) #Choosing the first entry in case there are multiple for whatever reason
                        if str("GO_pathway") not in entry_dict.keys():
                            GO_pathway_list.append("No GO pathway found")

        if count == 0: #If there's not an entry for that item in the genbank dictionary, which shouldn't happen but I guess it could
            results_list.append(item + " NOT FOUND IN GB")
            locus_tag_list.append("NOT IN GB FILE (probably a new gene)")
            original_locus_tag_list.append("NA")
            product_list.append("NA")
            protein_id_list.append("NA")
            GO_function_list.append("NA")
            GO_pathway_list.append("NA")

#We will just compare all the lists to make sure they are the same length and everything looks OK. If they aren't, you're goin to have to troubleshoot!
    if len(results_list) != len(old_locus_tags_list) != len(original_locus_tag_list)!= len(locus_tag_list) != len(log2FC_list)  != len(pval_list) != len(padj_list)  != len(protein_id_list) != len(product_list) != len(GO_function_list) != len(GO_pathway_list):
        print("WARNING: LENGTHS OF LISTS DO NOT MATCH!!!")
        exit()

#Time to write out the results! So exciting! I'm using pandas to merge the original file with the new annotation file
    df = pd.DataFrame(columns=('Gene', 'Gene_name_check', 'KEGG_locus_tag', 'New_locus_tag', 'Protein_ID', 'Product', 'GO_function', 'GO_pathway'))
    for i in range(len(locus_tag_list)):
        row = []
        row = [old_locus_tags_list[i], results_list[i], original_locus_tag_list[i], locus_tag_list[i], protein_id_list[i], product_list[i], GO_function_list[i], GO_pathway_list[i]]
        df.loc[len(df)] = row

    outputcsv = str(DEGS_csv.rstrip(".csv") + "_annotated.csv")
    original_file = pd.read_csv(DEGS_csv)
    merged = original_file.merge(df, on='Gene')    #This may throw an error if the gene names aren't identical...but in that case, you have a bigger problem!
    merged.to_csv(outputcsv, index=False)

#This little block just handles how many input variables you give it to make sure you have the correct number of arguments
if __name__ == '__main__':
    if len(sys.argv) == 3:
         annotate_DESeq2_results(sys.argv[1], sys.argv[2])
    else:
         print("Usage: DEGS csv, gb")
         sys.exit(0)

