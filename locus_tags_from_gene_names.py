#!/usr/bin/env python3

'''
locus_tags_from_gene_names.py - USE THIS TO CONVERT A LIST OF GENE NAMES/LOCUS TAGS TO ALL LOCUS TAGS
UPDATED: 8/27/20

This code will take in an input csv file of genome information (see below) and a input csv file with gene names
or mixed locus tags and gene names and output a csv file with the locus tags. If there is no locus tag avaliable
then the gene name will be provided.

    #input csv file structure by column (get this by exporting a .gff file, converting to csv, and then arranging the columns by hand)
    #1 - feature type
    #2 - start
    #3 - end
    #4 - strand
    #5 - gene
    #6 - locus tag
    #7 - product
    #8 - protein ID

'''

import csv
import sys

def annotate (inputcsv, locationscsv, outputcsv):
    gene_locus_dict = {}
    with open(locationscsv, 'r') as fh:
        fhcsv = csv.reader(fh, delimiter=',')
        field_names_list = next(fhcsv)

        totalentrycounter = 0
        genecounter = 0
        CDSandothercounter = 0

        uniquelocuslist = []

        for entry in fhcsv:
            gene = entry[4]
            gene = gene.lstrip()
            gene = gene.rstrip()

            locus = entry[5]
            locus = locus.rstrip()
            locus = locus.lstrip()

            uniquelocuslist.append(locus)

            totalentrycounter += 1

            if entry[0] == "Gene":
                locustag = []
                gene_entry = []
                locustag.append(locus)
                gene_entry.append(gene)
                gene_locus_dict[str(locustag)] = gene_entry

                genecounter += 1

            if entry[0] != "Gene":  # Adding back in all the gene entries that did not have gene records
                locustag = []
                gene_entry = []
                locustag.append(locus)
                gene_entry.append(gene)
                if str(locustag) not in gene_locus_dict.keys():
                    gene_locus_dict[str(locustag)] = gene_entry

                CDSandothercounter += 1


#Fills in the empty gene keys in locusgenedict, so that all values have entries
        for tag, gene in gene_locus_dict.items():
            if gene == ['']:
                gene_locus_dict[tag] = tag

    print(gene_locus_dict)


    ##KEEP ME IF YOU WANT A LOCUS TAG CONVERSION CSV FOR THE ENTIRE GENOME!
    ##Writes an output file of locusgenedict that has the locus tags in the first column and the gene name (or locus tag) in the second
    #        with open(outputlocustagcsv, 'w') as fh:
    #            writer = csv.writer(fh)
    #            for key, value in locusgenedict.items():
    #                writer.writerow([key, value])

    with open(inputcsv, 'r', encoding='utf-8-sig') as fh:
        fhcsv = csv.reader(fh, delimiter=',')

        inputcounter = 0
        outputcounter = 0

        for entry in fhcsv:

            inputcounter += 1
            entrycounter = 0

            entry = str(entry)
            entry = entry.rstrip("']")
            entry = entry.lstrip("['")
            entry = entry.replace(" ", "")

            for locustag, genename in gene_locus_dict.items():
                genename = str(genename)
                genename = genename.rstrip("']")
                genename = genename.lstrip("['")
                genename = genename.replace(" ", "")

                locustag = str(locustag)
                locustag = locustag.rstrip("']")
                locustag = locustag.lstrip("['")
                locustag = locustag.replace(" ", "")

                if genename == entry:
                    outputcounter += 1
                    entrycounter += 1
                    with open(outputcsv, 'a') as output:
                        writer = csv.writer(output)
                        writer.writerow([locustag])

            if entrycounter == 0:
                outputcounter += 1
                with open(outputcsv, 'a') as output:
                    writer = csv.writer(output)
                    writer.writerow([genename])

    if inputcounter != outputcounter:
        print("DIFFERENT NUMBER OF INPUT AND OUTPUT RESULTS!!! CHECK RESULTS CAREFULLY.")


if __name__ == '__main__':
    if len(sys.argv) == 4:
        annotate(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Useage: input list of gene names csv, input csv of genome information, output csv name")
