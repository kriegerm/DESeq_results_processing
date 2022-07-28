#!/usr/bin/env python3

'''
matchlocustags_withresults.py
UPDATED: 8/27/20

This specialized script will take in a shorter list of locus tags as a csv. It also requires a second input with locus tags in the first column and the rest of the 
data in the subsequent columns (I used this originally for DEseq results). It will go through that long data file and pull out on the results from your
shorter locus tag input. 

input file:
header row
locus_tag in the first column
All other data in the other columns
'''


import csv
import sys

def annotate (datacsv, inputcsv, outputcsv):

    datadict = {}

    with open(datacsv, 'r') as fh:
        fhcsv = csv.reader(fh, delimiter=',')
        field_names_list = next(fhcsv)
        count = 0

        for entry in fhcsv:
            count += 1
            locus = entry[0]
            locus = locus.rstrip()
            locus = locus.lstrip()
            locus = locus.replace(" ", "")

            datadict[str(locus)] = entry

    print(len(datadict))

    print("There are %d entries counted in the input data csv" %count)

    with open(inputcsv, 'r', encoding='utf-8-sig') as fh:
        fhcsv = csv.reader(fh, delimiter=',')

        for entry in fhcsv:
            entrycounter = 0

            entry = str(entry)
            entry = entry.rstrip("']")
            entry = entry.lstrip("['")
            entry = entry.replace(" ", "")

            for locustag, description in datadict.items():
                locustag = str(locustag)
                locustag = locustag.rstrip("']")
                locustag = locustag.lstrip("['")
                locustag = locustag.replace(" ", "")

                results = []
                writeline = []

                if locustag == entry:
                    entrycounter += 1
                    results.append(description)
                    writeline.append(locustag)
                    writeline.append(results)

                    writeline = str(writeline)
                    writeline = writeline.replace("'", "")
                    writeline = writeline.replace('"', "")
                    writeline = writeline.replace("]", "")
                    writeline = writeline.replace("[", "")
                    writeline = writeline.replace("/", "")
                    writeline = writeline.replace("\\", "")

                    with open(outputcsv, 'a') as output:
                        writer = csv.writer(output)
                        writer.writerow([writeline])


            if entrycounter == 0:
                with open(outputcsv, 'a') as output:
                    writer = csv.writer(output)
                    writer.writerow([writeline])

if __name__ == '__main__':
    if len(sys.argv) == 4:
        annotate(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("Useage: input csv with locus tags in first column and header row, input csv with list of locus tags, output csv name")
