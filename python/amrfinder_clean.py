#!/bin/python

import argparse
import os 
import re
import numpy
import time

def get_options():
        ## Just take in the list of fasta files and the length of the mfas to be output and a prefix for the files 
    

    purpose = ''' This is a scipt to take in a load of amrfinder tsvs and then output a single tsv with all consolidated resistances
    Usage:
    python amrfinder_clean.py --output <full_path_to_dir> --results [tsv files]'''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='fasta_list_creator.py')

    parser.add_argument('--results', required=True, help='amrfinder tsv files (required)', type=str, nargs="+")
    parser.add_argument('--output', required=True, help='prefix for outfile', type=str)
    
    args = parser.parse_args()

    return args


def main():
    start = time.perf_counter()
    inputs = get_options()

    fastas = inputs.results
    ## Create list of unique subclasses 
    ## open up the fasta list

    subclasses = []

    for res in fastas:
        with open(res, 'r') as file:
            for line in file:
                subclass = line.split("\t")[12] # This is the column index of the subclass 
                if subclass != "Subclass":
                    subclass2 = subclass + "-(" + line.split("\t")[11] + ")"
                    if subclass2 not in subclasses:
                        subclasses.append(subclass2)

    
    ## Now we have our list of subclasses 
    print("These are the %s subclasses present in the collection: " % len(subclasses))
    print("\n" + '\n'.join(subclasses) + "\n")

    ## Now we need to create the outdf 
    ## Have:
    ## id, abx1...abxn, num_res
    ## Lets do this line by line and append to the outfile 

    ## Set up the outfile
    outfile = inputs.output + ".tsv"

    if os.path.exists(outfile):
        print("Deleting already existing outfile %s" % outfile)
        os.remove(outfile)
    ## Set up the header line and a default line to be altered
    colnames = ["id"] + subclasses + ["num_resistant"]
    colline = '\t'.join(colnames) + "\n"
    test_line =  ["-"] + numpy.repeat(["No"], len(subclasses)).tolist() + [0]
    
    ## Going to append each line to the outfile by looping through the 
    ## files input and matching the resistance profiles to the header
    with open(outfile, 'a') as outer:
        ## Write the header line out
        outer.write(colline)
        for res in fastas:
            ## loop through the results tsvs now 
            if "GCF" in res:
                current_name = re.sub("_amrfinderplus.tsv","",res)
            else:
                current_name = re.sub("_amrfinderplus.tsv","_ukhsa",res)
            with open(res, 'r') as amr:
                ## Get the number of resistant classes of the file, none will just contain the header 
                classes = []
                for line in amr:
                    subclass = line.split("\t")[12] # This is the column index of the subclass 
                      
                    if subclass != "Subclass":
                        subclass2 = subclass + "-(" + line.split("\t")[11] + ")"
                        if subclass2 not in classes:
                            classes.append(subclass2)
                indexes = [colnames.index(x) for x in classes]
                current_line = test_line.copy()
                current_line[0] = current_name
                if len(indexes) > 0:
                    # Would be nice if there was a more elegant way to alter multiple indices at the same time,
                    # but can't find one at the moment. 
                    for x in indexes:
                        current_line[x] = "Yes"
                    current_line[len(current_line) - 1] = len(indexes)
                
                newline = '\t'.join(list(map(str,current_line))) + '\n'
                outer.write(newline)

    end = time.perf_counter()
    time_took = time.strftime('%H:%M:%S', time.gmtime((end - start)))
    print("Done, finished this in: %s" % time_took)

if __name__ == '__main__':
    main()

    




    