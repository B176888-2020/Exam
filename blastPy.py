#! /home/ninomoriaty/anaconda3/envs/PSBenv/bin/python
# Remember to change to /usr/bin/python3

################################## Modules and packages ##################################
# Packages and modules
## OS and IO
import os
import sys

## Regular Expression
import re

## String
import string

## Scientific Calculation
import pandas as pd
import Bio.Blast.Applications as bioblast
from Bio.Blast import NCBIXML

################################## Manual ##################################
manual = "\nNAME\n" \
         "    proTree: analyse the conservation level of protein family within taxonomy group\n" \
         "SYNOPSIS\n" \
         "    proTree [-options] [protein_families] [taxonomy_groups] [project_directory] [selection_modes] [selection]\n" \
         "DESCRIPTION\n" \
         "    proTree can analyse and visualise the conservation level of the protein family within taxonomy group.\n" \
         "    Functions in EMBOSS are also provdied to broaden the scale of motifs search and other usages.\n" \
         "\n" \
         "    Options:\n" \
         "\n" \
         "    -h, --help\n" \
         "        display this help and exit\n" \
         "    -v\n" \
         "        use interactive mode to conduct the analysis process. Users can provide inputs acccording to the requirements.\n" \
         "    -s\n" \
         "        use silent mode to conduct the analysis process if users can provide correct inputs.\n" \
         "\n" \
         "AUTHOR\n" \
         "    B176888-2020\n" \
         "REPORTINGS BUGS\n" \
         "    Please feel free to Report bugs to <https://github.com/B176888-2020/Exam/issues>\n" \
         "SEE ALSO\n" \
         "    Full program and documentation at <https://github.com/B176888-2020/Exam>."


################################## Functions ##################################
def how_to_blast(file, dir, type_iq, type_db, dir_iq, dir_db):
    dir_xml = dir + file + ".xml"
    dir_out = dir + file + ".out"
    if (type_iq == "nucl") and (type_db == "nucl") :
        res_blast5 = bioblast.NcbiblastnCommandline(query=dir_iq, db=dir_db, outfmt=5, out=dir_xml, evalue=0.001,)
        res_blast7 = bioblast.NcbiblastnCommandline(query=dir_iq, db=dir_db, outfmt=7, out=dir_out, evalue=0.001,)
    elif (type_iq == "prot") and (type_db == "nucl"):
        res_blast5 = bioblast.NcbitblastnCommandline(query=dir_iq, db=dir_db, outfmt=5, out=dir_xml, evalue=0.001)
        res_blast7 = bioblast.NcbitblastnCommandline(query=dir_iq, db=dir_db, outfmt=7, out=dir_out, evalue=0.001)
    elif (type_iq == "nucl") and (type_db == "prot"):
        res_blast5 = bioblast.NcbiblastxCommandline(query=dir_iq, db=dir_db, outfmt=5, out=dir_xml, evalue=0.001)
        res_blast7 = bioblast.NcbiblastxCommandline(query=dir_iq, db=dir_db, outfmt=7, out=dir_out, evalue=0.001)
    elif (type_iq == "prot") and (type_db == "prot"):
        res_blast5 = bioblast.NcbiblastpCommandline(query=dir_iq, db=dir_db, outfmt=5, out=dir_xml, evalue=0.001)
        res_blast7 = bioblast.NcbiblastpCommandline(query=dir_iq, db=dir_db, outfmt=7, out=dir_out, evalue=0.001)
    print(res_blast5)
    print(res_blast7)
    res_blast5()
    res_blast7()


def parse_blast(file, dir):
    E_VALUE_THRESH = 0.04
    dir_xml = dir + file + ".xml"
    res_xml = open(dir_xml, "r")
    blast_records = NCBIXML.parse(res_xml)
    dir_records = dir + file + ".txt"
    res_txt = open(dir_records, "w")
    for blast_record in blast_records:
        res_txt.write("\n################################## Query ##################################")
        res_txt.write("\nquery_sequence: " + str(blast_record.query))
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    res_txt.write("\n**** Hits ****")
                    res_txt.write("\nref_sequence: " + str(alignment.title))
                    res_txt.write("\nref_length: " + str(alignment.length))
                    res_txt.write("\ne value: " + str(hsp.expect))
                    res_txt.write("\n" + hsp.query[0:75] + "...")
                    res_txt.write("\n" + hsp.match[0:75] + "...")
                    res_txt.write("\n" + hsp.sbjct[0:75] + "...")
    res_txt.close()
    res_xml.close()


def main(type_blast, dir_inquiry, dir_database, mk_db, projectSpace):
    print("\n################################## Validate Inputs ##################################")
    # Check if there are any None objects
    while True:
        if (type_blast is None):
            print("\nInput Error: parameter type_blast is None.")
            type_blast = input("\nPlease enter proper value for type_blast or enter EXIT.")
            if (type_blast == "EXIT"):
                exit()
            else:
                continue
        elif (dir_inquiry is None):
            print("\nInput Error: dir_inquiry is None.")
            if (dir_inquiry == "EXIT"):
                exit()
            else:
                continue
        elif (dir_database is None):
            print("\nInput Error: dir_database is None.")
            if (dir_database == "EXIT"):
                exit()
            else:
                continue
        elif (mk_db is None):
            print("\nInput Error: parameter mk_db is None.")
            if (mk_db == "EXIT"):
                exit()
            else:
                continue
        elif (projectSpace is None):
            print("\nInput Error: projectSpace is None")
            if (projectSpace == "EXIT"):
                exit()
            else:
                continue
        else:
            break

    # Check if the shortcut of the data type is correct
    inq_name = re.split("/", dir_inquiry)[-1].split(".")[0]
    ref_name = re.split("/", dir_database)[-1].split(".")[0]
    list_types = re.split("-", type_blast)
    type_inquiry = list_types[0]
    type_database = list_types[1]
    while True:
        if (type_database != "prot") and (type_database != "nucl"):
            print("\nWarning: The dbtype of your database is wrong. " +
                  "You can correct it as follow or exit by press EXIT")
            type_database = input("\nPlease enter prot or nucl as the corrected dbtype or enter EXIT.")
            if (type_database == "EXIT"):
                exit()
            else:
                continue
        elif (type_inquiry != "prot") and (type_inquiry != "nucl"):
            print("\nWarning: The type of your inquiry sequence is wrong. " +
                  "You can correct it as follow or exit by press EXIT")
            type_inquiry = input("\nPlease enter prot or nucl as the corrected type of inquiry sequence or enter EXIT.")
            if (type_inquiry == "EXIT"):
                exit()
            else:
                continue
        else:
            break

    dirPro = projectSpace + inq_name + "_" + ref_name + "/"
    os.makedirs(dirPro, exist_ok=True)

    print("\n################################## Data Collection and Selection ##################################")
    # Get the directory of the data and create the folders to store the intermediate files as well as results.
    print("\n**** Collecting the data by esearch and efetch... ****")
    dir_data = dirPro + "data/"
    dir_results = dirPro + "sum_data/"
    os.makedirs(dir_data, exist_ok=True)
    os.makedirs(dir_results, exist_ok=True)

    # Make BLAST database
    if mk_db:
        print("\n**** Making BLAST database... ****")
        os.system("makeblastdb -in " + dir_database + " -dbtype " + type_database + " -out " + dir_data + ref_name)
        dir_database = dir_data + ref_name
    else:
        print("\n**** Skip the makeblastdb process and use existing BLAST database... ****")
        dir_database = dir_database

    print("\n################################## BLAST Analysis ##################################")
    # Do the BLAST analysis
    how_to_blast(inq_name, dir_results, type_inquiry, type_database, dir_inquiry, dir_database)
    parse_blast(inq_name, dir_results)
    output_xml = dir_results + inq_name + ".xml"
    output_out = dir_results + inq_name + ".out"
    print("\nDone. The output files:" + output_xml + " and " + output_out +
          " will be stored in the sum_data directory in " + str(dirPro))

    print("\nFinish: " + type_blast + " BLAST analysis with query: " + dir_inquiry + " and database:" + dir_database)
    print("\nThe outputs will be stored in the sum_data directory in " + str(dirPro))


################################## Main program ##################################
# Add edirect to python environemnt PATH in server
# os.environ["PATH"] += os.pathsep + "/localdisk/data/BPSM/Assignment2/"
# os.environ["PATH"] += os.pathsep + "/localdisk/home/$USER/edirect/"

# Some default values for the input arguments
projectSpace = "./"

# Get proper inputs from users
print("\n################################## Inputs ##################################")
dict_inputs = {}
if len(sys.argv) == 1:
    reply = input("\nInteractive Mode(Default): \n" +
                  "Please input the file directory to the "
                  "More detailed information about proTree usage can be checked by enter -h or --help. \n" +
                  "If you want to exit, please enter EXIT or press the Ctrl+C to exit the programme. \n")
    if reply == "EXIT":
        exit()
    elif reply == "-h" or reply == "--help":
        print(manual)
        exit()
    # TODO: Enter the information interactively
elif len(sys.argv) > 1:
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(manual)
        exit()
    elif sys.argv[1] == "-m":
        arg_file = sys.argv[2]
        if len(sys.argv) > 3:
            projectSpace = sys.argv[3]
    else:
        dict_inputs["type_blast"] = [sys.argv[1]]
        dict_inputs["dir_inquiry"] = [sys.argv[2]]
        dict_inputs["dir_database"] = [sys.argv[3]]
        dict_inputs["mk_db"] = [sys.argv[4]]
        # TODO: The parameters about the blast
        if len(sys.argv) > 5:
            projectSpace = sys.argv[5]

if sys.argv[1] == "-m":
    if len(sys.argv) > 3:
        print("InputError: When using the file to submit multiple inputs, " +
              "there should be only an argument provding the directory of the file containing the arguments list")
    else:
        names = ["type_blast", "dir_inquiry", "dir_database", "mk_db"]
        args_content = pd.read_csv(arg_file, names=names,
                                   delim_whitespace=True)
        for name in names:
            dict_inputs[name] = args_content[name]

counter = 0
for counter in list(range(len(dict_inputs["type_blast"]))):
    main(dict_inputs["type_blast"][counter], dict_inputs["dir_inquiry"][counter],
         dict_inputs["dir_database"][counter], dict_inputs["mk_db"][counter],
         projectSpace)