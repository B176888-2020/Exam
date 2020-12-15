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
         "    blastPy: conduct the blast analysis in different combination of query and databases,.\n" \
         "SYNOPSIS\n" \
         "    blastPy -m [The file containing multiple combinations of query and database directories and arguments] [project directory]\n" \
         "    blastPy [type of blast] [directory of query file] [directory of database file] [makeblastdb or not] [project directory]\n" \
         "DESCRIPTION\n" \
         "    blastPy can call blastn, blastp, blastx, tblastn to conduct different types of BLAST analysis.\n" \
         "    The parsed results of blastPy can also be found and checked in corresponding directories.\n" \
         "\n" \
         "    Options:\n" \
         "\n" \
         "    -h, --help\n" \
         "        display this help manual and exit.\n" \
         "    -m [The file containing multiple combinations of query and database directories and arguments]\n" \
         "        use multiple mode to conduct the BLAST analysis process for multiple queries and databases.\n" \
         "\n" \
         "AUTHOR\n" \
         "    B176888-2020\n" \
         "REPORTINGS BUGS\n" \
         "    Please feel free to report bugs to <https://github.com/B176888-2020/Exam/issues>\n" \
         "SEE ALSO\n" \
         "    Full program and documentation at <https://github.com/B176888-2020/Exam>."


################################## Functions ##################################
def how_to_blast(file, dir, type_iq, type_db, dir_iq, dir_db):
    dir_xml = dir + file + ".xml"
    dir_out = dir + file + ".out"
    if (type_iq == "nucl") and (type_db == "nucl"):
        res_blast5 = bioblast.NcbiblastnCommandline(query=dir_iq, db=dir_db, outfmt=5, out=dir_xml, evalue=0.001, )
        res_blast7 = bioblast.NcbiblastnCommandline(query=dir_iq, db=dir_db, outfmt=7, out=dir_out, evalue=0.001, )
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
        if type_blast is None:
            print("\nInputError: parameter type_blast is None.")
            type_blast = input("\nPlease enter proper value for type_blast or enter EXIT.\n")
            if type_blast.upper() == "EXIT":
                exit()
            else:
                continue
        elif dir_inquiry is None:
            print("\nInputError: dir_inquiry is None.")
            dir_inquiry = input("\nPlease enter proper value(dir of query file)  for dir_inquiry or enter EXIT.\n")
            if dir_inquiry.upper() == "EXIT":
                exit()
            else:
                continue
        elif dir_database is None:
            print("\nInputError: dir_database is None.")
            dir_database = input("\nPlease enter proper value(dir of database file) for dir_database or enter EXIT.\n")
            if dir_database == "EXIT":
                exit()
            else:
                continue
        elif mk_db is None:
            print("\nInputError: parameter mk_db is None.")
            mk_db = input("\nPlease enter proper value(True or False) for mk_db or enter EXIT.\n")
            if mk_db.upper() == "EXIT":
                exit()
            else:
                continue
        elif projectSpace is None:
            print("\nInputError: projectSpace is None")
            projectSpace = input("\nPlease enter proper value(dir of your project) for projectSpace or enter EXIT.\n")
            if projectSpace.upper() == "EXxIT":
                exit()
            else:
                continue
        else:
            break

    # Check if the file/directories are existing
    while True:
        if not os.path.exists(dir_inquiry):
            print("\nInputError: The directory of the query sequence is not found.")
            dir_inquiry = input("\nPlease enter a qualified directory for parameter dir_inquiry or enter EXIT.\n")
            if dir_inquiry.upper() == "EXIT":
                exit()
            else:
                continue
        elif (not os.path.exists(dir_database)) and (not os.path.exists("/".join(dir_database.split("/")[:-1]))):
            print("\nInputError: The directory of the database is not found.")
            dir_database = input("\nPlease enter a qualified directory for parameter dir_database or enter EXIT.\n")
            if dir_database.upper() == "EXIT":
                exit()
            else:
                continue
        elif not os.path.exists(projectSpace):
            print("\nInputError: The directory of the database is not found.")
            projectSpace = input("\nPlease enter a qualified directory for parameter projectSpace or enter EXIT.\n")
            if projectSpace.upper() == "EXIT":
                exit()
            else:
                continue
        else:
            break

    # Check if the file contents are empty
    if os.stat(dir_inquiry).st_size == 0:
        print("\nInputError: The query sequence file is empty. Please check if the query sequence file is correct.")
        exit()
    elif os.path.isfile(dir_database):
        if os.stat(dir_database).st_size == 0:
            print("\nInputError: The database file is empty. Please check if the database file is correct.")
            exit()
    elif len(os.listdir("/".join(dir_database.split("/")[:-1]))) == 0:
        print("\nInputError: The database directory is empty. Please check if the database file is correct.")
        exit()

    # Check if the mk_db is boolean type
    while True:
        if (mk_db != "True") and (mk_db != "False"):
            print("\nInputError: The mk_db should be a boolean string.")
            mk_db = input("\nPlease reassign True or False to parameter mk_db.\n")
            continue
        else:
            if mk_db == "True":
                mk_db = True
            elif mk_db == "False":
                mk_db = False
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
            type_database = input("\nPlease enter prot or nucl as the corrected dbtype or enter EXIT.\n")
            if type_database.upper() == "EXIT":
                exit()
            else:
                continue
        elif (type_inquiry != "prot") and (type_inquiry != "nucl"):
            print("\nWarning: The type of your inquiry sequence is wrong. " +
                  "You can correct it as follow or exit by press EXIT")
            type_inquiry = input(
                "\nPlease enter prot or nucl as the corrected type of inquiry sequence or enter EXIT.\n")
            if type_inquiry.upper() == "EXIT":
                exit()
            else:
                continue
        else:
            break

    dirPro = projectSpace + inq_name + "_" + ref_name + "/"
    os.makedirs(dirPro, exist_ok=True)

    # Check if the data is too much for further analysis
    Seq = open(dir_inquiry, "r")
    SeqContent = Seq.read()
    Seq.close()
    SeqCount = SeqContent.count(">")
    if (SeqCount > 1000):
        print("\nWarning: The number of query sequences is larger than 1000")
        while True:
            dec = input("\nDo you want to continue the BLAST analysis? Please enter YES or NO to make a decision.\n")
            if dec.upper() == "NO":
                exit()
            elif dec.upper() == "YES":
                break
            else:
                print("ResponseError: Sorry, please use YES or NO to make the decision.\n")
                continue

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
          " will be stored in " + str(dirPro) + "sum_data")

    print("\nFinish: " + type_blast + " BLAST analysis with query: " + dir_inquiry + " and database:" + dir_database)
    print("\nThe outputs will be stored in " + str(dirPro) + "sum_data")


################################## Main program ##################################
# Add edirect to python environemnt PATH in server
# os.environ["PATH"] += os.pathsep + "/localdisk/data/BPSM/Assignment2/"
# os.environ["PATH"] += os.pathsep + "/localdisk/home/$USER/edirect/"

# Some default values for the input arguments
projectSpace = "./"

# Get proper inputs from users
dict_inputs = {}
if len(sys.argv) == 1:
    reply = input("\nInputDefection: \n" +
                  "Please provide proper inputs of the query and database files and arguments based on the user manual or --help\n"
                  "Then you can conduct BLAST analysis of one(normal mode) or more(multiple mode activated by option -m) pairs of different query sequences and databases\n"
                  "More detailed information about proTree usage can be checked by enter -h or --help. \n" +
                  "If you want to exit, please enter EXIT or press the Ctrl+C to exit the programme. \n")
    if reply == "EXIT" or reply is None:
        exit()
    elif reply == "-h" or reply == "--help":
        print(manual)
    exit()
elif len(sys.argv) > 1:
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(manual)
        exit()
    elif sys.argv[1] == "-m":
        print("\n################################## Inputs ##################################")
        if len(sys.argv) < 3:
            print(
                "InputError: Input defects. Please check if there are any arguments are not provided as the --help manual")
            exit()
        arg_file = sys.argv[2]
        if len(sys.argv) > 3:
            projectSpace = sys.argv[3]
    else:
        print("\n################################## Inputs ##################################")
        if len(sys.argv) < 5:
            print(
                "InputError: Input defects. Please check if there are any arguments are not provided as the --help manual")
            exit()
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
