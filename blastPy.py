#! /home/ninomoriaty/anaconda3/envs/PSBenv/bin/python

################################## Modules and packages ##################################
# Packages and modules
## OS and IO
import os
import shutil
import subprocess
import sys

## Regular Expression
import re

## String
import string

## Dictionary
from collections import OrderedDict

## Scientific Calculation
import numpy as np
import matplotlib.pyplot as plt
import pandas
import scipy
import Bio

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

################################## Main program ##################################
# Add edirect to python environemnt PATH in server
# os.environ["PATH"] += os.pathsep + "/localdisk/data/BPSM/Assignment2/"
# os.environ["PATH"] += os.pathsep + "/localdisk/home/$USER/edirect/"

# Get proper inputs from users
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
elif len(sys.argv) > 1:
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(manual)
        exit()
    else:
        type_blast = sys.argv[1]
    # TODO: When the dataabse is not provided in the path
        dir_database = sys.argv[2]
        dir_inquiry = sys.argv[3]
    # TODO: if the database is existing
        mk_db = sys.argv[4]
    # TODO: The parameters about the blast



print("\n################################## Data Collection and Selection ##################################")
# Check if the shortcut of the data type is correct
filename = re.split("/", dir_database)[-1]
list_types = re.split("-", type_blast)
type_database = list_types[0]
type_inquiry = list_types[1]
while True:
    if (type_database != "prot") or (type_database != "nucl"):
        print("\nWarning: The dbtype of your database is wrong. " +
              "You can correct it as follow or exit by press EXIT")
        type_database = input("\nPlease enter prot or nucl as the corrected dbtype or enter EXIT.")
    if (type_inquiry != "prot") or (type_inquiry != "nucl"):
        print("\nWarning: The type of your inquiry sequence is wrong. " +
              "You can correct it as follow or exit by press EXIT")
        type_inquiry = input("\nPlease enter prot or nucl as the corrected type of inquiry sequence or enter EXIT.")
    if type_database == "EXIT" or type_inquiry == "EXIT":
        exit()
    else:
        break


# Get the primary protein sequence by `esearch` and `edirect`
print("\n**** Collecting the data by esearch and efetch... ****")
os.makedirs("./data", exist_ok=True)
os.makedirs("./sum_data", exist_ok=True)

# Make BLAST database
if mk_db:
    print("\n**** Making BLAST database... ****")
    os.system("makeblastdb -in " + dir_database + " --dbtype " + type_database + " -out ./data/" + filename)
else:
    print("\n**** Skip the makeblastdb process and use existing BLAST database... ****")

# Do the BLAST analysis
































