#!/usr/bin/env python

import sys
import os
import getopt
import pickle

# Variables that contains the inputs files 
OPT_INPUT_FILE=""

def Usage():
  print "\nreadCodonTableAsPickleStructure.py is a program that read an pickle file and print its content in the stdout. Only for internal use.\n"
  print "Usage:"
  print "\treadCodonTableAsPickleStructure.py -i [PICKLE FILE]\n"
  print "\nMandatory options:"
  print "\t-i, --input=FILE"
  print "\t\tPickle file (*.pkl) that its load and its content is printed in the stdout."
  print "\t-h, --help"
  print "\t\tShow the options of the program."
  print "\n"
  sys.exit(1)

# Function that read and parse the command line arguments.
def SetOptions(argv):
  if len(argv) == 0:
    Usage()
  options, remaining = getopt.getopt(argv, 'i:h', ['input=','help'])
  opt_flag = {'i': False}
  global OPT_INPUT_FILE
  for opt, argu in options:
    if opt in ('-i', '--input'):
      if not opt_flag['i']:
        if os.path.exists(argu):
          OPT_INPUT_FILE = argu
          opt_flag['i'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: File or path of the input file does not exist. ", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the input file. Option -i / --input was already set.\n"
        sys.exit(1)
    elif opt in ('-h', '--help'):
      Usage()
  if not opt_flag['i']:
    print >> sys.stderr , "[ERROR]: Input FASTA ffle is not defined. Option -i / --input.\n"
    sys.exit(1)

###########################################
################   MAIN   #################
###########################################

# Parse command line
SetOptions(sys.argv[1:])

# Load the pickle object.
fh_pkl_codon_table = open(OPT_INPUT_FILE, 'rb')
codon_table = pickle.load(fh_pkl_codon_table)
fh_pkl_codon_table.close()

print codon_table

