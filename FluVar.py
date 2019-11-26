#!/usr/bin/python

import sys
import os
import getopt
import subprocess
import tempfile

# Variables that contains the inputs files 
OPT_INPUT_FILE=""
OPT_LABEL_FILE=None
OPT_INPUT_FORMAT="T"
OPT_INPUT_TYPE="N"
OPT_STRAIN=""
OPT_OUTPUT_SUFFIX=None
OPT_PYMOL_SESSION=False
OPT_2D_PLOT=False
ALLOWED_STRAINS=["H1N1","H3N2"]
OPT_USER_FASTA_FILE=None

def Usage():
  print "\nfluvinput.py is a program that transformat input of variants positions in ... .\n"
  print "Usage:"
  print "\tfluvinput.py -i [TSV file] -s [H1N1|H3N2]\n"
  print "\nMandatory options:"
  print "\t-i, --input=FILE"
  print "\t\tTab delimited file that contains the nucleotide or protein positions to be highlighted."
  print "\t-s, --strain=H1N1|H3N2"
  print "\t\tThe strain in which the positions will be highlighted.Available options are H1N1 and H3N2."
  print "\t-o, --output-suffix=FILE"
  print "\t\tSuffix for the output files generated in the program."
  print "\t-t, --input-type=N|P"
  print "\t\tOption that indicates whether nucleotide [N] positions or amino acid [P] positions are indicaded. Nucleotides  [N] as default."

  print "\nOther options:"
  print "\t-f, --input-format=T|V"
  print "\t\tOption that indicates the type of input which could be TSV [T] or VCF [V]. By default TSV [V] is the input type."
  print "\t-l, --label-file=FILE"
  print "\t\tTab delimited file that contains the id conversion between user segments/protein id and the expected id for segments and proteins as:\n\n\t\t\tuser_id  segment_id[1-8]\tor\tuser_id  protein_id[PB1,PB2,PA,HA,NP,NA,M1,M2,NS1,NEP].\n"
  print "\t-r, --fasta-seq=FILE"
  print "\t\tUser FASTA file that is aligned to the reference sequence used in this software to proper pair positions in the used sequence to the reference."
  print "\t-p, --plot-2D=[T|F]"
  print "\t\tOption that indicate if a 2D plot of the variation across influenza segment or protein must be generated. By default is set as [F]alse."
  print "\t-c, --crystallographic=[T|F]"
  print "\t\tOption that indicate if a pymol session (*.pse) of the variations across crystallographic influenza protein structures is generated. By default is set as [F]alse."
  print "\t-h, --help"
  print "\t\tShow the options of the program."
  print "\n"
  sys.exit(1)

# Function that read and parse the command line arguments.
def SetOptions(argv):
  if len(argv) == 0:
    Usage()
  options, remaining = getopt.getopt(argv, 'i:s:f:t:l:r:o:p:c:h', ['input=','strain=','input-format=','input-type=','label-file=','fasta-seq=','output-suffix=','plot-2D','crystallographic','help'])
  opt_flag = {'i': False, 's':False,'f':False,'t':False,'l':False,'r':False,'o':False,'p':False,'c':False}
  global OPT_INPUT_FILE, OPT_INPUT_FORMAT, OPT_STRAIN, ALLOWED_STRAINS, OPT_INPUT_TYPE, OPT_LABEL_FILE, OPT_USER_FASTA_FILE, OPT_OUTPUT_SUFFIX, OPT_PYMOL_SESSION, OPT_2D_PLOT
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
    elif opt in ('-l', '--label-file'):
      if not opt_flag['l']:
        if os.path.exists(argu):
          OPT_LABEL_FILE = argu
          opt_flag['l'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: File or path of the label conversion file does not exist. ", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the label conversion file. Option -l / --label-file was already set.\n"
        sys.exit(1)
    elif opt in ('-r', '--fasta-seq'):
      if not opt_flag['r']:
        if os.path.exists(argu):
          OPT_USER_FASTA_FILE = argu
          opt_flag['r'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: File or path of the user FASTA file does not exist. ", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the user FASTA file. Option -r / --fasta-seq was already set.\n"
        sys.exit(1)
    elif opt in ('-p', '--plot-2D'):
      if not opt_flag['p']:
        if argu == "T":
          OPT_2D_PLOT = True
          opt_flag['p'] = True
        elif argu == "F":
          OPT_2D_PLOT = False
          opt_flag['p'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: Whether 2D plot is generated or not, only options [T]rue or [F]alse are allowed. Check option -p / --plot-2D. \n\tUnknown argument:", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine whether 2D plot is generated or not. Option -p / --plot-2D was already set.\n"
        sys.exit(1)
    elif opt in ('-c', '--crystallographic'):
      if not opt_flag['c']:
        if argu == "T":
          OPT_PYMOL_SESSION = True
          opt_flag['c'] = True
        elif argu == "F":
          OPT_PYMOL_SESSION = False
          opt_flag['c'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: Whether a pymol session (*.pse) of the variations across crystallographic influenza protein structures is generated or not, only allow options [T]rue or [F]alse. Check option -s / --structure. \n\tUnknown argument:", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine whether a pymol session (*.pse) of the variations across crystallographic influenza protein structures is generated or not. Option -s / --structure was already set.\n"
        sys.exit(1)
    elif opt in ('-s', '--strain'):
      if not opt_flag['s']:
        if argu.upper() in ALLOWED_STRAINS:
          OPT_STRAIN = argu.upper()
          opt_flag['s'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: Allowed strains are H1N1 and H3N2. The following strain is not supported: ", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the strain. Option -s / --strain was already set.\n"
        sys.exit(1)
    elif opt in ('-f', '--input-format'):
      if not opt_flag['f']:
        if argu.upper() in ["T","V"]:
          OPT_INPUT_FORMAT=argu.upper()
          opt_flag['f'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: Allowed input format are TSV [T] and VCF [V]. The following input type is not supported: ", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the input format. Option -f / --input-format was already set.\n"
        sys.exit(1)
    elif opt in ('-o', '--output-suffix'):
      if not opt_flag['o']:
        OPT_OUTPUT_SUFFIX = argu
        opt_flag['o'] = True
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the output suffix name. Option -o / --output-suffix was already set.\n"
        sys.exit(1)
    elif opt in ('-t', '--input-type'):
      if not opt_flag['t']:
        if argu.upper() in ["N","P"]:
          OPT_INPUT_TYPE=argu.upper()
          opt_flag['t'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: Allowed input types are Nucleotide [N] or Protein [P]. The following input type is not supported: ", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the input type. Option -t / --input-type was already set.\n"
        sys.exit(1)
    elif opt in ('-h', '--help'):
      Usage()
  if not opt_flag['i']:
    print >> sys.stderr , "[ERROR]: Input file is not defined. Option -i / --input.\n"
    sys.exit(1)
  if not opt_flag['s']:
    print >> sys.stderr , "[ERROR]: The strain is not defined. Option -s / --strain.\n"
    sys.exit(1)
  if not opt_flag['o']:
    print >> sys.stderr , "[ERROR]: Output suffix name is not defined. Option -o / --output-suffix.\n"
    sys.exit(1)
  if not opt_flag['t']:
    print >> sys.stderr , "[ERROR]: Input type was not defined either as [N]ucleotide or [P]rotein. Option -t / --input-type.\n"
    sys.exit(1)


###########################################
################   MAIN   #################
###########################################

# Parse command line
SetOptions(sys.argv[1:])

########### Append program and options for the execution of fluvinput sotfware ###########
cmd_fluvinput=['./fluvinput.py','-i',OPT_INPUT_FILE,'-s',OPT_STRAIN,'-t',OPT_INPUT_TYPE, '-f',OPT_INPUT_FORMAT]

if OPT_LABEL_FILE != None:
  cmd_fluvinput.append("-l")
  cmd_fluvinput.append(OPT_LABEL_FILE)

if OPT_USER_FASTA_FILE != None:
  cmd_fluvinput.append("-r")
  cmd_fluvinput.append(OPT_USER_FASTA_FILE)

# execute and wait for the program to finish
tmp_output_fluvinput = tempfile.NamedTemporaryFile()
#print "\n\ttemp: "+tmp_output_fluvinput.name+"\n"

print "Executing:\n\t"+" ".join(cmd_fluvinput)
proc = subprocess.Popen(cmd_fluvinput,stdout=tmp_output_fluvinput)
proc.wait()
exit_status = proc.wait()
if exit_status == 1: # the program in the subprocess exit with an error
  sys.exit(1)


########### Append program and options for the execution of fluvsearch sotfware ###########
cmd_fluvsearch=['./fluvsearch.py','-i',tmp_output_fluvinput.name,'-t',OPT_INPUT_TYPE, '-o',OPT_OUTPUT_SUFFIX]

if OPT_2D_PLOT:
  cmd_fluvsearch.append("-p")
  cmd_fluvsearch.append("T")

if OPT_PYMOL_SESSION:
  cmd_fluvsearch.append("-c")
  cmd_fluvsearch.append("T")

if OPT_USER_FASTA_FILE != None and OPT_INPUT_TYPE == "N":
  cmd_fluvsearch.append("-n")
  cmd_fluvsearch.append("UserSeqRefCodonTable.pkl")

# execute and wait for the program to finish
print "Executing:\n\t"+" ".join(cmd_fluvsearch)
proc = subprocess.Popen(cmd_fluvsearch)
exit_status = proc.wait()
if exit_status == 1: # the program in the subprocess exit with an error
  tmp_output_fluvinput.close() 
  sys.exit(1)

# close and erase tmp file
tmp_output_fluvinput.close()
