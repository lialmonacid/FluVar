#!/usr/bin/env python

import sys
import os
import getopt
from Bio import SeqIO
from Bio import Seq
import subprocess
import pickle

# Variables that contains the inputs files 
OPT_INPUT_FILE=""
OPT_CDS_COORD_FILE=""
OPT_STRAIN=None
ALLOWED_STRAINS=["H1N1","H3N2"]
H1N1_ref_file={'N':"./ref_seq/H1N1_ref.nt.fasta",'P':"./ref_seq/H1N1_ref.aa.fasta"}
H3N2_ref_file={'N':"./ref_seq/H3N2_ref.nt.fasta",'P':"./ref_seq/H3N2_ref.aa.fasta"}
H1N1_cds_coordinates_file="./CDS_coordenates_H1N1.txt"
H3N2_cds_coordinates_file="./CDS_coordenates_H3N2.txt"

def Usage():
  print "\ngenerateCodonTableStructureFromRef.py is a program that generate the codon table of a reference file. Only for internal use.\n"
  print "Usage:"
  print "\tgenerateCodonTableStructureFromRef.py -i [FASTA_FILE] -c [CDS_COORDINATES FILE]\n"
  print "\nMandatory options:"
  print "\t-i, --input=FILE"
  print "\t\tFASTA File used as input. It must be nuecleotides sequences with and specific id."
  print "\t-c, --coordinates=[FILE]"
  print "\t\tIt is a file that contain the coordenates of all protein cds in the FASTA file in a specific format."
  print "\t-s, --strain=H1N1|H3N2"
  print "\t\tThe influenza strain of fasta file. Available options are H1N1 and H3N2."
  print "\t-h, --help"
  print "\t\tShow the options of the program."
  print "\n"
  sys.exit(1)

# Function that read and parse the command line arguments.
def SetOptions(argv):
  if len(argv) == 0:
    Usage()
  options, remaining = getopt.getopt(argv, 'i:c:s:h', ['input=','coordinatesc=','strain=','help'])
  opt_flag = {'i': False,'c': False,'s':False}
  global ALLOWED_STRAINS, OPT_INPUT_FILE,OPT_CDS_COORD_FILE, OPT_STRAIN
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
    elif opt in ('-c', '--coordinates'):
      if not opt_flag['c']:
        if os.path.exists(argu):
          OPT_CDS_COORD_FILE = argu
          opt_flag['c'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: File or path of the cds coordinates file does not exist. ", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the cds coordinates file. Option -c / --coordinates was already set.\n"
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
    elif opt in ('-h', '--help'):
      Usage()
  if not opt_flag['i']:
    print >> sys.stderr , "[ERROR]: Input FASTA ffle is not defined. Option -i / --input.\n"
    sys.exit(1)
  if not opt_flag['c']:
    print >> sys.stderr , "[ERROR]: The cds coordinates file is not defined. Option -c / --coordinates.\n"
    sys.exit(1)
  if not opt_flag['s']:
    print >> sys.stderr , "[ERROR]: The influenza strain of the fasta file is not defined. Option -s / --strain.\n"
    sys.exit(1)

# Function that read a single fasta file and return the a dictionary
# Each sequence id follow the 1_HxNx notation where the first number is the segment.
def readFasta(fasta_file):
  single_dic={}
  handle = open(fasta_file, "rU")
  for record in SeqIO.parse(handle, "fasta"):
    record.id=record.id.split("_")[0] # renumbering the id to the segment number
    #record.id=record.id[0]dd
    single_dic[record.id]=record
  handle.close()
  return single_dic

# Function that read a fasta file which belong the reference nucleotide sequence of H1N1 and H3N2 human influenza genome.
def readFastaRefSeqs(fasta_file,strain):
  fasta_dic={}
  fasta_dic[strain]=readFasta(fasta_file)

  return fasta_dic

# Function that convert the start and end position string into numeric tupples
def getStartEndPos(region):
  fields=region.strip().split("-")
  start=int(fields[0])
  end=int(fields[1])
  return (start,end)

# Function that read the annotation file and load it in a dictionary of dictionaries.
def readAnnotationFile(annotation_file):
  features = {}
  for line in open(annotation_file,"r"):
    fields = line.strip().split()
    feature_name = fields[0] # should be a number between 1 and 8.
    feature_annotation = fields[1].upper() # should be a one of the influenza protein names [PB2, PB1, PA, HA, NP, NA, M1, M2, NS1 and NEP].
    feature_annotation_pos = getStartEndPos(fields[2])
    if len(fields) != 3:
      print >> sys.stderr , "\n[ERROR]: \n\t\t",line.strip(),"\n\nCDS coordinates file must be a tab-delimited file with three fields as shown below:\n\n\t\t4\tHA\t\33-1733\n"
      sys.exit(1)

    if feature_name not in ["1","2","3","4","5","6","7","8"]:
      print >> sys.stderr , "\n[ERROR]: \n\t\t",line.strip(),"\n\nCDS coordinates file in the first field must have an number between 1 and 8, depicting the influenza segments.\n"
      sys.exit(1)

    if feature_annotation not in ["PB2", "PB1", "PA", "HA", "NP", "NA", "M1", "M2", "NS1", "NEP"]:
      print >> sys.stderr , "\n[ERROR]: \n\t\t",line.strip(),"\n\nCDS coordinates file in the second field must have the one of the influenza protein names [PB2, PB1, PA, HA, NP, NA, M1, M2, NS1 or NEP].\n"
      sys.exit(1)

    if feature_name in features:
      if feature_annotation in features[feature_name]:
        features[feature_name][feature_annotation].append(feature_annotation_pos)
      else:
        features[feature_name][feature_annotation] = [feature_annotation_pos]
    else:
      features[feature_name]= { feature_annotation : [feature_annotation_pos] }

  return features

# Function that read the nucleotide coordinates of the CDS of 10 influenza proteins  
def readCdsCoordinates(cds_coordinates_file,strain):
  cds_coordinates_dic={}
  cds_coordinates_dic[strain] = readAnnotationFile(cds_coordinates_file)

  return cds_coordinates_dic

#
def lengthPreviousExons(exons,current_exon_index):
  length_previous_exons=0
  i=current_exon_index-1
  while i>= 0:
    exon_start = exons[i][0]
    exon_end = exons[i][1]
    if exon_end < exon_start:
      print >> sys.stderr , "\n[ERROR]: \n\t\t",str(exons[i]),"\n\nCoordinates in the annotation file must be the smallest coordinate and then the biggest, for example,  330-500.\n"
      sys.exit(1)
    length_previous_exons = length_previous_exons + ((exon_end-exon_start)+1)
    i = i - 1

  return length_previous_exons

#
def searchProteinAndCodonPos(cds_coordinates_dict, pos):
  proteins_pos=[]
  for protein in cds_coordinates_dict:
    for ele,exon in enumerate(cds_coordinates_dict[protein]):
      exon_start=exon[0]
      exon_end=exon[1]
      if pos >= exon_start and pos <= exon_end:
        cds_i = ( pos - exon_start ) + 1 + lengthPreviousExons(cds_coordinates_dict[protein],ele)
        cds_i_mode = cds_i%3
        putative_codon_num = int(cds_i)/int(3)
        if cds_i_mode != 0:
          codon_num=int(putative_codon_num)+1
        else:
          codon_num=int(putative_codon_num)
        proteins_pos.append(protein+":"+str(codon_num))

  return proteins_pos

def getMatchingSegmentToProtein(protein_name,strain,ref_cds_coordinates_dict):
  protein_name_copy = protein_name.upper()
  for segment in ref_cds_coordinates_dict[strain]:
    if protein_name_copy in ref_cds_coordinates_dict[strain][segment].keys():
      return segment
  
  print >> sys.stderr , "\n[ERROR]: Unknown protein \""+str(protein_name)+"\" so no matching segment was found. Contact the author because it is a major error.\n"
  sys.exit(1)


def getCodonTableInfo(codon_table_dict,ref_cds_coordinates_dict,ref_seqs_dict,proteins_pos_list,strain,nucleotide_pos):
  codon_table_dict_copied = codon_table_dict
  for protein_pos in proteins_pos_list:
    protein_fields=protein_pos.split(":")
    protein_name=protein_fields[0]
    segment=getMatchingSegmentToProtein(protein_name,strain,ref_cds_coordinates_dict)
    protein_codon_number=int(protein_fields[1])
    if protein_name not in codon_table_dict_copied[strain]: # if protein does not exist yet
      codon_table_dict_copied[strain][protein_name]={protein_codon_number:[None,str(nucleotide_pos),None,None,None]}
    else:
      if protein_codon_number not in codon_table_dict_copied[strain][protein_name]: # if protein codon number does not exist yet
        codon_table_dict_copied[strain][protein_name][protein_codon_number]=[None,str(nucleotide_pos),None,None,None]
      else:
        if codon_table_dict_copied[strain][protein_name][protein_codon_number][2] == None: # if the secod position of the codon has not been filled
          codon_table_dict_copied[strain][protein_name][protein_codon_number][2]=str(nucleotide_pos)
        else:
          if codon_table_dict_copied[strain][protein_name][protein_codon_number][3] == None: # if the third position of the codon has not been filled
            codon_table_dict_copied[strain][protein_name][protein_codon_number][3]=str(nucleotide_pos)
            codon_letter_1st = str(ref_seqs_dict[strain][segment].seq[int(codon_table_dict_copied[strain][protein_name][protein_codon_number][1])-1])
            codon_letter_2nd = str(ref_seqs_dict[strain][segment].seq[int(codon_table_dict_copied[strain][protein_name][protein_codon_number][2])-1])
            codon_letter_3rd = str(ref_seqs_dict[strain][segment].seq[int(codon_table_dict_copied[strain][protein_name][protein_codon_number][3])-1])
            codon = codon_letter_1st + codon_letter_2nd + codon_letter_3rd
            codon_table_dict_copied[strain][protein_name][protein_codon_number][0]=codon
            
            aa_code = Seq.translate(codon,to_stop=False,stop_symbol='*')
            codon_table_dict_copied[strain][protein_name][protein_codon_number][4] = aa_code
          else:
            print >> sys.stderr , "\n[ERROR]: The codon \""+str(protein_codon_number)+"\" is already set in the codon table as, "+codon_table_dict_copied[strain][protein_name][protein_codon_number]+". Contact the author because this is a major issue.\n"
            sys.exit(1)
  return codon_table_dict_copied

# 
def nucleotidePosToCodonPos(ref_cds_coordinates,ref_seqs):
  #print ref_cds_coordinates
  #print ref_seqs
  codon_table={}
  conversion_table={}
  for strain in ref_seqs: # H1N1 and H3N2
    codon_table[strain] = {}
    conversion_table[strain] = {}
    for segment in ref_seqs[strain]:
      conversion_table[strain][segment] = {}
      i=1
      segment_length = len(ref_seqs[strain][segment].seq)
      while i <= segment_length: # for each positions in the segment
        proteins_pos = searchProteinAndCodonPos(ref_cds_coordinates[strain][segment],i)
        conversion_table[strain][segment][i] = proteins_pos
        if len(proteins_pos)!=0: # for those nucleotide site that have protein positions
          codon_table=getCodonTableInfo(codon_table,ref_cds_coordinates,ref_seqs,proteins_pos,strain,i) # update codon table
        i = i +1
  
  #subset_codon_table=codon_table['H1N1']['M2'].items()
  #subset_codon_table.sort(key = lambda x: x[0])
  #print subset_codon_table
  return codon_table,conversion_table

###########################################
################   MAIN   #################
###########################################

# Parse command line
SetOptions(sys.argv[1:])

# Read the Reference Fasta file
ref_seqs=readFastaRefSeqs(OPT_INPUT_FILE,OPT_STRAIN)
ref_cds_coordinates=readCdsCoordinates(OPT_CDS_COORD_FILE,OPT_STRAIN)

codon_table,conversion_pos_table_AA=nucleotidePosToCodonPos(ref_cds_coordinates,ref_seqs)

# Save the object as pickle file for further use.
fh_pkl_codon_table = open(OPT_STRAIN + '_codonTable' + '.pkl', 'wb')
pickle.dump(codon_table , fh_pkl_codon_table , pickle.HIGHEST_PROTOCOL)
fh_pkl_codon_table.close()


