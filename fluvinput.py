#!/usr/bin/python

import sys
import os
import getopt
import vcf
import pickle
from Bio import SeqIO
from Bio import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

# Variables that contains the inputs files 
OPT_INPUT_FILE=""
OPT_LABEL_FILE=None
OPT_INPUT_FORMAT="T"
OPT_INPUT_TYPE="N"
OPT_STRAIN=""
OPT_OUTPUT_SUFFIX=None
OPT_USER_FASTA_FILE=None
ALLOWED_STRAINS=["H1N1","H3N2"]
H1N1_ref_file={'N':"./ref_seq/H1N1_ref.nt.fasta",'P':"./ref_seq/H1N1_ref.aa.fasta"}
H3N2_ref_file={'N':"./ref_seq/H3N2_ref.nt.fasta",'P':"./ref_seq/H3N2_ref.aa.fasta"}
H1N1_cds_coordinates_file="./ref_seq/CDS_coordenates_H1N1.txt"
H3N2_cds_coordinates_file="./ref_seq/CDS_coordenates_H3N2.txt"

def Usage():
  print "\nfluvinput.py is a program that transformat input of variants positions in ... .\n"
  print "Usage:"
  print "\tfluvinput.py -i [TSV file] -s [H1N1|H3N2]\n"
  print "\nMandatory options:"
  print "\t-i, --input=FILE"
  print "\t\tTab delimited file that contains the nucleotide or protein positions to be highlighted."
  print "\t-s, --strain=H1N1|H3N2"
  print "\t\tThe strain in which the positions will be highlighted.Available options are H1N1 and H3N2."
  print "\nOther options:"
  print "\t-f, --input-format=T|V"
  print "\t\tOption that indicates the type of input which could be TSV [T] or VCF [V]. By default TSV [V] is the input type."
  print "\t-l, --label-file=FILE"
  print "\t\tTab delimited file that contains the id conversion between user segments/protein id and the expected id for segments and proteins as:\n\n\t\t\tuser_id  segment_id[1-8]\tor\tuser_id  protein_id[PB1,PB2,PA,HA,NP,NA,M1,M2,NS1,NEP].\n"
  print "\t-t, --input-type=N|P"
  print "\t\tOption that indicates whether nucleotide [N] positions or amino acid [P] positions are indicaded. Nucleotides  [N] as default."
  print "\t-r, --fasta-seq=FILE"
  print "\t\tUser FASTA file that is aligned to the reference sequence used in this software to proper pair positions in the used sequence to the reference."
  print "\t-o, --output-suffix=FILE"
  print "\t\tSuffix for the output files generated in the program."
  print "\t-h, --help"
  print "\t\tShow the options of the program."
  print "\n"
  sys.exit(1)

# Function that read and parse the command line arguments.
def SetOptions(argv):
  if len(argv) == 0:
    Usage()
  options, remaining = getopt.getopt(argv, 'i:s:f:t:l:r:o:h', ['input=','strain=','input-format=','input-type=','label-file=','fasta-seq=','output-suffix=','help'])
  opt_flag = {'i': False, 's':False,'f':False,'t':False,'l':False,'r':False,'o':False}
  global OPT_INPUT_FILE, OPT_INPUT_FORMAT, OPT_STRAIN, ALLOWED_STRAINS, OPT_INPUT_TYPE, OPT_LABEL_FILE,OPT_USER_FASTA_FILE, OPT_OUTPUT_SUFFIX
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
    elif opt in ('-o', '--output-suffix'):
      if not opt_flag['o']:
        OPT_OUTPUT_SUFFIX = argu
        opt_flag['o'] = True
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the output suffix name. Option -o / --output-suffix was already set.\n"
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

def getFeaturesDataLength(strain,input_type):
  features={}
  if strain == "H1N1":
    if input_type == "N":
      features['1'] = 2341
      features['2'] = 2341
      features['3'] = 2233
      features['4'] = 1777
      features['5'] = 1565
      features['6'] = 1458
      features['7'] = 1027
      features['8'] = 890
    else: # P, proteins
      features['PB2'] = 760
      features['PB1'] = 758
      features['PA'] = 717
      features['HA'] = 567
      features['NP'] = 499
      features['NA'] = 470
      features['M1'] = 253
      features['M2'] = 98
      features['NS1'] = 220
      features['NEP'] = 122
      features['NS2'] = 122 # dummy option because NEP and NS2 are the same.
  else: # H3N2
    if input_type == "N":
      features['1'] = 2341
      features['2'] = 2341
      features['3'] = 2233
      features['4'] = 1762
      features['5'] = 1567
      features['6'] = 1466
      features['7'] = 1027
      features['8'] = 890
    else: # P, proteins
      features['PB2'] = 760
      features['PB1'] = 758
      features['PA'] = 717
      features['HA'] = 567
      features['NP'] = 499
      features['NA'] = 470
      features['M1'] = 253
      features['M2'] = 98
      features['NS1'] = 231
      features['NEP'] = 122
      features['NS2'] = 122 # dummy option because NEP and NS2 are the same.
  return features

# Take as input an string and it compared with the protein names and segments numbers and 
# it return the corresponding segment number where the protein is coded. In case it is
# a segment number, it return the same segment number.
def getSegmentAndPosition(feature_name,feature_pos,input_type,features_and_lengths_allowed):
  if feature_name in features_and_lengths_allowed:
    if feature_pos < 1 or feature_pos > features_and_lengths_allowed[feature_name]:
      print >> sys.stderr , "\n[ERROR]: ",feature_name,feature_pos,"position out of range. It should be bewteen 1 and ",str(features_and_lengths_allowed[feature_name]),"\n"
      sys.exit(1)
  else:
    print >> sys.stderr , "\n[ERROR]: The following segment/protein is not recognized: ",feature_name, ".\n For option -t ",input_type," it should be one of the following: "
    for key in features_and_lengths_allowed:
      print >> sys.stderr , str(key)
    sys.exit(1)

  return feature_name,int(feature_pos)

# Function that takes an single string as input and return a list with two element [position,["A","T"]].
def readPosition(position):
  position_data = []
  try:
    pos_as_int = int(position)
    string_section = []
  except ValueError as e: # in case that position is something like 350A, which cannot be converted to int
    i = 0
    numeric_section = ""
    string_section = []
    encounter_non_numeric_symbol = False
    while i < len(position):
      if position[i].isdigit() and not encounter_non_numeric_symbol:
        numeric_section = numeric_section + position[i]
      else:
        string_section.append(position[i].upper())
      i = i + 1
    pos_as_int = int(numeric_section)
  position_data.append(pos_as_int)
  position_data.append(string_section)

  return position_data

# Function that read the tab-delimited file and extract the segments and the postion
def readTSV(file_input,strain,input_type,id_conversion):
  variants_per_feature={}
  input_order=[]
  nucleotide_used_in_variation = {}
  features_and_lengths=getFeaturesDataLength(strain,input_type)
  for line in open(file_input,"r"):
    fields=line.strip().split()

    if len(fields)!=2:
      print >> sys.stderr , "\n[ERROR]: \n\t\t",line.strip(),"\n\nInput file must be a tab-delimited file with two fields as shown below:\n\n\t\tSegment_1\t1538\n"
      sys.exit(1)

    if id_conversion!=None:
      try:
        feature_name=id_conversion[fields[0]]
      except KeyError as e:
        if (input_type=="N" and (fields[0] in ["1","2","3","4","5","6","7","8"])) or (input_type=="P" and (fields[0].upper() in ["PB2","PB1","PA","HA","NP","NA","M1","M2","NS1","NEP"])):
          feature_name = fields[0].upper()
        else:
          print >> sys.stderr , '\n[ERROR]: ID %s missing in the ID conversion file. Check file in option -l / --label-file \n' % str(e)
          sys.exit(1)
    else:
      feature_name=fields[0]

    feature_pos , single_specific_var = readPosition(fields[1])
    segment , position = getSegmentAndPosition(feature_name,feature_pos,input_type,features_and_lengths) # perform validation

    input_order.append(feature_name + "|" + str(feature_pos))
    #nucleotide_used_in_variation[feature_name+"|"+strain+"|"+str(feature_pos)] = single_specific_var

    nucleotide_used_in_variation[feature_name+"|"+str(feature_pos)] = single_specific_var

    if segment in variants_per_feature:
      variants_per_feature[segment].append(position)
    else:
      variants_per_feature[segment]=[position]
  return input_order, variants_per_feature, nucleotide_used_in_variation 

# F
def readALTrecord(alt, position, segment):
  #alt could be: [A] or [G, T] or [None] or [G, GTACT]
  list_alt=[]
  for nucleotide in alt:
    if nucleotide != None:
      if len(nucleotide) == 1: # avoid insertions 
        if nucleotide.upper() in ["A","C", "G", "T"]:
          list_alt.append(nucleotide.upper())
        else:
          print >> sys.stderr , "\n[WARNING]: \n\tSkipping alternative nucleotide \""+nucleotide+"\" at position \""+ str(position)+"\" in segment \""+str(segment)+"\" because only A, C, T and G are allowed.\n"
      else:
        print >> sys.stderr , "\n[WARNING]: \n\tSkipping alternative indel \""+nucleotide+"\" at position \""+ str(position)+"\" in segment \""+str(segment)+"\".\n"

  return list_alt

# Function that read the tab-delimited file and extract the segments and the postion
def readVCF(file_input,strain,input_type,id_conversion):
  variants_per_feature={}
  input_order=[]
  nucleotide_used_in_variation = {}
  features_and_lengths=getFeaturesDataLength(strain,input_type)
  vcf_reader = vcf.Reader(open(file_input, 'r'))
  for record in vcf_reader:

    if id_conversion!=None:
      try:
        feature_name=id_conversion[record.CHROM]
      except KeyError as e:
        if (input_type=="N" and (record.CHROM in ["1","2","3","4","5","6","7","8"])) or (input_type=="P" and (record.CHROM.upper() in ["PB2","PB1","PA","HA","NP","NA","M1","M2","NS1","NEP"])):
          feature_name = record.CHROM.upper()
        else:
          print >> sys.stderr , '\n[ERROR]: ID %s missing in the ID conversion file. Check file in option -l / --label-file \n' % str(e)
          sys.exit(1)
    else:
      feature_name=record.CHROM

    segment, position = getSegmentAndPosition(feature_name,int(record.POS),input_type,features_and_lengths)

    single_specific_var = readALTrecord(record.ALT)

    input_order.append(feature_name+"|"+str(record.POS))

    nucleotide_used_in_variation[feature_name+"|"+str(position)] = single_specific_var

    if segment in variants_per_feature:
      variants_per_feature[segment].append(position)
    else:
      variants_per_feature[segment]=[position]

  return input_order, variants_per_feature, nucleotide_used_in_variation


#
def strFormatingAltResi(pos,data_alt):
  # Data should be a list of strings
  if pos == "NO_REF":
    return str(pos)

  if data_alt == None:
    return str(pos)
  else:
    if len(data_alt) > 0:
      return str(pos)+"".join(data_alt)
    else:
      return str(pos)

# Function that print one variant per row as:
# Seg/protein  Position  Strain Position_corrected
def printFormattedOutput(input_order, variants, strain, position_convertion_table, output_suffix, nucleotides_in_variants):
  
  if output_suffix != None:
    f_h=open(output_suffix+".fluvinput.txt","w")

  for element in input_order:
    feature_name = element.split("|")[0]
    feture_pos = element.split("|")[1]
    original_pos = strFormatingAltResi(feture_pos,nucleotides_in_variants[element])
    if position_convertion_table != None:
      pos_relative_to_reference = strFormatingAltResi(position_convertion_table[feature_name][int(feture_pos)],nucleotides_in_variants[element])
    else:
      pos_relative_to_reference = original_pos

    output_string = str(feature_name)+"\t"+original_pos+"\t"+strain+"\t"+pos_relative_to_reference

    if output_suffix == None:
      print output_string
    else:
      f_h.write(output_string + "\n")
  
  if output_suffix != None:
    f_h.close()

# Function that read the id conversion file an load it in dictionary.
def readAndLoadLabelConversionFile(file_input,input_type):
  id_conversion={}
  for line in open(file_input,"r"):
    fields=line.strip().split()
    user_id=fields[0]
    expected_id=fields[1].upper()
    if input_type=="N": # Nucleotide
      if expected_id not in ["1","2","3","4","5","6","7","8"]:
        print >> sys.stderr , "\n[ERROR]: \t",line.strip(),"\n Id conversion has to pair \"" + user_id + "\" with a segment number [value of 1 to 8], which depict one of the influenza segments. Check file in option -l / --label-file.\n"
        sys.exit(1)
    else: # Protein
      if expected_id not in ["PB2","PB1","PA","HA","NP","NA","M1","M2","NS1","NEP"]:
        print >> sys.stderr , "\n[ERROR]: \t",line.strip(),"\n Id conversion has to pair \"" + user_id + "\" with a protein depicting one of the influenza proteins [PB2, PB1, PA, HA, NP, NA, M1, M2, NS1 and NEP]. Check file in option -l / --label-file.\n"
        sys.exit(1)

    id_conversion[user_id]=expected_id
  
  return id_conversion

# 
def readUserFasta(fasta_file,labels_conversion):
  single_dic={}
  handle = open(fasta_file, "rU")
  for record in SeqIO.parse(handle, "fasta"):

    fields = record.id.split("_")
    if (fields[0] in ["1","2","3","4","5","6","7","8"]) or (fields[0].upper() in ["PB2","PB1","PA","HA","NP","NA","M1","M2","NS1","NEP"]):
      record.id=fields[0].upper() # renumbering the id to the segment number
    else:
      fields = record.id.split("|")
      if (fields[0] in ["1","2","3","4","5","6","7","8"]) or (fields[0].upper() in ["PB2","PB1","PA","HA","NP","NA","M1","M2","NS1","NEP"]):
        record.id=fields[0].upper()
      else:
        if labels_conversion != None:
          try:
            putative_id = labels_conversion[record.id].upper()
            if (putative_id in ["1","2","3","4","5","6","7","8"]) or (putative_id in ["PB2","PB1","PA","HA","NP","NA","M1","M2","NS1","NEP"]):
              record.id=fields[0].upper()
            else:
              print >> sys.stderr , "\n[ERROR]: Unable to parse the sequences ID. It must be only the influenza segment number (1-8) or the protein name [PB2, PB1, PA, HA, NP, NA, M1, M2, NS1 and NEP]. Check data in \""+fasta_file+"\".\n"
              sys.exit(1)
          except KeyError as e:
            print >> sys.stderr , '\n[ERROR]: ID %s missing in the ID conversion file. Check IDs in the fasta reference file and in the ID conversion file.\n' % str(e)
            sys.exit(1)
        else:
          print >> sys.stderr , "\n[ERROR]: Unable to parse the sequences ID. It must be only the influenza segment number (1-8) or the protein name [PB2, PB1, PA, HA, NP, NA, M1, M2, NS1 and NEP]. Check data in \""+fasta_file+"\".\n"
          sys.exit(1)
    #record.id=record.id[0]
    single_dic[record.id]=record
  handle.close()
  return single_dic

# Function that read a single fasta file and return the a dictionary
# Each sequence id follow the 1_HxNx notation where the first number is the segment.
def readFasta(fasta_file):
  single_dic={}
  handle = open(fasta_file, "rU")
  for record in SeqIO.parse(handle, "fasta"):
    record.id=record.id.split("_")[0] # renumbering the id to the segment number
    #record.id=record.id[0]
    single_dic[record.id]=record
  handle.close()
  return single_dic

# Function that read a fasta file which belong the reference nucleotide sequence of H1N1 and H3N2 human influenza genome.
def readFastaRefSeqs(H1N1_ref_file, H3N2_ref_file):
  fasta_dic={}
  fasta_dic['H1N1']=readFasta(H1N1_ref_file)
  fasta_dic['H3N2']=readFasta(H3N2_ref_file)

  return fasta_dic

# function that align two set of sequence
def retrieveBestAlignmentOfTwoSequences(ref, seq):
  gap_open_penalty = -10
  gap_extend_penalty = -0.5

  aligments = pairwise2.align.globalds(ref,seq,matlist.blosum62,gap_open_penalty,gap_extend_penalty)
  ref_aln=aligments[0][0]
  seq_aln=aligments[0][1]
#  for protein_header in aligments_scores:
#      print >> sys.stderr, str(protein_header) + "\t" + str(aligments_scores[protein_header]) + "\t" + str(len(protein_frames[protein_header]['protein']))
#  if not found_highest_frame:
#      print >> sys.stderr, "[ERROR]: The best aligning frame could not be resolved"
#      sys.exit(1)
  return ref_aln,seq_aln

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
def searchProteinAndCodonPos(cds_coordinates_dict,pos):
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

#
def getCodonTableInfo(codon_table_dict,ref_cds_coordinates_dict,proteins_pos_list,strain,nucleotide_pos,nucleotide,segment):
  codon_table_dict_copied = codon_table_dict
  for protein_pos in proteins_pos_list:
    protein_fields=protein_pos.split(":")
    protein_name=protein_fields[0]
    protein_codon_number=int(protein_fields[1])
    if protein_name not in codon_table_dict_copied[strain]: # if protein does not exist yet
      codon_table_dict_copied[strain][protein_name]={protein_codon_number:[nucleotide,str(nucleotide_pos),None,None,None]}
    else:
      if protein_codon_number not in codon_table_dict_copied[strain][protein_name]: # if protein codon number does not exist yet
        codon_table_dict_copied[strain][protein_name][protein_codon_number]=[nucleotide,str(nucleotide_pos),None,None,None]
      else:
        if codon_table_dict_copied[strain][protein_name][protein_codon_number][2] == None: # if the secod position of the codon has not been filled
          codon_table_dict_copied[strain][protein_name][protein_codon_number][2]=str(nucleotide_pos)
          codon_table_dict_copied[strain][protein_name][protein_codon_number][0]=codon_table_dict_copied[strain][protein_name][protein_codon_number][0]+nucleotide
        else:
          if codon_table_dict_copied[strain][protein_name][protein_codon_number][3] == None: # if the third position of the codon has not been filled
            codon_table_dict_copied[strain][protein_name][protein_codon_number][3]=str(nucleotide_pos)
            codon_table_dict_copied[strain][protein_name][protein_codon_number][0]=codon_table_dict_copied[strain][protein_name][protein_codon_number][0]+nucleotide
            codon = codon_table_dict_copied[strain][protein_name][protein_codon_number][0]

            aa_code = Seq.translate(codon,to_stop=False,stop_symbol='*')
            codon_table_dict_copied[strain][protein_name][protein_codon_number][4] = aa_code
          else:
            print >> sys.stderr , "\n[ERROR]: The codon \""+str(protein_codon_number)+"\" is already set in the codon table as, "+codon_table_dict_copied[strain][protein_name][protein_codon_number]+". Contact the author because this is a major issue.\n"
            sys.exit(1)
  return codon_table_dict_copied

#
def saveCodonTableAsPickleFileStructure(data_object):
  fh_pkl_codon_table = open('UserSeqRefCodonTable' + '.pkl', 'wb')
  pickle.dump(data_object , fh_pkl_codon_table , pickle.HIGHEST_PROTOCOL)
  fh_pkl_codon_table.close()

# 
def findPositionRelationBetweenUserAndRefSeq(ref_seqs, user_seqs, segments_to_phase,ref_nt_cds_coordinates,nt_or_aa,strain):
  user_to_ref_pos={}
  codon_table={strain:{}}
  for segment in segments_to_phase:

    if segment not in user_seqs:
      print >> sys.stderr , "\n[ERROR]: Segment or protein \"" + str(segment) + "\" was not found in the USER Fasta file entered. Check file in option -r / --fasta-seq.\n"
      sys.exit(1)
    if segment not in ref_seqs:
      print >> sys.stderr , "\n[ERROR]: Major error because the reference file built-in the software does not contain the segment or protein \""+str(segment)+"\". Please contact the author of the sotfware.\n"
      sys.exit(1)

    ref_seq_aln,user_seq_aln = retrieveBestAlignmentOfTwoSequences(str(ref_seqs[segment].seq) , str(user_seqs[segment].seq))
    i=0
    ref_pos=1
    user_pos=1
    s_user_to_ref_pos={}
    while i<len(user_seq_aln):
      if user_seq_aln[i] != "-":
        if ref_seq_aln[i] != "-":
          s_user_to_ref_pos[user_pos]=ref_pos
          
          if nt_or_aa == "N":
            proteins_pos = searchProteinAndCodonPos(ref_nt_cds_coordinates[segment],ref_pos) # search if the position is within a CDS of the segment
            if len(proteins_pos)!=0: # for those nucleotide site that have protein positions
              codon_table = getCodonTableInfo(codon_table , ref_nt_cds_coordinates , proteins_pos , strain , user_pos , user_seq_aln[i] , segment) # update codon table
          
          user_pos=user_pos+1
          ref_pos=ref_pos+1
        else:
          #print ref_seq_aln
          #print user_seq_aln
          s_user_to_ref_pos[user_pos]="NO_REF"
          user_pos=user_pos+1
      else:
        if ref_seq_aln[i] != "-":
          ref_pos=ref_pos+1
      i=i+1
    user_to_ref_pos[segment]=s_user_to_ref_pos
  
  if nt_or_aa == "N": # Nucleotide
    saveCodonTableAsPickleFileStructure(codon_table)

  return user_to_ref_pos

# Function that convert the start and end position string into numeric tupples
def getStartEndPos(region):
  fields=region.strip().split("-")
  start=int(fields[0])
  end=int(fields[1])
  return (start,end)

# Function that read the annotation file and load it in a dictionary of dictionaries.
def readAnnotationFile(annotation_file):
  features={}
  for line in open(annotation_file,"r"):
    fields=line.strip().split()
    feature_name=fields[0]
    feature_annotation=fields[1]
    feature_annotation_pos=getStartEndPos(fields[2])
    if len(fields) != 3:
      print >> sys.stderr , "\n[ERROR]: \n\t\t",line.strip(),"\n\nAnnotation file must be a tab-delimited file with three fields as shown below:\n\n\t\tPB1\tNLS\t449-495\n"
      sys.exit(1)
    if feature_name in features:
      if feature_annotation in features[feature_name]:
        features[feature_name][feature_annotation].append(feature_annotation_pos)
      else:
        features[feature_name][feature_annotation]=[feature_annotation_pos]
    else:
      features[feature_name]= { feature_annotation : [feature_annotation_pos] }

  return features

# Function that read the nucleotide coordinates of the CDS of 10 influenza proteins  
def readCdsCoordinates(H1N1_cds_coordinates_file, H3N2_cds_coordinates_file):
  cds_coordinates_dic={}
  cds_coordinates_dic['H1N1'] = readAnnotationFile(H1N1_cds_coordinates_file)
  cds_coordinates_dic['H3N2'] = readAnnotationFile(H3N2_cds_coordinates_file)

  return cds_coordinates_dic

###########################################
################   MAIN   #################
###########################################

# Parse command line
SetOptions(sys.argv[1:])

# if the ID conversion file is set, the conversion table is loaded.
if OPT_LABEL_FILE != None:
  labels_conversion=readAndLoadLabelConversionFile(OPT_LABEL_FILE,OPT_INPUT_TYPE)
else:
  labels_conversion=None

# Read input file
if OPT_INPUT_FORMAT == "T": # For TSV file format. Default.
  input_order, variants, nucleotides_in_variants = readTSV(OPT_INPUT_FILE,OPT_STRAIN,OPT_INPUT_TYPE,labels_conversion)
elif OPT_INPUT_FORMAT == "V": # for CSV file format
  input_order, variants, nucleotides_in_variants = readVCF(OPT_INPUT_FILE,OPT_STRAIN,OPT_INPUT_TYPE,labels_conversion)

# Correct positions with a ref_strain
if OPT_USER_FASTA_FILE != None:
  ref_seqs = readFastaRefSeqs(H1N1_ref_file[OPT_INPUT_TYPE],H3N2_ref_file[OPT_INPUT_TYPE])
  ref_cds_coordinates = readCdsCoordinates(H1N1_cds_coordinates_file,H3N2_cds_coordinates_file) # as nucleotide. Check what happend if protein is used
  user_seq = readUserFasta(OPT_USER_FASTA_FILE,labels_conversion)
  position_convertion_table = findPositionRelationBetweenUserAndRefSeq(ref_seqs[OPT_STRAIN],user_seq,variants.keys(),ref_cds_coordinates[OPT_STRAIN],OPT_INPUT_TYPE,OPT_STRAIN)
else:
  position_convertion_table = None

# Print variants
printFormattedOutput(input_order,variants,OPT_STRAIN,position_convertion_table,OPT_OUTPUT_SUFFIX,nucleotides_in_variants)

