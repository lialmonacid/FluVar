#!/usr/bin/env python

import sys
import os
import getopt
import vcf
import pickle
from Bio import SeqIO
from Bio import Seq
import subprocess
from pymol import cmd

# Variables that contains the inputs files 
OPT_LABEL_FILE=None
OPT_INPUT_FORMAT="T"
OPT_INPUT_TYPE="N"
OPT_INPUT_FILE=""
OPT_STRAIN=""
OPT_2D_PLOT = False
OPT_PYMOL_SESSION = False
OPT_OUTPUT_SUFFIX=""
OPT_CODON_TABLE_FILE=None
ALLOWED_STRAINS=["H1N1","H3N2"]
H1N1_ref_file={'N':"./ref_seq/H1N1_ref.nt.fasta",'P':"./ref_seq/H1N1_ref.aa.fasta"}
H3N2_ref_file={'N':"./ref_seq/H3N2_ref.nt.fasta",'P':"./ref_seq/H3N2_ref.aa.fasta"}
H1N1_cds_coordinates_file="./ref_seq/CDS_coordenates_H1N1.txt"
H3N2_cds_coordinates_file="./ref_seq/CDS_coordenates_H3N2.txt"
codon_table_ref_files={'H1N1':"./ref_seq/H1N1_codonTable.pkl",'H3N2':"./ref_seq/H3N2_codonTable.pkl"}
R_scripts_folder="./R_scripts/"
Pymol_structure_folder="./Structures/"
annot_coordinates_folder="./Annotation_coordenates/"

annot_coordinates={
  'N':{
    'H3N2':annot_coordinates_folder+"Nucleotide_Annotation_H3N2.txt",
    'H1N1':annot_coordinates_folder+"Nucleotide_Annotation_H1N1.txt"
  },
  'P':{
    'H3N2':annot_coordinates_folder+"Protein_Annotation_H3N2.txt",
    'H1N1':annot_coordinates_folder+"Protein_Annotation_H1N1.txt"
  }
}

R_scripts={
  'N':{
    'H3N2':{
      '1':R_scripts_folder+"H3N2_1_nt.R",
      '2':R_scripts_folder+"H3N2_2_nt.R",
      '3':R_scripts_folder+"H3N2_3_nt.R",
      '4':R_scripts_folder+"H3N2_4_nt.R",
      '5':R_scripts_folder+"H3N2_5_nt.R",
      '6':R_scripts_folder+"H3N2_6_nt.R",
      '7':R_scripts_folder+"H3N2_7_nt.R",
      '8':R_scripts_folder+"H3N2_8_nt.R"
    },
    'H1N1':{
      '1':R_scripts_folder+"H1N1_1_nt.R",
      '2':R_scripts_folder+"H1N1_2_nt.R",
      '3':R_scripts_folder+"H1N1_3_nt.R",
      '4':R_scripts_folder+"H1N1_4_nt.R",
      '5':R_scripts_folder+"H1N1_5_nt.R",
      '6':R_scripts_folder+"H1N1_6_nt.R",
      '7':R_scripts_folder+"H1N1_7_nt.R",
      '8':R_scripts_folder+"H1N1_8_nt.R"
    }
  },
  'P':{
    'H3N2':{
      'PB2':R_scripts_folder+"H3N2_PB2_aa.R",
      'PB1':R_scripts_folder+"H3N2_PB1_aa.R",
      'PA':R_scripts_folder+"H3N2_PA_aa.R",
      'HA':R_scripts_folder+"H3N2_HA_aa.R",
      'NP':R_scripts_folder+"H3N2_NP_aa.R",
      'NA':R_scripts_folder+"H3N2_NA_aa.R",
      'M1':R_scripts_folder+"H3N2_M1_aa.R",
      'M2':R_scripts_folder+"H3N2_M2_aa.R",
      'NS1':R_scripts_folder+"H3N2_NS1_aa.R",
      'NEP':R_scripts_folder+"H3N2_NEP_aa.R"
    },
    'H1N1':{
      'PB2':R_scripts_folder+"H1N1_PB2_aa.R",
      'PB1':R_scripts_folder+"H1N1_PB1_aa.R",
      'PA':R_scripts_folder+"H1N1_PA_aa.R",
      'HA':R_scripts_folder+"H1N1_HA_aa.R",
      'NP':R_scripts_folder+"H1N1_NP_aa.R",
      'NA':R_scripts_folder+"H1N1_NA_aa.R",
      'M1':R_scripts_folder+"H1N1_M1_aa.R",
      'M2':R_scripts_folder+"H1N1_M2_aa.R",
      'NS1':R_scripts_folder+"H1N1_NS1_aa.R",
      'NEP':R_scripts_folder+"H1N1_NEP_aa.R"
    }
  }
}

pymol_structures={
  'H3N2':{
    'PB2':[Pymol_structure_folder+"PB2_5FMQ_3CW4_2JDQ.pse"],
    'PB1':[Pymol_structure_folder+"PB1_2ZTT.pse"],
    'PA':[Pymol_structure_folder+"PA_4AWH.pse",Pymol_structure_folder+"PA_2ZNL.pse"],
    'HA':[Pymol_structure_folder+"HA_H3N2_3HMG.pse",Pymol_structure_folder+"HA_H3N2_1HTM.pse"],
    'NP':[Pymol_structure_folder+"NP_3ZDP.pse"],
    'NA':[Pymol_structure_folder+"NA_H1N1_4B7R.pse"],
    'M1':[Pymol_structure_folder+"M1_4PUS.pse"],
    'M2':[Pymol_structure_folder+"M2_2RLF.pse"],
    'NS1':[Pymol_structure_folder+"NS1_4OPA_2ZKO_3L4Q_2RHK.pse"],
    'NEP':[Pymol_structure_folder+"NEP_1PD3.pse"]
  },
  'H1N1':{
    'PB2':[Pymol_structure_folder+"PB2_5FMQ_3CW4_2JDQ.pse"],
    'PB1':[Pymol_structure_folder+"PB1_2ZTT.pse"],
    'PA':[Pymol_structure_folder+"PA_4AWH.pse",Pymol_structure_folder+"PA_2ZNL.pse"],
    'HA':[Pymol_structure_folder+"HA_H1N1_3UBQ.pse"],
    'NP':[Pymol_structure_folder+"NP_3ZDP.pse"],
    'NA':[Pymol_structure_folder+"NA_H1N1_4B7R.pse"],
    'M1':[Pymol_structure_folder+"M1_4PUS.pse"],
    'M2':[Pymol_structure_folder+"M2_2RLF.pse"],
    'NS1':[Pymol_structure_folder+"NS1_4OPA_2ZKO_3L4Q_2RHK.pse"],
    'NEP':[Pymol_structure_folder+"NEP_1PD3.pse"]
  }
}
def Usage():
  print "\nfluvsearch.py is a program that compare the annotation with a set of position across the influenza genome.\n"
  print "Usage:"
  print "\tfluvsearch.py -a [Annotation file] \n"
  print "\nMandatory options:"
  print "\t-i, --input-pos=FILE"
  print "\t\tTab delimited file with the positions that are going to be merged with the annotation."
  print "\t-o, --output-suffix=FILE"
  print "\t\tSuffix for the output files generated in the program."
  print "\t-t, --input-type=N|P"
  print "\t\tOption that indicates whether nucleotide [N] positions or amino acid [P] positions are indicaded. Nucleotides [N] as default."
  print "\nOther options:"
  print "\t-p, --plot-2D=[T|F]"
  print "\t\tOption that indicate if a 2D plot of the variation across influenza segment or protein must be generated. By default is set as [F]alse."
  print "\t-c, --crystallographic=[T|F]"
  print "\t\tOption that indicate if a pymol session (*.pse) of the variations across crystallographic influenza protein structures is generated. By default is set as [F]alse."
  print "\t-n, --nucleotide-refcodons=[FILE]"
  print "\t\tCodons file that are present in the nucleotide sequence generated with fluvinput.py programs (file is in pickle format). This option should be use when a influenza genome used to search the variants is different to the default one."
  print "\t-h, --help"
  print "\t\tShow the options of the program."
  print "\n"
  sys.exit(1)

# Function that read and parse the command line arguments.
def SetOptions(argv):
  if len(argv) == 0:
    Usage()
  options, remaining = getopt.getopt(argv, 'i:p:c:o:t:n:h', ['input-pos=','plot-2D=','crystallographic=','output-suffix=','input-type=','nucleotide-refcodons=','help'])
  opt_flag = {'i': False,'p': False,'c': False, 'o':False,'t':False,'n':False}
  global OPT_INPUT_FORMAT, OPT_STRAIN, ALLOWED_STRAINS, OPT_INPUT_TYPE, OPT_LABEL_FILE, OPT_INPUT_FILE, OPT_2D_PLOT, OPT_PYMOL_SESSION, OPT_OUTPUT_SUFFIX, OPT_CODON_TABLE_FILE
  for opt, argu in options:
    if opt in ('-i', '--input-pos'):
      if not opt_flag['i']:
        if os.path.exists(argu):
          OPT_INPUT_FILE = argu
          opt_flag['i'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: File or path of the input file does not exist. ", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the input file. Option -i / --input-pos was already set.\n"
        sys.exit(1)
    elif opt in ('-n', '--nucleotide-refcodons'):
      if not opt_flag['n']:
        if os.path.exists(argu):
          OPT_CODON_TABLE_FILE = argu
          opt_flag['n'] = True
        else:
          print >> sys.stderr , "\n[ERROR]: File or path of the codon table file (*.pkl) does not exist. ", argu, "\n"
          sys.exit(1)
      else:
        print >> sys.stderr , "\n[ERROR]: Trying to redefine the codon table file. Option -n / --nucleotide-refcodons was already set.\n"
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
    elif opt in ('-h', '--help'):
      Usage()
  if not opt_flag['i']:
    print >> sys.stderr , "[ERROR]: Input file is not defined. Option -i / --input-pos.\n"
    sys.exit(1)
  if not opt_flag['o']:
    print >> sys.stderr , "[ERROR]: Output suffix name is not defined. Option -o / --output-suffix.\n"
    sys.exit(1)
  if not opt_flag['t']:
    print >> sys.stderr , "[ERROR]: Input type was not defined either as [N]ucleotide or [P]rotein. Option -t / --input-type.\n"
    sys.exit(1)

# Function that takes an single string as input and return a list with two element [position,["A","T"]].
def readPosition(position):

  if position == "NO_REF":
    return ["NO_REF", []]

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
    #print position+"\t"+numeric_section
    pos_as_int = int(numeric_section)
  position_data.append(pos_as_int)
  position_data.append(string_section)
  
  return position_data

# Function that read the tab-delimited file and extract the segments and the postion
def readTSV(file_input):
  variants_per_feature={}
  input_order=[]
  nucleotide_used_in_variation = {}
  corrected_variants_relation = {}
  for line in open(file_input,"r"):
    fields=line.strip().split()

    if len(fields)!=4:
      print >> sys.stderr , "\n[ERROR]: \n\t\t",line.strip(),"\n\nInput file must be a tab-delimited file with four fields as exampled below:\n\n\t\tPB2\t1538\tH1N1\t\tor\t\t4\t1180\tH1N1 \n"
      sys.exit(1)

    feature_name=fields[0]
    feature_pos , single_specific_var = readPosition(fields[3]) # modified position so is relative to the buil-in reference
    feature_strain=fields[2]

    feature_pos_original , single_specific_var_original = readPosition(fields[1]) # original position to be searched

    input_order.append(feature_name+"|"+fields[1]+"|"+feature_strain+"|"+str(feature_pos))
    if feature_pos != "NO_REF":
      nucleotide_used_in_variation[feature_name+"|"+feature_strain+"|"+str(feature_pos)] = single_specific_var
      
      corrected_variants_relation[feature_name+"|"+feature_strain+"|"+str(feature_pos)] = feature_pos_original

      if feature_strain in variants_per_feature:
        if feature_name in variants_per_feature[feature_strain]:
          variants_per_feature[feature_strain][feature_name].append(feature_pos)
        else:
          variants_per_feature[feature_strain][feature_name]=[feature_pos]
      else:
        variants_per_feature[feature_strain]= { feature_name : [feature_pos] }

  return input_order,variants_per_feature,nucleotide_used_in_variation , corrected_variants_relation

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

# Function that read the annotation file and load it in a dictionary of dictionaries.
def readAnnotation(annotation_files):
  annotations={}
  for input_type in annotation_files: # N and P.
    annotations[input_type]={}
    for strain in annotation_files[input_type]: # H1N1 and H3N2
      annotations[input_type][strain]=readAnnotationFile(annotation_files[input_type][strain])

  return annotations

# 
def getAnnotationForPos(pos,features_dict):
  list_features_of_pos = []
  for feature in features_dict:
    for start,end in features_dict[feature]:
      if pos>=start and pos<=end:
        list_features_of_pos.append(feature)
        break
  if len(list_features_of_pos)!=0:
    return "; ".join(list_features_of_pos)
  else:
    return None

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

# Function that read the nucleotide coordinates of the CDS of 10 influenza proteins  
def readCdsCoordinates(H1N1_cds_coordinates_file, H3N2_cds_coordinates_file):
  cds_coordinates_dic={}
  cds_coordinates_dic['H1N1'] = readAnnotationFile(H1N1_cds_coordinates_file)
  cds_coordinates_dic['H3N2'] = readAnnotationFile(H3N2_cds_coordinates_file)

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

# 
def nucleotidePosToCodonPos(ref_cds_coordinates,ref_seqs):
  conversion_table={}
  for strain in ref_seqs: # H1N1 and H3N2
    conversion_table[strain] = {}
    for segment in ref_seqs[strain]:
      conversion_table[strain][segment] = {}
      i=1
      segment_length = len(ref_seqs[strain][segment].seq)
      while i <= segment_length: # for each positions in the segment
        proteins_pos = searchProteinAndCodonPos(ref_cds_coordinates[strain][segment],i)
        conversion_table[strain][segment][i] = proteins_pos
        i = i +1
  
  return conversion_table

#
def highlightsVariationsInPymolSession(pdb_file,list_positions,pymol_session_output):
  
  if not os.path.exists(pdb_file):
    print >> sys.stderr , "\n[ERROR]: File or path to the PDB file does not exist: ", pdb_file, "\n"
    sys.exit(1)

  print "reinitialize()"
  cmd.reinitialize()
  print "load(\"" + pdb_file + "\")"
  cmd.load(pdb_file)
  for pos in list_positions:
    selection_name = "var_" + str(pos)
    selection="(chain A and resi " + str(pos) + ")"
    cmd.select(selection_name , selection)
    count_atoms=cmd.count_atoms(selection_name)
    if count_atoms != 0:
      cmd.color("red", selection_name)
    else:
      cmd.delete(selection_name)
    cmd.deselect()
  print "save(\""+ pymol_session_output +"\")"
  cmd.save(pymol_session_output)

# Function that assing the corresponding amino numering to each position. This function might be usefull in the case of incomplete aa sequences, otherwise is a dummy function. 
def aminoacidPosToCodonPos(ref_seqs):
  conversion_table={}
  for strain in ref_seqs: # H1N1 and H3N2
    conversion_table[strain] = {}
    for protein in ref_seqs[strain]:
      conversion_table[strain][protein] = {}
      i=1
      protein_length = len(ref_seqs[strain][protein].seq)
      while i <= protein_length: # for each positions in the segment
        proteins_pos = i
        conversion_table[strain][protein][i] = [protein+":"+str(proteins_pos)]
        i = i +1
  return conversion_table

#
def getCodonInfo(aa_list,codon_table,strain):
  data=[]
  for element in aa_list:
    fields=element.split(":")
    feature_name = fields[0]
    feature_aa_pos = int(fields[1])
    data.append(codon_table[strain][feature_name][feature_aa_pos])
  return data

#
def strFormatingCodonInfo(data):
  # Data should be a list of list
  if data == None:
    return ""
  else:
    data_as_str=[]
    # Data is a list of list
    for element in data:
      element_as_str = ":".join(element)
      data_as_str.append(element_as_str)
    return "; ".join(data_as_str)

#
def strFormatingAltResi(data):
  # Data should be a list of strings
  if data == None:
    return ""
  else:
    if len(data) > 0:
      return "; ".join(data)
    else:
      return ""

#
def strFormatingCorrectedPos(data, alt_nucleotides):
  # data should be as follow:
  # data = 350 # integer
  # alt_nucleotides = ["A","C"] # list of single characters
  if len(alt_nucleotides) != 0:
    pos_formatted = str(data) + "".join(alt_nucleotides)
  else:
    pos_formatted = str(data)
  
  return pos_formatted

#
def readCodonTableFromPickleFile(pickle_file):
  global codon_table_ref_files
  if pickle_file != None:
    fh_pkl_codon_table = open(pickle_file, 'rb')
    codon_table = pickle.load(fh_pkl_codon_table)
    fh_pkl_codon_table.close()
  else: # use default file from the built-in reference genomes
    codon_table = {}
    for strain in codon_table_ref_files:
      fh_pkl_codon_table = open(codon_table_ref_files[strain], 'rb')
      ctable = pickle.load(fh_pkl_codon_table)
      fh_pkl_codon_table.close()
      codon_table.update(ctable)
  return codon_table

#
def getAltResidues(nucleotides_to_change_for, segment, position_in_seg, aa_list, codon_data, corrected_variants_pos_relation):
# nucleotides_to_change_for=['T', 'A'],  segment='5',  position_in_seg=350,  strain='H1N1',  aa_list=['NP:102'],  codon_data=[['GGA', '349', '350', '351', 'G']] ,corrected_variants_pos_relation=350
# nucleotides_to_change_for=['A'],  segment='2',  position_in_seg=350,  strain='H1N1',  aa_list=['PB1:96'],  [['GAA', '232', '233', '234', 'E']], corrected_variants_pos_relation=232
  
  original_position_in_seg = int(corrected_variants_pos_relation)

  alt_resi = []
  for i,prot_and_pos in enumerate(aa_list):
    prot , aa_pos = prot_and_pos.split(":")
    codon, codon_1, codon_2, codon_3, aa = codon_data[i]
    for nucleotide in nucleotides_to_change_for:
      if original_position_in_seg == int(codon_3):
        modif_codon = codon[0] + codon[1] + nucleotide
      elif original_position_in_seg == int(codon_2):
        modif_codon = codon[0] + nucleotide + codon[2]
      elif original_position_in_seg == int(codon_1):
        modif_codon = nucleotide + codon[1] + codon[2]
      else:
        print >> sys.stderr , "\n[ERROR]: Unable to find position \""+ str(original_position_in_seg)+"\" of segment \""+str(segment)+"\"in the codon data "+str(prot_and_pos)+"->:", codon_data, "\nPlease contact the author because is a major issue.\n"
        sys.exit(1)

      new_aa = Seq.translate(modif_codon, to_stop=False, stop_symbol='*')
      mut= aa + aa_pos + new_aa

      alt_resi.append(prot + ":" + str(original_position_in_seg) + nucleotide + ":" + mut)
  
  return alt_resi

#
def areValidNucleotides(nucleotides_to_change_for, input_type):
  allow_nucleotides = ["A", "C", "T", "G"]
  if input_type == "N":
    for key in nucleotides_to_change_for:
      segment, strain, pos_nt = key.split("|")
      for nucleotide in nucleotides_to_change_for[key]:
        if nucleotide not in allow_nucleotides:
          print >> sys.stderr , "\n[ERROR]: The provided alternative nucleotide \""+ str(nucleotide)+"\" of segment \""+str(segment)+"\" at position "+str(pos_nt)+" is not allowed. Only A, C, T and G nucleotides are permitted. Check the input files.\n"
          sys.exit(1)
          return False
    return True
  else: # P. Does not check for proteins aminoacids 
    return True


###########################################
################   MAIN   #################
###########################################

# Parse command line
SetOptions(sys.argv[1:])

# Read the Annotation file.
annotation=readAnnotation(annot_coordinates)

# Read the Reference Fasta file
ref_seqs=readFastaRefSeqs(H1N1_ref_file[OPT_INPUT_TYPE],H3N2_ref_file[OPT_INPUT_TYPE])
codon_table=readCodonTableFromPickleFile(OPT_CODON_TABLE_FILE)
if OPT_INPUT_TYPE == "N":
  ref_cds_coordinates=readCdsCoordinates(H1N1_cds_coordinates_file,H3N2_cds_coordinates_file)
  conversion_pos_table_AA=nucleotidePosToCodonPos(ref_cds_coordinates,ref_seqs)
else: # P
  conversion_pos_table_AA=aminoacidPosToCodonPos(ref_seqs)
  codon_table=None

# Read the Variation file.
input_order , variants , nucleotides_in_variants , corrected_variants_relation = readTSV(OPT_INPUT_FILE)

areValidNucleotides(nucleotides_in_variants,OPT_INPUT_TYPE)

#print input_order
#print variants
#print nucleotides_in_variants
#print corrected_variants_relation
#sys.exit(1)

# Comparing positions with annotation and printing
output={}
aa_positions={}
nt_positions={}
for key in variants: # strain in key
  aa_positions[key] = {}
  nt_positions[key] = {}
  for feature in variants[key]: # protein or segments in feature
    nt_positions[key][feature]=list(set(variants[key][feature]))
    nt_positions[key][feature].sort()
    for pos in variants[key][feature]:
      annotation_intersection=getAnnotationForPos(pos,annotation[OPT_INPUT_TYPE][key][feature])
      aa_list = conversion_pos_table_AA[key][feature][pos]

      if OPT_INPUT_TYPE == "N":
        codon_info = getCodonInfo(aa_list, codon_table, key)
        alt_resi = getAltResidues(nucleotides_in_variants[str(feature)+"|"+str(key)+"|"+str(pos)], feature, pos, aa_list, codon_info, corrected_variants_relation[str(feature)+"|"+str(key)+"|"+str(pos)])
      else:
        codon_info = None
        alt_resi = None

      if annotation_intersection != None:
        output[feature+"|"+str(pos)+"|"+key] = [feature , pos , aa_list ,codon_info , alt_resi, key , annotation_intersection]
      else:
        output[feature+"|"+str(pos)+"|"+key] = [feature , pos , aa_list ,codon_info , alt_resi, key , ""]

      for aa_pos in aa_list:
        protein,pos = aa_pos.split(":")
        if protein in aa_positions[key]:
          if int(pos) not in aa_positions[key][protein]:
            aa_positions[key][protein].append(int(pos))
            aa_positions[key][protein].sort()
        else:
          aa_positions[key][protein] = [int(pos)]

# Printing the output in the same order as the input.
f_var = open(OPT_OUTPUT_SUFFIX+"_var_annotated.txt","w")
for element in input_order:
  element_fields=element.split("|")
  actual_key = element_fields[0]+"|"+element_fields[3]+"|"+element_fields[2]
  
  if element_fields[3] != "NO_REF":
    s_feature_name = str(output[actual_key][0])
    s_original_pos = element_fields[1]
    s_pos_relative_to_builtin_ref = strFormatingCorrectedPos(output[actual_key][1],nucleotides_in_variants[element_fields[0]+"|"+element_fields[2]+"|"+element_fields[3]])
    s_protein_and_pos_list = "; ".join(output[actual_key][2])
    s_codon_data = strFormatingCodonInfo(output[actual_key][3])
    s_alt_resi = strFormatingAltResi(output[actual_key][4])
    s_strain = output[actual_key][5]
    s_annot = output[actual_key][6]
  else:
    s_feature_name = str(element_fields[0])
    s_original_pos = str(element_fields[1])
    s_pos_relative_to_builtin_ref = str(element_fields[3])
    s_strain = str(element_fields[2])
    s_protein_and_pos_list = s_codon_data = s_alt_resi = s_annot = ""

  output_line = s_feature_name+"\t"+s_original_pos+"\t"+s_pos_relative_to_builtin_ref+"\t"+s_protein_and_pos_list+"\t"+s_codon_data+"\t"+s_alt_resi+"\t"+s_strain+"\t"+s_annot

  #output_line = feature_name+"\t"+element_fields[1]+"\t"+strFormatingCorrectedPos(output[actual_key][1],nucleotides_in_variants[element_fields[0]+"|"+element_fields[2]+"|"+element_fields[3]])+"\t"+"; ".join(output[actual_key][2])+"\t"+strFormatingCodonInfo(output[actual_key][3])+"\t"+strFormatingAltResi(output[actual_key][4])+"\t"+output[actual_key][5]+"\t"+output[actual_key][6]
  f_var.write(output_line+"\n")
f_var.close()

# Section to generate the 2D plots in R software. 
if OPT_2D_PLOT:
  # run R to generate image
  if OPT_INPUT_TYPE == "N": # nucleotides
    for strain in nt_positions:
      for segment in nt_positions[strain]:
        # write the tmp file with the sites of variations to be highlight 
        var_file_name=strain+"."+str(segment)+".nt.txt"
        f = open(var_file_name,"w")
        pos_as_string=";".join(map(str,nt_positions[strain][segment]))
        f.write(pos_as_string+"\n")
        f.close()

        # Run the R script to produce the image with the variatons
        # Execution command example: Rscript --vanilla H1N1_4_nt.R INPUT_FILE OUTPUT_FILE
        print "Rscript --vanilla " + R_scripts[OPT_INPUT_TYPE][strain][segment] + " " + var_file_name
        proc = subprocess.Popen(['Rscript' , '--vanilla' , R_scripts[OPT_INPUT_TYPE][strain][segment] , var_file_name ],stdout=subprocess.PIPE)
        proc.wait()

        # Erase tmp file
  else:
    for strain in nt_positions:
      for segment in nt_positions[strain]: # in segment is the protein name
        # write the tmp file with the sites of variations to be highlight
        var_file_name=strain+"."+str(segment)+".aa.txt"
        f = open(var_file_name,"w")
        pos_as_string=";".join(map(str,nt_positions[strain][segment]))
        f.write(pos_as_string+"\n")
        f.close()

        # Run the R script to produce the image with the variatons
        # Execution command example: Rscript --vanilla H1N1_4_nt.R INPUT_FILE OUTPUT_FILE
        print "Rscript --vanilla " + R_scripts[OPT_INPUT_TYPE][strain][segment] + " " + var_file_name
        proc = subprocess.Popen(['Rscript' , '--vanilla' , R_scripts[OPT_INPUT_TYPE][strain][segment] , var_file_name ],stdout=subprocess.PIPE)
        proc.wait()

# Section to generate a pymol session to highlight the variants
if OPT_PYMOL_SESSION:
  for strain in aa_positions:
    for protein in aa_positions[strain]:
      for pdb_file in pymol_structures[strain][protein]:
        pymol_session_output = os.path.splitext(os.path.basename(pdb_file))[0]+"_var.pse"
        highlightsVariationsInPymolSession(pdb_file , aa_positions[strain][protein] , pymol_session_output)


