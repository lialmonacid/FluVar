# Rscript --vanilla H3N2_7_nt.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=1027
pot<-ProteinAxisTrack(littleTicks = TRUE)

############ Position to Highlight ############
# Read selected position to highlight
positions_sel <- read.table(file = args[1] , sep = ";", header = FALSE)
positions_sel <- as.vector(t(positions_sel[1,]))

# Create data frame with selected positions across all the sequence
data_var <- as.data.frame( matrix(rep(c(0),times=prot.length), ncol=1) , stringsAsFactors=FALSE )
colnames(data_var)<-c("SelectedVariants")

# Mark position to highlight in the final figure.
data_var[positions_sel,]<-1

################################## Variants track ############################################
highlighted_variants <- DTrack(data = data_var$SelectedVariants, 
                               start = as.integer(rownames(data_var)), 
                               width=0, 
                               name = "Variations",
                               type="h")
displayPars(highlighted_variants) <- list(showAxis = FALSE)

################################## Annotation track ############################################

# Gene structure
Track_1<-ATrack(start=c(1,26,26,740,1008), end=c(25,784,51,1007,1027), 
                name ="Gene structure",
                id=c("5-NCR","M1-CDS","M2-CDS","M2-CDS","3-NCR"),
                feature=c("5-NCR","M1-CDS","M2-CDS","M2-CDS","3-NCR"),
                just.group = "above")
group(Track_1) <-c("5-NCR","M1-CDS","M2-CDS","M2-CDS","3-NCR")

# M1 Features
Track_2<-ATrack(start=c(326,467,26,518), end=c(340,517,253,781),
                name ="M1 features",
                id=c("NLS","RNA binding","RNP interacting","RNP interacting"),
                feature=c("NLS","RNA binding","RNP interacting","RNP interacting"),
                just.group = "above")
group(Track_2) <-c("NLS","RNA binding","RNP interacting","RNP interacting")

# M2 architecture
Track_3<-ATrack(start=c(26,740,786,846), end=c(51,785,845,899),
                name ="M2 architecture",
                id=c("Ectodomain","Ectodomain","Transmembrane","Amphipathic helix"),
                feature=c("Ectodomain","Ectodomain","Transmembrane","Amphipathic helix"),
                just.group = "above")
group(Track_3) <-c("Ectodomain","Ectodomain","Transmembrane","Amphipathic helix")

# M2 interactions
Track_4<-ATrack(start=c(924,984,852), end=c(932,995,878),
                name ="M2 interactions",
                id=c("M1","LC3","Caveolin"),
                feature=c("M1","LC3","Caveolin"),
                just.group = "above")

# M2 features
Track_5<-ATrack(start=c(849,867,873,849,867,879,876,882,891,876,882,894), end=c(851,869,875,851,869,881,878,884,893,878,884,896),
                name ="M2 CRAC motif",
                id=c("I","I","I",
                     "II","II","II",
                     "III","III","III",
                     "IV","IV","IV"),
                feature=c("I","I","I",
                          "II","II","II",
                          "III","III","III",
                          "IV","IV","IV"),
                just.group = "above")
group(Track_5) <-c("I","I","I",
                   "II","II","II",
                   "III","III","III",
                   "IV","IV","IV")

# Structure coverage
Structure_coverage<-ATrack(start=c(26,780), end=c(499,893),
                           name ="Structure coverage",
                           id=c("4PUS-M1","2RLF-M2"),
                           feature=c("4PUS-M1","2RLF-M2"),
                           just.group = "above",
                           "4PUS-M1"="#00cc99",
                           "2RLF-M2"="#008060")
group(Structure_coverage) <-c("4PUS-M1","2RLF-M2")

pdf("M_annotation.pdf",width = 9,height=6.2)
plotTracks(trackList=c(Track_1, Track_2, Track_3, Track_4, Track_5, highlighted_variants, pot, Structure_coverage),
           sizes = c(1,1,1,0.8,1,1,1,0.5),
           from=-10, to=prot.length+10,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           "5-NCR"="#000000",
           "M1-CDS"="#009999",
           "M2-CDS"="#cc6699",
           "3-NCR"="#000000",
           "RNA binding" = "#d9d9d9",
           "NLS"="#000000",
           "RNP interacting"="#3f0080",
           "Ectodomain"="#993399",
           "Transmembrane"="#c49a6c",
           "Amphipathic helix"="#d9d9d9",
           "M1"="#ed1c24",
           "LC3"="#993399",
           "Caveolin"="#993399",
           "I"="#00cc66",
           "II"="#d7df23",
           "III"="#9933ff",
           "IV"="#993399"
           )
dev.off()
