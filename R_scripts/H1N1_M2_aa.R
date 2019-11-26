# Rscript --vanilla H1N1_M2_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=98
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

# M2 architecture
Track_3<-ATrack(start=c(1,25,45), end=c(24,44,62),
                name ="M2 architecture",
                id=c("Ectodomain","Transmembrane","Amphipathic helix"),
                feature=c("Ectodomain","Transmembrane","Amphipathic helix"),
                just.group = "above")
group(Track_3) <-c("Ectodomain","Transmembrane","Amphipathic helix")

# M2 interactions
Track_4<-ATrack(start=c(71,91,47), end=c(73,94,55),
                name ="M2 interactions",
                id=c("M1","LC3","Caveolin"),
                feature=c("M1","LC3","Caveolin"),
                just.group = "above")


# M2 features
Track_5<-ATrack(start=c(46,52,54,46,52,56,55,57,60,55,57,61), end=c(46,52,54,46,52,56,55,57,60,55,57,61),
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
Structure_coverage<-ATrack(start=c(23), end=c(60),
                           name ="Structure coverage",
                           id=c("2RLF-M2"),
                           feature=c("2RLF-M2"),
                           just.group = "above",
                           "2RLF-M2"="#008060")
group(Structure_coverage) <-c("2RLF-M2")

pdf("M2_annotation.pdf",width = 9,height=6.2)
plotTracks(trackList=c(Track_3, Track_4, Track_5, highlighted_variants, pot, Structure_coverage),
           sizes = c(1,1,1,1,1,0.5),
           from=1, to=prot.length,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
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
