# Rscript --vanilla H3N2_NA_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=470
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

# Antigenic sites
# Track_2<-ATrack(start=c(441,540,633,207,894,678,261,309,855), end=c(515,557,671,239,917,734,272,359,872),
#                 name ="Antigenic sites",
#                 id=c("A","B","B","C","C","D","E","E","E"),
#                 feature=c("A","B","B","C","C","D","E","E","E"))
# group(Track_2) <-c("A","B","B","C","C","D","E","E","E")

# Features
Track_3<-ATrack(start=c(119,152,292), end=c(119,152,292),
                name ="Features",
                id=c("Catalytic residues","Catalytic residues","Catalytic residues"),
                feature=c("Catalytic residues","Catalytic residues","Catalytic residues"))
group(Track_3) <-c("Catalytic residues","Catalytic residues","Catalytic residues")

# Structure coverage
Structure_coverage<-ATrack(start=c(83), end=c(469),
                           name ="Structure coverage",
                           id=c("4B7R"),
                           feature=c("4B7R"),
                           "4B7R"="#00cc99")
group(Structure_coverage) <-c("4B7R")

pdf("NA_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_3, highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,1,0.5),
           from=-5, to=prot.length+5,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "Catalytic residues"="#ff0000")
dev.off()
