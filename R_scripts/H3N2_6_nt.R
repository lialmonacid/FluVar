# Rscript --vanilla H3N2_6_nt.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=1466
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
Track_1<-ATrack(start=c(1,20,1430), end=c(19,1429,1466), 
                name ="Gene structure",
                id=c("5-NCR","NA-CDS","3-NCR"),
                feature=c("5-NCR","NA-CDS","3-NCR"))


# Antigenic sites
# Track_2<-ATrack(start=c(441,540,633,207,894,678,261,309,855), end=c(515,557,671,239,917,734,272,359,872),
#                 name ="Antigenic sites",
#                 id=c("A","B","B","C","C","D","E","E","E"),
#                 feature=c("A","B","B","C","C","D","E","E","E"))
# group(Track_2) <-c("A","B","B","C","C","D","E","E","E")

# Features
Track_3<-ATrack(start=c(374,473,893), end=c(376,475,895),
                name ="Features",
                id=c("Catalytic residues","Catalytic residues","Catalytic residues"),
                feature=c("Catalytic residues","Catalytic residues","Catalytic residues"))
group(Track_3) <-c("Catalytic residues","Catalytic residues","Catalytic residues")

# Structure coverage
Structure_coverage<-ATrack(start=c(266), end=c(1426),
                           name ="Structure coverage",
                           id=c("4B7R"),
                           feature=c("4B7R"),
                           "4B7R"="#00cc99")
group(Structure_coverage) <-c("4B7R")

pdf("NA_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_1, Track_3, highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,1,0.5),
           from=-20, to=prot.length+15,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "5-NCR"="#000000",
           "NA-CDS"="#ff0000",
           "3-NCR"="#000000",
           "Catalytic residues"="#ff0000")
dev.off()
