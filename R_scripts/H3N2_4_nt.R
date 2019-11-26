# Rscript --vanilla H3N2_4_nt.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=1762
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
Track_1<-ATrack(start=c(1,30,1731,78,1065), end=c(29,1730,1762,1064,1727), 
                name ="Gene structure",
                id=c("5-NCR","HA-CDS","3-NCR","HA1","HA2"),
                feature=c("5-NCR","HA-CDS","3-NCR","HA1","HA2"))


# Antigenic sites
Track_2<-ATrack(start=c(441,540,633,207,894,678,261,309,855), end=c(515,557,671,239,917,734,272,359,872),
                name ="Antigenic sites",
                id=c("A","B","B","C","C","D","E","E","E"),
                feature=c("A","B","B","C","C","D","E","E","E"))
group(Track_2) <-c("A","B","B","C","C","D","E","E","E")

# Features
Track_3<-ATrack(start=c(468,636,729,1062,1617), end=c(479,662,752,1112,1691),
                name ="Features",
                id=c("Receptor site","Receptor site","Receptor site","Fusion peptide","Transmembrane"),
                feature=c("Receptor site","Receptor site","Receptor site","Fusion peptide","Transmembrane"))
group(Track_3) <-c("Receptor site","Receptor site","Receptor site","Fusion peptide","Transmembrane")

# Structure coverage
Structure_coverage<-ATrack(start=c(78,1065,78,1176), end=c(1061,1589,158,1589),
                           name ="Structure coverage",
                           id=c("3HMG","3HMG","1HTM","1HTM"),
                           feature=c("3HMG","3HMG","1HTM","1HTM"),
                           "3HMG"="#00cc99",
                           "1HTM"="#008060")
group(Structure_coverage) <-c("3HMG","3HMG","1HTM","1HTM")

pdf("HA_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_1, Track_2, Track_3,highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,0.5,1,0.5),
           from=-20, to=prot.length+20,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "5-NCR"="#000000",
           "HA-CDS"="#ff0000",
           "3-NCR"="#000000",
           "HA1"="#ffd1b3",
           "HA2"="#00cc99",
           "Receptor site"="#fa8072",
           "Fusion peptide"="#de1e08",
           "Transmembrane"="#fdd3ce",
           "A"="#009900",
           "B"="#990099",
           "C"="#6699ff",
           "D"="#ff9933",
           "E"="#9966ff")
dev.off()
