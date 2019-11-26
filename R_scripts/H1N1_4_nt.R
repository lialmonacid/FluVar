# Rscript --vanilla H1N1_4_nt.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=1777
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
Track_1<-ATrack(start=c(1,33,1734,84,1065), end=c(32,1733,1777,1064,1730), 
                name ="Gene structure",
                id=c("5-NCR","HA-CDS","3-NCR","HA1","HA2"),
                feature=c("5-NCR","HA-CDS","3-NCR","HA1","HA2"))


# Antigenic sites
Track_2<-ATrack(start=c(453,540,558,633,579,690,786,492,744,291), end=c(458,554,575,668,593,698,794,509,749,308),
                name ="Antigenic sites",
                id=c("Sa","Sa","Sa","Sb","Ca1","Ca1","Ca1","Ca2","Ca2","Cb"),
                feature=c("Sa","Sa","Sa","Sb","Ca1","Ca1","Ca1","Ca2","Ca2","Cb"))
group(Track_2) <-c("Sa","Sa","Sa","Sb","Ca1","Ca1","Ca1","Ca2","Ca2","Cb")

# Features
Track_3<-ATrack(start=c(471,639,732,1065,1617), end=c(482,665,755,1115,1700),
                name ="Features",
                id=c("Receptor site","Receptor site","Receptor site","Fusion peptide","Transmembrane"),
                feature=c("Receptor site","Receptor site","Receptor site","Fusion peptide","Transmembrane"))
group(Track_3) <-c("Receptor site","Receptor site","Receptor site","Fusion peptide","Transmembrane")

# Structure coverage
Structure_coverage<-ATrack(start=c(84,1065), end=c(1052,1580), 
                           name ="Structure coverage",
                           id=c("3UBQ","3UBQ"),
                           feature=c("3UBQ","3UBQ"),
                           "3UBQ"="#00cc99")
group(Structure_coverage) <-c("3UBQ","3UBQ")

pdf("HA_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_1, Track_2, Track_3, highlighted_variants, pot, Structure_coverage),
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
           "Sa"="#009900",
           "Sb"="#990099",
           "Ca1"="#6699ff",
           "Ca2"="#ff9933",
           "Cb"="#9966ff")
dev.off()
