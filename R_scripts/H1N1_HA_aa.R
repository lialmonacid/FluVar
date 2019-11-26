# Rscript --vanilla H1N1_4_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=567
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
Track_1<-ATrack(start=c(18,345), end=c(344,566), 
                name ="Protein structure",
                id=c("HA1","HA2"),
                feature=c("HA1","HA2"))

# Antigenic sites
Track_2<-ATrack(start=c(141,170,176,201,183,220,252,154,238,87), end=c(142,174,181,212,187,222,254,159,239,92),
                name ="Antigenic sites",
                id=c("Sa","Sa","Sa","Sb","Ca1","Ca1","Ca1","Ca2","Ca2","Cb"),
                feature=c("Sa","Sa","Sa","Sb","Ca1","Ca1","Ca1","Ca2","Ca2","Cb"))
group(Track_2) <-c("Sa","Sa","Sa","Sb","Ca1","Ca1","Ca1","Ca2","Ca2","Cb")

# Features
Track_3<-ATrack(start=c(147,203,234,345,529), end=c(150,211,241,361,556),
                name ="Features",
                id=c("Receptor site","Receptor site","Receptor site","Fusion peptide","Transmembrane"),
                feature=c("Receptor site","Receptor site","Receptor site","Fusion peptide","Transmembrane"))
group(Track_3) <-c("Receptor site","Receptor site","Receptor site","Fusion peptide","Transmembrane")

# Structure coverage
Structure_coverage<-ATrack(start=c(18,345), end=c(340,516), 
                           name ="Structure coverage",
                           id=c("3UBQ","3UBQ"),
                           feature=c("3UBQ","3UBQ"),
                           "3UBQ"="#00cc99")
group(Structure_coverage) <-c("3UBQ","3UBQ")

pdf("HA_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_1, Track_2, Track_3, highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,0.5,1,0.5),
           from=-5, to=prot.length+5,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
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
