# Rscript --vanilla H3N2_M1_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=253
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

# M1 Features
Track_2<-ATrack(start=c(101,148,1,165), end=c(105,164,76,252),
                name ="M1 features",
                id=c("NLS","RNA binding","RNP interacting","RNP interacting"),
                feature=c("NLS","RNA binding","RNP interacting","RNP interacting"),
                just.group = "above")
group(Track_2) <-c("NLS","RNA binding","RNP interacting","RNP interacting")

# Structure coverage
Structure_coverage<-ATrack(start=c(1), end=c(158),
                           name ="Structure coverage",
                           id=c("4PUS-M1"),
                           feature=c("4PUS-M1"),
                           just.group = "above",
                           "4PUS-M1"="#00cc99")
group(Structure_coverage) <-c("4PUS-M1")

pdf("M1_annotation.pdf",width = 9,height=4.2)
plotTracks(trackList=c(Track_2, highlighted_variants, pot, Structure_coverage),
           sizes = c(1,1,1,0.5),
           from=1, to=prot.length,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           "RNA binding" = "#d9d9d9",
           "NLS"="#000000",
           "RNP interacting"="#3f0080"
           )
dev.off()
