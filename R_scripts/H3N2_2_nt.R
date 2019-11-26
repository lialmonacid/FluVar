# Rscript --vanilla H3N2_2_nt.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=2341
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

################# Variants track ######################
highlighted_variants <- DTrack(data = data_var$SelectedVariants, 
                               start = as.integer(rownames(data_var)), 
                               width=0, 
                               name = "Variations",
                               type="h")
displayPars(highlighted_variants) <- list(showAxis = FALSE)

################# Annotation track ######################

# Gene structure
Track_1<-ATrack(start=c(1,25,2299), end=c(24,2298,2341), 
                name ="Gene structure",
                id=c("5-NCR","PB1-CDS","3-NCR"),
                feature=c("5-NCR","PB1-CDS","3-NCR"))


# PPI
Track_2<-ATrack(start=c(25,1540), end=c(168,2001),
                name ="Protein interactions",
                id=c("PA","PB2"),
                feature=c("PA","PB2"))


# Features
Track_3<-ATrack(start=c(583,655,25,721,823,913,1210), end=c(585,657,441,771,2295,960,1473), 
                name ="Features",
                id=c("NLS","NLS","RNA binding","RNA binding","RNA binding","Polymerase activity","Polymerase activity"),
                feature=c("NLS","NLS","RNA binding","RNA binding","RNA binding","Polymerase activity","Polymerase activity"))
group(Track_3) <-c("NLS","NLS","RNA binding","RNA binding","RNA binding","Polymerase activity","Polymerase activity")

# Structure coverage
Structure_coverage<-ATrack(start=c(2077), end=c(2295), 
                           name ="Structure coverage",
                           id=c("2ZTT"),
                           feature=c("2ZTT"),
                           "2ZTT"="#00cc99")
group(Structure_coverage) <-c("2ZTT")

pdf("PB1_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_1, Track_2, Track_3,highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,0.5,1,0.5),
           from=-40, to=prot.length+40,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "5-NCR"="#000000",
           "PB2-CDS"="#ff0000",
           "3-NCR"="#000000",
           "RNA binding" = "#d9d9d9",
           "Polymerase activity" = "#ff4dff",
           "NLS"="#000000",
           "PA"="#ffff66",
           "PB2"="#3333ff")
dev.off()
