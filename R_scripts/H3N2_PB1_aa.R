# Rscript --vanilla H3N2_PB1_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=758
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

# PPI
Track_2<-ATrack(start=c(1,506), end=c(48,659),
                name ="Protein interactions",
                id=c("PA","PB2"),
                feature=c("PA","PB2"))


# Features
Track_3<-ATrack(start=c(187,211,1,233,267,297,396), end=c(187,211,139,249,757,312,483), 
                name ="Features",
                id=c("NLS","NLS","RNA binding","RNA binding","RNA binding","Polymerase activity","Polymerase activity"),
                feature=c("NLS","NLS","RNA binding","RNA binding","RNA binding","Polymerase activity","Polymerase activity"))
group(Track_3) <-c("NLS","NLS","RNA binding","RNA binding","RNA binding","Polymerase activity","Polymerase activity")

# Structure coverage
Structure_coverage<-ATrack(start=c(685), end=c(757), 
                           name ="Structure coverage",
                           id=c("2ZTT"),
                           feature=c("2ZTT"),
                           "2ZTT"="#00cc99")
group(Structure_coverage) <-c("2ZTT")

pdf("PB1_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_2, Track_3,highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,1,0.5),
           from=-5, to=prot.length+5,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "RNA binding" = "#d9d9d9",
           "Polymerase activity" = "#ff4dff",
           "NLS"="#000000",
           "PA"="#ffff66",
           "PB2"="#3333ff")
dev.off()
