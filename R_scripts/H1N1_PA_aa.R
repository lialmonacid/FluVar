# Rscript --vanilla H1N1_PA_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=717
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

# PPI
Track_2<-ATrack(start=c(163,601,493,557), end=c(178,692,512,574),
                name ="Protein interactions",
                id=c("cRNA binding","PB1","hCLE","hCLE"),
                feature=c("cRNA binding","PB1","hCLE","hCLE"))
group(Track_2) <-c("cRNA binding","PB1","hCLE","hCLE")

# Features
Track_3<-ATrack(start=c(124,154,186), end=c(139,154,247), 
                name ="Features",
                id=c("NLS","NLS","NLS"),
                feature=c("NLS","NLS","NLS"))
group(Track_3) <-c("NLS","NLS","NLS")

# Structure coverage
Structure_coverage<-ATrack(start=c(1,257,354,398,558), end=c(197,348,371,549,716), 
                           name ="Structure coverage",
                           id=c("4AWH","2ZNL","2ZNL","2ZNL","2ZNL"),
                           feature=c("4AWH","2ZNL","2ZNL","2ZNL","2ZNL"),
                           "4AWH"="#00cc99",
                           "2ZNL"="#008060")
group(Structure_coverage) <-c("4AWH","2ZNL","2ZNL","2ZNL","2ZNL")

pdf("PA_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_2, Track_3, highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,1,0.5),
           from=-5, to=prot.length+5,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "cRNA binding"="#d9d9d9",
           "NLS"="#000000",
           "PB1"="#ff0000",
           "hCLE"="#993399")
dev.off()
