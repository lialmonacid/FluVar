# Rscript --vanilla H1N1_3_nt.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=2233
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
Track_1<-ATrack(start=c(1,25,2176), end=c(24,2175,2233), 
                name ="Gene structure",
                id=c("5-NCR","PA-CDS","3-NCR"),
                feature=c("5-NCR","PA-CDS","3-NCR"))


# PPI
Track_2<-ATrack(start=c(511,1825,1501,1693), end=c(558,2100,1560,1746),
                name ="Protein interactions",
                id=c("cRNA binding","PB1","hCLE","hCLE"),
                feature=c("cRNA binding","PB1","hCLE","hCLE"))
group(Track_2) <-c("cRNA binding","PB1","hCLE","hCLE")

# Features
Track_3<-ATrack(start=c(394,484,580), end=c(441,486,765), 
                name ="Features",
                id=c("NLS","NLS","NLS"),
                feature=c("NLS","NLS","NLS"))
group(Track_3) <-c("NLS","NLS","NLS")

# Structure coverage
Structure_coverage<-ATrack(start=c(25,793,1084,1216,1696), end=c(615,1068,1137,1671,2172), 
                           name ="Structure coverage",
                           id=c("4AWH","2ZNL","2ZNL","2ZNL","2ZNL"),
                           feature=c("4AWH","2ZNL","2ZNL","2ZNL","2ZNL"),
                           "4AWH"="#00cc99",
                           "2ZNL"="#008060")
group(Structure_coverage) <-c("4AWH","2ZNL","2ZNL","2ZNL","2ZNL")

pdf("PA_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_1, Track_2, Track_3,highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,0.5,1,0.5),
           from=-35, to=prot.length+20,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "5-NCR"="#000000",
           "PA-CDS"="#ff0000",
           "3-NCR"="#000000",
           "cRNA binding"="#d9d9d9",
           "NLS"="#000000",
           "PB1"="#ff0000",
           "hCLE"="#993399")
dev.off()
