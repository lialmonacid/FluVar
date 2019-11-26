# Rscript --vanilla H3N2_NS1_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=231
pot<-ProteinAxisTrack(littleTicks = TRUE)

############ Position to Highlight ############
# Read selected position to highlight

positions_sel <- read.table(file = args[1] , sep = ";", header = FALSE)
#positions_sel <- read.table(file = "/Users/lalmon/Documents/PROGRAMS/FluVar/H3N2.8.nt.txt" , sep = ";", header = FALSE)
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

# NS1 Protein interaction
Track_2<-ATrack(start=c(1,1,1,81,147,123,81), end=c(81,81,81,113,188,127,113),
                name ="NS1 interaction",
                id=c("PABP1","RIG-I","E1B-AP5","CPSF","CPSF","PKR","eIF4G1"),
                feature=c("PABP1","RIG-I","E1B-AP5","CPSF","CPSF","PKR","eIF4G1"),
                groupAnnotation="id",
                just.group = "above")
#group(Track_2) <-c("PABP1","RIG-I","E1B-AP5","CPSF","CPSF","PKR","eIF4G1")

# NS1 Features
Track_3<-ATrack(start=c(34,138,1), end=c(38,147,81),
                name ="NS1 features",
                id=c("NLS","NES","dsRNA"),
                feature=c("NLS","NES","dsRNA"),
                groupAnnotation="id",
                just.group = "above")

# Structure coverage
Structure_coverage<-ATrack(start=c(4,85,1,83,85), end=c(79,201,70,202,203),
                           name ="Structure coverage",
                           id=c("4OPA-NS1","4OPA-NS1","2ZKO-NS1","3L4Q-NS1","2RHK-NS1"),
                           feature=c("4OPA-NS1","4OPA-NS1","2ZKO-NS1","3L4Q-NS1","2RHK-NS1"),
                           just.group = "above",
                           groupAnnotation="id",
                           "4OPA-NS1"="#00cc99",
                           "2ZKO-NS1"="#008060",
                           "3L4Q-NS1"="#b3ffec",
                           "2RHK-NS1"="#00664d")
group(Structure_coverage) <-c("4OPA-NS1","4OPA-NS1","2ZKO-NS1","3L4Q-NS1","2RHK-NS1")

pdf("NS1_annotation.pdf",width = 9,height=6.5)
plotTracks(trackList=c(Track_2, Track_3, highlighted_variants, pot, Structure_coverage),
           #sizes = c(0.3,0.3,1,0.5,1),
           from=-1, to=prot.length+1,
           showFeatureId = FALSE,
           "PABP1"="#993399",
           "RIG-I"="#993399",
           "E1B-AP5"="#993399",
           "eIF4G1"="#993399",
           "CPSF"="#993399",
           "PKR"="#993399",
           "NLS"="#000000",
           "NES"="#000000",
           "dsRNA"="#d9d9d9"
)
dev.off()
