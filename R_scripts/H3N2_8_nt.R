# Rscript --vanilla H3N2_8_nt.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=890
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

# Gene structure
Track_1<-ATrack(start=c(1,27,27,529,865), end=c(26,719,56,864,890),
                name ="Gene structure",
                id=c("5-NCR","NS1-CDS","NEP-CDS","NEP-CDS","3-NCR"),
                feature=c("5-NCR","NS1-CDS","NEP-CDS","NEP-CDS","3-NCR"),
                just.group = "above",
                mergeGroups=TRUE,
                groupAnnotation="id"
                )
group(Track_1)<-c("5-NCR","NS1-CDS","NEP-CDS","NEP-CDS","3-NCR")

# NS1 Protein interaction
Track_2<-ATrack(start=c(27,27,27,267,465,393,267), end=c(269,269,269,365,590,407,365),
                name ="NS1 interaction",
                id=c("PABP1","RIG-I","E1B-AP5","CPSF","CPSF","PKR","eIF4G1"),
                feature=c("PABP1","RIG-I","E1B-AP5","CPSF","CPSF","PKR","eIF4G1"),
                groupAnnotation="id",
                just.group = "above")
#group(Track_2) <-c("PABP1","RIG-I","E1B-AP5","CPSF","CPSF","PKR","eIF4G1")

# NS1 Features
Track_3<-ATrack(start=c(126,438,27), end=c(140,467,269),
                name ="NS1 features",
                id=c("NLS","NES","dsRNA"),
                feature=c("NLS","NES","dsRNA"),
                groupAnnotation="id",
                just.group = "above")

# NEP Domains
Track_4<-ATrack(start=c(27,529,658), end=c(56,571,861),
                name ="NEP domains",
                id=c("Protease sensitive N-terminal","Protease sensitive N-terminal","Protease resistant C-terminal"),
                feature=c("Protease sensitive N-terminal","Protease sensitive N-terminal","Protease resistant C-terminal"),
                mergeGroups=TRUE,
                groupAnnotation="id",
                just.group = "above")
group(Track_4)<-c("Protease sensitive N-terminal","Protease sensitive N-terminal","Protease resistant C-terminal")

# NEP Features
Track_5<-ATrack(start=c(688,778,532,565,730), end=c(753,843,561,573,732),
                name ="NEP features",
                id=c("C1 alpha-helix","C2 alpha-helix","NES","Phosphorylation serine-rich motif","M1 interacting"),
                feature=c("C1 alpha-helix","C2 alpha-helix","NES","Phosphorylation serine-rich motif","M1 interacting"),
                group=c("C1 alpha-helix","C2 alpha-helix","NES","Phosphorylation serine-rich motif","M1 interacting"),
                groupAnnotation="id",
                just.group = "above")

# Structure coverage
Structure_coverage<-ATrack(start=c(36,279,27,273,279,685), end=c(263,629,236,632,635,846),
                           name ="Structure coverage",
                           id=c("4OPA-NS1","4OPA-NS1","2ZKO-NS1","3L4Q-NS1","2RHK-NS1","1PD3-NEP"),
                           feature=c("4OPA-NS1","4OPA-NS1","2ZKO-NS1","3L4Q-NS1","2RHK-NS1","1PD3-NEP"),
                           groupAnnotation="id",
                           just.group = "above",
                           "4OPA-NS1"="#00cc99",
                           "2ZKO-NS1"="#008060",
                           "3L4Q-NS1"="#b3ffec",
                           "2RHK-NS1"="#00664d",
                           "1PD3-NEP"="#ccfff2")
group(Structure_coverage) <-c("4OPA-NS1","4OPA-NS1","2ZKO-NS1","3L4Q-NS1","2RHK-NS1","1PD3-NEP")

pdf("NS_annotation.pdf",width = 9,height=6.5)
plotTracks(trackList=c(Track_1, Track_2, Track_3, Track_4, Track_5, highlighted_variants, pot, Structure_coverage),
           #sizes = c(0.3,0.3,1,0.5,1),
           from=-7, to=prot.length+7,
           showFeatureId = FALSE,
           "5-NCR"="#000000",
           "NS1-CDS"="#009999",
           "NEP-CDS"="#cc6699",
           "3-NCR"="#000000",
           "PABP1"="#993399",
           "RIG-I"="#993399",
           "E1B-AP5"="#993399",
           "eIF4G1"="#993399",
           "CPSF"="#993399",
           "PKR"="#993399",
           "NLS"="#000000",
           "NES"="#000000",
           "dsRNA"="#d9d9d9",
           "Protease sensitive N-terminal"="#993399",
           "Protease resistant C-terminal"="#C49A6C",
           "C1 alpha-helix"="#DA1C5C",
           "C2 alpha-helix"="#2E3192",
           "Phosphorylation serine-rich motif"="#D7DF23",
           "M1 interacting"="#ED1C24"
)
dev.off()
