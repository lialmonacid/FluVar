# Rscript --vanilla H3N2_NEP_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=122
pot<-ProteinAxisTrack(littleTicks = TRUE)

############ Position to Highlight ############
# Read selected position to highlight

positions_sel <- read.table(file = args[1] , sep = ";", header = FALSE)
#positions_sel <- read.table(file = "/Users/lalmon/Documents/PROGRAMS/FluVar/H3N2.NEP.aa.txt" , sep = ";", header = FALSE)
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

# NEP Domains
Track_4<-ATrack(start=c(1,54), end=c(53,121),
                name ="NEP domains",
                id=c("Protease sensitive N-terminal","Protease resistant C-terminal"),
                feature=c("Protease sensitive N-terminal","Protease resistant C-terminal"),
                mergeGroups=TRUE,
                groupAnnotation="id",
                just.group = "above")
group(Track_4)<-c("Protease sensitive N-terminal","Protease resistant C-terminal")

# NEP Features
Track_5<-ATrack(start=c(64,94,12,23,78), end=c(85,115,21,25,78),
                name ="NEP features",
                id=c("C1 alpha-helix","C2 alpha-helix","NES","Phosphorylation serine-rich motif","M1 interacting"),
                feature=c("C1 alpha-helix","C2 alpha-helix","NES","Phosphorylation serine-rich motif","M1 interacting"),
                group=c("C1 alpha-helix","C2 alpha-helix","NES","Phosphorylation serine-rich motif","M1 interacting"),
                groupAnnotation="id",
                just.group = "above")

# Structure coverage
Structure_coverage<-ATrack(start=c(63), end=c(116),
                           name ="Structure coverage",
                           id=c("1PD3-NEP"),
                           feature=c("1PD3-NEP"),
                           just.group = "above",
                           groupAnnotation="id",
                           "1PD3-NEP"="#ccfff2")
group(Structure_coverage) <-c("1PD3-NEP")

pdf("NEP_annotation.pdf",width = 8.8,height=4.2)
plotTracks(trackList=c(Track_4, Track_5, highlighted_variants, pot, Structure_coverage),
           sizes = c(1,1,1,1,1),
           from=1, to=prot.length,
           showFeatureId = FALSE,
           "NES"="#000000",
           "Protease sensitive N-terminal"="#993399",
           "Protease resistant C-terminal"="#C49A6C",
           "C1 alpha-helix"="#DA1C5C",
           "C2 alpha-helix"="#2E3192",
           "Phosphorylation serine-rich motif"="#D7DF23",
           "M1 interacting"="#ED1C24"
)
dev.off()
