# Rscript --vanilla H3N2_NP_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=499
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
Track_2<-ATrack(start=c(189,371,1,255,1,338,405), end=c(358,498,161,498,20,338,405),
                name ="Protein interactions",
                id=c("NP","NP","PB2","PB2","BAT1","Actin","Actin"),
                feature=c("NP","NP","PB2","PB2","BAT1","Actin","Actin"))
group(Track_2) <-c("NP","NP","PB2","PB2","BAT1","Actin","Actin")

# Features
Track_3<-ATrack(start=c(1,198,336,1), end=c(13,216,345,180),
                name ="Features",
                id=c("NLS","NLS","NAS","RNA binding"),
                feature=c("NLS","NLS","NAS","RNA binding"))
group(Track_3) <-c("NLS","NLS","NAS","RNA binding")

# Structure coverage
Structure_coverage<-ATrack(start=c(21,405), end=c(393,498),
                           name ="Structure coverage",
                           id=c("3ZDP","3ZDP"),
                           feature=c("3ZDP","3ZDP"),
                           "3ZDP"="#00cc99")
group(Structure_coverage) <-c("3ZDP","3ZDP")

pdf("NP_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_2, Track_3, highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,1,0.5),
           from=-5, to=prot.length+5,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "NLS"="#000000",
           "NAS"="#cc9966",
           "RNA binding"="#d9d9d9",
           "NP"="#3366ff",
           "PB2"="#3333ff",
           "BAT1"="#993399",
           "Actin"="#ffffff")
dev.off()
