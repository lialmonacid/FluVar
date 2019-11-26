# Rscript --vanilla H1N1_5_nt.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=1565
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
Track_1<-ATrack(start=c(1,46,1543), end=c(45,1542,1565),
                name ="Gene structure",
                id=c("5-NCR","NP-CDS","3-NCR"),
                feature=c("5-NCR","NP-CDS","3-NCR"))


# PPI
Track_2<-ATrack(start=c(610,1156,46,808,46,1057,1258), end=c(1119,1539,528,1539,105,1059,1260),
                name ="Protein interactions",
                id=c("NP","NP","PB2","PB2","BAT1","Actin","Actin"),
                feature=c("NP","NP","PB2","PB2","BAT1","Actin","Actin"))
group(Track_2) <-c("NP","NP","PB2","PB2","BAT1","Actin","Actin")

# Features
Track_3<-ATrack(start=c(46,637,1051,46), end=c(84,693,1080,585),
                name ="Features",
                id=c("NLS","NLS","NAS","RNA binding"),
                feature=c("NLS","NLS","NAS","RNA binding"))
group(Track_3) <-c("NLS","NLS","NAS","RNA binding")

# Structure coverage
Structure_coverage<-ATrack(start=c(106,1258), end=c(1224,1539),
                           name ="Structure coverage",
                           id=c("3ZDP","3ZDP"),
                           feature=c("3ZDP","3ZDP"),
                           "3ZDP"="#00cc99")
group(Structure_coverage) <-c("3ZDP","3ZDP")

pdf("NP_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_1, Track_2, Track_3, highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,0.5,1,0.5),
           from=-20, to=prot.length+20,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "5-NCR"="#000000",
           "NP-CDS"="#ff0000",
           "3-NCR"="#000000",
           "NLS"="#000000",
           "NAS"="#cc9966",
           "RNA binding"="#d9d9d9",
           "NP"="#3366ff",
           "PB2"="#3333ff",
           "BAT1"="#993399",
           "Actin"="#ffffff")
dev.off()
