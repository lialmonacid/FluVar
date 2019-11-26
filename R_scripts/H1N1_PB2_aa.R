# Rscript --vanilla H1N1_PB2_aa.R FILE 
args <- commandArgs(trailingOnly = TRUE)

library(Pviz)

################## Axis ####################
prot.length=760
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

highlighted_variants <- DTrack(data = data_var$SelectedVariants, 
                 start = as.integer(rownames(data_var)), 
                 width=0, 
                 name = "Variations",
                 type="h")
displayPars(highlighted_variants) <- list(showAxis = FALSE)

################# Annotation track ######################

# PPI
Track_2<-ATrack(start=c(1,580,1,580,678), end=c(124,759,269,683,759),
                             name ="Protein interactions",
                             id=c("PB1","PB1","NP","NP","ImportinAlpha5"),
                             feature=c("PB1","PB1","NP","NP","ImportinAlpha5"))
group(Track_2) <-c("PB1","PB1","NP","NP","ImportinAlpha5")

# Features
Track_3<-ATrack(start=c(242,363,404,538,449,736,753), end=c(282,363,404,577,495,736,753), 
                      name ="Features",
                      id=c("Cap-binding","Cap-binding","Cap-binding","Cap-binding","NLS","NLS","NLS"),
                      feature=c("Cap-binding","Cap-binding","Cap-binding","Cap-binding","NLS","NLS","NLS"))
group(Track_3) <-c("Cap-binding","Cap-binding","Cap-binding","Cap-binding","NLS","NLS","NLS")

# Structure coverage
Structure_coverage<-ATrack(start=c(252,535,686,750), end=c(676,742,742,757), 
                           name ="Structure coverage",
                           id=c("5FMQ","3CW4","2JDQ","2JDQ"),
                           feature=c("5FMQ","3CW4","2JDQ","2JDQ"),
                           "5FMQ"="#00cc99",
                           "3CW4"="#008060",
                           "2JDQ"="#b3ffec")
group(Structure_coverage) <-c("5FMQ","3CW4","2JDQ","2JDQ")

pdf("PB2_annotation.pdf",width = 9,height=5.2)
plotTracks(trackList=c(Track_2, Track_3,highlighted_variants, pot, Structure_coverage),
           sizes = c(0.5,0.5,0.5,1,0.5),
           from=-5, to=prot.length+5,
           featureAnnotation="feature",
           groupAnnotation="feature",
           showFeatureId = FALSE,
           just.group = "above",
           "Cap-binding"="#ff9933",
           "NLS"="#000000",
           "PB1"="#77b300",
           "NP"="#3366ff")
dev.off()
