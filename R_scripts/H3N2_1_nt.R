# Rscript --vanilla 1.nt.R FILE 
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

highlighted_variants <- DTrack(data = data_var$SelectedVariants, 
                 start = as.integer(rownames(data_var)), 
                 width=0, 
                 name = "Variations",
                 type="h")
displayPars(highlighted_variants) <- list(showAxis = FALSE)

################# Annotation track ######################

# Gene structure
Track_1<-ATrack(start=c(1,28,2308), end=c(27,2307,2341), 
                      name ="Gene structure",
                      id=c("5-NCR","PB2-CDS","3-NCR"),
                      feature=c("5-NCR","PB2-CDS","3-NCR"))
group(Track_1) <-c("5-NCR","PB2-CDS","3-NCR")

# PPI
Track_2<-ATrack(start=c(28,1765,28,1765,2059), end=c(399,2304,834,2076,2304),
                             name ="Protein interactions",
                             id=c("PB1","PB1","NP","NP","ImportinAlpha5"),
                             feature=c("PB1","PB1","NP","NP","ImportinAlpha5"))
group(Track_2) <-c("PB1","PB1","NP","NP","ImportinAlpha5")

# Features
Track_3<-ATrack(start=c(751,1114,1237,1639,1372,2233,2284), end=c(873,1116,1239,1758,1512,2235,2286), 
                      name ="Features",
                      id=c("Cap-binding","Cap-binding","Cap-binding","Cap-binding","NLS","NLS","NLS"),
                      feature=c("Cap-binding","Cap-binding","Cap-binding","Cap-binding","NLS","NLS","NLS"))
group(Track_3) <-c("Cap-binding","Cap-binding","Cap-binding","Cap-binding","NLS","NLS","NLS")

# Structure coverage
Structure_coverage<-ATrack(start=c(781,1630,2083,2275), end=c(2055,2253,2253,2298), 
                           name ="Structure coverage",
                           id=c("5FMQ","3CW4","2JDQ","2JDQ"),
                           feature=c("5FMQ","3CW4","2JDQ","2JDQ"),
                           just.group = "above",
                           "5FMQ"="#00cc99",
                           "3CW4"="#008060",
                           "2JDQ"="#b3ffec")
group(Structure_coverage) <-c("5FMQ","3CW4","2JDQ","2JDQ")

pdf("PB2_annotation.pdf",width = 9,height=5.2)
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
           "Cap-binding"="#ff9933",
           "NLS"="#000000",
           "PB1"="#77b300",
           "NP"="#3366ff")
dev.off()
