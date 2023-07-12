
 
#pe inverted results
in.pe.csv <- read.csv("~/Desktop/InvertedMethylationResults_PE.xls.csv")
tr.chromosomes <- in.pe.csv$ChrClean
startvalue <- in.pe.csv$ER_Start
endvalue <- in.pe.csv$Extraend
pevalues <- in.pe.csv$BW_Log2Rto
tracks = BioCircosBarTrack("bartrack", chromosomes= tr.chromosomes, starts= startvalue,ends= endvalue, values = pevalues, maxRadius = 0.9, minRadius = 0.7, color = "#40B9D4")

#pe normal results
pe.csv <- read.csv("~/Desktop/normal methylation data pe.xlsx - NormalMethylationResults_PE.xls.csv")
tr.chromosomes <- pe.csv$ChrClean
startvalue <- pe.csv$ER_Start
endvalue <- pe.csv$Extraend
pevalues <- pe.csv$BW_Log2Rto
tracks = tracks + BioCircosBarTrack("bartracktwo", chromosomes= tr.chromosomes, starts= startvalue,ends= endvalue, values = pevalues, maxRadius = 0.9, minRadius = 0.7, color = "orange")
# Add background
tracks = tracks + BioCircosBackgroundTrack("bars_background",maxRadius = 0.9, minRadius = 0.7, colors = "#2222EE")

#iugr inverted results
in.iugr.csv <- read.csv("~/Desktop/InvertedMethylationResults_IUGR.csv")
tr.chromosomes3 <- in.iugr.csv$ChrClean
startvalue3 <- in.iugr.csv$ER_Start
endvalue3 <- in.iugr.csv$ExtraEnd
pevalues3 <- in.iugr.csv$BW_Log2Rto
tracks = tracks + BioCircosBarTrack("bartrackthree", chromosomes= tr.chromosomes3, starts= startvalue3,ends= endvalue3, values = pevalues3, maxRadius = 0.65, minRadius = 0.45, color = "#40B9D4")

#iugr normal results
iugr.csv <- read.csv("~/Desktop/Normaliugr.csv")
tr.chromosomes3 <- iugr.csv$ChrClean
startvalue3 <- iugr.csv$ER_Start
endvalue3 <- iugr.csv$ExtraEnd
pevalues3 <- iugr.csv$BW_Log2Rto
tracks = tracks + BioCircosBarTrack("bartrackfour", chromosomes= tr.chromosomes3, starts= startvalue3,ends= endvalue3, values = pevalues3, maxRadius = 0.65, minRadius = 0.45, color = "orange")

# Add background
tracks = tracks + BioCircosBackgroundTrack("bars_background",maxRadius = 0.65, minRadius = 0.45, colors = "#2222EE")

#Add a Title
tracks = tracks + BioCircosTextTrack("testText", 'PE & IUGR Differentially Methylated Regions', weight = "lighter", 
                                     x = - .98, y = - 1.34)

BioCircos(tracks, genomeFillColor = "Spectral", chrPad = 0.02, displayGenomeBorder = FALSE, yChr =  FALSE, genomeTicksDisplay = FALSE,  genomeLabelTextSize = 15)




