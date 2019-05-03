#MetSkye : targeted signal extraction
# Author : Dinesh Kumar Barupal dinkumar@ucdavis.edu 12112018

install.packages("pacman")
install.packages("mzR")
pacman::p_load(mzR)

## Copy the mzML files and the target list into a directory. The path to the data location should not contain spaces.

# Initial parameters

directory <- "Complete path to the data"
rtDelta <- 0.15
mzDelta <- 0.008

## Read the data
dir1 <- getwd()
setwd(directory)

mzrtlist <- read.delim("LCMS_MZRT_file.txt",header = T, stringsAsFactors = F)

## Save mzR R objects.
filelist <- dir()
filelist <- filelist[grep(".mzML", filelist)]
for( j in 1:length(filelist)) {
  xfile <- mzR::openMSfile(filelist[j])
  peakTable <- header(xfile) # this gets table of details for each spectra
  spectraList <- spectra(xfile) # this gets the spectra values.
  save(peakTable, file=gsub(".mzML","_peaktable.RData",filelist[j]))
  save(spectraList, file = gsub(".mzML","_spectra.RData",filelist[j]))
}

## Extract Peak Heights from each file.

filelist <- dir(directory)
filelist <- filelist[grep("_peaktable.RData", filelist)]

for( i in 1:length(filelist)) {
  targetdf <- mzrtlist
  targetdf$sampleintensity  <- 0
  targetdf$deltart <- 0
  targetdf$ApexRT <- 0
  load(file=gsub("_peaktable.RData","_peaktable.RData",filelist[i]))
  load(file = gsub("_peaktable.RData","_spectra.RData",filelist[i]))
  rtvec <- round(peakTable$retentionTime[which(peakTable$msLevel==1)]/60,2)
  spectra.sb <- spectraList[which(peakTable$msLevel==1)]
  spectra.sb <- lapply(1:length(spectra.sb), function(x) {cbind(spectra.sb[[x]],rtvec[x])})
  spectra.sb <- data.frame(do.call(rbind, spectra.sb))
  names(spectra.sb)  <- c("mz","intensity","RT")
  spectra.sb <- spectra.sb[round(spectra.sb$mz,1)%in%round(targetdf$MZ,1),]
  for(j in 1:nrow(targetdf)) {
    mzindex <- which(spectra.sb$RT < targetdf$RT[j]+ rtDelta & spectra.sb$RT > targetdf$RT[j] - rtDelta & spectra.sb$mz < targetdf$MZ[j]+mzDelta & spectra.sb$mz > targetdf$MZ[j]-mzDelta)
    if(length(mzindex)!=0) {
      intenvec <- spectra.sb$intensity[mzindex]
      rtvec <- spectra.sb$RT[mzindex]
      targetdf$ApexRT[j] <- rtvec[which.max(intenvec)]
      targetdf$sampleintensity[j] <- max(intenvec)
    }
  }
  con1 <- file(gsub("_peaktable.RData","_peakHeights.txt",filelist[i]),"w")
  write.table(targetdf, con1, col.names = T, row.names = F, quote = F, sep="\t")
  close(con1)
}

## Merge all the files and export the data table.

filelist <- dir(directory)
filelist.input <- filelist[grep("_peaktable.RData", filelist)]
peakHeightFiles <- gsub("_peaktable.RData","_peakHeights.txt",filelist.input)

peakHeightFiles.list <- lapply(peakHeightFiles, function(x){
  read.delim(x, header = T, stringsAsFactors = F)
})
peakHeight_db <- read.delim(peakHeightFiles[1], header = T, stringsAsFactors = F)
con1 <- file("master_dataset_peakheight.txt","w")
writeLines(paste(c(names(peakHeight_db),gsub("_peaktable.RData","",filelist.input)),collapse = "\t"), con1)
for(i in 1:nrow(peakHeight_db)) {
  intenvec <- sapply(peakHeightFiles.list, function(x){x$sampleintensity[i]})
  writeLines(paste(c(as.character(peakHeight_db[i,]),round(intenvec, digits = 0)), collapse = "\t"), con1)
}
close(con1)
closeAllConnections()
setwd(dir1)
