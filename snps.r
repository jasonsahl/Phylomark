library(Biostrings)
Args <- commandArgs();
alignment <- read.DNAStringSet(Args[4], "fasta")
var<-consensusMatrix(alignment)
y<-as.data.frame(var)
write.table(y, Args[5], quote = FALSE, sep = "\t")
masim <- sapply(names(alignment), function(x) sapply(names(alignment), function(y) pid(PairwiseAlignedXStringSet(alignment[[x]], alignment[[y]]), type="PID2")))
madist <- 1-(masim/100)
disttree <- hclust(as.dist(madist))
consCol <- function(myalign=alignment) {
        myalign <- as.matrix(myalign)
        consV <- sapply(seq(along=myalign[1,]), function(x) max(table(myalign[,x][!myalign[,x]
        %in% "-"]))/length(myalign[,1]))
        return(consV)
}
consV <- consCol(myalign=alignment)

slidingWindow <- function(mydata=consV, win=c(1, 100), start=1, end=length(consV)) {
        mydata <- mydata[start:end]
        windex <- t(sapply(0:length(mydata), function(x) win+x))
        mywind <- sapply(seq(along=mydata), function(x) mean(mydata[windex[x,1]:windex[x,2]],
        na.rm=TRUE))
        mywind <- mywind[1:(length(mywind)-win[2])]
        return(mywind)
}
slideV <- slidingWindow(mydata=consV, win=c(1, 20)) # Returns conservation vector for window size 20.
par(mar=c(2,4,4,9))
pdf(Args[6])
plot(slideV, type="l", lwd=2, col="red", main="Sliding Window Conservation Analysis", sub="Window Size: 20", xlab="Alignment Position", ylab="Conservation: max=1.0")
plot(as.dendrogram(disttree), edgePar=list(col=4, lwd=5), horiz=T)
dev.off()
