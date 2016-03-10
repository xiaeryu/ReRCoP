args <- commandArgs(TRUE)

## args[1] : removal matrix
## args[2] : output plot, in pdf format

## Read and format data
data <- read.table(args[1], header=T, row.names=1)
vec1 <- c(0,0,seq(ncol(data)-3,0,-1))
vec2 <- c(0,0,seq(ncol(data)-2.2,0.8,-1))
input <- rbind(vec1,vec2,data)

## Draw a plot
pdf(args[2])
par(mar=c(4,8,4,4))

# Frame
plot(c(1,1), xlim=c(0,max(input[,2])), ylim=c(0,ncol(input)-2), type='n', xlab='Position(bp)', ylab='', yaxt="n")

# Draw absent genes
for(i in seq(3, nrow(input))){
    for(j in seq(3, ncol(input))){
        rect(input[i,1], input[1,j], input[i,2], input[2,j], col=ifelse(input[i,j]==1, 'gray', 'white'), border=ifelse(input[i,j]==1, 'gray', 'white'))
    }
}

# Add text
for(j in seq(3,ncol(input))){
    mtext(text=colnames(input[j]), las=1, side=2, line=1, at=(input[1,j]+input[2,j])/2)
}

dev.off()
