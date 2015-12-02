source('analysis/precision_recall.R')
library(ggplot2)
library(reshape)

depth.data <- function() {
    files <- c( "array2_indel.csv", "array2_snv_mnv.csv", "array2_sv.csv", "array3_indel.csv", "array3_snv_mnv.csv", "array3_sv.csv", "array4_indel.csv", "array4_snv_mnv.csv", "array4_sv.csv" )
    arrays <- c( "array2", "array2", "array2", "array3", "array3", "array3", "array4", "array4", "array4" )
    type <- c( "indel", "snv_mnv", "sv", "indel", "snv_mnv", "sv", "indel", "snv_mnv", "sv" )
    callers <- list( indel_callers, snv_callers, sv_callers, indel_callers, snv_callers, sv_callers, indel_callers, snv_callers, sv_callers)


    l <- list()

    for ( i in 1:length(files) ) {
        df <- ingest_csv(files[i], callers[[i]], keep.lowdepth=TRUE)
        df <- df[,c("val_tdepth","val_ndepth")]
        colnames(df) <- c("Tumour", "Normal")
        df$array <- arrays[i]
        df$variant <- type[i]

        l[[i]] <- df
    }

    data <- do.call(rbind, l)
    data$variant <- factor(data$variant, levels=c("snv_mnv","indel","sv"))
    return(data)
}

depth.plot <- function() {
    data <- depth.data()
    counts <- aggregate(data$Tumour, by=list(data$array,data$variant), FUN=length)
    colnames(counts) <- c("array","variant","ncalls")

    faildata <- data[data$Tumour < 2 & data$Normal < 2,]
    failcounts <- aggregate(faildata$Tumour, by=list(faildata$array,faildata$variant), FUN=length)
    counts$failcounts <- failcounts[,3]

    lowdata <- data[data$Tumour > 1 & data$Normal > 1 & (data$Tumour < 25 | data$Normal < 25),]
    lowcounts <- aggregate(lowdata$Tumour, by=list(lowdata$array,lowdata$variant), FUN=length)
    counts$lowcounts <- lowcounts[,3]

    counts$goodcounts <- counts$ncalls - counts$failcounts - counts$lowcounts
    counts$labels <- sapply(1:nrow(counts), function(x) paste(counts$failcounts[x],'/',counts$lowcounts[x],'/',counts$goodcounts[x],sep=''))
    counts$labels[1] <- paste("fail/low/good = ",counts$labels[1],sep='')
    counts$x <- 100
    counts$y <- 100000

    line <- counts[,c("array","variant")]
    line$x <- 25
    line <- rbind(line,line)
    n <- nrow(counts)
    line$y <- c(rep(1,n), rep(10000,n))

    data <- melt(data, id=c("array","variant"))
    colnames(data) <- c("array","variant","sample","depth")
    ggplot(data) + geom_density(aes(x=depth, color=sample, y=..count..)) + facet_grid(array~variant) + 
            geom_text(data=counts, aes(x=x,y=y, label=labels), size=3) + geom_line(data=line,aes(x=x,y=y),color="red",size=1.25) +
            scale_x_log10() + scale_y_log10(limits=c(1,100000)) 
    ggsave("plots/results/depth_plot.pdf", width=10, height=10)
}

depth.plot()
