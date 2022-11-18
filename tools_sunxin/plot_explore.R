in_arg = commandArgs(TRUE)

a = read.table(file = paste('explore_', in_arg, sep = ''), header = TRUE, 
               sep = '\t', stringsAsFactors = FALSE)

b = a

library(ggplot2)
library(Cairo)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}


CairoPDF(file = paste('plot_', in_arg, sep = ''), height = 8, width = 25)

for (i in 1:ncol(b))
{
    b[,i] = as.numeric(b[,i])
    sb = as.matrix(summary(as.numeric(b[,i])))
    ptitle = paste(paste(rownames(sb), collapse = '    '), '\n',
                   paste(sb[,1], collapse = '    '), sep = '')
    
    pl  = ggplot(data = b) + geom_density(aes(x=b[,i])) +
        xlab(colnames(b)[i]) + theme_bw() + ggtitle(ptitle)
    pll = pl + xlim(0, sb["3rd Qu.",])
    plll = pl + xlim(0, sb["Mean",] + sb["3rd Qu.",])
    print(multiplot(pl, pll, plll, cols =3))
}

dev.off()

