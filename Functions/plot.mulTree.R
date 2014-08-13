##########################
#Plots the results of a mulTree analysis
##########################
#Plots a boxplots of the fixed and random terms of the summarized multi tree MCMCglmm
#v0.1
##########################
#SYNTAX :
#<mulTree.mcmc> a mcmc chain written by the mulTree function. Can be either a unique file or a chain name referring to multiple files. Use read.mulTree() to properly load the chains




#<...> any additional argument to be passed to plot() function
#<horizontal> whether to plot the boxplots horizontally or not (default=TRUE)
#<add.table> a data.frame with the row names corresponding to the names of the "hdr.mcmc" to add to the plot (default=NULL). If add.table is not NULL, horizontal is set to TRUE.
#<LaTeX> whether to print the latex code
##########################
#----
#guillert(at)tcd.ie - 12/08/2014
##########################
#Requirements:
#-R 3
#-R package "MCMCglmm"
#-R package "coda"
#-R package "xtable" (optional)
##########################


plot.mulTree<-function(mulTree.mcmc, horizontal=FALSE, add.table=FALSE, LaTeX=FALSE ,     probs = c(95, 75, 50), xlabels = NULL, ylabels= NULL, type = "boxes", clr = gray((9:1)/10), scl = 1, xspc = 0.5, prn = FALSE, leg = FALSE, ct = "mode",ylims=NULL)
{
#HEADER
    require(MCMCglmm)
    require(coda)

#DATA
    #mulTree.mcmc
    if (class(mulTree.mcmc) != "table.mulTree") {
        dat<-table.mulTree(mulTree.mcmc)
    }





    #horizontal
    if(class(horizontal) != 'logical'){
        stop('"horizontal" must be logical.')
    }

    #add.table
    if(is.null(add.table)) {
        add.table=FALSE
    } else {
        #Checking if add.table is a data.frame object with the same column names as in hdr.mcmc
        if(class(add.table) != 'data.frame') {
            stop(as.character(substitute(add.table))," must be a \"data.frame\" object.", call.=FALSE) 
        } else {
            if(any(sort(names(hdr.mcmc)) != sort(row.names(add.table)))){
                stop(as.character(substitute(add.table))," must have the same row names as the names in ", as.character(substitute(hdr.mcmc)), call.=FALSE) 
            } else {
                table<-add.table
                add.table=TRUE
                horizontal=TRUE
            }
        }
    }

    #LaTeX
    if(class(LaTeX) != 'logical'){
        stop('"LaTeX" must be logical.')
    } else {
        if(LaTeX == TRUE) {
            require(xtable)
        }
    }

#FUNCTION

    #Density Plot function (from densityplot.R by Andrew Jackson - a.jackson@tcd.ie)
    FUN.densityplot <- function (dat, probs = c(95, 75, 50), xlabels = NULL, ylabels= NULL, type = "boxes", clr = gray((9:1)/10), scl = 1, xspc = 0.5, prn = FALSE, leg = FALSE, ct = "mode",ylims=NULL)
    {

    #dev.new()

    n <- ncol(dat)

    # test
    #probs = c(95, 75, 50)
    #xlabels = NULL
    #type = "boxes"
    #clr = gray((9:1)/10)
    #scl = 1
    #xspc = 0.5
    #prn = FALSE
    #leg = FALSE
    #ct = "mode"

    if (is.null(ylabels)){ylabels="value"}
        
    # Set up the plot
    if (is.null(ylims)){ylims<-c(min(dat) - 0.1*min(dat), max(dat) + 0.1*(max(dat)))}

    plot(1,1, xlab = "Source", ylab = ylabels, main = paste("","", sep = ""),
                xlim = c(1 - xspc, n + xspc),
                ylim = ylims, type = "n",
                xaxt = "n")
            if (is.null(xlabels)) {
                axis(side = 1, at = 1:n,
                    labels = (as.character(names(dat))))
            } else {
                axis(side = 1, at = 1:n,
                    labels = (xlabels))
            }



    clrs <- rep(clr, 5)
    for (j in 1:n) {
            temp <- hdr(dat[, j], probs, h = bw.nrd0(dat[,j]))
            line_widths <- seq(2, 20, by = 4) * scl
            bwd <- c(0.1, 0.15, 0.2, 0.25, 0.3) * scl
            if (prn == TRUE) {
                cat(paste("Probability values for Column", j, "\n"))
            }
            for (k in 1:length(probs)) {
                temp2 <- temp$hdr[k, ]
                if (type == "boxes") {
                    polygon(c(j - bwd[k], j - bwd[k], j + bwd[k], j + bwd[k]),
                      c(min(temp2[!is.na(temp2)]), max(temp2[!is.na(temp2)]), 
                      max(temp2[!is.na(temp2)]), min(temp2[!is.na(temp2)])),
                      col = clrs[k])
                    if (ct == "mode") {points(j,temp$mode,pch=19)}
                    if (ct == "mean") {points(j,mean(dat[,j]),pch=19)}
                    if (ct == "median") {points(j,median(dat[,j]),pch=19)}

                }
                if (prn == TRUE) {
                    cat(paste("\t", probs[k], "% lower =", format(max(min(temp2[!is.na(temp2)]),
                      0), digits = 2, scientific = FALSE), "upper =",
                      format(min(max(temp2[!is.na(temp2)]), 1), digits = 2,
                        scientific = FALSE), "\n"))
                }
            } # close the loop across probs
        } # close the loop across teh columns in dat
    }


#PLOTTING THE MCMCglmm RESULTS

    FUN.densityplot(dat)
#OUTPUT

#End
}