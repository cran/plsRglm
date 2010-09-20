plots.confints.bootpls = function (ic_bootobject, indices = NULL, legendpos = "topleft", 
    prednames = TRUE, articlestyle = TRUE, xaxisticks=TRUE, ...) 
{
    nr <- nrow(ic_bootobject)
    if (is.null(indices)) {
        indices <- 1:nr
    }
    plotpos <- (1:nr)[1:length(indices)]
    if (articlestyle) {
        oldpar <- par()
        par(mar = c(2, 2, 1, 1) + 0.1, mgp = c(2, 1, 0))
    }
    plot(c(1, 1), xlab = "", ylab = "", type = "n", xlim = c(1, 
        length(indices) + 0.5), ylim = c(min(ic_bootobject), 
        max(ic_bootobject)), xaxt = "n", ...)
    if (attr(ic_bootobject, "typeBCa")) {
        if (prednames) {
            if(xaxisticks){
              axis(1, at = plotpos + 0.225, labels = rownames(ic_bootobject)[indices])
            }
            else
            {
              axis(1, at = plotpos + 0.225, labels = rownames(ic_bootobject)[indices],lwd.ticks=0)
            }
        }
        else {
            if(xaxisticks){
              axis(1, at = plotpos + 0.225, labels = paste("x", (1:nr)[indices], sep = ""))
            }
            else
            {
              axis(1, at = plotpos + 0.225, labels = paste("x", (1:nr)[indices], sep = ""),lwd.ticks=0)
            }
        }
        arrows(plotpos, ic_bootobject[indices, 1], plotpos, ic_bootobject[indices, 
            2], lend = "butt", lwd = 2, lty = 2, col = "darkgreen", 
            code = 3, angle = 90, length = 0.1)
        arrows(plotpos + 0.15, ic_bootobject[indices, 3], plotpos + 
            0.15, ic_bootobject[indices, 4], lend = "butt", lwd = 2, 
            lty = 4, col = "blue", code = 3, angle = 90, length = 0.1)
        arrows(plotpos + 0.3, ic_bootobject[indices, 5], plotpos + 
            0.3, ic_bootobject[indices, 6], lend = "butt", lwd = 2, 
            lty = 5, col = "red", code = 3, angle = 90, length = 0.1)
        arrows(plotpos + 0.45, ic_bootobject[indices, 7], plotpos + 
            0.45, ic_bootobject[indices, 8], lend = "butt", lwd = 2, 
            lty = 1, col = "black", code = 3, angle = 90, length = 0.1)
        abline(h = 0, lty = 3, lwd = 2)
        legend(legendpos, legend = c("Normal", "Basic", "Percentile", 
            "BCa"), lty = c(2, 4, 5, 1), col = c("darkgreen", 
            "blue", "red", "black"), lwd = 2)
    }
    else {
        if (prednames) {
          if(xaxisticks){
            axis(1, at = plotpos + 0.15, labels = rownames(ic_bootobject)[indices])
          }
          else
          {
            axis(1, at = plotpos + 0.15, labels = rownames(ic_bootobject)[indices],lwd.ticks=0)          
          }
        }
        else {
          if(xaxisticks){
            axis(1, at = plotpos + 0.15, labels = paste("x",(1:nr)[indices], sep = ""))
            }
          else
          {
            axis(1, at = plotpos + 0.15, labels = paste("x",(1:nr)[indices], sep = ""),lwd.ticks=0)        
          }
        }
        arrows(plotpos, ic_bootobject[indices, 1], plotpos, ic_bootobject[indices, 
            2], lend = "butt", lwd = 2, lty = 2, col = "darkgreen", 
            code = 3, angle = 90, length = 0.1)
        arrows(plotpos + 0.15, ic_bootobject[indices, 3], plotpos + 
            0.15, ic_bootobject[indices, 4], lend = "butt", lwd = 2, 
            lty = 4, col = "blue", code = 3, angle = 90, length = 0.1)
        arrows(plotpos + 0.3, ic_bootobject[indices, 5], plotpos + 
            0.3, ic_bootobject[indices, 6], lend = "butt", lwd = 2, 
            lty = 5, col = "red", code = 3, angle = 90, length = 0.1)
        abline(h = 0, lty = 3, lwd = 2)
        legend(legendpos, legend = c("Normal", "Basic", "Percentile"), 
            lty = c(2, 4, 5), col = c("darkgreen", "blue", "red"), 
            lwd = 2)
    }
    if (articlestyle) {
        par(oldpar)
    }
}
