# Circle Chord Plot for Tree Representation of Brain Connectomes
# This code plots circle chord plot for the Tree Representations of 
# Brain Structural Connectivity via Persistent Homology by D. Li, P. Nguyen, 
# Z. Zhang and D. Dunson
#
# Author: Phuc Nguyen, Aug. 2020

# Level 4
chord_sslobe <- function(ribcol="green", gap.after=1, labeled=TRUE,
                         h.ratio=0.5, cex=0.6, A=NULL, corcol=NULL, track.col="lightgray") {
  midid <- (length(L5)/2)
  factors <- c(L5[1:midid], "blank1",
               L5[(midid + 1):length(L5)], "blank2")
  ribsize <- SSLOBE_RIBSIZE(gap.after = gap.after)
  xlim <- rbind(cbind(rep(0, midid), ribsize[1:midid]), 
                c(0,10),
                cbind(rep(0, midid), ribsize[(midid+1):length(L5)]),
                c(-1,10))
  circos.par(start.degree=80, gap.after=gap.after, cell.padding=c(0.02, 0, 0.02, 0))
  circos.initialize(factors=factors,
                    xlim=xlim)
  circos.trackPlotRegion(factors = factors, 
                         ylim = c(0, 1),
                         track.height = 0.02,
                         track.margin = c(0, 0.3),
                         cell.padding = rep(0.001, 4),
                         bg.border = track.col,
                         bg.col = track.col,
                         panel.fun = function(x, y) {
                           name <- get.cell.meta.data("sector.index")
                           xlim <- get.cell.meta.data("xlim")
                           ylim <- get.cell.meta.data("ylim")
                           if (length(grep("blank", name)) > 0) {
                             # Plot white rectangles to create space at top and bottom.
                             circos.rect(xleft=xlim[1], # removing -1 hear so the gray border extends to all chords
                                         ybottom=ylim[1]-1, 
                                         xright=xlim[2]+1, 
                                         ytop=ylim[1]+2, 
                                         col="white", border="white")
                           } else {
                             if (labeled) {
                               #text direction (dd) and adjusmtents (aa)
                               theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                               dd <- ifelse(theta < 90 || theta > 270, 
                                            "clockwise", "reverse.clockwise")
                               aa = c(1, 0.1)
                               if(theta < 90 || theta > 270)  aa = c(0, 0.1)
                               circos.text(x=mean(xlim), 
                                           y=ylim[2]+1, 
                                           labels=name,
                                           facing=dd, adj = aa, #text angle
                                           cex=cex # text size
                               )
                             }
                           }
                         })
  # Plotting only some of the top connections
  l4id <- unlist(sapply(L4, function(n) which(NODES == n)))
  Agiven <- !is.null(A)
  if (!Agiven) {
    # If a matrix A is not given, use the mean brain matrix
    # and plotting a few representative connectomes
    A <- get_chord_mat(NODES, MEAN_BRAIN, l4id)
    links <- get_top_links(A)
  } else {
    # When A is supplied, we want the value to be mean because its mean correlations
    # and plot only the significant connect ions that don't have zero correlations
    # with the traits
    Amax <- max(A)
    Amin <- min(A)
    A <- get_chord_mat(NODES, A, l4id, fun="mean")
    links <- get_sig_links(A)
  }
  
  if (nrow(links) > 0) {
    ribbonsize <- xlim[,2]
    names(ribbonsize) <- factors
    links[["point1"]] <- sapply(links$from, function(n) ribbonsize[[n]])
    links[["point2"]] <- sapply(links$to, function(n) ribbonsize[[n]])
    if (!Agiven) {
      if (length(ribcol) > 1) {
        parents12 <- sapply(links$from, function(n) get_parent(n))
        links[["col"]] <- sapply(parents12, function(n) ribcol[[n]])
      } else {
        links[["col"]] <- ribcol
      }
    } else {
      stopifnot(!is.null(corcol))
      ribcol <- cut(c(links$value, Amax, Amin), # max(links$value), -max(links$value),min(links$value), -min(links$value)),
                    breaks=N_CORCOL-1, labels = F)[1:nrow(links)]
      #links[["col"]] <- corcol[11+10*round(links$value, 1)]
      links[["col"]] <- corcol[ribcol]
    }
    # Make a black border
    # for (i in 1:nrow(links)) {
    #   circos.link(sector.index1=links$from[i], point1=c(-0.3,links$point1[i]+0.3),
    #               sector.index2=links$to[i], point2=c(-0.3,links$point2[i]+0.3),
    #               col="black",
    #               h.ratio=h.ratio)
    # }
    for (i in 1:nrow(links)) {
      if (links$col[i] != "blank") {
        circos.link(sector.index1=links$from[i], point1=c(0,links$point1[i]),
                    sector.index2=links$to[i], point2=c(0,links$point2[i]),
                    col=links$col[i],
                    h.ratio=h.ratio)
      }
    }
  }
  circos.clear()
}
