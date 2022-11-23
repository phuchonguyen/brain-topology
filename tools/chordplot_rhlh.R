# Circle Chord Plot for Tree Representation of Brain Connectomes
# This code plots circle chord plot for the Tree Representations of 
# Brain Structural Connectivity via Persistent Homology by D. Li, P. Nguyen, 
# Z. Zhang and D. Dunson
#
# Author: Phuc Nguyen, Aug. 2020

# Level 1 has 2 regions: left and right hemisphere
chord_rhlh <- function(ribcol="blue", gap.after=3, labeled=TRUE, A=NULL, corcol=NULL, 
                       track.col="lightgray"
                       ) {
  circos.par(start.degree = 80)
  circos.initialize(factors = c("rh", "blanktop", "lh", "blankbot"), 
                    # Setting xlim to 101 so that the superimposed plot gray boarder covers all chords
                    # If plotting rh-lh alone, change to c(0, 100)
                    xlim = rbind(c(0,100), c(0,10), c(0,101), c(0,10)))
  circos.trackPlotRegion(factors = c("rh", "blanktop", "lh", "blankbot"), 
                         ylim = c(0, 1),
                         track.height = 0.02,
                         track.margin = c(0, 0.3),
                         cell.padding = rep(0.002, 4),
                         bg.border = track.col,
                         bg.col = track.col,
                         panel.fun = function(x, y) {
                           name <- get.cell.meta.data("sector.index")
                           xlim <- get.cell.meta.data("xlim")
                           ylim <- get.cell.meta.data("ylim")
                           if (length(grep("blank", name)) > 0) {
                             # Plot white rectangles to create space at top and bottom.
                             circos.rect(xleft=xlim[1]-1, 
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
                                           y=ylim[2]+5, 
                                           labels=name,
                                           facing=dd, adj=aa, #text angle
                                           cex=1 # text size
                               )}
                           }
                         }
  )
  if (is.null(A)) {
    if (ribcol != "blank") {
      circos.link(sector.index1="rh", point1=c(0,100),
                  sector.index2="lh", point2=c(0,100),
                  col=ribcol) 
    }
  } else {
    value <- A["lh-insula", "rh-insula"]
    if (is.null(value)) stop("rh-insula not in A")
    if (value > 0) {
      stopifnot(!is.null(corcol))
      Acol <- A
      Acol[] <- cut(c(as.vector(A)), #, max(A), -max(A), min(A), -min(A)),
                    breaks=N_CORCOL-1, labels = F)[1:length(A)]
      ribcol <- corcol[Acol["lh-insula", "rh-insula"]]
      #ribcol <- corcol[11+10*round(value,1)]
      # circos.link(sector.index1="lh", point1=c(-0.3,100.3),
      #             sector.index2="rh", point2=c(-0.3,100.3),
      #             col="black")
      circos.link(sector.index1="rh", point1=c(0,100),
                  sector.index2="lh", point2=c(0,100),
                  col=ribcol) 
    }
  }
  circos.clear()
}