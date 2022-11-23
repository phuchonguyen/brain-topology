# Circle Chord Plot for Tree Representation of Brain Connectomes
# This code plots circle chord plot for the Tree Representations of 
# Brain Structural Connectivity via Persistent Homology by D. Li, P. Nguyen, 
# Z. Zhang and D. Dunson
#
# Author: Phuc Nguyen, Aug. 2020

# Level 2 has 5 lobes each hemisphere.
chord_lobe <- function(ribcol="blue", gap.after=3, labeled=TRUE,
                       h.ratio=0.5, A=NULL, corcol=NULL, img.labeled=TRUE,
                       track.col="lightgray"
                       ) {
  factors <- c(L3_BAC_FRO_RIG,
               "blank1",
               L3_FRO_BAC_LEF,
               "blank2"
  )
  # size of each regions based on # of sub-regions
  # xlim <- rbind(cbind(rep(0, 6), c(10, 6, 4, 10, 28, 42)), 
  #               c(0,10),
  #               cbind(rep(0, 6), rev(c(10, 6, 4, 10, 28, 42))),
  #               c(0,10))
  xlim <- rbind(cbind(rep(0, length(LOBE_RIBSIZE_BAC_FRO)), LOBE_RIBSIZE_BAC_FRO), 
                c(0,10),
                cbind(rep(0, length(LOBE_RIBSIZE_BAC_FRO)), rev(LOBE_RIBSIZE_BAC_FRO)),
                c(0,10))
  circos.par(start.degree=80, gap.after=gap.after)
  circos.initialize(factors=factors,
                    xlim=xlim)
  circos.trackPlotRegion(factors = factors, 
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
                                           y=ylim[2]+1, 
                                           labels=name,
                                           facing=dd, adj = aa, #text angle
                                           cex=0.6 # text size
                               ) 
                             }
                             if (img.labeled) {
                               xplot=get.cell.meta.data("xplot", sector.index = name)
                               x=.85*cos(aspace::as_radians((xplot[1]+xplot[2])/2))
                               y=.85*sin(aspace::as_radians((xplot[1]+xplot[2])/2))
                               if (name == "rh-occipitallobe") {
                                 y = 0.02353117
                                 x = x
                               }
                               readPNG(file.path("figures/paper/ggseg-dk/", 
                                                 REGION_IMG_MAP[[name]])) %>% 
                                 rasterImage(x-0.1, y-0.1, 
                                             x+0.13, y+0.05)
                             }
                           }
                         }
  )
  # Plotting only some of the top connections
  Agiven <- !is.null(A)
  if (!Agiven) {
    # If a matrix A is not given, a few representative connectomes
    #A <- get_chord_mat(NODES, MEAN_BRAIN, L2ID)
    links <- get_lobe_links()
  } else {
    # When A is supplied, we want the value to be mean because its mean correlations
    # and plot only the significant connections that don't have zero correlations
    # with the traits
    Amax <- max(A)
    Amin <- min(A)
    A <- get_chord_mat(NODES, A, L2ID, fun="mean")
    links <- get_sig_links(A)
  }
  
  if (nrow(links)==0) {
    circos.clear()
    return()
  }
  
  ribbonsize <- c(LOBE_RIBSIZE_BAC_FRO, 10, rev(LOBE_RIBSIZE_BAC_FRO), 10)
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
                  breaks = N_CORCOL-1, labels = F)[1:nrow(links)]
    #links[["col"]] <- corcol[11+10*round(links$value, 1)]
    links[["col"]] <- corcol[ribcol]
  }
  
  for (i in 1:nrow(links)) {
    if (links$col[i] != "blank") {
      circos.link(sector.index1=links$from[i], point1=c(0,links$point1[i]+1), # +1 to fix weird gap in hemis bundle
                  sector.index2=links$to[i], point2=c(0,links$point2[i]+1),
                  col=links$col[i],
                  h=0.4, h2=0.35,
                  # h.ratio=h.ratio
      )
    }
  }
  circos.clear()
}