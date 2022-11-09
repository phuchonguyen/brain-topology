chord_adjmat <- function(gap.after=1.5, labeled=TRUE,
                         h.ratio=0.5, cex=0.6, A=NULL, topn=50,
                         palette="redwhiteblue", title="", legend=TRUE,
                         legend.x=1, legend.y=0.3, img.labeled=FALSE,
                         lwd=1) {
  # Color palette
  if (palette =="redwhiteblue") {
    corcolfunc <- colorRampPalette(c("#0072B2",
                                     "#D3D3D3", #"#FFFFFF",
                                     "#D55E00"))
  } else if (palette =="viridis") {
    corcolfunc <- viridis::viridis
  } else if (palette == "turbo") {
    corcolfunc <- viridis::turbo
  } else stop("palette not implemented")
  corcol <- corcolfunc(N_CORCOL)
  
  # Chord plot
  # Plot the ggseg images
  ggseg_factors <- c(L3_BAC_FRO_RIG,
                     "blank1",
                     L3_FRO_BAC_LEF,
                     "blank2"
  )
  ggseg_xlim <- rbind(cbind(rep(0, length(LOBE_RIBSIZE_BAC_FRO)), LOBE_RIBSIZE_BAC_FRO), 
                      c(0,10),
                      cbind(rep(0, length(LOBE_RIBSIZE_BAC_FRO)), rev(LOBE_RIBSIZE_BAC_FRO)),
                      c(0,10))
  circos.par("canvas.xlim" = c(-1, 1), 
             "canvas.ylim" = c(-1, 1),
             start.degree=80, gap.after=3)
  circos.initialize(factors=ggseg_factors,
                    xlim=ggseg_xlim)
  circos.trackPlotRegion(factors = ggseg_factors, 
                         ylim = c(0, 1),
                         track.height = 0.02,
                         track.margin = c(0, 0.3),
                         cell.padding = rep(0.002, 4),
                         bg.border = "white",
                         bg.col = "white",
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
  circos.clear()
  par(new = TRUE)
  # Plot the tracks
  midid <- (length(L6)/2)
  factors <- c(L6[1:midid], "blank1",
               L6[(midid + 1):length(L6)], "blank2")
  ribsize <- LEAVES_RIBSIZE()
  xlim <- rbind(cbind(rep(0, midid), ribsize[1:midid]), 
                c(0,10),
                cbind(rep(0, midid), ribsize[(midid+1):length(L6)]),
                c(0,10))
  circos.par(
    "canvas.xlim" = c(-1, 1), 
    "canvas.ylim" = c(-1, 1),
    start.degree=80, gap.after=1, cell.padding=c(0.002, 0, 0.002, 0))
  circos.initialize(factors=factors,
                    xlim=xlim)
  circos.trackPlotRegion(factors = factors, 
                         ylim = c(0, 1),
                         track.height = 0.02,
                         track.margin = c(0, 0.3),
                         cell.padding = rep(0.001, 4),
                         bg.border = "lightgray",
                         bg.col = "lightgray",
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
                             if (labeled & !img.labeled) {
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
  Agiven <- !is.null(A)
  if (Agiven) {
    links <- get_sig_links(A, topn=topn)
  } else {
    A <- MEAN_BRAIN
    links <- get_top_links(A)
  }
  
  if (nrow(links)==0) {
    circos.clear()
    return()
  }
  
  ribbonsize <- xlim[,2]
  names(ribbonsize) <- factors
  links[["point1"]] <- sapply(links$from, function(n) floor(ribbonsize[[n]]/2))
  links[["point2"]] <- sapply(links$to, function(n) floor(ribbonsize[[n]]/2))
  
  if (Agiven) {
    stopifnot(!is.null(corcol))
    ribcol <- cut(c(links$value, max(A), min(A)), # max(links$value), -max(links$value),min(links$value), -min(links$value)), 
                  breaks=N_CORCOL-1, labels = F)[1:nrow(links)]
    #links[["col"]] <- corcol[11+10*round(links$value, 1)]
    links[["col"]] <- corcol[ribcol]
  } else {
    links[["col"]] <- links$value
  }
  
  for (i in 1:nrow(links)) {
    from <- links$from[i]
    to <- links$to[i]
    fromid <- match(from, factors)
    toid <- match(to, factors)
    h.ratio <- 0.5 + max(0, (34-abs(fromid-toid))*0.3/34)  # So links between nearby nodes are shorter
    circos.link(sector.index1=from, point1=links$point1[i],
                sector.index2=to, point2=links$point2[i],
                col=links$col[i],
                # lwd=abs(scale(links$value))*lwd, # make linewidth proportional to correlation
                h.ratio=h.ratio)
  }
  circos.clear()
  if (legend) {
    legend(x=legend.x, y=legend.y, 
           #title="correlation",
           legend = c(formatC(max(A), format = "e", digits = 2), 
                      rep("", N_CORCOL-2), 
                      formatC(min(A), format = "e", digits = 2)), #c(-mx, -mx/2, 0, mx/2, mx),#seq(-1,1,0.5), 
           col = rev(corcolfunc(N_CORCOL)), 
           pch = 19, 
           bty = "n", 
           pt.cex = 2,  # size of points 
           cex = 0.75,  # size of text
           text.col = "black", 
           y.intersp = 0.2,
           horiz = F)
  }
  title(main=title)
}