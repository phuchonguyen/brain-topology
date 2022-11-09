chord_tree <- function(ribcol="blue", gap.after=1, labeled=TRUE, label.choice=c(F,F,F,T),
                       h.ratio=0.5, cex=0.6, A=NULL, tree=NULL, legend=TRUE,
                       title="", palette="redwhiteblue", track.col="lightgray",
                       legend.x=1, legend.y=0.3, img.labeled=FALSE) {
  if (palette =="redwhiteblue") {
    corcolfunc <- colorRampPalette(c("#0072B2",
                                     "#D3D3D3", #"#FFFFFF",
                                     "#D55E00"))
  } else if (palette =="viridis") {
    corcolfunc <- viridis::viridis
  } else stop("palette not implemented")
  corcol <- corcolfunc(N_CORCOL)
  if (img.labeled) {
    label.choice=c(F,F,F,F)
    track.col = rgb(255, 255, 255,alpha=0,maxColorValue=255) # clear white
  }
  
  chord_rhlh(ribcol=ribcol[1], 
             labeled = label.choice[1], A=A,
             corcol = corcol,
             track.col = track.col
             )
  par(new=TRUE)
  chord_lobe(ribcol=ribcol, labeled = label.choice[2],
             h.ratio=0.4, 
             A=A, 
             corcol = corcol,
             gap.after = 1.5,
             img.labeled = img.labeled)  # gap.after to get rid of weird gaps in hemis bundle
  par(new=TRUE)
  chord_sublobe(ribcol=ribcol, labeled = label.choice[3],
                h.ratio=0.65,
                A=A,
                corcol = corcol,
                gap.after = 1.5,
                track.col = track.col
                )
  par(new=TRUE)
  chord_sslobe(ribcol=ribcol, labeled = label.choice[4],
               h.ratio=0.8,
               A=A,
               corcol = corcol,
               gap.after = 1.5,
               track.col = track.col
               )
  if (legend) {
    if (!is.null(A)) {
      mx <- max(A)
      mn <- min(A)
    } else if (!is.null(tree)) {
      mx <- max(tree)
      mn <- min(tree)
    } else {
      mx <- max(MEAN_BRAIN)
      mn <- min(MEAN_BRAIN)
    }
    if (palette == "redwhiteblue") {
      leg <- c(formatC(max(abs(mn), abs(mx)), digits=1, format="e"), 
                  rep("", N_CORCOL-2),
                  formatC(-max(abs(mn), abs(mx)), digits=1, format="e"))
    } else if (palette =="viridis") {
      leg <- c(round(mx), rep("", nblank-1), round((mx+mn)/2), rep("", nblank), round(mn))
      #leg <- round(seq(round(mx), round(mn), length.out = 6))  # Hacky for figure 2
    }
    
    legend(x=legend.x, 
           y=legend.y, 
           #title="correlation",
           legend = leg, #c(-mx, -mx/2, 0, mx/2, mx),#seq(-1,1,0.5), 
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
