# Plot saved blocked-sampler results; no resampling. Use up to 500 evenly spaced
# retained draws for the displayed surfaces, all draws for parameter summaries.
plot_blocked <- function(case,out) {
  d <- unclass(posterior::as_draws_matrix(case$draws))
  take <- unique(round(seq(1,nrow(d),length.out=min(500,nrow(d)))))
  ds <- d[take,,drop=FALSE]
  transform <- fourier_operator(case$settings$n,case$settings$m,
    round(case$grid$u*case$settings$n),2*pi*case$grid$f)
  hdraw <- sapply(seq_len(nrow(ds)),function(i) Mod(transform(
    ds[i,"A"]*chirp_time(case$settings$n,ds[i,c("f0","fdot")])) )^2)
  hs <- apply(hdraw,1,median)
  cs <- ds[,sprintf("c[%d]",1:ncol(case$B)),drop=FALSE]
  psd <- exp(apply(cs %*% t(case$B),2,median))
  u <- sort(unique(case$grid$u)); f <- sort(unique(case$grid$f))
  png(file.path(out,"joint.png"),1500,1100,res=150,bg="white")
  par(oma=c(1,0,0,0))
  layout(matrix(1:6,2,byrow=TRUE),widths=c(1,1,0.18))
  palette <- hcl.colors(100,"viridis")
  for(row in 1:2) {
    values <- if(row==1) list(Mod(case$h)^2,hs) else list(case$S,psd)
    limits <- range(unlist(values))
    for(j in 1:2) {
      par(mar=c(4,4,3,0.5))
      image(u,f,matrix(values[[j]],length(u)),col=palette,zlim=limits,
            xlab="Rescaled time u",ylab="Frequency [Hz]",useRaster=FALSE,
            main=c("Injected moving Fourier power","Posterior median signal power","True LS2 PSD","Posterior median PSD")[(row-1)*2+j])
    }
    par(mar=c(4,0.5,3,3))
    zz <- seq(limits[1],limits[2],length.out=100)
    image(1:2,zz,matrix(rep(zz,each=2),2),col=palette,axes=FALSE,xlab="",ylab="")
    axis(4,las=1); box()
  }
  dev.off()
  png(file.path(out,"signal_diagnostics.png"),1200,1000,res=150,bg="white")
  par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(2,0,0,0))
  colors <- c("#2378A8","#D17520","#407D50","#9A529A")
  labels <- c(A="A",f0="f0 [Hz]",fdot="fdot [Hz/s]")
  for(v in c("A","f0","fdot")) {
    x <- matrix(case$draws[,,v],ncol=4)
    ix <- unique(round(seq(1,nrow(x),length.out=min(2000,nrow(x)))))
    matplot(ix,x[ix,],type="l",lty=1:4,col=colors,xlab="Retained iteration",ylab=labels[v])
    abline(h=case$truth[v],lty=2)
    hist(d[,v],breaks=40,col="grey80",border="white",main=paste(v,"posterior"),xlab=labels[v])
    abline(v=case$truth[v],col="red3",lwd=2)
    abline(v=quantile(d[,v],c(.05,.95)),lty=3)
  }
  mtext(paste("LS2 + chirp; dashed/red = injection; dotted = 90% posterior interval. Gate:",case$ok),
        1,outer=TRUE,line=.4,cex=.8)
  dev.off()
  # MH has acceptance rates and chain diagnostics, not NUTS divergence/E-BFMI fields.
  png(file.path(out,"psd_diagnostics.png"),1200,900,res=150,bg="white")
  par(mfrow=c(2,2),mar=c(3,4,2,1),oma=c(2,0,0,0))
  for(v in c("phi[1]","phi[2]")) {
    x <- matrix(case$draws[,,v],ncol=4)
    ix <- unique(round(seq(1,nrow(x),length.out=min(2000,nrow(x)))))
    matplot(ix,x[ix,],type="l",lty=1:4,col=colors,xlab="Retained iteration",ylab=v)
    ac <- sapply(1:4,function(j) as.numeric(acf(x[,j],lag.max=100,plot=FALSE)$acf))
    matplot(0:100,ac,type="l",lty=1:4,col=colors,xlab="Lag",ylab="ACF",main=v)
  }
  mtext(sprintf("All parameters: max R-hat %.3f; min bulk/tail ESS %.0f",
                max(case$summary$rhat),min(case$summary$ess_bulk,case$summary$ess_tail)),1,outer=TRUE,line=.5)
  dev.off()
}
