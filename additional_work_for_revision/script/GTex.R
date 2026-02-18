rm(list = ls())
library(fsusieR)
#load("C:/Document/Serieux/Travail/Data_analysis_and_papers/fsusie-experiments/additional_work_for_revision/data_analysis/GTeX/result_GaoPTK2B.RData")
plot(out$res)
#plot(out$res$lBF[[5]])


  
Y=out$Y 
X=out$X 
base_pair=out$pos
 print(dim(Y))
plot( base_pair, Y[1,], main=" Observed base pair count for individual 1")

 print( dim( Y))
pos = out$pos
print(pos[1:10])
bin_size=pos[2]-pos[1]
print( bin_size)
plot(pos, Y[1,], main=" Observed binned count for individual 1")
  size_factor=rowSums(Y)/ mean(rowSums(Y))
Y_cor =   (  Y /as.vector(size_factor))

plot(pos, log1p(Y[1,] ))
 k=3
x=X[,out$res$cs[[k]][1]]

 uni_res= univariate_functional_regression(Y=log1p(Y_cor),
                                          X=as.matrix(X[,out$res$cs[[k]][1]],
                                                      ncol=1),
                                          method="smash",
                                          alpha = 0.1
                                          )
plot(pos, uni_res$effect_estimate, 
     ylab = "estimated effect",
     xlab="base pair",
     type="l",
     lwd=2,
     col="royalblue",
     ylim=c (min(uni_res$cred_band),max(uni_res$cred_band)))
lines( pos,uni_res$cred_band[1,], 
       lty=2,
       lwd=1.4,
       col="royalblue")
lines(pos, uni_res$cred_band[2,], 
        lty=2,
       lwd=1.4,
       col="royalblue")
abline(h=0)
 
 
library(AnnotationHub)
library(org.Hs.eg.db)
library(GenomicRanges)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
pos0 = pos[1]
pos1= pos[length(pos)]+bin_size
positions=pos
chr=out$chr
 

effect=effect=rbind(uni_res$effect_estimate,
                    uni_res$cred_band,
                    rep(0,length(uni_res$effect_estimate)))
group_colors <- c("black" ,"royalblue","royalblue","royalblue" )
group_lwd= c(1,2,1,1)

read_counts <-  rbind( colMeans(log1p(Y_cor[x == 1,])),
                       colMeans(log1p(Y_cor[x == 0,])) 
)
obs_effect=rbind (read_counts[2,]-read_counts[1,],
                  rep(0, length(read_counts[2,] )) ) 

ylim= c( min( c(effect,obs_effect)), max(c(  effect,obs_effect)))

group_cred= c(1:3,0)  
effect_track <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            ylim=ylim,
            data = effect,
            genome = "hg38",
            groups= group_cred,
            name = "Effect SNP 261",
            type = "l",   
            lwd= group_lwd,  
            col = group_colors,
            track.margin = 0.05,
            cex.title = 0.6,
            cex.axis = 0.6,
            col.axis = "black",
            col.title = "black",
            fontface = "plain",
            background.title = "white",
            fontface.title = 1, 
            legend=FALSE)
 

#### Create the gene track
 
# Create a "gene region" track.
gene_track <- GeneRegionTrack(txdb,
                              genome = "hg38",
                              chromosome = chr,
                              pos0 = pos0,
                              pos1 = pos1,
                              name = "",
                              showId = TRUE,
                              geneSymbol = TRUE,
                              col.axis = "black",
                              col.title = "black",
                              transcriptAnnotation = "symbol",
                              rotation.title = 0,
                              cex.title = 0.6,
                              col = "salmon",
                              fill = "salmon",
                              background.title = "white")
 

#### Create the observed data point track
 
read_counts <-  rbind( colMeans(log1p(Y_cor[x == 1,])),
                       colMeans(log1p(Y_cor[x == 0,])) 
)

groups <- c("SNP   = 0",
            "SNP   = 1" )
geno_colors <- c("turquoise", "navyblue" )
lab_y = "avg. log1p count" 
data_track <- DataTrack(range = GRanges(seqnames = chr,
                                        ranges = IRanges(start = positions,
                                                         end = positions + 1)),
                        data = read_counts,
                        genome = "hg38",
                        groups = groups,
                        name = lab_y  , 
                        type = "p" , #"p",#type = "l",
                        col = geno_colors  ,
                        track.margin = 0.05,
                        cex.title = 0.6,
                        cex.axis = 0.6,
                        col.axis = "black",
                        col.title = "black",
                        fontface = "plain",
                        background.title = "white",
                        fontface.title = 1,
                        cex.legend = 0.6,
                        cex=0.2)
 obs_diff_track <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            data =  read_counts[1,]- read_counts[2,],
            genome = "hg38",
            ylim=ylim,             
            name = "observed difference \n between genotype",
            type = "p", 
            lwd= 2,
            col = "blue",
            track.margin = 0.05,
            cex.title = 0.6,
            cex.axis = 0.6,
            col.axis = "black",
            col.title = "black",
            fontface = "plain",
            background.title = "white",
            fontface.title = 1, 
            legend=FALSE)
obs_diff_track2 <-
  DataTrack(range = GRanges(seqnames = chr,
                            ranges = IRanges(start = positions,
                                             end = positions + 1)),
            data =  rep(0, 1024), 
            genome = "hg38",
            ylim= ylim, 
            name ="Observed difference ",
            type = c("l"  ),
            col = c("black" ),
            track.margin = 0.05,
            cex.title =  0.6,
            cex.axis =  0.6,
            col.axis = "black",
            col.title = "black",
            fontface = "plain",
            background.title = "white",
            fontface.title = 1,
            legend = FALSE )
obs_diff_track =OverlayTrack(trackList = list(obs_diff_track,
                                              obs_diff_track2),background.title = "white")


  
tracks <- c( 
  effect_track,
  obs_diff_track,
  data_track,
  gene_track 
)

plotTracks(tracks, from = pos0, to = pos1, sizes = c(  3,3, 3, 1 ))
 
   


print(pos[1:10])
bin_size= 
  print( bin_size)
plot(pos, Y[1,], main=" Observed binned count for individual 1")
size_factor=rowSums(Y)/ mean(rowSums(Y))
Y_cor =   (  Y /as.vector(size_factor))

plot(pos, log1p(Y[1,] ))

uni_res= univariate_functional_regression(Y=log1p(Y_cor),
                                          X=as.matrix(X[,261],
                                                      ncol=1),
                                          method="smash",
                                          alpha = 0.1
)



flr_nxT <- function(Y_nxT, X, t = NULL, k = 25) {
  stopifnot(
    is.matrix(Y_nxT),
    length(X) == nrow(Y_nxT)
  )
  
  n <- nrow(Y_nxT)
  T <- ncol(Y_nxT)
  
  if (is.null(t)) {
    t <- seq_len(T)
  }
  
  stopifnot(length(t) == T)
  
  library(mgcv)
  
  df <- data.frame(
    y  = as.vector(t(Y_nxT)),
    t  = rep(t, each = n),
    X  = rep(X, times = T),
    id = rep(seq_len(n), times = T)
  )
  
  fit <- bam(
    y ~
      s(t, k = k) +
      s(t, by = X, k = k),
    data = df,
    discrete = TRUE
  )
  
  return(fit)
}



fit <- flr_nxT(
  Y_nxT = log1p(Y_cor),     # n x T matrix
  X     = X[,123] # length n
)
summary(fit)$s.table[2,4]

pv=c()

for ( k in 1:ncol(X)){
  fit <- flr_nxT(
    Y_nxT = log1p(Y_cor),     # n x T matrix
    X     = X[,k] # length n
  )
  pv=c(pv, summary(fit)$s.table[2,4])
  
  print(k)
}





plot(-log10(pv+1e-12))
image(cor(X))




X_perm= X[sample(1:nrow(X)),]
pv0=c()

for ( k in 1:ncol(X)){
  fit <- flr_nxT(
    Y_nxT = log1p(Y_cor),     # n x T matrix
    X     = X_perm[,k] # length n
  )
  pv0=c(pv0, summary(fit)$s.table[2,4])
  
  print(k)
}





plot(-log10(pv0+1e-12))
image(cor(X))





pvlm =c()
pvlm0 =c()

for ( k in 1:ncol(X)){
  for ( j in 1: ncol(Y)){
    fit <-  lm(Y[,j]~X[,k ]) # n x T matrix
    pvlm=c(pvlm, summary(fit)$coefficients[2,4])
    fit <-  lm(Y[,j]~X_perm[,k ]) # n x T matrix
    pvlm0=c(pvlm0, summary(fit)$coefficients[2,4])
 
   
  }
 
 print(k)
 
}





plot(-log10(pvlm+1e-12))

plot(-log10( pvlm0+1e-12))



image(cor(X))




X_perm= X[sample(1:nrow(X)),]
pv0=c()

for ( k in 1:ncol(X)){
  fit <- flr_nxT(
    Y_nxT = log1p(Y_cor),     # n x T matrix
    X     = X_perm[,k] # length n
  )
  pv0=c(pv0, summary(fit)$s.table[2,4])
  
  print(k)
}




