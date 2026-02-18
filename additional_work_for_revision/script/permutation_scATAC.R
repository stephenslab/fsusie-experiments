
library(mvf.susie.alpha)
library(fsusieR)

set.seed(1)
path= paste0("/project2/mstephens/wdenault/anjing_data/data/" )
lf=list.files(path)
tt0=list()
tt1= list()
tt2=list()

tt3=list()
tt4=list()
tt_param=list()
filter_ncell =0

ncells <- t(read.delim(paste0(path,
                              "kellis_snatac.combined_genome_ncell.tss6.filtered_84samples.6celltype.txt"),
                       row.names=1))
ncells
rownames(ncells) <- gsub("\\.", "-", rownames(ncells))


all_count=t(read.delim(paste0(path,"updated.kellis_snatac.combined_genome_reads.tss6.filtered_84samples.6celltype.txt"),
                       row.names=1))
rownames(all_count) <- gsub("\\.", "-", rownames(all_count))


counts <- readRDS(paste0( path,
                          "xiong_atact_seq_multi_trait.169regions.51.2kb_expanded.binned_coverage_count.reordered_1mb_maf005.rds")
)


o =1
for ( null_w in c( 0.1, 0.5))
{
  for ( thres_c in  c(0.01,0.05,0.1, 0.5)){


    for (k in 1:10){

      reg_num=sample(1: length(counts))[1]
      tcount= counts[[reg_num]]

      X= counts[[reg_num]]$filtered_geno

      str(tcount)
      scaling_list=list()
      Y_f=list()
      k=1
      par(mfrow=c(3,2))



      for ( k in 1:(length(tcount)-1)){


        tt= tcount[[k]]$pheno
        common_rows <- intersect(rownames(tt), rownames(ncells))

        # Reorder `ncells` to match the order of `tt`
        ncells_ordered <- ncells[common_rows, , drop = FALSE]

        common_rows <- intersect(rownames(tt), rownames(all_count))

        # Reorder `ncells` to match the order of `tt`
        all_count_ordered <- all_count[common_rows, , drop = FALSE]

        average_total_read= all_count_ordered[,k]/ncells_ordered[,k] # average total number of read per cell in an individual

        scaling= average_total_read/mean(average_total_read, na.rm = TRUE)
        temp_Y_f =  tt*(1/ncells_ordered[,k] )
        scaling_list[[k]]=scaling
        Y_f[[k]]=log1p(temp_Y_f/scaling)

      }

      for ( i in 1: nrow(X)){
        for ( k in 1:length(Y_f)){

          if (  is.na(ncells[i,k]) )
            Y_f[[k]][i,]=NA
        }
      }

      if( filter_ncell>0){


        for ( i in 1: nrow(X)){
          for ( k in 1:length(Y_f)){

            if (ncells[i,k]<filter_ncell || is.na(ncells[i,k]) )
              Y_f[[k]][i,]=NA
          }
        }

      }
      pos=list(pos1=1:1024,
               pos2=1:1024,
               pos3=1:1024,
               pos4=1:1024,
               pos5=1:1024,
               pos6=1:1024
      )

      for ( k in 1:ncol(X)){
        if ( sum(is.na(X[,k]))>0){
          X[which(is.na(X[,k])),k]= median(X[-which(is.na(X[,k])),k])
        }
      }
      # 5% MAF filtering
      if( length(which (apply(X,2,mean)/2 <0.1))>0){
        X=X[,-which (apply(X,2,mean)/2 <0.1) ]
      }
      # there is actually a lot of flipped allele
      if( length(which (apply(X,2,mean)/2 >0.9))>0){
        X=X[,-which (apply(X,2,mean)/2 >0.9) ]
      }
      Y=list(Y_u=NULL,
             Y_f=Y_f)




      X_perm=as.matrix(   X[ sample(1:nrow(X)),])
      tt_m= multfsusie(Y=Y, X=  X_perm, pos=pos, L=3, cor_small = FALSE,
                       post_processing = "none" ,
                       quantile_trans = TRUE,
                       thresh_lowcount=thres_c)

      tt_msmall= multfsusie(Y=Y, X=  X_perm, pos=pos, L=3, cor_small = TRUE,
                            nullweight =null_w,
                            post_processing = "none" ,
                            quantile_trans = TRUE,
                            thresh_lowcount=thres_c)


      id=sample(1:6)[1]
      tt= fsusieR::susiF(Y=log1p(Y$Y_f[[id]]), X=X_perm, L=2,
                         nullweight =null_w,
                         post_processing = "none" ,
                         quantile_trans = TRUE,
                         thresh_lowcount=thres_c
      )


      tt_small= fsusieR::susiF(Y=log1p(Y$Y_f[[id]]), X=X_perm, L=2,
                               nullweight =null_w,
                               post_processing = "none" ,
                               quantile_trans = TRUE,
                               cor_small = TRUE,
                               thresh_lowcount=thres_c
      )
      lol=susieR::susie(y=log1p(apply(out$Y/out$size_factor_local,1, sum)), X=X_perm)
      tt0[[o]]=tt_small$cs
      tt1[[o]]=tt$cs

      tt2[[o]]=lol$sets
      tt3[[o]]=tt_m$cs
      tt4[[o]]=tt_msmall$cs
      tt_param[[o]]= c(null_w, thres_c)
      out= list(mfsusie=tt3 ,
                mfsusie_corsmall= tt4 ,

                fsusie=tt1,
                fsusie_corsmall=tt0,
                susie=tt2,
                param= tt_param)

      save(out , file="/home/wdenault/fsusi_simu/correlated_sim/permutation_scATAC.RData")


    }

  }
}
