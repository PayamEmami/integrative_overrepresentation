

#' Extracts KEGG database entries given user define data
#'
#' This functions extracts KEGG database. The main use is to keep track of number of genes and compounds in the KEGG data.
#' 
#'
#' @param references A named list of data.frames. See details
#' @param pathway_column Which column in the data.frames shows the pathway ids.
#' @param hit_column Which column in the data.frames shows the hits.
#' @param name_column Which column in the data.frames shows the name of the compounds.
#' @param ncores Number of cores to use. See details.
#' @examples
#'metabolomics_data <- data.frame(name=c("metabo1","metabo2","metabo3"),
#'                                pathway=c("map00930,map01100,map01120,map01220","map01502","map01502,map02010,map04978,map0523"),
#'                                hit=c(FALSE,FALSE,TRUE))
#'
#'transcript_data <- data.frame(name=c("gene1","gene2","gene3","gene4"),
#'                              pathway=c("map00930,map01100","map01502,map00930","map01502,map02010,map04978","map01502,map02010"),
#'                              hit=c(FALSE,FALSE,TRUE,TRUE))
#'
#'mirna_data <- data.frame(name=c("m1","m2","m3","m4"),
#'                         pathway=c("map00930","map01502,map00930","map01502,map02010","map01502,map02010"),
#'                         hit=c(FALSE,FALSE,TRUE,TRUE))
#'references <- list(transcripts=transcript_data,metabolites=metabolomics_data,mirna=mirna_data)
#'database <- create_database(references)
#'
#' @return
#' A list of pathways with the entries information.
#' 
#' @details 
#' When number of cores is too high, KEGG database might block the connection. It is better to use a reasonable number. 
#' Some number from 1 to 5 should generally be ok. If you realize there are too many warnings, decrease the number.
#' In practice some of pathway maps will not match. So you need to check them to make sure the warning is not because of missing pathway
#' 
#' @import KEGGREST
#' @import future.apply
#' @import foreach
#' @import doFuture
#' 
#' 
create_database <- function(references,
                pattern="map\\d+",
                pathway_column="pathway",
                hit_column = "hit",
                name_column="name",
                ncores=1){
  
  
pre_ids<-unique(as.vector(unlist(  sapply(references,function(omics){
  omics[,pathway_column]
}))))
map_ids<-unique(unlist(regmatches(pre_ids, gregexpr(pattern, pre_ids))))
library(KEGGREST)
library(future.apply)
library(foreach)
library(doFuture)
registerDoFuture()
plan(tweak(multisession,workers=ncores))
database <- foreach(pathway_id=map_ids)%dopar%
  {
    # Fetch details of the pathway
    
    pathway<-list()
    result = tryCatch({
      pathway_detail <- keggGet(pathway_id)
      ko_id <- pathway_detail[[1]]$KO_PATHWAY
      pathway_detail <- keggGet(ko_id)
      pathway$KOID <- ko_id
      pathway$MAPID <- pathway_id
      pathway$NAME <- pathway_detail[[1]]$NAME
      pathway$DESCRIPTION <- paste(pathway_detail[[1]]$DESCRIPTION,collapse = " ")
      pathway$ORTHOLOGY_LENGTH <- length(names(pathway_detail[[1]]$ORTHOLOGY))
      pathway$ORTHOLOGY <- names(pathway_detail[[1]]$ORTHOLOGY)
      pathway$COMPOUND_LENGTH <- length(names(pathway_detail[[1]]$COMPOUND))
      pathway$COMPOUND <- names(pathway_detail[[1]]$COMPOUND)
      
    }, error = function(e) {
      
      cat("pathway with id:",pathway_id,"has given error!!", "The error is:",e$message,"\n")
    })
    
    pathway
  }

return(database)
}


#' Performs integrative over-representation analysis using a set of pre-defined KEGG pathway
#'
#' This functions performs integrative over-representation analysis using a set of pre-defined KEGG pathway
#' 
#'
#' @param references A named list of data.frames. See details
#' @param omics_type A vector of character strings, the same length as `references`. Each element of this vector can be either ORTHOLOGY or COMPOUND. See details.
#' @param database Is a special list that must be created using `create_database` function.
#' @param test_type Which test to perform. fisher or hypergeometric. default: hypergeometric
#' @param background_type Use the local database or global KEGG to perform the test. See details. Default: local
#' @param weight_background_type Use the local database or global KEGG to weighting different omics. See details. Default: global
#' @param merged_analysis If true, pathway weighting is not perform instead the features from all modalities will be merged. 
#' @param weight_type Can be either "unweighted","overall","pathway". See details. Default: overall
#' @param pattern Pattern for extracting KEGG IDs. Don't change unless you know what you are doing!
#' @param pathway_column Which column in the data.frames shows the pathway ids.
#' @param hit_column Which column in the data.frames shows the hits.
#' @param name_column Which column in the data.frames shows the name of the compounds.
#' @export
#' @examples
#'metabolomics_data <- data.frame(name=c("metabo1","metabo2","metabo3"),
#'                                pathway=c("map00930,map01100,map01120,map01220","map01502","map01502,map02010,map04978,map0523"),
#'                                hit=c(FALSE,FALSE,TRUE))
#'
#'transcript_data <- data.frame(name=c("gene1","gene2","gene3","gene4"),
#'                              pathway=c("map00930,map01100","map01502,map00930","map01502,map02010,map04978","map01502,map02010"),
#'                              hit=c(FALSE,FALSE,TRUE,TRUE))
#'
#'mirna_data <- data.frame(name=c("m1","m2","m3","m4"),
#'                         pathway=c("map00930","map01502,map00930","map01502,map02010","map01502,map02010"),
#'                         hit=c(FALSE,FALSE,TRUE,TRUE))
#'references <- list(transcripts=transcript_data,metabolites=metabolomics_data,mirna=mirna_data)

#'omics_type <-c("ORTHOLOGY","COMPOUND","ORTHOLOGY")
#'database <- create_database(references)
#'do_pathway_overrepresentation(references = references,omics_type = omics_type,database_input  = database)
#'
#' @return
#' A list with pathway information for each omics and also the overlapping pathways
#'
#' @details
#' Each data.frame must contain three columns pathway, hit and name. 
#' The name columns must be unique name of the feature (for example ID of the compound of gene name etc). 
#' The pathway column must be KEGG map id (for example map00010). Multiple pathways for each feature can be sperated by comma or space.
#' The hit column must be of logical type. The TRUE value in the hit column means this feature is selected and FALSE means it's not selected for example by t-test or similar.
#' So in summary, the "references" parameter is a list of such data.frames where each element of the list is a data modality (e.g omics).
#' All data.frames must contain all the quantified data, meaning it must including all significat and not significant data.
#' `omics_type` must be a vector of with the same length as `references`. That is basically indicating which dataset is of which type.
#' Each database can be either ORTHOLOGY or COMPOUND. ORTHOLOGY refers to genes and COMPOUND refers to metabolites.
#' 
#' `background_type` is used to set which background to use. If set to local, the background is calculated from the given references.
#' Otherwise the background will be total entries for ORTHOLOGY or COMPOUND in the KEGG database.
#' 
#' In the integrative analysis, the fused p-values are weighted by some value. This value can be equal weights or different weights depending on
#' how many entries are part of a particular pathway. Setting `weight_background_type` to global, use the KEGG compounds to compute the weights. Setting 
#' this parameter to local, will force the function to calculate the weights based on the `references`
#' 
#' `weight_type` can be either "unweighted","overall","pathway". Given pathway A for omics O1 and O2, the unweighted option gives each omics a weight of 0.5.
#' Using the `overall` AND `weight_background_type` set to global, the function will calculate the proprtion of compounds and genes in the databas and use them as weights.
#' The `pathway` option AND `weight_background_type` set to global will instruct the function to calculate these propotions for each pathways separately.
#' 
#' If `weight_background_type` set to `local`, the function will use the reference for weight estimation in this way:
#' If `weight_type` is `overall` then weight is the fraction of number of features in each omics relative to the total number of features in all the omics
#' If `weight_type` is `pathway` the weight is the fraction of mapped reference features to each pathway. 
#' 
#' In general the better option would be to set `weight_background_type` to `global`.
#' 
#' 
#' 
do_pathway_overrepresentation <-function(references,
                                         omics_type,
                                         database_input=database,
                                         test_type=c("hypergeometric","fisher")[1],
                                         background_type=c("local","global")[1],
                                         weight_background_type=c("local","global")[2],
                                         merged_analysis = FALSE,
                                         weight_type = c("unweighted","overall","pathway")[2],
                                         pattern="map\\d+",
                                         pathway_column="pathway",
                                         hit_column = "hit",
                                         name_column="name"

                                         )
{
  
  ## define required functions
  
  performWeightedZtest <- function(p, weights=c(1,1)){
    zp <- (qnorm(p, lower.tail = FALSE) %*% weights)/sqrt(sum(weights^2));
    res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE));
    res;
  }
  
  
  
  GetFisherPvalue <- function(numSigMembers, numSigAll, numMembers, numAllMembers){
    z <- cbind(numSigMembers, numSigAll-numSigMembers, numMembers-numSigMembers, numAllMembers-numMembers-numSigAll+numSigMembers);
    z <- lapply(split(z, 1:nrow(z)), matrix, ncol=2);
    z <- lapply(z, fisher.test, alternative = 'greater');
    p.values <- as.numeric(unlist(lapply(z, "[[", "p.value"), use.names=FALSE));
    return(p.values);
  }
  
  
  ## extract IDs
pathway_ids<-lapply(references,function(x){
  outp <- regmatches(x[,pathway_column], gregexpr(pattern, x[,pathway_column]))
  names(outp)<-x[,name_column]
  outp
})

pathway_hits<-lapply(references,function(x){
  
  outp <-regmatches(x[x[,hit_column]==TRUE,pathway_column], 
             gregexpr(pattern, x[x[,hit_column]==TRUE,pathway_column]))
  
  names(outp)<-x[x[,hit_column]==TRUE,name_column]
  outp

})

## create global database


overal_pathways_ids<-unlist(pathway_ids,recursive = F)
overal_pathway_hits<-unlist(pathway_hits,recursive = F)
overal_omics_type<-unlist(omics_type,recursive = F)
overal_references<-do.call(rbind,references)


list_of_ids <- unique(unlist(overal_pathways_ids))
tota_number_of_features <- nrow(overal_references)
y <-overal_pathway_hits
x<-overal_pathways_ids
overal_database <- lapply(list_of_ids,function(id)
{
  ## mapped features to each pathways
  total_number = sum(sapply(x,function(x1){any(id%in%x1)}))
  
  ## mapped and significant to each pathway
  total_hits = sum(sapply(y,function(y1){any(id%in%y1)}))
  
  ## find the current compound in the database
  hit_pathway_index <- 0
  for(db in 1:length(database_input))
  {
    
    
    if(length(database_input[[db]])>0 && (database_input[[db]]$MAPID==id || database_input[[db]]$KOID==id))
    {
      hit_pathway_index <- db
      break
    }
    
  }
  
  list(id=id,total_number = total_number,
       total_hits = total_hits,
       ## total number of features in the dataset
       tota_number_of_features = tota_number_of_features,
       ## total number of hits in each dataset
       total_number_of_hits = sum(overal_references[,hit_column]),
       ## total number of compounds
       compounds_in_kegg = ifelse(hit_pathway_index!=0,
                                  length(unique(database_input[[hit_pathway_index]]$COMPOUND)),-1),
       ## total number of genes in the database. 
       genes_in_kegg = ifelse(hit_pathway_index!=0,
                              length(unique(database_input[[hit_pathway_index]]$ORTHOLOGY)),-1)
  )
  
})




local_database <- mapply(function(x,y,z,ref) {

  list_of_ids <- unique(unlist(x))
  tota_number_of_features <- nrow(ref)
  out_list <- lapply(list_of_ids,function(id){

    ## mapped features to each pathways
    total_number = sum(sapply(x,function(x1){any(id%in%x1)}))
    
    ## mapped and significant to each pathway
    total_hits = sum(sapply(y,function(y1){any(id%in%y1)}))
    
    ## find the current compound in the database
    hit_pathway_index <- 0
    for(db in 1:length(database_input))
    {
      
 
    if(length(database_input[[db]])>0 && (database_input[[db]]$MAPID==id || database_input[[db]]$KOID==id))
    {
      hit_pathway_index <- db
      break
    }
      
    }
    
    list(id=id,total_number = total_number,
         total_hits = total_hits,
         ## total number of features in the dataset
         tota_number_of_features = tota_number_of_features,
         ## total number of hits in each dataset
         total_number_of_hits = sum(ref[,hit_column]),
         ## total number of compounds
         compounds_in_kegg = ifelse(hit_pathway_index!=0,
                                    length(unique(database_input[[hit_pathway_index]]$COMPOUND)),-1),
         ## total number of genes in the database. 
         genes_in_kegg = ifelse(hit_pathway_index!=0,
                                length(unique(database_input[[hit_pathway_index]]$ORTHOLOGY)),-1)
         )
  })
  
}, x=pathway_ids, y=pathway_hits,z=as.list(omics_type),ref=references)








## count total number of unique genes in the database
total_orthology <- length(unique(unlist(sapply(database_input,function(x){x$ORTHOLOGY}))))
total_compound <- length(unique(unlist(sapply(database_input,function(x){x$COMPOUND}))))


overal_database <- list(overal=overal_database)

if(merged_analysis)
{
  
  enrchiment_analysis<-lapply(overal_database,function(omics){
    
    
    lapply(omics,function(current_path){
      hit.num <- current_path$total_hits
      q.size <- current_path$total_number_of_hits
      skip_test <- FALSE
      if(background_type=="local")
      {
        set.num <- current_path$total_number
        uniq.count <- current_path$tota_number_of_features
      }else{
      
          set.num <- current_path$genes_in_kegg + current_path$compounds_in_kegg
          uniq.count <- total_orthology + total_compound
          if(set.num<1)
          {
            skip_test <- TRUE
          }
 
        
      }
      
      res.mat<-matrix(NA, nrow=1, ncol=6);
      colnames(res.mat)<-c("Total", "Expected", "Hits", "Raw p","notes","pathway_id");
      
      res.mat[,1]<-set.num;
      res.mat[,2]<-q.size*(set.num/uniq.count);
      res.mat[,3]<-hit.num;
      
      if(!skip_test)
      {
        if(test_type=="hypergeometric"){
          res.mat[,4] <- phyper(hit.num-1, set.num, uniq.count-set.num, q.size, lower.tail=F);
        }else if(test_type == "fisher"){
          res.mat[,4] <- GetFisherPvalue(hit.num, q.size, set.num, uniq.count);
        }
      }else{
        res.mat[,5]<-"zero total number of features in the pathway"
      }
      
      res.mat[,6] <-current_path$id
      
      current_path$results <- res.mat
      current_path
      
    })
    
  })
  
  
  added_res<-lapply(enrchiment_analysis,function(omics){
    
    results_out<-do.call(rbind,lapply(omics,function(x){x$results}))
    results_out
  })
  

  
}else{
  enrchiment_analysis<-mapply(function(omics,database_type){
    
    
    lapply(omics,function(current_path){
      hit.num <- current_path$total_hits
      q.size <- current_path$total_number_of_hits
      skip_test <- FALSE
      if(background_type=="local")
      {
        set.num <- current_path$total_number
        uniq.count <- current_path$tota_number_of_features
      }else{
        if(database_type=="ORTHOLOGY")
        {
          set.num <- current_path$genes_in_kegg
          uniq.count <- total_orthology
          if(set.num<1)
          {
            skip_test <- TRUE
          }
          
        }else if(database_type=="COMPOUND"){
          
          set.num <- current_path$compounds_in_kegg
          uniq.count <- total_orthology
          if(set.num<1)
          {
            skip_test <-TRUE
          }
        }else{
          stop("Wrong omics_type!")
        }
        
      }
      
      res.mat<-matrix(NA, nrow=1, ncol=6);
      colnames(res.mat)<-c("Total", "Expected", "Hits", "Raw p","notes","pathway_id");
      
      res.mat[,1]<-set.num;
      res.mat[,2]<-q.size*(set.num/uniq.count);
      res.mat[,3]<-hit.num;
      
      if(!skip_test)
      {
        if(test_type=="hypergeometric"){
          res.mat[,4] <- phyper(hit.num-1, set.num, uniq.count-set.num, q.size, lower.tail=F);
        }else if(test_type == "fisher"){
          res.mat[,4] <- GetFisherPvalue(hit.num, q.size, set.num, uniq.count);
        }
      }else{
        res.mat[,5]<-"zero total number of features in the pathway"
      }
      res.mat[,6] <-current_path$id
      
      
      current_path$results <- res.mat
      current_path
      
    })
    
  },omics=local_database,database_type=as.list(omics_type))
  
  ### here we have to combine the p-values from multiple datasets
  
  #extract name of the pathways in all the datasets only if number of hits are higher than zero!
  
  list_of_ids_stats <- lapply(enrchiment_analysis,function(omics){
    
na.omit(    sapply(omics,function(pathway){
  if(pathway$total_hits>0)
    pathway$id
  else
    NA
}))
  })
  
  
  
  list_of_pathways_output <- c()
  # find intesection of different omics pathways
  overlapping_pathways <- Reduce(intersect, list_of_ids_stats)
  if(length(overlapping_pathways)>0)
  {
    
    # find index of intersection
    overlapping_pathways_index <- lapply(enrchiment_analysis,function(omics){
      match(overlapping_pathways,
      sapply(omics,function(pathway){
        pathway$id
      }))
      
    })
    
    for(ol in 1:length(overlapping_pathways))
    {
      to_be_merged <-mapply(function(omics,match_index){

        list(omics[[match_index[ol]]])
        
      },omics=enrchiment_analysis,match_index=overlapping_pathways_index)
      
      
      ## calculate weights
      w.vec<-c()
      if(weight_background_type=="local")
      {
        
        if(weight_type == "unweighted"){ # unweighted
          w.vec <- c(w.vec,0.5);
        }else if(weight_type == "overall"){ # overall
          total_count<-sum(sapply(pathway_ids,length))
          w.vec<-as.vector(sapply(pathway_ids,length)/total_count)
        }else{ # pathway level
          
          w.vec<-as.vector(unlist( mapply(function(omics,match_index){
           
           list(omics[[match_index[ol]]]$total_number)
           
         },omics=enrchiment_analysis,match_index=overlapping_pathways_index)))
          
          w.vec<-w.vec/sum(w.vec)
          
        }
       
       
        
      }else{
        
        total_orthology <- length(unique(unlist(sapply(database_input,function(x){x$ORTHOLOGY}))))
        total_compound <- length(unique(unlist(sapply(database_input,function(x){x$COMPOUND}))))
        total.count<-total_orthology+total_compound
        ow.m <- total_compound/total.count
        ow.g <- total_orthology/total.count
        
        
        # for pathway based weights
        cmpd.counts<-to_be_merged[[1]]$compounds_in_kegg
        gene.counts<-to_be_merged[[1]]$genes_in_kegg
        
        path.uniq.lens<- gene.counts+cmpd.counts
        pw.m <- cmpd.counts/path.uniq.lens;
        pw.g <- gene.counts/path.uniq.lens;
        
        w.vec<-c()
        
        for(otp in omics_type)
        {
          if(otp=="ORTHOLOGY")
          {
            if(weight_type == "unweighted"){ # unweighted
              w.vec <- c(w.vec,0.5);
            }else if(weight_type == "overall"){ # overall
              w.vec <- c(w.vec,ow.g);
            }else{ # pathway level
              w.vec <- c(w.vec,pw.g);
            }
            
          }
          else if(otp=="COMPOUND")
          {
            if(weight_type == "unweighted"){ # unweighted
              w.vec <- c(w.vec,0.5);
            }else if(weight_type == "overall"){ # overall
              w.vec <- c(w.vec,ow.g);
            }else{ # pathway level
              w.vec <- c(w.vec,pw.g);
            }
          }
          
        }
      }
      
      
      
      ## integrate p-values
      extraced_pvalues <- as.numeric(as.vector(sapply(to_be_merged, function(x){x$results[,"Raw p"]})))
      int_pvalue <- as.vector(performWeightedZtest(extraced_pvalues,w.vec)$p)
      
      ## merge dataframes into one
      ## empty data.frame
      tmp_output<-to_be_merged[[1]]$results
      tmp_output[]<-""
      
      merged_names<-names(to_be_merged)[1]
      tmp_output[]<-paste(merged_names,":",to_be_merged[[merged_names]]$results[1,], sep = "")
      for(merged_names in names(to_be_merged)[-1])
      {
        tmp_output[]<-paste(tmp_output,paste(merged_names,":",to_be_merged[[merged_names]]$results[1,], sep = ""),sep=",")
      }
      tmp_output[,"Raw p"]<-int_pvalue
      tmp_output[,"pathway_id"]<-overlapping_pathways[ol]
      
      list_of_pathways_output<-rbind(list_of_pathways_output,tmp_output)
    }
  
  }
  
  added_res<-lapply(enrchiment_analysis,function(omics){
    
    results_out<-do.call(rbind,lapply(omics,function(x){x$results[,,drop=FALSE]}))
  
    ## clean from overlaps
    
    results_out<-results_out[!results_out[,"pathway_id"]%in%overlapping_pathways,,drop=FALSE]
    results_out
    })

  added_res$overlapping<-list_of_pathways_output

}

## add the feature names

for(i in 1:length(added_res))
{

  ## add additional columns depending on the input
  keep_names <- colnames(added_res[[i]])
  for(x in names(references))
  {
    keep_names<-colnames(added_res[[i]])
    added_res[[i]]<-cbind(added_res[[i]],NA)
    colnames(added_res[[i]])<-c(keep_names,x)
    added_res[[i]][,x]
    for(j in 1:nrow(added_res[[i]]))
    {
      id <- as.character(added_res[[i]][j,"pathway_id"])
      
      selected <- sapply(regmatches(references[[x]][,pathway_column], gregexpr(pattern, references[[x]][,pathway_column])),function(xy){
        id%in%xy
      })
      if(any(selected))
      {
        added_res[[i]][j,x]<-paste(references[[x]][selected,name_column],collapse = ",")
      }
    
    }
    
  }
}

return(added_res)
}




