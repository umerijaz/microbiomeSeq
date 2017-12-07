
#establishing pairwise co-occurence patterns.
co_occur_pairs<-function(dataset, rhos){
 
  final.results<-data.frame()

  trts<-as.vector(unique(dataset$trt))
 
  for(t in 1:length(trts)){
   
    dataset_trt<-subset(dataset, dataset$trt==trts[t])
    dataset_trt_no0<-subset(dataset_trt, dataset_trt$ab1 > 0 & dataset_trt$ab2 > 0)
    
    dataset_trt_no0$pairs<-paste(dataset_trt_no0$taxa1,dataset_trt_no0$taxa2)
    
    for(r in length(rhos)){
     
      if(rhos[r] < 0){temp<-subset(dataset_trt_no0, dataset_trt_no0$rho <= rhos[r])}
      if(rhos[r] > 0){temp<-subset(dataset_trt_no0, dataset_trt_no0$rho >= rhos[r])}
      if(dim(temp)[1]>1){
        
        temp.graph<-igraph::simplify(igraph::graph.edgelist(as.matrix(temp[,c(2,3)]),directed=FALSE))
        edge_list<-data.frame(igraph::get.edgelist(temp.graph,names=TRUE))
        
        edge_list$pairs<-paste(edge_list$X1,edge_list$X2)
        edge_list_pvals<-merge(edge_list,dataset_trt_no0,by="pairs",sort=FALSE  )
        
        edge_list_pvals$rho_cut<-rhos[r]
        edge_list_pvals$trt<-trts[t]
        
        edge_list_pvals$qval<-p.adjust(edge_list_pvals$p.value,method="fdr",n=length(edge_list_pvals$p.value))
        as.matrix(names(edge_list_pvals))
        
        final.results<-rbind(final.results,edge_list_pvals[,-c(2:3)])	}
    }
  }
  return(final.results)
}

comm_stats<-function(dataset, rhos){
  
  final.results<-data.frame()
  
  trts<-as.vector(unique(dataset$trt))
  
  for(t in 1:length(trts)){
    
    dataset_trt<-subset(dataset, trt==trts[t])
    
    for(r in 1:length(rhos)){
      
      if(rhos[r] < 0){temp<-subset(dataset_trt, rho <= rhos[r])}
      if(rhos[r] > 0){temp<-subset(dataset_trt, rho >= rhos[r])}
      if(dim(temp)[1]>1){
        
        temp.graph<-igraph::simplify(igraph::graph.edgelist(as.matrix(temp[,c(2,3)]),directed=FALSE))
        temp_comm<-igraph::edge.betweenness.community(temp.graph, directed=FALSE)
        
        member_data<-cbind(row.names(as.matrix(membership(temp_comm))),as.matrix(membership(temp_comm)))
        row.names(member_data)<-NULL
        
        
        rho_cut<-rep(rhos[r],dim(member_data)[1])
        trt<-rep(trts[t],dim(member_data)[1])
        stats<-cbind(trt,rho_cut,member_data)
        colnames(stats)<-c("trt","rho_cut","taxon","module")
        
        
        final.results<-rbind(final.results,stats)	}
    }
  }
  return(final.results)
}