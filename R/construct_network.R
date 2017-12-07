
#producing visualisation of co-occurance analysis results.
#Three items to be plotted.
#1 . The network showing sub communities present in the community.
#2. The plot of Betweeness Vs Eigen Value centrality.
#3. plot of the four network statistal measures
#unrequired plots may be turned off via arguments.

construct_network <- function(edge_lists, comm_results, plotNetStats=F, plotNetwork=F,plotBetweennessEeigenvalue=F, scale.vertex.size=1, scale.edge.width=1, ...){
  
  colours <- microbiomeseq_cols()
  
  rho_cuts <- as.numeric(levels(comm_results$rho_cut)) 
  
  positive_rho_cuts <- rho_cuts[rho_cuts>0]
  
  g1 <- NULL
  
  for (i in unique(edge_lists$trt)){
  
      for (j in positive_rho_cuts){
    
        g1<-igraph::graph.edgelist(as.matrix(subset(edge_lists,trt==i)[,3:4]),directed=FALSE)
        
        E(g1)$color<-ifelse(subset(edge_lists,trt==i)[,"rho"]<0,"red","blue")
        
        E(g1)$width <- as.numeric(subset(edge_lists,trt==i)[,"rho"])*scale.edge.width
        
        V(g1)$size<-igraph::degree(g1)*scale.vertex.size
        
        #Give colours to subcommunities
        V(g1)$color<-as.character(sapply(V(g1)$name,function(x) ifelse(x %in% as.character(subset(comm_results,trt==i & rho_cut==j)[,"taxon"]),
                      colours[as.numeric(as.character(subset(comm_results,trt==i & rho_cut==j & taxon==x)[,"module"]))],"white")))
        
        V(g1)$module<-as.character(sapply(V(g1)$name,function(x) if(x %in% as.character(subset(comm_results,trt==i & rho_cut==j)[,"taxon"])){
                                                              return(subset(comm_results,trt==i & rho_cut==j & taxon==x)[,"module"])}))
        
        #plot network
        if(plotNetwork){
          pdf(paste("network_",paste(unique(edge_lists$trt),collapse="_vs_"),"_",i,"_",j,".pdf",sep=""),width=30,height=30)
          par(mai=c(0,0,1,0))       #this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
          plot(g1, layout=layout.fruchterman.reingold, main=i, vertex.label.dist=0.05,  vertex.frame.color='blue', 	
               vertex.label.color='black', vertex.label.font=2,vertex.label=V(g1)$name, ...)
          dev.off()
          
        }
        
        #plot betweenness - Eigen value plots
        if(plotBetweennessEeigenvalue){
          p <- plot_betweenness_eigenvalue(g1)
          pdf(paste("Betweeness_Eigenvalue_",paste(unique(edge_lists$trt),collapse="_vs_"),"_",i,"_",j,"_",".pdf",sep=""))
          print(p)
          dev.off()
        }
       
      }
  }
  out <- list(edgelists=edge_lists, comm_results=comm_results, graph=g1)
  return(out)
}


#produce a plot of betweenness Vs eigen value centrality.
plot_betweenness_eigenvalue <- function(g1){
  
  df<-data.frame(Betweenness=igraph::betweenness(g1),Eigenvalue=igraph::evcent(g1)$vector,Subcommunity=V(g1)$module, Taxa=V(g1)$name,Degree=igraph::degree(g1),Closeness=igraph::closeness(g1))
  df_community_mean<-aggregate(df[,c("Eigenvalue","Betweenness")],by=list(Subcommunity=df$Subcommunity),FUN="median")
  df_community_summary<-aggregate(df[,c("Eigenvalue","Betweenness")],by=list(Subcommunity=df$Subcommunity),FUN="summary")
  colnames(df_community_summary$Eigenvalue)<-paste("Eigenvalue",colnames(df_community_summary$Eigenvalue))
  colnames(df_community_summary$Betweenness)<-paste("Betweenness",colnames(df_community_summary$Betweenness))
  df_community_mean<-data.frame(df_community_mean,df_community_summary$Eigenvalue,df_community_summary$Betweenness)
  df_community_mean<-df_community_mean[df_community_mean$Subcommunity!="white",]
  colnames(df_community_mean)<-gsub("_$","",gsub("\\.","_",colnames(df_community_mean)))
  
  p<-ggplot(df)
  p<-p+geom_point(aes(x=Eigenvalue, y=Betweenness,size=Degree,fill=Subcommunity), colour="blue", pch=21,alpha=0.2)
  p<-p+scale_fill_manual(values=unique(get.vertex.attribute(g1, "color")))
  p<-p+theme_bw()
  p<-p+xlab("Eigenvalue Centrality")+ylab("Betweenness Centrality")
  p<-p+geom_point(aes(x=Eigenvalue,y=Betweenness,fill=Subcommunity),data=df_community_mean,colour="white",pch=23,size=7)

  if(dim(df_community_mean)[1]>0){
    p<-p+geom_errorbar(aes(x=Eigenvalue,y=Betweenness,ymin =Betweenness_1st_Qu,ymax = Betweenness_3rd_Qu,colour=Subcommunity),width=0.01,data=df_community_mean)
    p<-p+geom_errorbarh(aes(x=Eigenvalue,y=Betweenness,xmin =Eigenvalue_1st_Qu,xmax = Eigenvalue_3rd_Qu,colour=Subcommunity),height=0.01,data=df_community_mean)
    p<-p+scale_colour_manual(values=levels(df_community_mean$Subcommunity),guide = FALSE)
  }
  return(p)
}


