biplot<-function(phy,color=NULL,shape=NULL,top=10,pointsize=5,alpha=0.7,taxa="Phylum",ellipse=FALSE,biplot=TRUE,show=TRUE){
  ###c("edgernorm", "varstab", "randomsubsample", "proportion", "relative", "log-relative", "scale")
  #plog<-normalise_data(phy,norm.method = "relative")
  plog=transform_sample_counts(phy,function(x)x/sum(x))
  #plog=transform_sample_counts(phy,function(x)log(x+1))
  pord<- ordinate(plog,  method = "MDS", distance = "bray")
  p<-plot_ordination(plog, pord, color = color,type="samples",shape=shape) +geom_point(size=pointsize,alpha=alpha,aes_string(shape=shape))+theme_light(base_size = 15)+
    scale_color_brewer(type="qual", palette="Set1")
  p1<-plot_ordination(plog, pord,type="taxa",shape=taxa) +geom_point(size=pointsize,alpha=alpha,aes_string(shape=shape))+theme_light(base_size = 15)+
    scale_color_brewer(type="qual", palette="Set2")
  pp<-p1$data
  ll=gsub('.*;','',gsub(';NA','',apply(pp[3:ncol(pp)],1,function(x)paste(x,collapse =";"))))
  if(show==FALSE){
    lx=rownames(pp)
  }else{
    lx<-paste(rownames(pp),ll,sep="\n")
  }
  pp$labels=lx
  pp$dist=p1$data[,1]^2+p1$data[,2]^2
  pp=pp[order(pp$dist,decreasing = T),]
  pp=pp[1:top,]
  arrowhead = arrow(length = unit(0.02, "npc"))
  p2<-p+geom_segment(aes(xend=1.3*Axis.1,yend=1.3*Axis.2,x=0,y=0),size=0.5,color="darkgray",arrow=arrowhead,data=pp)+ 
    geom_text_repel(aes(x=1.3*Axis.1,y=1.3*Axis.2,label=labels),color="black",data=pp,show.legend = FALSE)+scale_color_brewer(type="qual", palette="Set1")
  if(ellipse==TRUE){
    p2<-p2+stat_ellipse()
  }
  if(biplot==FALSE){
    p3<-p
  }else{
    p3<-p2
  }
  p3$layers<-p3$layers[-1]
  p3
}
