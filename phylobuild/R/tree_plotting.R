# Tree plotting and branch-length helpers used by the tree-building step.
# Copied verbatim from the farm my_functions/ library: plot_tree.R and plot_tree_examples.R
# (derived from Nick Williams' tree plotting code), phytools_scripts.R and
# Pop_size_estimation_functions.R. Only functions reachable from the tree-building step are kept.

# ---- from plot_tree.R ----
set_cedge=function(parent,tree){
  for(child in get_node_children(parent,tree)){
    #cat("\nsetting parent - child ",parent,child,"\n")
    child.idx=which(tree$edge[,2]==child)
    parent.idx=which(tree$edge[,2]==parent)
    if(length(parent.idx)==0){
      pedge=0
    }else{
      pedge=tree$cedge[parent.idx]
    }
    #cat("before:",length(which(tree$cedge>0)),"\n")
    tree$cedge[child.idx]=pedge+tree$edge.length[child.idx];
    #cat("after:",length(which(tree$cedge>0)),"\n")
    #cat(parent,child,":",tree$cedge)
    tree=set_cedge(child,tree)
  }
  tree
}

get_height=function(tree,node){
  #if(is.null(tree$height)){
  #  tree$height=rep(NA,1:length(tree$edge.length))
  #}
  N=length(tree$tip.label)
  if(node<=N){
    return(node)#tree$height[which(tree$edge[,2]==node)]=node
  }else{
    children=get_node_children(node,tree)
    return(mean(sapply(children,function(child) get_height(tree,child))))
  }
}

set_height=function(tree){
  tree$height_end=sapply(tree$edge[,2],function(i) get_height(tree,i))
  tree$height_start=tree$height[match(tree$edge[,1],tree$edge[,2])]
  N=length(tree$tip.label)
  root=N+1
  idx=which(tree$edge[,1]==root)
  tree$height_start[idx]=mean(tree$height_end[idx])
  tree
}

elbow=function(x0,x1,y0,y1,...){ 
  arrows(x0=x0,y0=y0,y1=y1,length=0,...)
  arrows(x0=x0,x1=x1,y0=y1,length=0,...)
}

elbowv=function(x0,x1,y0,y1,...){ 
  #browser()
  arrows(x0=x0,x1=x1,y0=y0,length=0,...)
  arrows(x0=x1,y0=y0,y1=y1,length=0,...)
  
}

set_tree_coords=function(atree){
  ##get the cumulative distance from the root.
  tree=atree
  tree$cedge=rep(0,length(tree$edge.length))
  N=length(tree$tip.label)
  root=N+1
  tree=set_cedge(root,tree)
  tt=set_height(tree)
  atree$coords=data.frame(a0=tt$cedge-tt$edge.length,a1=tt$cedge,
                          b0=tt$height_start,b1=tt$height_end,stringsAsFactors = FALSE)
  if(!is.null(atree$color)){
    atree$coords$color=atree$color
  }
  atree
}

plot_tree=function(tree,direction="down",cex.label=5,offset=0,plot_axis=T,title=NULL,b_do_not_plot=FALSE,lwd=1,bars=NULL,default_edge_color="darkgrey",ymax=NULL,cex.terminal.dots=0,vspace.reserve=0,cex.axis=1){
  par(mar=c(1, 1, 1, 3) + 0.1)
  #browser()
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords
  
  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ymax=max(coords$a1)*1.05
    }
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  if(is.null(bars)){
    ymin=0-ymax*0.05-vspace.reserve*ymax
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(ymin,ymax),xlab="",ylab="")
  }else{
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=max(coords$a1)
    ##browser()
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    if(cex.label>0){
    text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,pos = 1,srt=270,col=tipcol)
    }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=c("darkgrey", "blueviolet","deeppink")[Y_loss], cex=cex.terminal.dots,pch=15)
    }
  }
  tree$direction=direction
  tree$top=top
  #scale =10
  #(the small steps are for trees scaled to a small age, e.g. foetuses run with age 1)
  scales=c(0.01,0.02,0.05,0.1,0.2,0.5,2,5,10,100,200,500,1000,2000,5000)
  scale=scales[max(c(1,which(ymax/4>scales)))]
  #scale=scales[max(which(ymax/2>=scales))]
  #browser()
  cat("scale=",scale,"\n")
  if(plot_axis){
    axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2, cex.axis = cex.axis, lwd = 1, lwd.ticks = 1, col = "black") 
  }
  if(!is.null(title)) {
    text(x = 0,y=ymax,pos = 4,labels = title)
  }
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("Mutation number",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,border = "black",lwd=0.25,col="darkred")
    
  }
  tree$ymax=ymax
  tree$vspace.reserve=vspace.reserve
  tree
}

get_all_node_children=function(node,tree){
  children=tree$edge[which(tree$edge[,1]==node),2]
  offspring=children
  for(child in children){
    offspring=c(offspring,get_all_node_children(child,tree))
  }
  offspring
}

get_node_children=function(node,tree){
  tree$edge[which(tree$edge[,1]==node),2]
}

get_samples_in_clade=function(node,tree){
  if(node<=length(tree$tip.label)){
    return(tree$tip.label[node])
  }
  tree$tip.label[intersect(get_all_node_children(node,tree),1:length(tree$tip.label))]
}

get_y_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$top-tree$coords[idx,c("a1","a0")])
}

get_x_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$coords[idx,c("b1","b0")])
}

# ---- from plot_tree_examples.R ----
get_idx_for_node=function(details,node){
  which(details$node==node)
}

get_edge_info=function(tree,details,node){
  y=get_y_range(tree,node)
  x=get_x_range(tree,node)
  idx=get_idx_for_node(details,node)
  samples=get_samples_in_clade(node,tree)
  list(yb=y[1],yt=y[2],x=x[1],xm=x[2],idx.in.details=idx,samples=samples)
}

add_annotation=function(tree,details=NULL,matrices=NULL,annot_function,...){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(tree,details,matrices,tree$edge[i,2],...))
}

add_simple_labels=function(tree,##<< enhanced phylo returned from plot_tree
                           details,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                           matrices,##<< a list of matrices parallel to details with columns named by tip labels
                           node,##<< Node (see details)
                           query.field,##<< Name of column in details to query against
                           query.allowed.df,##<< Values of query field which should be annotated. data.frame value,col,pch columns.
                           label.field,##<< Name of column in details specifying the label text.
                           cex.label=1,
                           b.add.label=TRUE,
                           b.add.marker=TRUE,
                           ... ##<< paremeters for points (not color)
){
  info=get_edge_info(tree,details,node)
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>50){
    stop("too many variants to annotate")
  }
  query.value=details[[query.field]][idx]
  idx.match=match(query.value,query.allowed.df$value)
  cols=query.allowed.df$col[idx.match]
  pch=query.allowed.df$pch[idx.match]
  
  vlabels=details[[label.field]][idx]
  ## spread out
  N=length(idx)
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  if(N>0){
    yd=info$yt-info$yb
    if(N==1){
      y=0.5*(info$yb+info$yt)
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    if(b.add.marker){
      points(rep(info$x,N),y,col=cols,pch=pch,...)
    }
    if(b.add.label){
      text(rep(info$x,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}

add_vaf=function(tree,##<< enhanced phylo returned from plot_tree
                 details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                 matrices,##<< a list of matrices parallel to details with columns named by tip labels
                 node,
                 samples=NULL,
                 b.plot.bars=TRUE,
                 lwd.rect=1,
                 min.depth=1,
                 vc.field,
                 vc.df,
                 filter.on=NULL,
                 ...
){
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  #browser()
  ##Can do additional filter e.g. missense 
  #if(!is.null(filter.on)){
  #info$idx=info$idx,3)#info$idx.in.details[which(details$VC=="missense")]
  # }
  
  if(length(info$idx)==0){
    return(NULL)
  }
  if(b.plot.bars){
    plotF=plotBars
  }else{
    plotF=plotPie
  }
  if(is.null(samples)){
    samples=info$samples
  }
  cat(length(info$idx),"\n")
  if(length(samples)>1){
    if(length(info$idx)>1){
      df=data.frame(mtr=rowSums(matrices$mtr[info$idx,samples],na.rm = TRUE),
                    dep=rowSums(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }else{
      df=data.frame(mtr=sum(matrices$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }
  }else{
    df=data.frame(mtr=matrices$mtr[info$idx,samples],
                  dep=matrices$dep[info$idx,samples],stringsAsFactors = FALSE)
  }
  df=cbind(df,details[info$idx,])
  df=df[which(df$dep>=min.depth),]
  df$vaf=df$mtr/df$dep
  df=df[which(!is.na(df$vaf)),]
  N=dim(df)[1]
  if(N==0){
    return(df)
  }
  
  
  df=df[order(df$vaf),]
  yd=info$yt-info$yb
  
  
  if(N==1){
    y=0.5*(info$yb+info$yt)
    width=yd
  }else{
    y1=seq(info$yb,info$yt,length.out = N+2)
    #Don't use the ends..
    y=y1[2:(N+1)]
    width=y[2]-y[1]
  }
  
  if(!b.plot.bars){
    r=0.8  ##r>0.5 will cause likely overlap problems
  }else{
    r=0.4
  }
  #arrows(x0=info$x-w,x1=info$x-w+df$vaf*2*w,y0=y,lend=1,length=0,col="black")
  #arrows(x0=info$x-w+df$vaf*2*w,x1=info$x+w,y0=y,lend=1,length=0,col="grey")
  for(i in 1:N){
    vaf=min(df$vaf[i],0.999)
    if(is.na(vaf)){
      plotF(x=info$x,y = y[i],radius=r,col=c("lightgray","lightgray"),prop=c(0,1),border="lightgray",width=width)
    }else{
      plotF(x=info$x,y = y[i],radius = r,col=c("black","white"),prop = c(vaf,1-vaf),width = width)
    }
  }
  if( !b.plot.bars){
    
    return(df)
  }
  ##Now check to see if we need to highlight
  ##Test if VAF is significantly > 0.05 or significantly < 0.45
  ##Can also do a binomial test...
  MTR=sum(df$mtr)
  DEP=sum(df$dep)
  min.mean.vaf=0.45
  z1=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
  z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
  max.p.value=max(z1$p.value,z2$p.value)
  txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
  if(z1$p.value>0.05 & z2$p.value<0.05) {
    border.color="green"
  } else if(max.p.value<0.05){
    if(max.p.value<0.05/dim(tree$edge)[1]){
      border.color="red"
    }else{
      border.color="blue"
    }
  }else{
    border.color="darkgrey"
  }
  
  
  rect(xleft=info$x-r,xright=info$x+r,ybottom=y[1]-width/2,ytop=y[N]+width/2,border=border.color,lwd=lwd.rect)
  if(border.color!="darkgrey"){
    text(txt,x=info$x,y=y[1]+0.3*(y[N]-y[1]),col="black",cex=0.6)
  }
  arrows(x0=info$x,y0=info$yb,y1=info$yt,lwd=0.5,col="black",length=0,lend=2)
  
  #df
  
}

plotBars=function(x,y,radius,col,prop,border="black",width=1){
  #cat(prop,"\n")
  if(width<2){
    arrows(x0 = x-radius,y0=y,x1=x-radius+2*radius*prop[1],col="darkgrey",lend=2,length=0)
    arrows(x0 = x-radius+2*radius*prop[1],y0=y,x1=x-radius+2*radius,col=rgb(0.98,0.98,0.98),lend=2,length=0)
  }else{
    rect(xleft = x-radius,xright =x-radius+2*radius*prop[1],ybottom = y-width/2,ytop=y+width/2,border = NA,col="darkgrey")
    rect(xleft =  x-radius+2*radius*prop[1],xright =x+radius,ybottom = y-width/2,ytop=y+width/2,border = NA,col=rgb(0.98,0.98,0.98))
  }
  1
}

plotPie=function(x,y,radius,col,prop,llwd=0.5,border="black",width=NA){
  lims=par("usr")
  as=dev.size()
  asr=as[1]/as[2]
  yscale=asr*(lims[4]-lims[3])/(lims[2]-lims[1])
  prop=prop/sum(prop)
  cutpoint=c(0,cumsum(prop)*2*pi)
  
  N=2*pi/0.05
  n=ceiling(N*diff(cutpoint)/(2*pi))
  d=diff(cutpoint)/n
  if(length(prop)>1){
    for(i in 2:length(cutpoint)){
      polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              border=border,col=col[i-1],lwd=llwd)
    }
  }else{
    i=2
    polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            border=border,col=col[1],lwd=llwd)
  }
  yscale
}

plot_tree_vaf=function(tree,details,matrices,samples=NULL,b.plot.bars =TRUE,filter.on=NULL){
  res=add_annotation(tree,details,matrices,
                     function(tree,details,matrices,node){
                       add_vaf(tree,details,matrices,node,samples=samples,b.plot.bars = b.plot.bars,filter.on=filter.on)})
}

boxtext <- function(x, y, labels = NA, col.text = NULL, col.bg = NA, 
                    border.bg = NA, adj = NULL, pos = NULL, offset = 0.5, 
                    padding = c(0.5, 0.5), cex = 1, font = graphics::par('font')){
  
  ## The Character expansion factro to be used:
  theCex <- graphics::par('cex')*cex
  
  ## Is y provided:
  if (missing(y)) y <- x
  
  ## Recycle coords if necessary:    
  if (length(x) != length(y)){
    lx <- length(x)
    ly <- length(y)
    if (lx > ly){
      y <- rep(y, ceiling(lx/ly))[1:lx]           
    } else {
      x <- rep(x, ceiling(ly/lx))[1:ly]
    }       
  }
  
  ## Width and height of text
  textHeight <- graphics::strheight(labels, cex = theCex, font = font)
  textWidth <- graphics::strwidth(labels, cex = theCex, font = font)
  
  ## Width of one character:
  charWidth <- graphics::strwidth("e", cex = theCex, font = font)
  
  ## Is 'adj' of length 1 or 2?
  if (!is.null(adj)){
    if (length(adj == 1)){
      adj <- c(adj[1], 0.5)            
    }        
  } else {
    adj <- c(0.5, 0.5)
  }
  
  ## Is 'pos' specified?
  if (!is.null(pos)){
    if (pos == 1){
      adj <- c(0.5, 1)
      offsetVec <- c(0, -offset*charWidth)
    } else if (pos == 2){
      adj <- c(1, 0.5)
      offsetVec <- c(-offset*charWidth, 0)
    } else if (pos == 3){
      adj <- c(0.5, 0)
      offsetVec <- c(0, offset*charWidth)
    } else if (pos == 4){
      adj <- c(0, 0.5)
      offsetVec <- c(offset*charWidth, 0)
    } else {
      stop('Invalid argument pos')
    }       
  } else {
    offsetVec <- c(0, 0)
  }
  
  ## Padding for boxes:
  if (length(padding) == 1){
    padding <- c(padding[1], padding[1])
  }
  
  ## Midpoints for text:
  xMid <- x + (-adj[1] + 1/2)*textWidth + offsetVec[1]
  yMid <- y + (-adj[2] + 1/2)*textHeight + offsetVec[2]
  
  ## Draw rectangles:
  rectWidth <- textWidth + 2*padding[1]*charWidth
  rectHeight <- textHeight + 2*padding[2]*charWidth    
  graphics::rect(xleft = xMid - rectWidth/2, 
                 ybottom = yMid - rectHeight/2, 
                 xright = xMid + rectWidth/2, 
                 ytop = yMid + rectHeight/2,
                 col = col.bg, border = border.bg)
  
  ## Place the text:
  graphics::text(xMid, yMid, labels, col = col.text, cex = theCex, font = font, 
                 adj = c(0.5, 0.5))    
  
  ## Return value:
  if (length(xMid) == 1){
    invisible(c(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                yMid + rectHeight/2))
  } else {
    invisible(cbind(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                    yMid + rectHeight/2))
  }    
}

add_simple_labels_line=function(tree,##<< enhanced phylo returned from plot_tree
                           details,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                           matrices,##<< a list of matrices parallel to details with columns named by tip labels
                           node,##<< Node (see details)
                           query.field,##<< Name of column in details to query against
                           query.allowed.df,##<< Values of query field which should be annotated. data.frame value,col,pch columns.
                           label.field,##<< Name of column in details specifying the label text.
                           cex.label=2,
                           b.add.label=TRUE,
                           b.add.marker=TRUE,
                           lty=1,
                           lwd=3,
                           ... ##<< paremeters for points (not color)
){
  info=get_edge_info(tree,details,node)
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>1){
    print("Some branches have multiple variants, plotting only the first variant of each branch")
    query.value=details[[query.field]][idx]
    vlabels=paste0(details[[label.field]][idx],collapse="\n")
  } else {
    query.value=details[[query.field]][idx]
    idx.match=match(query.value,query.allowed.df$value)
    cols=query.allowed.df$col[idx.match]
    vlabels=details[[label.field]][idx]
  }

  ## spread out
  N=ifelse(length(idx)>0,1,0)
  ##Vertical offset so that labels sit slightly above the markers.
  if(N>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col="red",lend=1,lwd=lwd,lty=lty,...)
    if(b.add.label){
      #text(rep(info$x,N),y=info$yb+0.5*(info$yt-info$yb),labels = vlabels,pos = 2,offset = 0.25,cex=cex.label)
      boxtext(info$x-1,info$yb+runif(1,min=0.35*(info$yt-info$yb),max=0.65*(info$yt-info$yb)),col.bg="white",border.bg="black",padding = c(0.5, 5),labels = vlabels,pos=2,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}

plot_tree_labels=function(tree,details,
                          query.field="VC",
                          type="label",
                          query.allowed.df=data.frame(value=c("nonsense","frameshift"),
                                                      col=c("red","black"),pch=c(17,18)
                                                      ),
                          label.field="GENE",
                          cex.label=1,
                          lty=1,
                          lwd=1){
  if(type=="label") {res=add_annotation(tree,
                     details,list(),
                     function(tree,details,matrices,node){
                       add_simple_labels(tree,details,matrices,node,
                                         query.field =query.field,
                                         query.allowed.df = query.allowed.df,
                                         label.field = label.field,
                                         cex.label =cex.label)})
  with(query.allowed.df,legend("topleft",legend=value,col=col,pch=pch))
  }
  if(type=="line") {res=add_annotation(tree,
                                         details,list(),
                                         function(tree,details,matrices,node){
                                           add_simple_labels_line(tree,details,matrices,node,
                                                                  query.field =query.field,
                                                                  query.allowed.df = query.allowed.df,
                                                                  label.field = label.field,
                                                                  cex.label=cex.label,
                                                                  lty=lty,
                                                                  lwd=lwd)})}
}

plot_sharing_info=function(tree,details,matrices,node,donor_ID,recip_ID,sharing_cols=c("black","dark green","red"),...){  #sharing_cols is a vector of colours for "shared", "donor only" and "recipient only" branches.
  info=get_edge_info(tree,details,node=node)
  n_donor=sum(grepl(donor_ID,info$samples))
  n_recip=sum(grepl(recip_ID,info$samples))
  sharing_info=ifelse(n_donor>0&n_recip>0,"shared",ifelse(n_donor>0,"donor","recipient"))
  names(sharing_cols)=c("shared","donor","recipient")
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=sharing_cols[sharing_info],lend=1,lwd=ifelse(sharing_info=="shared",1.5,1),...)
  }
}

plot_MAV_mut=function(tree,
                      details,
                      matrices,
                      node,
                      lesion_node=NA,
                      mut1,
                      mut2=NULL,
                      colours=c("dark gray","red","blue"),
                      cex=0.4,
                      #show_pval=FALSE,
                      ...) {
  #Define the col.scale from the colours vector
  require(dichromat)
  mut_colfunc = colorRampPalette(colours[-1])
  mut_colscale = mut_colfunc(11)
  
  #Get the tips within the lesion node
  if(!is.na(lesion_node)) {
    expected_samples=getTips(tree,lesion_node)
  }
  
  #Get the vaf
  info=get_edge_info(tree,details,node=node)
  samples=info$samples
  mut1_reads=sum(matrices$NV[mut1,samples])
  if(is.null(mut2)){mut2_reads=0}else{mut2_reads=sum(matrices$NV[mut2,samples])}
  variant_reads=mut1_reads+mut2_reads
  total_reads=variant_reads+(sum(matrices$NR[mut1,samples]) - mut1_reads)
  wt_reads=total_reads-variant_reads
  
  if((mut1_reads+mut2_reads)!=0) {
    if(is.null(mut2)) {
      base_col=colours[2]
    } else {
      base_col=mut_colscale[1+round(10*mut1_reads/(mut1_reads+mut2_reads),digits=0)] #How red or blue should the "mut" element of the colour be
    }
  	final_colfunc=colorRampPalette(c("light gray",base_col))
  	final_colscale=final_colfunc(101)
  	branch_col=final_colscale[1+round(100*(mut1_reads+mut2_reads)/total_reads)] #Now how "concentrated" should the mut colour be
  } else {
  	branch_col="light gray"
  }  
  
  #Plot the branches using the colour scale
  if(length(tree$edge.length[tree$edge[,2]==node])>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=branch_col,lend=1,...)
  }
  
  #Print the total depth for that branch at the node
  if(!is.na(lesion_node)){
    if(node %in% which(tree$tip.label %in% expected_samples)) {
      text(info$x,info$yb,paste0(variant_reads,"/",total_reads),srt=90,cex = cex,col="black",font=2)
    }
  }

  #Print the mutation name
  if(node==1) { #Do for node==1 so that only prints the name once
    if(is.null(mut2)) {
      text(1,1,pos=4,mut1)
    } else {
      text(1,1,pos=4,paste0(mut1,"/",strsplit(x=mut2,split="-")[[1]][4]))
    }
  }
}

highlight_groups=function(tree,group1,group2,cols=c("#17698E","#17A258")) {
  any_group=c(group1,group2) #make combined list of samples assigned to at least one group
  
  #Iterate through all the nodes & work if (1) all the tips are in group1 (2) all are in group2 (3) mixture of both
  col_vec=sapply(tree$edge[,2],function(node) {
    node_daughters=getTips(tree,node)
    if(any(node_daughters%in%any_group)) {
      node_daughters<-node_daughters[node_daughters%in%any_group] #only include those tips that are assigned to a specific group
      if(all(node_daughters%in%group1)) {
        col<-cols[1]
      } else if(all(node_daughters%in%group2)) {
        col<-cols[2]
      } else {
        col<-"black" #If has tips in both groups, colour black
      }
    } else {
      col<-"lightgrey"
    }
  })
  return(col_vec)
}

# ---- from phytools_scripts.R ----
getTips = function(tree,node) {
  require(ape)
  if(node <= length(tree$tip.label)) {
    daughters <- tree$tip.label[node]
  } else {
    daughters <- extract.clade(tree, node = node)$tip.label
  }
  return(daughters)
}

correct_edge_length = function(node, tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE, get_edge_from_tree=FALSE) {
  daughters <- getTips(tree = tree, node = node)
  #correct SNVs on edge, or set to 0 if want an indel only tree
  if(include_SNVs == TRUE) {
    if(get_edge_from_tree) {
    	nSNV=tree$edge.length[tree$edge[,2]==node]
    	} else {
    		nSNV = sum(details$node == node & details$Mut_type == "SNV")
    		}   		
    all_sens_SNVs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"SNV_sensitivity"]
    branch_SNV_sens = 1 - prod(1-all_sens_SNVs)  
    new_nSNV = nSNV/branch_SNV_sens
  } else {
    new_nSNV <- 0
    }
  #correct INDELs on edge, or set to 0 if want an SNV only tree
  if(include_indels == TRUE) {
    nINDEL = sum(details$node == node & details$Mut_type == "INDEL")
    all_sens_INDELs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"INDEL_sensitivity"]
    branch_INDEL_sens = 1 - prod(1-all_sens_INDELs)
    new_nINDEL = nINDEL/branch_INDEL_sens
  } else {
    new_nINDEL <- 0
  }
  new_edge_length = new_nSNV + new_nINDEL
  return(new_edge_length)
}

get_subset_tree = function(tree, details, v.field = "Mut_type", value = "SNV") {
  get_new_edge_length = function(node, tree, details,v.field,value) {
    sum(details$node == node & details[v.field] == value)
  }
  tree_subset = tree
  tree_subset$edge.length = sapply(tree$edge[,2], get_new_edge_length, tree = tree, details = details,v.field = v.field,value=value)
  return(tree_subset)
}

get_corrected_tree = function(tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE,get_edge_from_tree=FALSE) {
  tree_c = tree
  tree_c$edge.length = sapply(tree$edge[,2], correct_edge_length, tree = tree, details = details, sensitivity_df = sensitivity_df, include_indels = include_indels, include_SNVs=include_SNVs,get_edge_from_tree=get_edge_from_tree)
  return(tree_c)
}

get_mut_burden = function(tree) {
  mut_burden = nodeHeights(tree)[tree$edge[,2] %in% 1:length(tree$tip.label),2]
  return(mut_burden)
}

# ---- from Pop_size_estimation_functions.R ----
find.distance <- function(tree, from, to) {
  path <- nodepath(tree, from, to)
  res <- 0
  for (j in 2:length(path)) {
    index <- which(tree$edge[,2] == path[j] & tree$edge[,1] == path[j-1], arr.ind = TRUE)
    res <- res + tree$edge.length[index]
  }
  return(res)
}

length.normalise <- function(orig.tree, new.tree, curr.node, remaining.stick) {
    curr.node.children <- unlist(Descendants(orig.tree, curr.node, "children"))

    for (j in curr.node.children) {
      index <- which(orig.tree$edge[,2] == j & orig.tree$edge[,1] == curr.node, arr.ind = TRUE)

      if (j <= length(orig.tree$tip.label)) { #PhyloBuild fix: j is a node number (was compared with tip labels)
        new.tree$edge.length[index] <- remaining.stick
      } else {
        curr.node.tips <- unlist(Descendants(orig.tree, j, "tips"))
        curr.dist <- find.distance(orig.tree, curr.node, j)
        if (curr.dist == 0) {curr.dist <- 0.01} # So that no edge lengths are zero
        desc.lengths <- sapply(curr.node.tips, FUN = find.distance, tree = orig.tree, from = curr.node)
        new.tree$edge.length[index] <- if (mean(desc.lengths) > 0) remaining.stick * curr.dist / mean(desc.lengths) else 0 #PhyloBuild fix: avoid 0/0
        shorter.stick <- remaining.stick - new.tree$edge.length[index]

        # Call as recursive function
        new.tree <- length.normalise(orig.tree, new.tree, j, shorter.stick)
      }
    }
    return(new.tree)
}

make.ultrametric.tree <- function(tree) {
    root.number <- length(tree$tip.label) + 1
    ultra.tree <- length.normalise(tree, tree, root.number, 1)
    return(ultra.tree)
}

