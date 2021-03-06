

#' A wrapper that swaps the channel names with marker names 
#' 
#' It is useful to prepare multiple GatingSets for merging when the marker (instead of channel) names are consistent accross batches.
#' It basically updates the dimension names of the underling flow data as well as the gate definitions stored in GatingSets. 
#' It's important to be aware that the gates in the original GatingSet will be updated even a new GatingSet object will be returend as the result. 
#' 
#' @return a new GatingSet
#' @param gs GatingSet objects
#' @export 
swapChannelMarker <- function(gs){
  fs <- getData(gs)
  
  for(sn in sampleNames(fs)){
    fr <- fs@frames[[sn]]
    fs@frames[[sn]] <- swapChannelMarker_flowframe(fr, use.exprs = FALSE)
  }
  
  newColNames <- flowCore::colnames(fs@frames[[sn]])
  flowCore::colnames(fs) <- newColNames #assuming the order of colnames between fr and fs were consistent
  flowData(gs) <- fs
  
  fr <- fs@frames[[sn]]
  pd <- pData(parameters(fr))
  pd <- pd[!is.na(pd$desc), 2:1]
  colnames(pd) <- c("old", "new")
  
  updateChannels(gs, pd)
  
}

#' Preprocesses a Cytotrol flowFrame object
#'
#' Our goal here is to use swap the marker names and the channel names within a
#' \code{flowFrame} object to ensure that the \code{flowFrame} objects across
#' centers can be merged into a single \code{flowSet}.
#'
#'
#' @param fr the \code{flowFrame} object to preprocess
#' @param use.exprs logical indicate that if the data matrix will be updated as well.
#' @return the updated \code{flowFrame} object containing only the markers of
#' interest
#' @export 
swapChannelMarker_flowframe <- function(fr, use.exprs = TRUE) {
  
  
  fr_rownames <- rownames(pData(parameters(fr)))
  
  # Preprocesses each of the columns in the flow_frame
  for (j in seq_len(length(flowCore::colnames(fr)))) {
    
    marker_idx <- paste0(fr_rownames[j], "S")
    channel_idx <- paste0(fr_rownames[j], "N")
    
    marker <- description(fr)[[marker_idx]]
    channel <- description(fr)[[channel_idx]]
    
    # In the case the marker name is given, we swap the marker and channel
    # names.
    if (!is.null(marker)) {
      # Converts the marker names to a common name
      marker <- as.vector(marker)
      
      # Updates the channel with the marker
      description(fr)[[channel_idx]] <- marker 
      pData(parameters(fr))[j, "name"] <- marker
      
      # Updates the marker information in the flow_frame with the channel
      description(fr)[[marker_idx]] <- channel
      pData(parameters(fr))[j, "desc"] <- channel
    }
  }
  
  if(use.exprs)
    colnames(exprs(fr)) <- flowCore::colnames(fr)
  
  # Subset to markers of interest
  fr
}




#' a uility function to match fcs files based on the channels used by one example FCS file 
#' 
#' It uses \link{readFCSPar} to read parameters from FCS header to select target files, 
#' thus be used as a prefilter before \code{read.flowSet} or \link{read.ncdfFlowSet} call. 
#' 
#' @param x \code{character} vector giving the list of fcs files to match 
#' @param pattern \code{character} the example FCS file that contains the channels of interest
#' @return a \code{character} vector of fcs files that has the identical channels with \code{subset}
#' @export 
#' @examples 
#' \dontrun{
#' grep.FCS(pattern = bcells[2],  x = c(bcells,tcells))
#'  #return TRUE  TRUE FALSE FALSE
#' }
grep.FCS <- function(pattern, x){
#      browser()
      targetChnls <- readFCSPar(pattern)
      unname(sapply(x, function(thisFile){
                          thisChnls <- readFCSPar(thisFile)
                          setequal(thisChnls, targetChnls)
                          
                        })
                )
      
      
    }


#' a wrapper for \code{save_gslist}
#' @return a copy of original gslist with modified cdf path when cdf == "move"
#' otherwise, it behaves the same as \code{save_gslist}
#' @export 
save_gslist_labkey <- function(gslist, path, cdf, ...){
  
  save_gslist(gslist, path, cdf = cdf, ...)
  
  if(cdf == "move"){
      
      newListOfGS <- lapply(gslist, function(thisGS){
      cdfName <- basename(flowData(thisGS)@file)
      newFullPath <- file.path(path, thisGS@guid, cdfName)
      flowData(thisGS)@file <-  newFullPath
      thisGS
      }, level = 1)
      
      GatingSetList(newListOfGS)
      
      }
}

#' plot by prarent index
#' 
#' This API is mainly used for labkey module. It takes a parent index instead of the actual gate index.
#' When there is no gate associated with the x,y channel specified by user, it simply plots the \code{xyplot} 
#' or \code{densityplot} without the gate. 
#' 
#' @param x \code{character} x channel
#' @param y \code{character} y channel, if \code{NULL},then try to do \code{densityplot}
#' @export 
#' @importFrom BiocGenerics colnames

plotGate_labkey <- function(G,parentID,x,y,smooth=FALSE,cond=NULL,xlab=NULL,ylab=NULL, overlay = NULL, overlay.symbol = NULL, key = NULL, ...){
  #get all childrens
  
  cids <- getChildren(G[[1]], parentID, showHidden = FALSE, path = "auto")
  if(length(cids)>0)
  {
    #try to match to projections
#		browser()
    isMatched<-lapply(cids,function(cid){
          g<-getGate(G[[1]],cid)
          if(class(g)!="booleanFilter") 
          {
            prj<-parameters(g)
            if(length(prj)==1)#1d gate
            {
              return (prj%in%c(x,y))
              
            }else
            {
              #2d gate but y is absent
              if(is.null(y))
                return (FALSE)
              #try to match x,y to 2d gate
              revPrj<-rev(prj)
              if((x==prj[1]&&y==prj[2])||(x==revPrj[1]&&y==revPrj[2]))
                return (TRUE)
              else
                return (FALSE)	
            }
          }else
            return (FALSE)
        })
    
    ind<-which(unlist(isMatched))
    if(length(ind)>0)
      isPlotGate<-TRUE
    else
      isPlotGate<-FALSE
  }else
    isPlotGate<-FALSE
#  browser()
  formula1 <- flowWorkspace:::mkformula(c(y,x),isChar=TRUE)
#  formula1<-paste("`",y,"`~`",x,"`",sep="")
  if(!is.null(cond))
    formula1<-paste(formula1,cond,sep="|")
  formula1 <- as.formula(formula1)
	
  type <- ifelse(is.null(y), "densityplot","xyplot")
  if(isPlotGate)
    plotGate(G,cids[ind],formula=formula1,smooth=smooth,xlab=xlab,ylab=ylab, type = type, overlay = overlay, ...)
  else
  {
    fs<-getData(G,parentID)
    axisObject <- flowWorkspace:::.formatAxis(x=G[[1]],data=fs[[1]],xParam=x,yParam=y,...)
    if(is.null(xlab)){
      xlab <- axisObject$xlab
    }
    if(is.null(ylab)){
      ylab <- axisObject$ylab
    }
    if(type == "xyplot"){
#      browser()
      if(!is.null(overlay)){
        if(is.null(overlay.symbol)){
#            browser()
          # set symbol color automatically if not given
          nOverlay <- length(overlay)
          overlay.fill <- RColorBrewer::brewer.pal(9,"Set1")[1:nOverlay]
          names(overlay.fill) <- overlay
          overlay.symbol <- sapply(overlay.fill, function(col)list(fill = col), simplify = FALSE)
        }
        #set legend automatically if it is not given
        if(is.null(key)){
          
          key = list(text = list(names(overlay.symbol), cex = 0.6)
              , points = list(col = sapply(overlay.symbol, "[[", "fill") 
                  , pch = 19
                  , cex = 0.5) 
              , columns = length(overlay.symbol)
              , between = 0.3
              , space = "bottom"
              , padding.text = 5)
        }
#        browser()
      overlay <- sapply(overlay, function(thisOverlay)getData(G,thisOverlay)[,c(y,x)])
    }
      xyplot(formula1
          ,fs
          ,smooth=smooth
          ,xlab=xlab
          ,ylab=ylab
          ,scales=axisObject$scales
          , overlay = overlay
          , overlay.symbol = overlay.symbol
          , key = key
          ,...
              )  
    }else{
      densityplot(formula1
          ,fs
          ,xlab=xlab
          ,scales=axisObject$scales
          ,...)
    }
    
  }
  
}



##merge gs 
#' @param force \code{logical} if TRUE, drop any channels if neccessry, 
#'                          otherwise, be conservative by only dropping unused channels
.mergeGS <- function(this_gslist, force = FALSE){
  
        
        .Defunct("dropRedundantChannels and dropRedundantNodes")
        
        if(length(this_gslist) > 1){
          #find the common colnames
          col_list <- lapply(this_gslist,function(this_gs)colnames(flowData(this_gs)))
          global_colnames <- Reduce(intersect, col_list)
          
          if(is.null(global_colnames))
            stop("Can't merge!no common channels.")
          
          this_gslist <- lapply(this_gslist,function(this_gs){
#                    browser()
                this_fs <- getData(this_gs)
                
                if(force){
                  toDrop <- setdiff(colnames(this_fs), global_colnames)
                  if(length(toDrop) >0)
                    message("drop ", toDrop)
                  flowData(this_gs) <- this_fs[,global_colnames]
                }else
                {
                  #drop the unused marker from fs                    
                  this_fs_colnames <- colnames(this_fs)
                  this_fr <- this_fs[[1]]
                  this_pd <- pData(parameters(this_fr))
                  within_common_chnnl <- this_fs_colnames%in%global_colnames
                  non_na_channel <- unname(!is.na(this_pd[,"desc"]))
                  to_include <- grepl(pattern="[FS]SC|[Tt]ime",this_pd[,"name"])
                  to_include <- to_include |  non_na_channel | within_common_chnnl
                  
                  if(length(which(to_include)) != nrow(this_pd)){
                    #drop channels from colnames of flowFrame
                    message("drop empty channel:",this_pd[!to_include,1])
                    fr_colnames <- colnames(this_fr)
                    fr_colnames <- fr_colnames[to_include]
                    #update the colnames of flowSet accordingly                                   
                    this_fs_colnames <- this_fs_colnames[match(fr_colnames,this_fs_colnames)]
                  }
                  
                  
                  if(!setequal(global_colnames,this_fs_colnames))
                    stop("merge failed!These channels are not common across data sets:\n"
						  , paste0(setdiff(this_fs_colnames,global_colnames))	 
                          , "\n If they are non-NA channels, try force = TRUE to drop them."
                          )
                }
                
                #reorder colnames of other gs by global_colnames
                flowData(this_gs) <- this_fs[,global_colnames] 
                
                
                this_gs
              })
        }
        
        GatingSetList(this_gslist)
      
}

#' this is the wrapper for labkey where only one gslist to be returned
#merge_gs_labkey <- function(x,...){
#  gs_groups <- .groupByTree(x, drop = TRUE)
# 
#  if(length(gs_groups)>1){
#      stop("Can't merge because multiple gating trees are present!")
#  }else{
#    res <- .mergeGS(gs_groups)
#    res[[1]]
#  }
#
#}

#TODO: to deprecate
#' cluster/merge GatingSets based on the gating tree structures.
#' 
#' merge GatingSets based on the gating tree structures.
#' 
#' @details Group the individual GatingSets by their gating schemes.It is done by comparing the node list returned by \code{\link{getNodes}},which assumes
#' they follow the same population naming conventions.
#' Meanwhile the unused channels are automatically dropped to make sure the flow data has identical data structure within each group.
#' In order to further merge multiple GatingSet objects into one, use \code{\link{rbind2}}.     
#' 
#' @param x A \code{list} of \code{GatingSet}s . 
#' @return A \code{\link{GatingSetList}} that contains multiple GatingSets each of which share the same gating and data structure.
#' @author Mike Jiang \email{wjiang2@@fhcrc.org}
#' @seealso \code{\link{rbind2}},\code{\link{GatingSetList}}
#' @examples \dontrun{
#' 	#load gatingsets from disk
#' 	#gs_toMerge is the path that stores multiple archived gatingsets
#' 	gs_list<-lapply(list.files("flowIncubator/output/gs_toMerge",full=T),function(this_folder){
#'       flowWorkspace:::load_gs(this_folder)
#'     })
#'     
#' 	gs_list <- merge_gs(gs_list)
#' 	gs_list
#' 	
#' }
#' @export 
merge_gs<-function(x,...){
#      browser()
  .Defunct("dropRedundantChannels and dropRedundantNodes")
      message("Grouping by Gating tree...")
      node_seq <-unlist(lapply(x,function(this_gs){
                this_gh <- this_gs[[1]]
                this_nodes <- getNodes(this_gh, showHidden = TRUE)
                paste(this_nodes,collapse = "")
                
              }))
      gs_groups <- split(x,node_seq)
        
      
      #start to merge      
      lapply(1:length(gs_groups),function(i){
#            browser()
            this_group <- gs_groups[[i]]

            #drop the unused marker from fs
            if(length(this_group) > 1){
              this_group <- lapply(this_group,function(this_gs){
#                    browser()
                    this_fs <- getData(this_gs)
                    
                    this_pd <- pData(parameters(this_fs[[1]]))
                    non_na_channel <- unname(!is.na(this_pd[,"desc"]))
                    to_include <- grepl(pattern="[FS]SC|[Tt]ime",this_pd[,"name"])
                    to_include <- to_include |  non_na_channel
                    if(length(which(to_include)) != nrow(this_pd)){
                    
                      message("drop empty channel:",this_pd[!to_include,1])
                      
                      flowData(this_gs) <- this_fs[,to_include]
                    
                    }
                    this_gs
                  })
            }
            
            GatingSetList(this_group)
          })
        
      
    }

plotGateM <- function(formula,gs, parent, children, ...){
  
  
#' get parent data
  fsParent <- lapply(nodes, function(node){
        fs <- getData(gs, parent)
        fs <- fs@data[[1]]
        fs <-as.flowSet(fs)
        pData(fs)[,"subset"] <- node
        sn <- sampleNames(fs)
        sampleNames(fs) <- paste0(sn, node)
        fs
      })
  fsParent <- ncdfFlowList(fsParent)
  
  #' get overlay data
  fslist <- lapply(children, function(node){
        fs <- getData(gs, node)
        fs <- fs@data[[1]]
        fs <-as.flowSet(fs)
        pData(fs)[,"subset"] <- node
        sn <- sampleNames(fs)
        sampleNames(fs) <- paste0(sn, node)
        fs
      })
  fslist <- ncdfFlowList(fslist)
  
  xyplot(formula, fsParent, overlay = fslist, ...)
  
}    
