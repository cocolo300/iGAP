library(ape)
library(ggdendro)
library(ggtree)
library(patchwork)
library(RColorBrewer)
library(RSDA)
library(ggESDA)
library(seriation)
library(dataSDA)

#' @name iGAP
#' @title Symbolic interval data with Generalized Association Plots
#' @description Let interval data be visualized as a heatmap, and the correlations heatmap
#' between interval data variables and the distances heatmap between concepts can be collaged to
#' explore the hidden structure of interval data.
#' @import ggplot2 ggthemes ape ggdendro ggtree patchwork RColorBrewer RSDA ggESDA seriation dataSDA
#' @param data A iGAP object. It can also be either RSDA object or
#' classical data frame.
#' @param standardize_condition Boolean variables, which mean the color
#' present by standardize condition (if TRUE) or matrix condition (if FALSE).
#' @param cor_up For the upper-side association correlation matrix using
#' methods (centers, BD, BG)
#' @param distance_right For the right side association distance matrix using
#' methods (Gowda.Diday, Ichino, Minkowski, Hausdorff)
#' @param normalize_right A logical value indicating whether normalize the data in the ichino or hausdorff method.
#' @param SpanNormalize_right A logical value indicating whether.
#' @param euclidea_right A logical value indicating whether use the euclidean distance.
#' @param method_up For the upper-side association distance matrix use sorting
#' methods (R2E, GW, HC, OLO, MDS, isoMDS, isomap, monoMDS, metaMDS, Sammon,
#' ARSA, GSA, SGD, QAP, TSP, Identity, Random,  Reverse, SPIN, VAT, ward.D,
#' ward.D2, single, complete", average , mcquitty , median, centroid)
#' @param method_right For the right side association distance matrix use sorting
#' methods (R2E, GW, HC, OLO, MDS, isoMDS, isomap, monoMDS, metaMDS, Sammon,
#' ARSA, GSA, SGD, QAP, TSP, Identity, Random, Reverse, SPIN, VAT, ward.D,
#'  ward.D2, single, complete", average , mcquitty , median, centroid)
#' @param Rcolor Color scale to use for interval data.
#' @param Rcolorr Color scale to use for the right side association distance matrix.
#' @param Rcoloru Color scale to use for the upper-side association correlation matrix.
#' @return Return a ggplot2 object.
#' @usage function(data = NULL,
#'          standardize_condition = TRUE,
#'          cor_up = "centers",
#'          distance_right = "Hausdorff",
#'          normalize_right = TRUE,
#'          SpanNormalize_right = TRUE,
#'          euclidea_right = TRUE,
#'          method_up = NULL,
#'          method_right = NULL,
#'          Rcolor = "Spectral",
#'          Rcolorr = "Spectral",
#'          Rcoloru = "RdBu")
#' @examples
#' p <- iGAP(facedata)
#' p
#'
#' test <- classic.to.sym(x = ex1_db2so, concept = c(state, sex), variables = c(county, group, age))
#' p <- iGAP(test, standardize_condition = F,
#'           cor_up = c("BG"),
#'           distance_right = "Hausdorff",
#'           normalize_right = T,
#'           SpanNormalize_right = T,
#'           euclidea_right = T,
#'           method_up = "R2E",
#'           method_right = "R2E",
#'           Rcolor = "Spectral",
#'           Rcolorr = "Spectral",
#'           Rcoloru = "RdBu")
#' p
#'
#' p <- iGAP(test, standardize_condition = F,
#'           cor_up = c("BG"),
#'           distance_right = "Hausdorff",
#'           normalize_right = T,
#'           SpanNormalize_right = T,
#'           euclidea_right = T,
#'           method_up = "single",
#'           method_right = "single",
#'           Rcolor = "Spectral",
#'           Rcolorr = "Spectral",
#'           Rcoloru = "RdBu")
#' p
#'
#' @export
iGAP <- function(data = NULL,standardize_condition=T,cor_up = "BG",distance_right = "Hausdorff",
                 normalize_right = TRUE,
                 SpanNormalize_right = TRUE,
                 euclidea_right = TRUE,
                 method_up = NULL,
                 method_right = NULL,
                 Rcolor="Spectral",
                 Rcolorr="Spectral",
                 Rcoloru="RdBu") {


  tree <- c("ward.D", "ward.D2", "single", "complete", "average" , "mcquitty" , "median", "centroid")
  seriation <- c("R2E","GW","HC","OLO","MDS","isoMDS","isomap","monoMDS","metaMDS","Sammon","ARSA","GSA","SGD","QAP","TSP","Identity","Random"
                 ,"Reverse","SPIN","VAT")

  ggSymData <- testData(data)
  iData <- ggSymData$intervalData
  myHeatMapNamesr <- rownames(iData)
  myHeatMapNamesc <- colnames(iData)
  p<-dim(data)[2]
  n<-dim(data)[1]

  qq1 <- sym.dist.interval(
    data,
    method = distance_right,
    normalize = normalize_right,
    SpanNormalize = SpanNormalize_right,
    euclidea = euclidea_right
  )

  A <- cor(data, method = cor_up)
  AA <- as.dist(A)

  if (!is.null(method_up) && !is.null(method_right)) {
    if (method_up %in% seriation && method_right %in% seriation) {

      o_up<-seriate((1-AA),method=method_up)
      oo_up<-get_order(o_up)
      one<-data[,oo_up]
    suppressWarnings({
      rownames(one) <- c(myHeatMapNamesr)
      })
      #右排序
      o_right<-seriate(qq1,method=method_right)
      oo_right<-get_order(o_right)
      two<-one[oo_right,]
     suppressWarnings({
      rownames(two) <- c(myHeatMapNamesr[oo_right])
      })
      #plot
      y<-ggInterval_indexImage999(two,standardize_condition=standardize_condition,Rcolor=Rcolor)

      y1<-Plot_right(two,Rcolorr=Rcolorr,
                     distance_right = distance_right,
                     normalize_right = normalize_right,
                     SpanNormalize_right = SpanNormalize_right,
                     euclidea_right = euclidea_right)

      y2<-Plot_up(two,Rcoloru=Rcoloru,cor_up = cor_up)

      if(n>2*p){
        layout <- "
    11#####
    11#####
    2233333
    2233333
    2233333
    2233333
    2233333
  "
      } else{
        layout <- "
    111###
    111###
    111###
    222333
    222333
    222333
  "
      }

      combined_plot <- y2+y+y1+guide_area()+plot_layout(design = layout,guides = 'collect')&
        theme(legend.position='bottom')

    } else if (method_up %in% tree && method_right %in% tree) {
      #上樹
      phup <- (1-AA) %>%
        hclust(method = method_up) %>%
        as.phylo.hclust() %>%
        ggtree() + scale_x_reverse()
      dendro <- (1-AA) %>%
        hclust(method = method_up)
      #倒
      three<-data[,rev(dendro$order)]
     suppressWarnings({
      rownames(three) <- c(myHeatMapNamesr)
      })
      #右樹
      qq2 <- sym.dist.interval(
        three,
        method = distance_right,
        normalize = normalize_right,
        SpanNormalize = SpanNormalize_right,
        euclidea = euclidea_right
      )
      phright <- qq2 %>%
        hclust(method = method_right) %>%
        as.phylo.hclust() %>%
        ggtree() + scale_x_reverse()
     suppressWarnings({
      myHeatMapNamesr1<-rownames(three)
      })
      dendro1 <- qq2 %>%
        hclust(method = method_right)
      #倒
      four<-three[rev(dendro1$order),]
    suppressWarnings({
      rownames(four) <- c(myHeatMapNamesr1[rev(dendro1$order)])
      })
      #plot
      y<-ggInterval_indexImage999(four,standardize_condition=standardize_condition,Rcolor=Rcolor)

      y1<-Plot_right(four,Rcolorr=Rcolorr,
                     distance_right = distance_right,
                     normalize_right = normalize_right,
                     SpanNormalize_right = SpanNormalize_right,
                     euclidea_right = euclidea_right)

      y2<-Plot_up(four,Rcoloru=Rcoloru,cor_up = cor_up)

      if(n>2*p){
        layout <- "
    112#####
    112#####
    33444445
    33444445
    33444445
    33444445
    33444445
  "
      }  else{
        layout <- "
    1112###
    1112###
    1112###
    3334445
    3334445
    3334445
  "
      }
      combined_plot <- y2+phup+y+y1+phright+guide_area()+plot_layout(design = layout,guides = 'collect')&
        theme(legend.position='bottom')
    } else if (method_up %in% seriation && method_right %in% tree) {
      o_up<-seriate((1-AA),method=method_up)
      oo_up<-get_order(o_up)
      one<-facedata[,oo_up]
     suppressWarnings({
      rownames(one) <- c(myHeatMapNamesr)
     })
      #右樹
      qq2 <- sym.dist.interval(
        one,
        method = distance_right,
        normalize = normalize_right,
        SpanNormalize = SpanNormalize_right,
        euclidea = euclidea_right
      )
      phright <- qq2 %>%
        hclust(method = method_right) %>%
        as.phylo.hclust() %>%
        ggtree() + scale_x_reverse()
      myHeatMapNamesr1<-rownames(one)
      dendro1 <- qq2 %>%
        hclust(method = method_right)
      #倒
      four<-one[rev(dendro1$order),]
    suppressWarnings({
      rownames(four) <- c(myHeatMapNamesr1[rev(dendro1$order)])
      })
      #plot
      y<-ggInterval_indexImage999(four,standardize_condition=standardize_condition,Rcolor=Rcolor)

      y1<-Plot_right(four,Rcolorr=Rcolorr,
                     distance_right = distance_right,
                     normalize_right = normalize_right,
                     SpanNormalize_right = SpanNormalize_right,
                     euclidea_right = euclidea_right)

      y2<-Plot_up(four,Rcoloru=Rcoloru,cor_up = cor_up)

      if(n>2*p){
        layout <- "
    112#####
    112#####
    33444445
    33444445
    33444445
    33444445
    33444445
  "
      }  else{
        layout <- "
    1112###
    1112###
    1112###
    3334445
    3334445
    3334445
  "
      }
      combined_plot <- y2+plot_spacer()+y+y1+phright+guide_area()+plot_layout(design = layout,guides = 'collect')&
        theme(legend.position='bottom')

    } else if (method_up %in% tree && method_right %in% seriation) {
      #上樹
      phup <- (1-AA) %>%
        hclust(method = method_up) %>%
        as.phylo.hclust() %>%
        ggtree() + scale_x_reverse()
      dendro <- (1-AA) %>%
        hclust(method = method_up)
      #倒
      three<-data[,rev(dendro$order)]
    suppressWarnings({
      rownames(three) <- c(myHeatMapNamesr)
      })
      #右排序
      o_right<-seriate(qq1,method=method_right)
      oo_right<-get_order(o_right)
      two<-three[oo_right,]
      suppressWarnings({
      rownames(two) <- c(myHeatMapNamesr[oo_right])
      })
      #plot
      y<-ggInterval_indexImage999(two,standardize_condition=standardize_condition,Rcolor=Rcolor)

      y1<-Plot_right(two,Rcolorr=Rcolorr,
                     distance_right = distance_right,
                     normalize_right = normalize_right,
                     SpanNormalize_right = SpanNormalize_right,
                     euclidea_right = euclidea_right)

      y2<-Plot_up(two,Rcoloru=Rcoloru,cor_up = cor_up)

      if(n>2*p){
        layout <- "
    112#####
    112#####
    33444445
    33444445
    33444445
    33444445
    33444445
  "
      }  else{
        layout <- "
    1112###
    1112###
    1112###
    3334445
    3334445
    3334445
  "
      }
      combined_plot <- y2+phup+y+y1+plot_spacer()+guide_area()+plot_layout(design = layout,guides = 'collect')&
        theme(legend.position='bottom')
    }
  } else if (!is.null(method_up)) {
    if (method_up %in% seriation) {
      o_up<-seriate((1-AA),method=method_up)
      oo_up<-get_order(o_up)
      one<-data[,oo_up]
    suppressWarnings({
      rownames(one) <- c(myHeatMapNamesr)
      })
      #plot
      y<-ggInterval_indexImage999(one,standardize_condition=standardize_condition,Rcolor=Rcolor)

      y1<-Plot_right(one,Rcolorr=Rcolorr,
                     distance_right = distance_right,
                     normalize_right = normalize_right,
                     SpanNormalize_right = SpanNormalize_right,
                     euclidea_right = euclidea_right)

      y2<-Plot_up(one,Rcoloru=Rcoloru,cor_up = cor_up)

      if(n>2*p){
        layout <- "
    11#####
    11#####
    2233333
    2233333
    2233333
    2233333
    2233333
  "
      } else{
        layout <- "
    111###
    111###
    111###
    222333
    222333
    222333
  "
      }
      combined_plot <- y2+y+y1+guide_area()+plot_layout(design = layout,guides = 'collect')&
        theme(legend.position='bottom')

    } else if (method_up %in% tree) {
      #上樹
      phup <- (1-AA) %>%
        hclust(method = method_up) %>%
        as.phylo.hclust() %>%
        ggtree() + scale_x_reverse()
      dendro <- (1-AA) %>%
        hclust(method = method_up)
      #倒
      three<-data[,rev(dendro$order)]
     suppressWarnings({
      rownames(three) <- c(myHeatMapNamesr)
      })
      #plot
      y<-ggInterval_indexImage999(three,standardize_condition=standardize_condition,Rcolor=Rcolor)

      y1<-Plot_right(three,Rcolorr=Rcolorr,
                     distance_right = distance_right,
                     normalize_right = normalize_right,
                     SpanNormalize_right = SpanNormalize_right,
                     euclidea_right = euclidea_right)

      y2<-Plot_up(three,Rcoloru=Rcoloru,cor_up = cor_up)

      if(n>2*p){
        layout <- "
    112#####
    112#####
    33444445
    33444445
    33444445
    33444445
    33444445
  "
      }  else{
        layout <- "
    1112###
    1112###
    1112###
    3334445
    3334445
    3334445
  "
      }
      combined_plot <- y2+phup+y+y1+plot_spacer()+guide_area()+plot_layout(design = layout,guides = 'collect')&
        theme(legend.position='bottom')
    }
  } else if (!is.null(method_right)) {
    if (method_right %in% seriation) {
      #右排序
      o_right<-seriate(qq1,method=method_right)
      oo_right<-get_order(o_right)
      two<-data[oo_right,]
    suppressWarnings({
      rownames(two) <- c(myHeatMapNamesr[oo_right])
      })
      #plot
      y<-ggInterval_indexImage999(two,standardize_condition=standardize_condition,Rcolor=Rcolor)

      y1<-Plot_right(two,Rcolorr=Rcolorr,
                     distance_right = distance_right,
                     normalize_right = normalize_right,
                     SpanNormalize_right = SpanNormalize_right,
                     euclidea_right = euclidea_right)

      y2<-Plot_up(two,Rcoloru=Rcoloru,cor_up = cor_up)

      if(n>2*p){
        layout <- "
    11#####
    11#####
    2233333
    2233333
    2233333
    2233333
    2233333
  "
      } else{
        layout <- "
    111###
    111###
    111###
    222333
    222333
    222333
  "
      }

      combined_plot <- y2+y+y1+guide_area()+plot_layout(design = layout,guides = 'collect')&
        theme(legend.position='bottom')
    } else if (method_right %in% tree) {
      #右樹
      phright <- qq1 %>%
        hclust(method = method_right) %>%
        as.phylo.hclust() %>%
        ggtree() + scale_x_reverse()
      dendro1 <- qq1 %>%
        hclust(method = method_right)
      #倒
      four<-data[rev(dendro1$order),]
      suppressWarnings({
      rownames(four) <- c(myHeatMapNamesr[rev(dendro1$order)])
      })
      #plot
      y<-ggInterval_indexImage999(four,standardize_condition=standardize_condition,Rcolor=Rcolor)

      y1<-Plot_right(four,Rcolorr=Rcolorr,
                     distance_right = distance_right,
                     normalize_right = normalize_right,
                     SpanNormalize_right = SpanNormalize_right,
                     euclidea_right = euclidea_right)

      y2<-Plot_up(four,Rcoloru=Rcoloru,cor_up = cor_up)

      if(n>2*p){
        layout <- "
    112#####
    112#####
    33444445
    33444445
    33444445
    33444445
    33444445
  "
      }  else{
        layout <- "
    1112###
    1112###
    1112###
    3334445
    3334445
    3334445
  "
      }

      combined_plot <- y2+plot_spacer()+y+y1+phright+guide_area()+plot_layout(design = layout,guides = 'collect')&
        theme(legend.position='bottom')
    }
  } else {
    #plot
    y<-ggInterval_indexImage999(data,standardize_condition=standardize_condition,Rcolor=Rcolor)

    y1<-Plot_right(data,Rcolorr=Rcolorr,
                   distance_right = distance_right,
                   normalize_right = normalize_right,
                   SpanNormalize_right = SpanNormalize_right,
                   euclidea_right = euclidea_right)

    y2<-Plot_up(data,Rcoloru=Rcoloru,cor_up = cor_up)

  if(n>2*p){
    layout <- "
    11#####
    11#####
    2233333
    2233333
    2233333
    2233333
    2233333
  "
  }else{
    layout <- "
    111###
    111###
    111###
    222333
    222333
    222333
  "
    }


    combined_plot <- y2+y+y1+guide_area()+plot_layout(design = layout,guides = 'collect')&
      theme(legend.position='bottom')
  }


  return(combined_plot)
}


scale_sym_table <- function(d, n, p){
  temp1 <- sapply(1:p, FUN = function(x) unlist(data.frame(d[[x]])))
  temp2 <- apply(temp1, 2, scale)
  newd <- data.frame(temp2[1:n, ], temp2[(n+1):(n*2), ])
  myd <- classic2sym(newd, groupby = "customize",
                     minData = temp2[1:n, ],
                     maxData = temp2[(n+1):(n*2), ])
  colnames(myd$intervalData) <- colnames(d)
  return(myd)
}
#test whether data can be used for ggplot
testData <- function(data){
  if("ggESDA" %in% class(data)){ # if ggESDA class?
    return(data)
  }else{
    if(("symbolic_tbl" %in% class(data))){#if RSDA class?
      return(RSDA2sym(data))
    }else{
      warning("Automatically transform a classical data to symbolic data")
      return(classic2sym(data))
    }
  }
}


####################interval data build
buildPlotData <- function(aa, bb,iData,adjustStrip,ggSymData, n, datasetMax){
  datasetMax<-max(ggSymData$statisticsDF$max)
  n<-dim(iData)[1]
  p<-dim(iData)[2]

  #  aa <- c()
  #  bb <- c()
  #  for(i in colnames(iData)){
  #    a<-iData[[i]]$min
  #    b<-iData[[i]]$max
  # aa <- c(aa,a)
  #  bb <- c(bb,b)
  # }

  aa <- sapply(iData, function(x) x$min)
  bb <- sapply(iData, function(x) x$max)

  adjustStrip <- min(bb-aa)/50
  datasetMax<-ceiling(datasetMax)
  datasetMax <- ifelse(datasetMax >200, 200, datasetMax)
  adjustCoef<-ifelse(datasetMax*n>2000,1,ceiling(2000/(datasetMax*n)))

  myAdjust <- datasetMax*adjustCoef

  #  yresult <- c()

  #  for (i in 0:(p - 1)) {
  #    for (j in 1:n) {
  #      for (k in 1:myAdjust) {
  #        yresult <- c(yresult, i * myAdjust + k)
  #      }
  #   }
  #  }
  i_vals <- rep(seq(0, (p - 1) * myAdjust, by = myAdjust), each = n * myAdjust)
  k_vals <- rep(seq(1, myAdjust), times = p * n)
  yresult <- i_vals + k_vals

  #I咬換ARRAY帶進去地I個ARRAY
  d2 <- data.frame(x = rep(n:1,p,each=datasetMax*adjustCoef)-0.5,
                   xend = rep((n+0.5):1.5,p,each=datasetMax*adjustCoef),
                   y = yresult)
  val2<-mapply(aa,bb,FUN=function(x,y) sort(runif(datasetMax*adjustCoef,x,y)))
  val<-matrix(val2,ncol=1)
  d<-data.frame(d2,value=val)

  return(d)
}

###################################interval plot
suppressWarnings({
  ggInterval_indexImage999<-function(data = NULL,standardize_condition=TRUE,Rcolor_interval="Spectral"){

    ggSymData <- testData(data)
    iData <- ggSymData$intervalData
    myHeatMapNames <- rownames(iData)
    reversed_myHeatMapNamesr <- myHeatMapNames[length(myHeatMapNames):1]


    p<-dim(iData)[2]
    n<-dim(iData)[1]


    datasetMin<-min(ggSymData$statisticsDF$min)
    datasetMax<-max(ggSymData$statisticsDF$max)

    datasetMax_size<-ceiling(datasetMax)
    adjustCoef<-ifelse(datasetMax_size*n>2000,1,ceiling(2000/(datasetMax_size*n)))
    size_size<-(-0.0009615385*datasetMax*adjustCoef*p)+1.7053846153846153
    size_size <- ifelse(size_size < 0, 0, size_size)
    with(data,{
      #add heatmap


      #get numerical data
      numericData <- unlist(lapply(data.frame(iData[1:dim(iData)[2]]) ,FUN = is.sym.interval))
      iData <- iData[,which(numericData)]

      #scale
      if(standardize_condition){
        iData <- scale_sym_table(iData, n, p)$intervalData
      }

      d <- buildPlotData(aa, bb,iData, adjustStrip, ggSymData,n, datasetMax)
      #whether column condition
      if(!standardize_condition){
        #midp <- allDataMean
        NAME <- "Matrix Condition"



        base <- ggplot(d, aes(x = x, xend = xend, y = y, yend = y,color=value)) +
          geom_segment(size = size_size) +
          scale_colour_distiller(palette = Rcolor_interval,
                                 name = 'Value',
                                 breaks=c(min(d$value),max(d$value)),
                                 labels = c(round(min(d$value)),round(max(d$value))))+
          coord_flip()+ scale_x_continuous(breaks=c(1:n),labels = reversed_myHeatMapNamesr,expand=c(0,0))+
          labs(x="Concepts",y="")+
          theme(panel.grid=element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank()) +scale_y_continuous(breaks = NULL,expand=c(0,0)) +xlab(NULL) + ylab(NULL)+xlab(NULL) + ylab(NULL)
        return(base)


      }else{# column condition
        NAME <- "Column Condition"


        p <- ggplot(d, aes(x = x, xend = xend, y = y, yend = y, color = value))

        p <- p+geom_segment(size = size_size)+
          scale_colour_distiller(palette = Rcolor_interval,
                                 breaks=c(min(d$value),max(d$value)),
                                 labels = c(round(min(d$value), digits = 3),round(max(d$value), digits = 3)))+
          scale_x_continuous(breaks=c(1:n),labels = reversed_myHeatMapNamesr,expand=c(0,0))+
          coord_flip()+
          guides(fill=none,alpha=none)+
          labs(x = "Concepts",y="")+
          theme(
            panel.grid=element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())+
          scale_y_continuous(breaks = NULL,expand=c(0,0)) +xlab(NULL) + ylab(NULL)
        return(p)
      }

    })
  }
})
###########################distance data build
buildPlotData_right <- function(data = NULL,distance_right = "Hausdorff", normalize_right = TRUE,
                                SpanNormalize_right = TRUE,
                                euclidea_right = TRUE){
  n<-dim(data)[1]
  qq1<-sym.dist.interval(
    data,
    method = distance_right,
    normalize = normalize_right,
    SpanNormalize = SpanNormalize_right,
    euclidea = euclidea_right)
  m <- as.matrix(qq1)
  df<-data.frame(m,row.names=rownames(data))
  df$names <- rownames(data)
  mm1 <- gather(df, 1:n, key="condition", value='distance')
  d2 <- data.frame(x = rep(n:1,n), y = rep(1:n, each=n))
  d<-data.frame(mm1,d2)
  return(d)
}


#####################################cor data build
buildPlotData_up <- function(data = NULL,cor_up = c("centers")){
  p<-dim(data)[2]
  mm<-cor(data,method = cor_up)
  mm$names <- colnames(data)
  mm1 <- gather(mm, 1:p, key="condition", value='correlation')
  d2 <- data.frame(x = rep(p:1,p), y = rep(1:p, each=p))
  d<-data.frame(mm1,d2)
  return(d)
}
####################################distance plot

Plot_right <- function(data = NULL,Rcolorr="Spectral",distance_right = "Hausdorff", normalize_right = TRUE,
                       SpanNormalize_right = TRUE,
                       euclidea_right = TRUE){
  r<-buildPlotData_right(data,
                         distance_right = distance_right,
                         normalize_right = normalize_right,
                         SpanNormalize_right = SpanNormalize_right,
                         euclidea_right = euclidea_right)
  rr<-ggplot(r, aes(x, y)) +
    geom_tile(aes(fill = distance))+scale_fill_distiller(palette=Rcolorr)+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+coord_flip()+
    theme(panel.grid=element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank()
          ,axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()
          ,axis.text.y = element_blank(),
          axis.title= element_blank()) +xlab(NULL) + ylab(NULL)
  return(rr)

}

####################################cor plot
Plot_up <- function(data = NULL,Rcoloru="RdBu",cor_up = c("centers")){
  p<-dim(data)[2]
  u<-buildPlotData_up(data,cor_up = cor_up)
  uu<-ggplot(u, aes(x, y)) +
    geom_tile(aes(fill = correlation))+
    scale_fill_distiller(palette=Rcoloru, limits = c(-1.000001, 1.000001))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(breaks = 1:p, label=colnames(data),expand=c(0,0))+coord_flip()+
    theme(panel.grid=element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.line.x = element_blank()) +xlab(NULL) + ylab(NULL)
  return(uu)

}


