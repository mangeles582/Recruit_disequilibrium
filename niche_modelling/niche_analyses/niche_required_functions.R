
#### occurrence filtering required functions 

ReduceSpatialClustering = function(data, minimum.distance){
  
  #count rows
  row<-1
  
  
  #repite la operación hasta que se cumple la condición de salida
  repeat{
    
    #contenido de la fila (para no tirar de toda la tabla en todas las operaciones)
    f<-data[row, ]
    
    #genera los l????mites de la cuadr????cula de búsqueda
    ymax<-f$lat + minimum.distance
    ymin<-f$lat - minimum.distance
    xmax<-f$lon + minimum.distance
    xmin<-f$lon - minimum.distance
    
    #selecciona de la tabla los datos con coordenadas dentro del rectángulo que no tienen las mismas coordenadas que la fila con la que estamos trabajando, y las elimina de la tabla
    data<-data[!((data$lat <= ymax) & (data$lat >= ymin) & (data$lon <= xmax) & (data$lon >= xmin) & (data$lat != f$lat | data$lon != f$lon)), ]
    
    #estima de filas por procesar
    print(paste("Processed rows: ", row, " out of ", nrow(data), sep=""))
    
    #suma 1 al contador de la fila
    row<-row+1
    
    #condición de salida cuando llega a la última fila
    if(row>=nrow(data))break
  }
  
  return(data)
  
}#reduce point density function
ecospat.mantel.correlogram1 <- function (dfvar, colxy, n, colvar, max, nclass, nperm) 
{ 
  require(ecodist)
  envnorm <- data.frame(t((t(dfvar[, colvar]) - apply(dfvar[, 
                                                            colvar], 2, mean))/apply(dfvar[, colvar], 2, sd)))
  row.rand <- sample(1:nrow(dfvar), n, replace = TRUE)
  envdist <- dist(envnorm[row.rand, ])
  geodist <- dist(dfvar[row.rand, colxy])
  b <- seq(from = min(geodist), to = max, length.out = nclass)
  crlg <- mgram(envdist, geodist, breaks = b, nperm = nperm)
  plot(crlg)
  abline(h = 0)
  return(crlg)
}

#### plot occurrences

plot_map_occ <- function( shapefile, x_ext=c(-35, 70), 
                          y_ext=c(18,70), df1=NULL,df2=NULL,df3=NULL, 
                          title=NULL, point_size=0.8, title_size=NULL, 
                          axis_title_size=NULL, axis_text_size=NULL){
                           
  gg1 <- ggplot()+
    geom_polygon(data=shapefile, aes(x=long, y = lat, group = group), 
                 fill='grey', size=.2, color='grey68')+
    coord_sf(crs =4326, xlim = x_ext,ylim = y_ext, expand = FALSE)+
    ggtitle(title)+
    scale_y_continuous(breaks=seq(y_ext[1],y_ext[2],(y_ext[2]-y_ext[1])/5))+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5, size= title_size),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(color="black", size=axis_title_size),
          axis.text=element_text(color="black" , size=axis_text_size),
          text=element_text(family="serif"))+
    NULL
  
  if (!is.null(df1)){
    gg1 <- gg1 +
      geom_point(data = df1, aes(lon, lat), col="black", shape=19, size=point_size)
    }
  if (!is.null(df2)){
    gg1 <- gg1 +
      geom_point(data = df2, aes(lon, lat), col="red", shape=1, size=point_size)
  }
  if (!is.null(df3)){
    gg1 <- gg1 +
      geom_point(data = df3, aes(lon, lat), col="green", shape=1, size=point_size)
  }
  return(gg1)
}

#### niche analysis required functions


#' Centre of Gravity or Mass calculations for spatial data
#' 
#' \code{COGravity} calculates the Centre of Gravity (or also known as Centre
#' of Mass) for point or raster spatial data.\cr \cr \bold{Note:} NA data is
#' automatically ommitted from analysis.
#' 
#' For raster-based data, if \code{wt} is missing, the values of the ascii are
#' assumed to be the weights; otherwise, the values are assumed to be the
#' \code{z} values.
#' 
#' @param x a vector of e.g., longitudes or eastings, or a raster of class
#' 'asc', 'RasterLayer' or 'SpatialGridDataFrame'.
#' @param y a vector of e.g., latitudes or northings.
#' @param z a vector of e.g., elevations.
#' @param wt a vector or raster of class 'asc', 'RasterLayer' or
#' 'SpatialGridDataFrame' representing weights for data.
#' @return Returns a named vector of data representing the Centre of Gravity in
#' x, y & z dimensions (depending on data supplied).
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' 
#' 
#' #create some points
#' x = seq(154,110,length=25)
#' y = seq(-10,-54,length=25)
#' z = seq(100,200,length=25)
#' wt = runif(25) #random weights
#' #calculate the Centre of Gravity for these points
#' COGravity(x,y,z,wt)
#' 
#' #create a simple objects of class 'asc'
#' x = as.asc(matrix(1:50,nr=50,nc=50))
#' wt = as.asc(matrix(runif(50),nr=50,nc=50))
#' 
#' #calculate COG with weighting defined in x 
#' COGravity(x)
#' #calculate COG with weighting defined in wt (values in x are assumed elevation (z)) 
#' COGravity(x,wt=wt)
#' 
#' 
#' @export COGravity
COGravity <- function(x,y=NULL,z=NULL,wt=NULL) {
  #check if raster from sp or raster package and convert if necessary
  if (any(class(x) %in% 'RasterLayer')) x = asc.from.raster(x)
  if (any(class(x) == 'SpatialGridDataFrame')) x = asc.from.sp(x)
  #check if x is vector or matrix
  if (is.vector(x)) { #if the data is a vector...do calculations
    if (is.null(wt)) {	#if no weighting supplied, calculate means & standard deviations
      out = c(COGx=mean(x,na.rm=TRUE),COGx.sd=sd(x,na.rm=TRUE))
      if (!is.null(y)) out = c(out,COGy=mean(y,na.rm=TRUE),COGy.sd=sd(y,na.rm=TRUE))
      if (!is.null(z)) out = c(out,COGz=mean(z,na.rm=TRUE),COGz.sd=sd(z,na.rm=TRUE))
    } else { #if weighting supplied, calculate weighted means and variances to get COG
      out = c(COGx=wt.mean(x,wt),COGx.sd=wt.sd(x,wt))
      if (!is.null(y)) out = c(out,COGy=wt.mean(y,wt),COGy.sd=wt.sd(y,wt))
      if (!is.null(z)) out = c(out,COGz=wt.mean(z,wt),COGz.sd=wt.sd(z,wt))
    }
  } else if (any(class(x) == 'asc')) { #if x is of class 'asc'
    if (is.null(wt)) { #if wt is null then assume that values in x are the weights
      pos = as.data.frame(which(is.finite(x),arr.ind=TRUE))
      pos$x = getXYcoords(x)$x[pos$row]
      pos$y = getXYcoords(x)$y[pos$col]
      pos$wt = x[cbind(pos$row,pos$col)]
      out = c(COGx=wt.mean(pos$x,pos$wt),COGx.sd=wt.sd(pos$x,pos$wt),COGy=wt.mean(pos$y,pos$wt),COGy.sd=wt.sd(pos$y,pos$wt))
    } else { #if wt is supplied, it must be of the same dim as x and then the values of x are assumed to be your z
      if (!all(dim(x)==dim(wt))) stop('the grids for x & weights must be of the same dimensions')
      pos = as.data.frame(which(is.finite(x),arr.ind=TRUE))
      pos$x = getXYcoords(x)$x[pos$row]
      pos$y = getXYcoords(x)$y[pos$col]
      pos$z = x[cbind(pos$row,pos$col)]
      pos$wt = wt[cbind(pos$row,pos$col)]
      out = c(COGx=wt.mean(pos$x,pos$wt),COGx.sd=wt.sd(pos$x,pos$wt),COGy=wt.mean(pos$y,pos$wt),COGy.sd=wt.sd(pos$y,pos$wt),COGz=wt.mean(pos$z,pos$wt),COGz.sd=wt.sd(pos$z,pos$wt))
    }
    
  }
  # return the output	
  return(out)
}



#' Weighted mean, variance and standard deviation calculations
#' 
#' \code{wt.mean} calculates the mean given a weighting of the values. \cr \cr
#' \code{wt.var} is the unbiased variance of the weighted mean calculation
#' using equations of GNU Scentific Library
#' (\url{http://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.htmland}.\cr\cr
#' \code{wt.sd} is the standard deviation of the weighted mean calculated as
#' the sqrt of \code{wt.var}. \cr \cr \bold{Note:} NA data is automatically
#' ommitted from analysis.
#' 
#' 
#' @param x is a vector of numerical data.
#' @param wt is a vector of equal length to \code{x} representing the weights.)
#' @return returns a single value from analysis requested.
#' @author Jeremy VanDerWal \email{jjvanderwal@@gmail.com}
#' @examples
#' #define simple data
#' x = 1:25 # set of numbers
#' wt = runif(25) #some arbitrary weights
#' 
#' #display means & variances (unweighted and then weighted)
#' mean(x); wt.mean(x,wt)
#' var(x); wt.var(x,wt)
#' sd(x); wt.sd(x,wt)
#' 
#' 
#' @export 
wt.mean <- function(x,wt) {
  s = which(is.finite(x*wt)); wt = wt[s]; x = x[s] #remove NA info
  return( sum(wt * x)/sum(wt) ) #return the mean
}

#' @rdname wt.mean
#' @export
wt.var <- function(x,wt) {
  s = which(is.finite(x + wt)); wt = wt[s]; x = x[s] #remove NA info
  xbar = wt.mean(x,wt) #get the weighted mean
  return( sum(wt *(x-xbar)^2)*(sum(wt)/(sum(wt)^2-sum(wt^2))) ) #return the variance
} 

#' @rdname wt.mean
#' @export
wt.sd <- function(x,wt) { 
  return( sqrt(wt.var(x,wt)) ) #return the standard deviation
} 




meteo_sim <- function(length = 5000, mean_t = NULL, sd_t = NULL, 
                      mean_p = NULL, sd_p= NULL, names = 'specie'){
  if (any(length(mean_t) != length(sd_t),
          length(mean_p) != length(sd_t),
          length(mean_p) != length(sd_p))){stop('vector length must be equal')}
  n <- length(mean_t)
  res <- list(NULL)
  for ( i in 1:n ){
    name <- paste0(names,'_',i)
    x <- rnorm(length, mean= mean_t[i], sd = sd_t[i])
    y <- rnorm(length, mean= mean_p[i], sd = sd_p[i])
    res[[i]] <- data_frame(specie = name,temp=x,precip=y)
  }
  # res2 <- bind_rows(res)
  return (res)
}

GMINdistance <- function(origin_coord_x, origin_coord_y, end_coord_x, end_coord_y){
  end.df <- data.frame(cbind(end_coord_x, end_coord_y))
  dist.df <- data.frame(cbind(origin_coord_x, origin_coord_y))
  names(dist.df) <- c("x", "y")
  for( i in 1: nrow(end.df)){
    for (k in 1: nrow(dist.df)){
      dist.df[k, "dist"] <- sqrt((dist.df[k,"x"] - end_coord_x[i])^2+ (dist.df[k, "y"] - end_coord_y [i])^2)
    }
    min.dist <- min(dist.df$dist)
    end.df[i, "min.dist"] <- min.dist
    end.df[i, "closest.x"] <- unique(dist.df[dist.df$dist==min.dist, "x"])
    end.df[i, "closest.y"] <- unique(dist.df[dist.df$dist==min.dist, "y"])
    
  }
  
  result <- end.df[, c("min.dist", "closest.x", "closest.y")]
  #names(result) <- c("x_coord", "y_coord", "min_dist")# to show names of each element of the vector
  return(result)
  
}# corregir
# end_coord is the coordinate that iteratively will measures the distance with origin_coord


GMAXdistance <- function(origin_coord_x, origin_coord_y, end_coord_x, end_coord_y){
  end.df <- data.frame(cbind(end_coord_x, end_coord_y))
  dist.df <- data.frame(cbind(origin_coord_x, origin_coord_y))
  names(dist.df) <- c("x", "y")
  for( i in 1: nrow(end.df)){
    for (k in 1: nrow(dist.df)){
      dist.df[k, "dist"] <- sqrt((dist.df[k,"x"] - end_coord_x[i])^2+ (dist.df[k, "y"] - end_coord_y [i])^2)
    }
    max.dist <- max(dist.df$dist)
    end.df[i, "max.dist"] <- max.dist
    end.df[i, "farthest.x"] <- unique(dist.df[dist.df$dist==max.dist, "x"])
    end.df[i, "farthest.y"] <- unique(dist.df[dist.df$dist==max.dist, "y"])
    
  }
  
  result <- end.df[, c("max.dist", "farthest.x", "farthest.y")]
  #names(result) <- c("x_coord", "y_coord", "min_dist")# to show names of each element of the vector
  return(result)
  
}


GMINdistance.univariate <- function(origin_coord_x, end_coord_x){
  end.df <- data.frame(end_coord_x)
  dist.df <- data.frame(origin_coord_x)
  names(dist.df) <- c("x")
  min.df <- data.frame()
  
  for(z in 1:nrow(end.df)){
    n_ref <- end.df[z, 1]
    dist.df[, 1+z] <- abs(dist.df[, 1]-n_ref)
    min.df0 <- dist.df[dist.df[,1+z]==min(dist.df[,1+z]),]
    min.df0 <- distinct(min.df0)
    names(min.df0) <- c("x", "dist")
    min.df <- rbind(min.df, min.df0)
  }
  
  minimum.dist <- min(min.df$dist)
  result <- min.df%>%
    dplyr::filter(dist==minimum.dist)%>%
    rename(closest.x=x,
           min.dist=dist)
  
  return(result)
  
}# corregir
# end_coord is the coordinate that iteratively will measures the distance with origin_coord


GMAXdistance.univariate <- function(origin_coord_x, end_coord_x){
  end.df <- data.frame(end_coord_x)
  dist.df <- data.frame(origin_coord_x)
  names(dist.df) <- c("x")
  max.df <- data.frame()
  
  for(z in 1:nrow(end.df)){
    n_ref <- end.df[z, 1]
    dist.df[, 1+z] <- abs(dist.df[, 1]-n_ref)
    max.df0 <- dist.df[dist.df[,1+z]==max(dist.df[,1+z]),]
    max.df0 <- distinct(max.df0)
    names(max.df0) <- c("x", "dist")
    max.df <- rbind(max.df, max.df0)
  }
  
  maximum.dist <- max(max.df$dist)
  result <- max.df%>%
    dplyr::filter(dist==maximum.dist)%>%
    rename(furthest.x=x,
           max.dist=dist)
  
  return(result)
  
}


intra.limit.univariate0 <- function(table = df, step= NULL){# table format df[, c(x, density)]
  
  if(nrow(table)>3){
  
  i0 <- min(table[, 1])
  max_row <- nrow(table)-1
  table_sub <- table[2:max_row, ]
  table_sub[1, "continuity"] <- if(abs((table_sub[1, 1]-i0)-step)<step/3){ 
    TRUE}else{
      FALSE
    }
  
  for (i in 2:nrow(table_sub)){
    
    table_sub[i, "continuity"] <- if(abs((table_sub[i, 1]-table_sub[i-1, 1])-step)<step/3){ 
      TRUE}else{
        FALSE
      }
    
    if(table_sub[i, "continuity"]==FALSE){
      table_sub[i-1, "continuity"] <- FALSE
    }else{
      TRUE }
  
    }
  
  }else{table_sub <- table
        table_sub[, "continuity"] <- TRUE
  }
 
  return(table_sub)  

}
  



intra.limit.univariate <- function(table = df, step= NULL){# table format df[, c(x, density)]
  
  if(nrow(table)>3){
    
    i0 <- min(table[, 1])
    max_row <- nrow(table)-1
    table_sub <- table[2:max_row, ]
    
    if(abs((table_sub[1, 1]-i0)-step)<step/3){ 
      table_sub[1, "continuity"] <- TRUE}else{
        table_sub[1, "continuity"] <-FALSE
      }
    
    if(table_sub[1, "continuity"]==F){ 
      table_sub[1, "position"] <- "end"}else{
        
      }
    
    for (i in 2:nrow(table_sub)){
      
       if(abs((table_sub[i, 1]-table_sub[i-1, 1])-step)<step/3){ 
         table_sub[i, "continuity"] <-TRUE}else{
           table_sub[i, "continuity"] <- FALSE
        }
    }
    
    for (i in 2:nrow(table_sub)){
      if(table_sub[i, "continuity"]==FALSE){
        table_sub[i, "position"] <- "start"
       }else{ 
        table_sub[i, "position"] <- NA }
      
    }
      
    for (i in 2:nrow(table_sub)){
      if(table_sub[i, "continuity"]==FALSE){
        table_sub[i-1, "continuity"] <- FALSE
        table_sub[i-1, "position"] <- if(table_sub[i-1, "position"]=="start"){
          "peak"}else{"end"}
      
        }
      
   
      }
    
  }else{table_sub <- table
  table_sub[, "continuity"] <- TRUE
  table_sub[, "position"] <- NA
  }
  
  return(table_sub)  
  
}


line_ecuation <-  function(x1, x2, y1, y2){
  m <- (y2-y1)/(x2-x1)
  n <-  y1 - (m*x1)
  result <- c(m, n)
  names(result) <- c("m","n")
  return(result)
}

raster_probability <- function (dudi_pca, rasters, niche_raster, ...) 
{
  z <- getValues(rasters[[1]])
  rasters <- if(length(z[which(!is.na(z))])>4000000){
    aggregate(rasters, fact=5, fun=mean)
  } else {
    rasters
  }
  raster_values <- data.frame(raster::getValues(rasters))
  raster_li <- suprow(dudi_pca, raster_values)$li
  projected_suit <- raster::extract(niche_raster, raster_li[, c(1:2)])
  raster_out = rasters[[1]]
  raster::values(raster_out) <- projected_suit
  return(raster_out)
}


## mode for density data.frame

mode_value <- function (density_df, density_col, ...) 
  {z <- which(names(density_df)%in%density_col[1])
   mode_value <- density_df[which.max(density_df[,z]),]
   return(mode_value)
  }


#' Adds a new scale to a plot
#'
#' Creates a new scale "slot". Geoms added to a plot after this function will
#' use a new scale definition.
#'
#' @param new_aes A string with the name of the aesthetic for which a new scale
#' swill be created.
#'
#' @details
#' `new_scale_color()`, `new_scale_colour()` and `new_scale_fill()` are just
#' aliases to `new_scale("color")`, etc...
#'
#' @examples
#' library(ggplot2)
#'
#' # Equivalent to melt(volcano), but we don't want to depend on reshape2
#' topography <- expand.grid(x = 1:nrow(volcano),
#'                           y = 1:ncol(volcano))
#' topography$z <- c(volcano)
#'
#' # point measurements of something at a few locations
#' measurements <- data.frame(x = runif(30, 1, 80),
#'                            y = runif(30, 1, 60),
#'                            thing = rnorm(30))
#'
#' ggplot(mapping = aes(x, y)) +
#'   geom_contour(data = topography, aes(z = z, color = stat(level))) +
#'   # Color scale for topography
#'   scale_color_viridis_c(option = "D") +
#'   # geoms below will use another color scale
#'   new_scale_color() +
#'   geom_point(data = measurements, size = 3, aes(color = thing)) +
#'   # Color scale applied to geoms added after new_scale_color()
#'   scale_color_viridis_c(option = "A")
#'
#' @export
new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

#' @export
#' @rdname new_scale
new_scale_fill <- function() {
  new_scale("fill")
}

#' @export
#' @rdname new_scale
new_scale_color <- function() {
  new_scale("colour")
}

#' @export
#' @rdname new_scale
new_scale_colour <- function() {
  new_scale("colour")
}

#' @export
#' @importFrom ggplot2 ggplot_add
ggplot_add.new_aes <- function(object, plot, object_name) {
  plot$layers <- bump_aes_layers(plot$layers, new_aes = object)
  plot$scales$scales <- bump_aes_scales(plot$scales$scales, new_aes = object)
  plot$labels <- bump_aes_labels(plot$labels, new_aes = object)
  plot
}

bump_aes_layers <- function(layers, new_aes) {
  lapply(layers, bump_aes_layer, new_aes = new_aes)
  
}

bump_aes_layer <- function(layer, new_aes) {
  original_aes <- new_aes
  
  old_aes <- names(layer$mapping)[remove_new(names(layer$mapping)) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  old_geom <- layer$geom
  
  old_setup <- old_geom$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup(data, params)
  }
  
  new_geom <- ggplot2::ggproto(paste0("New", class(old_geom)[1]), old_geom,
                               handle_na = new_setup)
  
  new_geom$default_aes <- change_name(new_geom$default_aes, old_aes, new_aes)
  new_geom$non_missing_aes <- change_name(new_geom$non_missing_aes, old_aes, new_aes)
  new_geom$required_aes <- change_name(new_geom$required_aes, old_aes, new_aes)
  new_geom$optional_aes <- change_name(new_geom$optional_aes, old_aes, new_aes)
  
  layer$geom <- new_geom
  
  old_stat <- layer$stat
  
  old_setup2 <- old_stat$handle_na
  new_setup <- function(self, data, params) {
    colnames(data)[colnames(data) %in% new_aes] <- original_aes
    old_setup2(data, params)
  }
  
  new_stat <- ggplot2::ggproto(paste0("New", class(old_stat)[1]), old_stat,
                               handle_na = new_setup)
  
  new_stat$default_aes <- change_name(new_stat$default_aes, old_aes, new_aes)
  new_stat$non_missing_aes <- change_name(new_stat$non_missing_aes, old_aes, new_aes)
  new_stat$required_aes <- change_name(new_stat$required_aes, old_aes, new_aes)
  new_stat$optional_aes <- change_name(new_stat$optional_aes, old_aes, new_aes)
  
  layer$stat <- new_stat
  
  layer$mapping <- change_name(layer$mapping, old_aes, new_aes)
  layer
}

bump_aes_scales <- function(scales, new_aes) {
  lapply(scales, bump_aes_scale, new_aes = new_aes)
}

bump_aes_scale <- function(scale, new_aes) {
  old_aes <- scale$aesthetics[remove_new(scale$aesthetics) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  scale$aesthetics[scale$aesthetics %in% old_aes] <- new_aes
  
  if (is.character(scale$guide)) {
    scale$guide <- match.fun(paste("guide_", scale$guide, sep = ""))()
  }
  scale$guide$available_aes[scale$guide$available_aes %in% old_aes] <- new_aes
  scale
}

bump_aes_labels <- function(labels, new_aes) {
  old_aes <-  names(labels)[remove_new(names(labels)) %in% new_aes]
  new_aes <- paste0(old_aes, "_new")
  
  names(labels)[names(labels) %in% old_aes] <- new_aes
  labels
}


change_name <- function(list, old, new) {
  UseMethod("change_name")
}

change_name.character <- function(list, old, new) {
  list[list %in% old] <- new
  list
}

change_name.default <- function(list, old, new) {
  nam <- names(list)
  nam[nam %in% old] <- new
  names(list) <- nam
  list
}

change_name.NULL <- function(list, old, new) {
  NULL
}


remove_new <- function(aes) {
  stringi::stri_replace_all(aes, "", regex = "(_new)*")
}



#standard error

se <- function(x) sd(x)/sqrt(length(x))

