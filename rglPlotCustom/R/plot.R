rglplot2        <- function(x, ...)
  UseMethod("rglplot", x)

rglplot.igraph <- function(x, ...) {

  require(rgl)
  
  graph <- x
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  create.edge <- function(v1, v2, r1, r2, ec, ew, am, as) {
    ## these could also be parameters:
    aw <- 0.005*3*as                      # arrow width
    al <- 0.005*4*as                      # arrow length    
    
    dist <- sqrt(sum((v2-v1)^2))   # distance of the centers

    if (am==0) {
      edge <- qmesh3d(c(-ew/2,-ew/2,dist,1, ew/2,-ew/2,dist,1, ew/2,ew/2,dist,1,
                        -ew/2,ew/2,dist,1,  -ew/2,-ew/2,0,1, ew/2,-ew/2,0,1,
                        ew/2,ew/2,0,1, -ew/2,ew/2,0,1),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8))
    } else if (am==1) {
      edge <- qmesh3d(c(-ew/2,-ew/2,dist,1, ew/2,-ew/2,dist,1,
                        ew/2,ew/2,dist,1, -ew/2,ew/2,dist,1,
                        -ew/2,-ew/2,al+r1,1, ew/2,-ew/2,al+r1,1,
                        ew/2,ew/2,al+r1,1, -ew/2,ew/2,al+r1,1,
                        -aw/2,-aw/2,al+r1,1, aw/2,-aw/2,al+r1,1,
                        aw/2,aw/2,al+r1,1, -aw/2,aw/2,al+r1,1, 0,0,r1,1),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8,
                        9,10,11,12, 9,12,13,13, 9,10,13,13, 10,11,13,13,
                        11,12,13,13))
    } else if (am==2) {
      box <- dist-r2-al
      edge <- qmesh3d(c(-ew/2,-ew/2,box,1, ew/2,-ew/2,box,1, ew/2,ew/2,box,1,
                        -ew/2,ew/2,box,1,  -ew/2,-ew/2,0,1, ew/2,-ew/2,0,1,
                        ew/2,ew/2,0,1, -ew/2,ew/2,0,1,
                        -aw/2,-aw/2,box,1, aw/2,-aw/2,box,1, aw/2,aw/2,box,1,
                        -aw/2,aw/2,box,1, 0,0,box+al,1),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8,
                        9,10,11,12, 9,12,13,13, 9,10,13,13, 10,11,13,13,
                        11,12,13,13))
    } else {
      edge <- qmesh3d(c(-ew/2,-ew/2,dist-al-r2,1, ew/2,-ew/2,dist-al-r2,1,
                        ew/2,ew/2,dist-al-r2,1, -ew/2,ew/2,dist-al-r2,1,
                        -ew/2,-ew/2,r1+al,1, ew/2,-ew/2,r1+al,1,
                        ew/2,ew/2,r1+al,1, -ew/2,ew/2,r1+al,1,
                        -aw/2,-aw/2,dist-al-r2,1, aw/2,-aw/2,dist-al-r2,1,
                        aw/2,aw/2,dist-al-r2,1, -aw/2,aw/2,dist-al-r2,1,
                        -aw/2,-aw/2,r1+al,1, aw/2,-aw/2,r1+al,1,
                        aw/2,aw/2,r1+al,1, -aw/2,aw/2,r1+al,1,
                        0,0,dist-r2,1, 0,0,r1,1),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8,
                        9,10,11,12, 9,12,17,17, 9,10,17,17, 10,11,17,17,
                        11,12,17,17,
                        13,14,15,16, 13,16,18,18, 13,14,18,18, 14,15,18,18,
                        15,16,18,18))
    }
      

    ## rotate and shift it to its position
    phi<- -atan2(v2[2]-v1[2],v1[1]-v2[1])-pi/2
    psi<- acos((v2[3]-v1[3])/dist)    
    rot1 <- rbind(c(1,0,0),c(0,cos(psi),sin(psi)), c(0,-sin(psi),cos(psi)))
    rot2 <- rbind(c(cos(phi),sin(phi),0),c(-sin(phi),cos(phi),0), c(0,0,1))
    rot <- rot1 %*% rot2
    edge <- transform3d(edge, rotationMatrix(matrix=rot))
    edge <- transform3d(edge, translationMatrix(v1[1], v1[2], v1[3]))

    ## we are ready 
    shade3d(edge, col=ec)
  }
  
  create.loop <- function(v, r, ec, ew, am, la, la2, as) {
    aw <- 0.005*3*as
    al <- 0.005*4*as
    wi <- aw*2                          # size of the loop
    wi2 <- wi+aw-ew                     # size including the arrow heads
    hi <- al*2+ew*2
    gap <- wi-2*ew

    if (am==0) {
      edge <- qmesh3d(c(-wi/2,-ew/2,0,1, -gap/2,-ew/2,0,1,
                        -gap/2,ew/2,0,1, -wi/2,ew/2,0,1,
                        -wi/2,-ew/2,hi-ew+r,1, -gap/2,-ew/2,hi-ew+r,1,
                        -gap/2,ew/2,hi-ew+r,1, -wi/2,ew/2,hi-ew+r,1,
                        wi/2,-ew/2,0,1, gap/2,-ew/2,0,1,
                        gap/2,ew/2,0,1, wi/2,ew/2,0,1,
                        wi/2,-ew/2,hi-ew+r,1, gap/2,-ew/2,hi-ew+r,1,
                        gap/2,ew/2,hi-ew+r,1, wi/2,ew/2,hi-ew+r,1,
                        -wi/2,-ew/2,hi+r,1, -wi/2,ew/2,hi+r,1,
                        wi/2,-ew/2,hi+r,1, wi/2,ew/2,hi+r,1
                        ),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7,
                        1,4,18,17,
                        9,10,11,12, 13,14,15,16, 9,10,14,13, 10,11,15,14,
                        11,12,16,15, 9,12,20,19,
                        5,13,19,17, 17,18,20,19, 8,16,20,18, 6,7,15,14
                        ))
    } else if (am==1 || am==2) {
      edge <- qmesh3d(c(-wi/2,-ew/2,r+al,1, -gap/2,-ew/2,r+al,1,
                        -gap/2,ew/2,r+al,1, -wi/2,ew/2,r+al,1,
                        -wi/2,-ew/2,hi-ew+r,1, -gap/2,-ew/2,hi-ew+r,1,
                        -gap/2,ew/2,hi-ew+r,1, -wi/2,ew/2,hi-ew+r,1,
                        wi/2,-ew/2,0,1, gap/2,-ew/2,0,1,
                        gap/2,ew/2,0,1, wi/2,ew/2,0,1,
                        wi/2,-ew/2,hi-ew+r,1, gap/2,-ew/2,hi-ew+r,1,
                        gap/2,ew/2,hi-ew+r,1, wi/2,ew/2,hi-ew+r,1,
                        -wi/2,-ew/2,hi+r,1, -wi/2,ew/2,hi+r,1,
                        wi/2,-ew/2,hi+r,1, wi/2,ew/2,hi+r,1,
                        # the arrow
                        -wi2/2,-aw/2,r+al,1, -wi2/2+aw,-aw/2,r+al,1,
                        -wi2/2+aw,aw/2,r+al,1, -wi2/2,aw/2,r+al,1,
                        -wi2/2+aw/2,0,r,1                   
                        ),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7,
                        1,4,18,17,
                        9,10,11,12, 13,14,15,16, 9,10,14,13, 10,11,15,14,
                        11,12,16,15, 9,12,20,19,
                        5,13,19,17, 17,18,20,19, 8,16,20,18, 6,7,15,14,
                        # the arrow
                        21,22,23,24, 21,22,25,25, 22,23,25,25, 23,24,25,25,
                        21,24,25,25
                        ))
    } else if (am==3) {
      edge <- qmesh3d(c(-wi/2,-ew/2,r+al,1, -gap/2,-ew/2,r+al,1,
                        -gap/2,ew/2,r+al,1, -wi/2,ew/2,r+al,1,
                        -wi/2,-ew/2,hi-ew+r,1, -gap/2,-ew/2,hi-ew+r,1,
                        -gap/2,ew/2,hi-ew+r,1, -wi/2,ew/2,hi-ew+r,1,
                        wi/2,-ew/2,r+al,1, gap/2,-ew/2,r+al,1,
                        gap/2,ew/2,r+al,1, wi/2,ew/2,r+al,1,
                        wi/2,-ew/2,hi-ew+r,1, gap/2,-ew/2,hi-ew+r,1,
                        gap/2,ew/2,hi-ew+r,1, wi/2,ew/2,hi-ew+r,1,
                        -wi/2,-ew/2,hi+r,1, -wi/2,ew/2,hi+r,1,
                        wi/2,-ew/2,hi+r,1, wi/2,ew/2,hi+r,1,
                        # the arrows
                        -wi2/2,-aw/2,r+al,1, -wi2/2+aw,-aw/2,r+al,1,
                        -wi2/2+aw,aw/2,r+al,1, -wi2/2,aw/2,r+al,1,
                        -wi2/2+aw/2,0,r,1,
                        wi2/2,-aw/2,r+al,1, wi2/2-aw,-aw/2,r+al,1,
                        wi2/2-aw,aw/2,r+al,1, wi2/2,aw/2,r+al,1,
                        wi2/2-aw/2,0,r,1                   
                        ),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7,
                        1,4,18,17,
                        9,10,11,12, 13,14,15,16, 9,10,14,13, 10,11,15,14,
                        11,12,16,15, 9,12,20,19,
                        5,13,19,17, 17,18,20,19, 8,16,20,18, 6,7,15,14,
                        # the arrows
                        21,22,23,24, 21,22,25,25, 22,23,25,25, 23,24,25,25,
                        21,24,25,25,
                        26,27,28,29, 26,27,30,30, 27,28,30,30, 28,29,30,30,
                        26,29,30,30
                        ))
    }

    # rotate and shift to its position
    rot1 <- rbind(c(1,0,0),c(0,cos(la2),sin(la2)), c(0,-sin(la2),cos(la2)))
    rot2 <- rbind(c(cos(la),sin(la),0),c(-sin(la),cos(la),0), c(0,0,1))
    rot <- rot1 %*% rot2
    edge <- transform3d(edge, rotationMatrix(matrix=rot))
    edge <- transform3d(edge, translationMatrix(v[1], v[2], v[3]))

    ## we are ready
    shade3d(edge, col=ec)
  }
  
  # Visual parameters
  params <- i.parse.plot.params(graph, list(...))
  labels <- params("vertex", "label")
  label.color <- params("vertex", "label.color")
  label.font <- params("vertex", "label.font")
  label.cex <- params("vertex", "label.cex")
  label.degree <- params("vertex", "label.degree")
  label.dist <- params("vertex", "label.dist")
  vertex.color <- params("vertex", "color")
  vertex.size <- (1/200) * params("vertex", "size")
  loop.angle <- params("edge", "loop.angle")
  loop.angle2 <- params("edge", "loop.angle2")

  edge.color <- params("edge", "color")
  edge.width <- (1/200) * params("edge", "width")
  edge.labels <- params("edge","label")
  arrow.mode <- params("edge","arrow.mode")
  arrow.size <- params("edge","arrow.size")
  
  layout <- params("plot", "layout")
  rescale <- params("plot", "rescale")

  # the new style parameters can't do this yet
  arrow.mode         <- i.get.arrow.mode(graph, arrow.mode)
  
  # norm layout to (-1, 1)
  if (ncol(layout)==2) { layout <- cbind(layout, 0) }
  if (rescale) {
    layout <- layout.norm(layout, -1, 1, -1, 1, -1, 1)
  }
  
  # add the edges, the loops are handled separately
  el <- get.edgelist(graph, names=FALSE)
  
  # It is faster this way
  par3d(skipRedraw=TRUE)

  # edges first
  for (i in seq(length=nrow(el))) {
    from <- el[i,1]
    to   <- el[i,2]
    v1 <- layout[from,]
    v2 <- layout[to,]
    am <- arrow.mode; if (length(am)>1) { am <- am[i] }
    ew <- edge.width; if (length(ew)>1) { ew <- ew[i] }
    ec <- edge.color; if (length(ec)>1) { ec <- ec[i] }
    r1 <- vertex.size; if (length(r1)>1) { r1 <- r1[from] }
    r2 <- vertex.size; if (length(r2)>1) { r2 <- r2[to] }

    if (from!=to) {
      create.edge(v1,v2,r1,r2,ec,ew,am,arrow.size)
    } else {
      la <- loop.angle; if (length(la)>1) { la <- la[i] }
      la2 <- loop.angle2; if (length(la2)>1) { la2 <- la2[i] }      
      create.loop(v1,r1,ec,ew,am,la,la2,arrow.size)
    }
    
  }
      
  # add the vertices
  if (length(vertex.size)==1) { vertex.size <- rep(vertex.size, nrow(layout)) }
  rgl.spheres(layout[,1], layout[,2], layout[,3], radius=vertex.size,
              col=vertex.color)

  # add the labels, 'l1' is a stupid workaround of a mysterious rgl bug
  labels[is.na(labels)] <- ""
  x <- layout[,1]+label.dist*cos(-label.degree)* 
    (vertex.size+6*10*log10(nchar(labels)+1))/200
  y <- layout[,2]+label.dist*sin(-label.degree)*
    (vertex.size+6*10*log10(nchar(labels)+1))/200
  z <- layout[,3]
  l1 <- labels[1]
  labels[1] <- ""
  rgl.texts(x,y,z, labels, col=label.color, cex=label.cex, adj=0)
  rgl.texts(c(0,x[1]), c(0,y[1]), c(0,z[1]),
            c("",l1), col=c(label.color[1],label.color[1]), cex=c(label.cex[1],label.cex[1]), adj=0)

  edge.labels[is.na(edge.labels)] <- ""
  if (any(edge.labels != "")) {
    x0 <- layout[,1][el[,1]]
    x1 <- layout[,1][el[,2]]
    y0 <- layout[,2][el[,1]]
    y1 <- layout[,2][el[,2]]
    z0 <- layout[,3][el[,1]]
    z1 <- layout[,4][el[,2]]
    rgl.texts((x0+x1)/2, (y0+y1)/2, (z0+z1)/2, edge.labels,
              col=label.color, cex=label.cex)
  }

  # draw everything
  par3d(skipRedraw=FALSE)
  
  invisible(NULL)
}