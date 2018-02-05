plot.scale <- function (col, val) {
  rect(14, 0, 15, 1, col="white")
  inc <- 1/(length(col)+2)
  st <- 0
  end <- 0
  for(i in 1:length(col)){
    st <- st+inc
    end <- st+inc
    rect(14.1, st, 14.9, end, col=col[i], border=F)
  }
  st <- st+inc
  text(14.5, 0, min(val))
  text(14.5, st, max(val))
}

deveg <- function(previous.veg.state, rate=1, k=1, step=1) {
  return(previous.veg.state-(step * rate * (1/previous.veg.state*k)))
  #return(previous.veg.state-(rate*(previous.veg.state*k)))
  #return(previous.veg.state-rate)
}

reveg <- function(previous.veg.state, jump=1) {
  return(previous.veg.state+jump)
}

flowering <- function (run, SYM.limit, veg.rate, deveg.rate, 
                       reveg.rate, IM.limit, FM.commit, print=F, jump=T) {
  
  if(print){
    plot(x=NULL, y=NULL, xlim=c(0, run), ylim=c(0, 8), xlab="Plastochrons", ylab="Vegetativeness")
    abline(h=SYM.limit, col='gray', lty=2)
    abline(h=IM.limit, col='gray')
    abline(h=FM.commit, col='gray')
    abline(v=c(0:9), col="gray", lty=3)
    polygon(c(-1, -1, run+1, run+1), c(FM.commit, IM.limit, IM.limit, FM.commit), col="#00000020")
    points(0, 7, pch=1)
    rect(xleft=0, ybottom=0, xright=2.5, ytop=3, col="white")
    text(x=0, y=2.5, labels=paste("deveg.rate:", deveg.rate), cex=0.75, pos=4)
    text(x=0, y=2, labels=paste("reveg.rate:", reveg.rate), cex=0.75, pos=4)
    text(x=0, y=1.5, labels=paste("SYM.limit:", SYM.limit), cex=0.75, pos=4)
    text(x=0, y=1, labels=paste("IM.limit:", IM.limit), cex=0.75, pos=4)
    text(x=0, y=0.5, labels=paste("FM.commit:", FM.commit), cex=0.75, pos=4)  	
  }
  
  
  wall <- F
  crossing <- F
  parent <- 0
  axis <- 1 
  inflorescence <- data.frame()
  inflorescence[1, 1] <- axis
  node <- 1
  step <- 0
  state <- 7
  inflorescence[node, 2] <- step
  inflorescence[node, 3] <- state
  inflorescence[node, 4] <- parent
  colnames(inflorescence) <- c("axis", "plastochron", "state", "parent")
  
  for(i in 1:run) { # steps = plastochrons
    axes <- max(inflorescence$axis)
    new.axis <- axes
    tmp <- data.frame()
    
    for(n in 1:axes) {
      if(n != 2) {
        previous.state <- inflorescence[inflorescence$axis == n & 
                                          inflorescence$plastochron == i-1, 3]
        if(length(previous.state) > 0 && previous.state > 0) {
          if(previous.state > IM.limit) {
            rate <- veg.rate
          } else {
            rate <- deveg.rate
          }
          
          current.state <- deveg(previous.state, rate=rate)
          
          if(current.state < IM.limit & previous.state > IM.limit){
            # state value at the crossing
            slope <-  (current.state - previous.state)
            inc.x <- (IM.limit - previous.state) / slope
            x.new <- (i-1) + inc.x  
            current.state <- deveg(IM.limit, rate=deveg.rate, step = (1 - inc.x))           
            crossing <- T
          }
          
          x0 <- i-1
          y0 <- previous.state
          x1 <- i
          y1 <- current.state
          
          if(i %% 1 == 0) {
            next.state <- reveg(current.state, jump=reveg.rate)
            
            if(current.state < SYM.limit & current.state > FM.commit) {
              parent   <- new.axis
              new.axis <- new.axis + 1
              if(print){
                arrows(x1, y1, x1, next.state, lty=2, length=0.1, col="grey")
              }
              
              type <- 5 # losange vide
              if(next.state > SYM.limit) {
                if(print){ 
                  arrows(x1, next.state, x1, 7, lty=2, length=0.1, col="grey")
                }
                if(jump) {
                  next.state <- 7
                }
                type <- 1 # rond vide
              }
              tmp <- data.frame(rbind(tmp, cbind(axis=new.axis,
                                                 plastochron=i,
                                                 state=next.state,
                                                 parent=parent)))
              
              if(print){
                points(x1, next.state, pch=type)
                #Sys.sleep(0.25)
              }
              
            } 
            if(current.state >= 0) {
              tmp <- data.frame(rbind(tmp,
                                      cbind(axis=n,
                                            plastochron=i,
                                            state=current.state,
                                            parent=0)))
              if(print){
                if(crossing == F){
                  segments(x0, y0, x1, y1)
                }else{
                  segments(x0, y0, x.new, IM.limit)
                  segments(x.new, IM.limit, x1, y1)
                  crossing <- F
                }
                if(current.state > IM.limit) {
                  points(x1, y1, pch=25, bg="green")
                }
              }
            } else {
              tmp <- data.frame(rbind(tmp,
                                      cbind(axis=n,
                                            plastochron=i,
                                            state=current.state,
                                            parent=0)))
              if(print){ 
                slope <- (current.state - previous.state)
                inc.x <- (0 - previous.state) / slope
                x.new <- (i-1) + inc.x  
                segments(x0, y0, x.new, 0)
                points(x.new, 0, pch=8, col="red")
              }
            }
          }
        }
      }
    }  
    inflorescence <- rbind(inflorescence, tmp)
  }
  return(inflorescence)
}

veg.rate      <- 5

run           <- 9
IM.limit      <- 4
FM.commit     <- 2
SYM.limit     <- 5
jump          <- F

deveg.rate    <- 9
reveg.rate    <- 0.2

inflorescence <- flowering(run, SYM.limit, veg.rate, deveg.rate, 
                           reveg.rate, IM.limit, FM.commit, print=T, jump=jump)

#length(inflorescence$state[inflorescence$state < 0])

meta.sim <- data.frame()
meta.node <- 1

for(IM.limit in seq(from=2, to=6, by=0.5)){  
   for(FM.commit in seq(from=0, to=IM.limit, by=0.5)){
  
  SYM.limit <- IM.limit+1
  
  sim <- data.frame()
  node <- 1
  
  for(deveg.rate in seq(from=0, to=20, by=1)){
    for(reveg.rate in seq(from=0, to=5, by=0.5)){
      
      inflorescence <- flowering(run, SYM.limit, veg.rate, deveg.rate, 
                                 reveg.rate, IM.limit, FM.commit, jump=jump)
      
      # Branching
      branching <- 0
      sympodial <- 0
      im <- 0
      sympodial.node <- 99
      flower <- 0
      flo.count <- 0
      
      for(i in 3:max(inflorescence$axis)) {
        tmp <- inflorescence[inflorescence$axis == i, ]
        # Branching
        b <- max(length(tmp$state[tmp$state < SYM.limit & tmp$state > FM.commit]))
        if(b > branching) {
          branching <- b
        }      
      }
      
      sympodial <- length(inflorescence$parent[inflorescence$state > SYM.limit & 
                                                 inflorescence$plastochron > 3 & inflorescence$parent > 0])
      sympodial.node <- min(inflorescence$parent[inflorescence$state > SYM.limit & 
                                                   inflorescence$plastochron > 3 & inflorescence$parent > 0])
      im <- length(inflorescence$parent[inflorescence$state <= SYM.limit & 
                                          inflorescence$plastochron > 3 & inflorescence$parent > 0])
      if(length(inflorescence$state[inflorescence$state < 0]) > 0){
        flower <- 1
        flo.count <- length(inflorescence$state[inflorescence$state < 0])
        if(flo.count == 1) flower <- 2
      }
      
      sim[node, 1] <- deveg.rate
      sim[node, 2] <- reveg.rate
      sim[node, 3] <- branching
      sim[node, 4] <- sympodial
      sim[node, 5] <- sympodial.node
      sim[node, 6] <- flower
      sim[node, 7] <- im
      sim[node, 8] <- flo.count
      
      node <- node + 1
    }
  }
  
  colnames(sim) <- c("deveg", "reveg", "branching", "sympodial", "sympodial.node", "flower", "im", "flo.count")
  sim$flo.count <- as.integer(sim$flo.count)
  
  sim$branched <- 1
  sim$branched[sim$branching >= 3] <- 3
  sim$branched[sim$branching >= 6] <- 4
  sim$branched[sim$branching < 3 & sim$branching > 1] <- 2
  br <- sort(unique(sim$branching))
  col.br <- rainbow(length(br)+1)
  
  col.2 <- rainbow(2)
  
  fl <- sort(unique(sim$flo.count), decreasing=T)
  fl.order <- order(fl)
  for(i in 1:length(fl)) sim$flo.count[sim$flo.count == fl[i]] <- fl.order[i]
  col.flo <- rainbow(length(fl.order)+1)
  
  sim$symped <- 0
  sim$symped[sim$sympodial >= 1] <- 1
  sp <- sort(unique(sim$sympodial))
  
  sim$sympodial.node[sim$sympodial.node == Inf] <- 0
  sim$sympodial.node[sim$sympodial.node == 99] <- 0
  sn <- sort(unique(sim$sympodial.node))
  col.symp <- rainbow(max(sn)+1)
  
  ## Ratio
  
  sim$ratio <- round(sim$sympodial / sim$im)
  sim$ratio[sim$ratio > 14] <- 14
  sim$ratio.class <- 0
  sim$ratio.class[sim$ratio >= 0 & sim$ratio < 0.8] <- 0
  sim$ratio.class[sim$ratio >= 0.8 & sim$ratio < 1.5] <- 1
  sim$ratio.class[sim$ratio >= 1.5] <- 2
  ratio <- sort(unique(sim$ratio), decreasing=T)
  ratio.order <- order(ratio)
  for(i in 1:length(ratio)) sim$ratio[sim$ratio == ratio[i]] <- ratio.order[i]
  col.ratio <- rainbow(length(ratio)+1)
  
  ## Different classes
  tmf        <- 1
  wt         <- 2
  s          <- 3
  sft        <- 8
  an         <- 4
  fa         <- 6
  s.sft      <- 7
  fa.sft      <- 5
  
  sim$type <- 0
  
  sim$type[sim$branched == 2 & sim$symped == 0 & sim$flower == 1 & sim$sympodial.node >= 0] <- wt
  sim$type[sim$branched == 3 & sim$symped == 0 & sim$flower == 1 & sim$sympodial.node >= 0] <- s
  sim$type[sim$branched == 2 & sim$symped == 1 & sim$flower == 1 & sim$sympodial.node >= 0] <- sft
  sim$type[sim$branched == 1 & sim$symped == 0 & sim$flower == 2 & sim$sympodial.node >= 0] <- tmf
  sim$type[sim$branched == 4 & sim$symped == 0 & sim$flower == 0 & sim$sympodial.node >= 0] <- an
  sim$type[sim$branched == 3 & sim$symped == 1 & sim$flower == 1 & sim$sympodial.node >= 0] <- s.sft
  sim$type[sim$branched == 4 & sim$symped == 1 & sim$flower == 0 & sim$sympodial.node >= 0] <- fa
  sim$type[sim$branched == 3 & sim$symped == 1 & sim$flower == 0 & sim$sympodial.node >= 0] <- fa.sft
  
  meta.sim[meta.node, 1] <- IM.limit
  meta.sim[meta.node, 2] <- FM.commit
  meta.sim[meta.node, 3] <- length(unique(sim$type)) 
  
  meta.node <- meta.node + 1
  
  if(length(unique(sim$type)) > 6){
  
    sim$name = ""
    sim$name[sim$branched == 2 & sim$symped == 0 & sim$flower == 1 & sim$sympodial.node >= 0] <- "wt"
    sim$name[sim$branched == 3 & sim$symped == 0 & sim$flower == 1 & sim$sympodial.node >= 0] <- "s"
    sim$name[sim$branched == 2 & sim$symped == 1 & sim$flower == 1 & sim$sympodial.node >= 0] <- "sft"
    sim$name[sim$branched == 1 & sim$symped == 0 & sim$flower == 2 & sim$sympodial.node >= 0] <- "tmf"
    sim$name[sim$branched == 4 & sim$symped == 0 & sim$flower == 0 & sim$sympodial.node >= 0] <- "an"
    sim$name[sim$branched == 3 & sim$symped == 1 & sim$flower == 1 & sim$sympodial.node >= 0] <-" s.sft"
    sim$name[sim$branched == 4 & sim$symped == 1 & sim$flower == 0 & sim$sympodial.node >= 0] <- "fa"
    sim$name[sim$branched == 3 & sim$symped == 1 & sim$flower == 0 & sim$sympodial.node >= 0] <- "fa.sft"
  
  
    x <- sort(unique(sim$deveg))
    y <- sort(unique(sim$reveg))
    color = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), space="Lab")
  
  ## Plot gradients
    png(file=paste("~/Desktop/flomato/flower/IM", IM.limit, "-FM", FM.commit , ".png", sep=""), width=600, height=500)
    m.flo <- matrix(data = sim$flo.count, nrow = length(unique(sim$deveg)), ncol = length(unique(sim$reveg)), byrow=T)
    filled.contour(x, y, m.flo, color = color, plot.axes = { axis(1, cex.axis=2.5); axis(2, cex.axis=2.5)}, key.axes = axis(4, cex.axis=3)) 
    dev.off()
  
    png(file=paste("~/Desktop/flomato/sympodial/IM", IM.limit, "-FM", FM.commit , ".png", sep=""), width=600, height=500)
    m.sympodial <- matrix(data = sim$sympodial, nrow = length(unique(sim$deveg)), ncol = length(unique(sim$reveg)), byrow=T)
    filled.contour(x, y, m.sympodial, color = color, plot.axes = { axis(1, cex.axis=2.5); axis(2, cex.axis=2.5)}, key.axes = axis(4, cex.axis=3)) 
    dev.off()
  
    png(file=paste("~/Desktop/flomato/branching/IM", IM.limit, "-FM", FM.commit , ".png", sep=""), width=600, height=500)
    m.branching <- matrix(data = sim$branching, nrow = length(unique(sim$deveg)), ncol = length(unique(sim$reveg)), byrow=T)
    filled.contour(x, y, m.branching, color = color, plot.axes = { axis(1, cex.axis=2.5); axis(2, cex.axis=2.5)}, 
                 key.axes = axis(4, cex.axis=3), main = paste("veg.rate = ",veg.rate)) 
    dev.off()
  
  ## Plot types
  
      png(file=paste("~/Desktop/flomato/morpho/IM", IM.limit, "-FM", FM.commit , ".png", sep=""), width=600, height=500)
      m.type <- matrix(data = sim$type, nrow = length(unique(sim$deveg)), ncol = length(unique(sim$reveg)), byrow=T)
      color = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"), space="Lab", interpolate="spline", bias=1)
      filled.contour(x, y, m.type, color = color, nlevels=15, plot.axes = {axis(1, cex.axis=2); axis(2, cex.axis=2)}, main = paste("veg.rate = ",leg)) 
      dev.off()
  
      png(file=paste("~/Desktop/flomato/morpho/name_IM", IM.limit, "-FM", FM.commit , ".png", sep=""), width=600, height=500)
      plot(1, 1, ylim=range(y), xlim=range(x), type='n')
      text(sim$deveg, sim$reveg, sim$name, cex=1, col=sim$type)
      dev.off()
    }
  }
}

colnames(meta.sim) <- c("IM", "FM", "morpho")

plot(meta.sim$IM, meta.sim$FM, col=meta.sim$morpho, type='n')
text(meta.sim$IM, meta.sim$FM, meta.sim$morpho)


