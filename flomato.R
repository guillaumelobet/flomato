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
  return(previous.veg.state-(step*rate*(1/previous.veg.state*k)))
  #return(previous.veg.state-(rate*(previous.veg.state*k)))
  #return(previous.veg.state-rate)
}

reveg <- function(previous.veg.state, jump=1) {
  return(previous.veg.state+jump)
}

flowering <- function (run, SYM.limit, veg.rate, deveg.rate, 
                       reveg.rate, IM.limit, FM.commit, print=F, jump=T, maxsim=10) {
  if(print){
    plot(x=NULL, y=NULL, xlim=c(0, 15), ylim=c(0, 8), xlab="Time [plastochrons]", ylab="Meristem vegetativeness[-]", cex.lab = 1.4, cex.axis=1.3)
   # polygon(c(-1, -1, 15+1, 15+1), c(SYM.limit, -1, -1, SYM.limit), col="#00000010", border="white")
    polygon(c(-1, -1, 15+1, 15+1), c(IM.limit, -1, -1, IM.limit), col="#00000010", border="white")
    polygon(c(-1, -1, 15+1, 15+1), c(FM.commit, -1, -1, FM.commit), col="#00000010", border="white")
   # abline(h=SYM.limit, col='darkgray', lwd=2)
    abline(h=IM.limit, col='darkgray', lwd=2)
    abline(h=FM.commit, col='darkgray', lwd=2)
    #text(15+0.4, SYM.limit+0.1, "Leaf repression", adj=c(1,0), col="#515151", cex=1.3)
    text(15+0.4, IM.limit+0.1, "Floral transition", adj=c(1,0), col="#515151", cex=1.3)
    text(15+0.4, FM.commit+0.1, "Floral commitment", adj=c(1,0), col="#515151", cex=1.3)
    points(0, 7, pch=1)
    rect(11, 0.3, 14.9, 2.3, col="white")
     text(x=14.5, y=1.5, labels=paste("dV = ", deveg.rate), cex=1.8, adj=c(1,0), col="#515151")
    text(x=14.5, y=0.5, labels=bquote(Delta*"V"==.(reveg.rate)), cex=1.8, adj=c(1,0), col="#515151")
  }
  
  node <- 1
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
  inflorescence[node, 5] <- node
  colnames(inflorescence) <- c("axis", "plastochron", "state", "parent", "node")
  
  for(i in 1:run) { # steps = plastochrons
    axes <- max(inflorescence$axis)
    new.axis <- axes
    tmp <- data.frame()
    
    for(n in 1:axes) {
      if(n != 2 & n <= maxsim) {
        previous.state <- inflorescence[inflorescence$axis == n & 
                                          inflorescence$plastochron == i-1, 3]
        previous.parent <- inflorescence[inflorescence$axis == n & 
                                           inflorescence$plastochron == i-1, 4]
        previous.node <- inflorescence[inflorescence$axis == n & 
                                         inflorescence$plastochron == i-1, 5]
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
              parent   <- n
              new.axis <- new.axis + 1
              if(print){
                if(n < maxsim-2) segments(x1, y1, x1, next.state, lty=3, lwd=2, col="#515151")
              }
              
              type <- 5 # losange vide
#               if(next.state > SYM.limit) {
# #                 if(print){ 
# #                   #segments(x1, next.state, x1, 7, lty=2, length=0.1, lwd=2, col="#515151")
# #                 }
#                 if(jump) {
#                   next.state <- 7
#                 }
#                 type <- 1 # rond vide
#               }
              if(n < maxsim-2){
                tmp <- data.frame(rbind(tmp, cbind(axis = new.axis,
                                                 plastochron = i,
                                                 state = next.state,
                                                 parent = parent,
                                                 node = 1)))
              }
              if(print){
                if(n < maxsim-2) points(x1, next.state, pch=type)
                if(current.state <= IM.limit){points(x1, current.state, pch=21, col="#515151", bg="#515151", lwd=1)}
                #Sys.sleep(0.25)
              }
              
            } 
            if(current.state >= 0) {
              if(current.state < IM.limit || i == 1){
                tmp <- data.frame(rbind(tmp,
                                        cbind(axis = n,
                                              plastochron = i,
                                              state = current.state,
                                              parent = previous.parent,
                                              node = previous.node + 1)))
              }
              
              if(print){
                if(crossing == F){
                  segments(x0, y0, x1, y1, lwd=2, col="#515151")
                }else{
                  segments(x0, y0, x.new, IM.limit, lwd=2, col="#515151")
                  segments(x.new, IM.limit, x1, y1, lwd=2, col="#515151")
                  crossing <- F
                }
                if(current.state > IM.limit) {
                  #points(x0, y0, pch=21, col="#1F9321", bg="#54B71E", lwd=2, cex=1.2)                  
                  points(x1, y1, pch=18,  col="#1F9321", lwd=2, cex=2)
                }
              }
            } else {
                tmp <- data.frame(rbind(tmp,
                                      cbind(axis = n,
                                            plastochron = i,
                                            state = current.state,
                                            parent = previous.parent,
                                            node = previous.node + 1)))
                if(print){ 
                slope <- (current.state - previous.state)
                inc.x <- (0 - previous.state) / slope
                x.new <- (i-1) + inc.x  
                segments(x0, y0, x.new, 0, lwd=2, col="#515151")
                points(x.new, 0, pch=21, col="red", bg="white", cex=1.6, lwd=1)
                points(x.new, 0, pch=8, col="red", lwd=1)
              }
            }
          }
        }
      }
    }  
    inflorescence <- rbind(inflorescence, tmp)
  }
tmp <- inflorescence[2, ]
points(tmp$plastochron, tmp$state, pch=21, col="#2C4D9E", bg="#1A6BBC", lwd=2, cex=1.2)

axis <- max(inflorescence$axis)
for(i in 1:axis){
  tmp <- inflorescence[inflorescence$axis == i,]
  if(nrow(tmp) > 0){
    tmp2 <- tmp[tmp$node==1,]
    points(tmp2$plastochron, tmp2$state, pch=21, col="#515151", bg="white", cex=1.2, lwd=2)
  } 
}
  return(inflorescence)
}



veg.rate      <- 5

run           <- 9
IM.limit      <- 6
FM.commit     <- 4.2
SYM.limit     <- IM.limit+1
jump          <- F


pdf(file="~/Desktop/simulation_single_mutants.pdf", width=10, height=12)
par(mfrow=c(4,2), mar=c(2,4,3,2), oma=c(4,4,0,4))
# Wild type
inflorescence <- flowering(run, SYM.limit,veg.rate, 8, 0.7, IM.limit, FM.commit, print=T, jump=jump)
mtext("wild-type", line=0.5, cex=1.2)

# s
inflorescence <- flowering(run, SYM.limit, veg.rate, 5, 0.7, IM.limit, FM.commit, print=T, jump=jump)
mtext("s", line=0.5, cex=1.2)

# an
inflorescence <- flowering(run, SYM.limit, veg.rate, 1, 0.7, IM.limit, FM.commit, print=T, jump=jump)
mtext("an", line=0.5, cex=1.2)

# fa
inflorescence <- flowering(run, SYM.limit, veg.rate, 1, 1, IM.limit, FM.commit, print=T, jump=jump)
mtext("fa", line=0.5, cex=1.2)


# sft
inflorescence <- flowering(run, SYM.limit, veg.rate, 10, 1.5, IM.limit, FM.commit, print=T, jump=jump)
mtext("sft", line=0.5, cex=1.2)

# j
inflorescence <- flowering(run, SYM.limit, veg.rate, 14, 1.8, IM.limit, FM.commit, print=T, jump=jump)
mtext("j", line=0.5, cex=1.2)
mtext("Time [plastochron]", side=1, line=3, cex=0.9)

# 
inflorescence <- flowering(run, SYM.limit, veg.rate, 18, 0.7, IM.limit, FM.commit, print=T, jump=jump)
mtext("tmf", line=0.5, cex=1.2)
mtext("Time [plastochron]", side=1, line=3, cex=0.9)

dev.off()




## Double mutants


pdf(file="~/Desktop/simulatio_double_mutants.pdf", width=10, height=12)
par(mfrow=c(4,2), mar=c(2,4,3,2), oma=c(4,4,0,4))

# Wild type
inflorescence <- flowering(run, SYM.limit,veg.rate, 8, 0.7, IM.limit, FM.commit, print=T, jump=jump)
mtext("wild-type", line=0.5, cex=1.2)

# j:fa
inflorescence <- flowering(run, SYM.limit, veg.rate, 1, 2.5, IM.limit, FM.commit, print=T, jump=jump)
mtext("j fa", line=0.5, cex=1.2)

# j:an
inflorescence <- flowering(run, SYM.limit, veg.rate, 1, 2, IM.limit, FM.commit, print=T, jump=jump)
mtext("j an", line=0.5, cex=1.2)

# j:s
inflorescence <- flowering(run, SYM.limit, veg.rate, 10, 1.8, IM.limit, FM.commit, print=T, jump=jump)
mtext("j s", line=0.5, cex=1.2)

# j:sft
inflorescence <- flowering(run, SYM.limit, veg.rate, 14, 3, IM.limit, FM.commit, print=T, jump=jump)
mtext("j sft", line=0.5, cex=1.2)

# sft:fa
inflorescence <- flowering(run, SYM.limit, veg.rate, 1, 2, IM.limit, FM.commit, print=T, jump=jump)
mtext("sft fa", line=0.5, cex=1.2)

# sft:an
inflorescence <- flowering(run, SYM.limit, veg.rate, 1, 1.5, IM.limit, FM.commit, print=T, jump=jump)
mtext("sft an", line=0.5, cex=1.2)
mtext("Time [plastochron]", side=1, line=3, cex=0.9)

# s:sft
inflorescence <- flowering(run, SYM.limit, veg.rate, 5, 1.5, IM.limit, FM.commit, print=T, jump=jump)
mtext("s sft", line=0.5, cex=1.2)
mtext("Time [plastochron]", side=1, line=3, cex=0.9)

dev.off()




