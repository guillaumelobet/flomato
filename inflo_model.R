deveg <- function(previous.veg.state, rate=1, k=1) {
  return(previous.veg.state-(rate*1/(previous.veg.state*k)))
}
reveg <- function(previous.veg.state, jump=1) {
  return(previous.veg.state+jump)
}
run           <- 10
veg.rate      <- 5

deveg.rate    <- 6
reveg.rate    <- 1

bump.interval <- 1
IM.limit      <- 4
FM.commit     <- 2
SYM.limit     <- 5

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

plot(x=NULL, y=NULL, xlim=c(0, run), ylim=c(0, 8), xlab="Plastochrons", ylab="Vegetativeness")
abline(h=SYM.limit, col='gray', lty=2)
abline(h=IM.limit, col='gray')
abline(h=FM.commit, col='gray')
points(0, 7, pch=1)
text(x=8, y=7.5, labels=paste("deveg.rate:", deveg.rate), cex=0.75, pos=4)
text(x=8, y=7, labels=paste("reveg.rate:", reveg.rate), cex=0.75, pos=4)
text(x=8, y=6.5, labels=paste("SYM.limit:", SYM.limit), cex=0.75, pos=4)
text(x=8, y=6, labels=paste("IM.limit:", IM.limit), cex=0.75, pos=4)
text(x=8, y=5.5, labels=paste("FM.commit:", FM.commit), cex=0.75, pos=4)

for(i in 1:run) { # steps = plastochrons
  print(i)
  axes <- max(inflorescence$axis)
  new.axis <- axes
  tmp <- data.frame()
  
  for(n in 1:axes) {
    if(n != 2) {
      previous.state <- inflorescence[inflorescence$axis == n & 
                                        inflorescence$plastochron == i-1, 3]
      if(length(previous.state) > 0) {
        if(previous.state > IM.limit) {
          rate <- veg.rate
        } else {
          rate <- deveg.rate
        }
        current.state <- deveg(previous.state, rate=rate)
        
        x0 <- i-1
        y0 <- previous.state
        x1 <- i
        y1 <- current.state
        
        if(i %% bump.interval == 0) {
          next.state <- reveg(current.state, jump=reveg.rate)
          if(current.state < SYM.limit & current.state > FM.commit) {
            parent   <- new.axis
            new.axis <- new.axis + 1
            arrows(x1, y1, x1, next.state, lty=2, length=0.1, col="grey")
            type <- 5 # losange vide
            if(next.state > SYM.limit) {
              arrows(x1, next.state, x1, 7, lty=2, length=0.1, col="grey")
              next.state <- 7
              type <- 1 # rond vide
            }
            tmp <- data.frame(rbind(tmp, cbind(axis=new.axis,
                                               plastochron=i,
                                               state=next.state,
                                               parent=parent)))
            
#             if(next.state > IM.limit) {
#               points(x1, next.state, pch=25, bg="green")
#             } else {
              points(x1, next.state, pch=type)
#             }
            segments(x0, y0, x1, y1)
            Sys.sleep(0.25)
            
          } 
          if(current.state >= 0) {
            tmp <- data.frame(rbind(tmp,
                                    cbind(axis=n,
                                          plastochron=i,
                                          state=current.state,
                                          parent=0)))
            segments(x0, y0, x1, y1)
            if(current.state > IM.limit) {
              points(x1, y1, pch=25, bg="green")
            }
          } else {
            #         tmp <- data.frame(rbind(tmp,
            #                                 cbind(axis=n,
            #                                       plastochron=i,
            #                                       state=0)))
            #segments(x0, y0, x1, 0)
            points(x0, y0, pch=8, col="red")
          }
        }
      }
    }
  }  
  inflorescence <- rbind(inflorescence, tmp)
}


# plot(x=NULL, y=NULL, xlim=c(0,i), ylim=c(0,8))
# abline(h=3.5, col='gray')
# abline(h=1.5, col='gray')
# for(j in 1:max(inflorescence$axis)) {
#   x <- inflorescence[inflorescence$axis == j, 2]
#   y <- inflorescence[inflorescence$axis == j, 3]
#   lines(x=x, y=y, type="o")
#   Sys.sleep(1)
# }