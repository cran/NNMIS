# basic function

NNMI <- function (y, xa=NULL, xb=NULL, time, event, MI = 10, NN = 5, w1 = 0.8, w2 = 0.2,
                  imputeCT = FALSE, NN.t=10, mfamily = NULL, datalevels = NULL, xc = NULL, Seed=NA)
{
  #require(survival)

  #------------------------------------------------------------------------------------------------------------------------
  if (mfamily != "gaussian") {
    y <- as.factor(y)
  }
  N <- length(time)
  xa <- as.matrix(xa); xb <- as.matrix(xb)
  p <- ncol(xa); q <- ncol(xb)
  #-------------------------------------------------------------------------------------------------------------------------

  # observation id
  ID <- seq(1:N)

  # observed survival time
  oT <- time

  # event indicator
  I.d <- event

  # log observed survival time
  lt <- log(oT)

  # create cumulative hazard function
  fit.s <- survival::survfit(survival::Surv(oT,I.d)~1)
  NA.H <- cumsum(fit.s$n.event/fit.s$n.risk)
  Ht <- rep(NA,N)
  for(i in 1:N){
    t.id <- max(fit.s$time[fit.s$time <= oT[i]])
    Ht[i] <- NA.H[fit.s$time == t.id]
  }

  # Missing indicator
  M1 <- rep(0,N)
  M1[which(is.na(y))] <- 1
  ID.1 <- ID[M1==1]            #all observation id with missing value
  ID.m <- sort(unique(ID.1))
  N.m <- length(ID.m)          #number of missing values

  # censering rate
  crate <- (1-mean(I.d))
  # missing rate
  m1rate <- mean(M1)

  # output preparation
  dat.NNMI <- dat.T.NNMI.t <- dat.Id.NNMI.t <- matrix(NA, nrow=N, ncol=MI)
  colnames(dat.NNMI) <- colnames(dat.T.NNMI.t) <- colnames(dat.Id.NNMI.t) <- paste("M",1:MI,sep="")

  # intercept
  X0 <- rep(1,N)

  # imputation target
  X1.o <- y
  for (i in 1:MI) {
    dat.NNMI[,i] <- X1.o
  }

  XXa <- as.matrix(cbind(X0,xa,Ht,I.d))
  XXb <- as.matrix(cbind(X0,xb,Ht,I.d))

  if(!is.na(Seed))
  set.seed(Seed)

  for(k in 1:MI){

    ID.B <- sample(seq(1:N), size=N, replace=T)          #bootstrap
    M1.B <- 1-M1[ID.B]
    xa.B <- xa[ID.B,]
    xb.B <- xb[ID.B,]
    X1.B <- X1.o[ID.B]
    Id.B <- I.d[ID.B]
    To.B <- oT[ID.B]
    lt.B <- lt[ID.B]

    fit.sb <- survival::survfit(survival::Surv(To.B,Id.B)~1)
    NA.H <- cumsum(fit.sb$n.event/fit.sb$n.risk)
    Ht.B <- rep(NA,N)

    for(i in 1:N){

      t.id<-max(fit.sb$time[fit.sb$time<=To.B[i]])

      Ht.B[i]<-NA.H[fit.sb$time==t.id]

    }

    XXa.B <- as.matrix(cbind(X0,xa.B,Ht.B,Id.B))
    XXb.B <- as.matrix(cbind(X0,xb.B,Ht.B,Id.B))

    #NNMI X1 missing probability model
    fit.mx1 <- stats::glm(M1.B ~ XXb.B-1, family = stats::binomial(link = 'logit'))
    if(sum(is.na(fit.mx1$coefficients))>0) stop('Predictors are linearly dependent.')

    ps.mx1 <- as.numeric(XXb.B %*% as.matrix(fit.mx1$coef))
    ps.mx1.c <- (ps.mx1-mean(ps.mx1))/stats::var(ps.mx1)^0.5


    #NNMI X1 missing outcome model
    if(mfamily != 'multinomial'){

      fit.x1.nnmi <- stats::glm(X1.B[M1.B==1] ~ XXa.B[M1.B==1,] - 1, family = mfamily)
      ps.x1 <- as.numeric(XXa.B %*% as.matrix(fit.x1.nnmi$coef))
      ps.x1.c <- (ps.x1-mean(ps.x1))/stats::var(ps.x1)^0.5   # standardized

    # impute missing values
      for(i in 1:N.m){ # by individual

        ID.i <- ID.m[i]

        if(M1[ID.i]==1){ # double check missingness

          #X1 NNMI
          ps.x1.c.i <- as.vector((XXa[ID.i,]%*%as.matrix(fit.x1.nnmi$coef)-mean(ps.x1))/stats::var(ps.x1)^0.5)
          ps.mx1.c.i <- as.vector((XXb[ID.i,]%*%as.matrix(fit.mx1$coef)-mean(ps.mx1))/stats::var(ps.mx1)^0.5)
          rs.x1 <- w1 * (ps.x1.c - ps.x1.c.i)^2 + w2 * (ps.mx1.c - ps.mx1.c.i)^2
          M1.s <- M1.B[order(rank(rs.x1))]
          X1.s <- X1.B[order(rank(rs.x1))]
          X1.s.o <- X1.s[M1.s == 1]
          X1.s.o.nn <- X1.s.o[1:NN]
          dat.NNMI[ID.i,k] <- sample(X1.s.o.nn, 1)

          }
      }

    }else{

      nn <- length(datalevels)-1

      fit.x1.nnmi <- ps.x1 <- ps.x1.c <- vector("list",nn)

      for(j in 1:nn){

        X1.B1 <- rep(0,length(X1.B))

        X1.B1[is.na(X1.B)] <- NA

        X1.B1[as.character(X1.B) == rev(datalevels)[j]] <- 1

        fit.x1.nnmi[[j]] <- stats::glm(X1.B1[M1.B==1] ~ XXa.B[M1.B==1,] - 1, family = 'binomial')
        ps.x1[[j]] <- as.numeric(XXa.B %*% as.matrix(fit.x1.nnmi[[j]]$coef))
        ps.x1.c[[j]] <- (ps.x1[[j]]-mean(ps.x1[[j]]))/stats::var(ps.x1[[j]])^0.5

      }
        # impute missing values
        for(i in 1:N.m){ # by individual

          ID.i <- ID.m[i]

          if(M1[ID.i]==1){ # double check missingness

            #X1 NNMI
            ps.x1.c.i <- rep(NA,nn)
            rs.w1 <- matrix(NA,ncol=nn,nrow=length(ps.mx1.c))

            for(l in 1:nn){

              ps.x1.c.i[l] <- (XXa[ID.i,]%*%as.matrix(fit.x1.nnmi[[l]]$coef)-mean(ps.x1[[l]]))/stats::var(ps.x1[[l]])^0.5

              rs.w1[,l] <- (ps.x1.c[[l]]-ps.x1.c.i[l])^2

            }

            ps.mx1.c.i <- as.vector((XXb[ID.i,]%*%as.matrix(fit.mx1$coef)-mean(ps.mx1))/stats::var(ps.mx1)^0.5)
            rs.x1 <- w1 * rowMeans(rs.w1) + w2 * (ps.mx1.c-ps.mx1.c.i)^2
            M1.s <- M1.B[order(rank(rs.x1))]
            X1.s <- X1.B[order(rank(rs.x1))]
            X1.s.o <- X1.s[M1.s == 1]
            X1.s.o.nn <- X1.s.o[1:NN]
            dat.NNMI[ID.i,k]<- sample(X1.s.o.nn, 1)

          }
        }
    }

    dat.NNMI <- as.data.frame(dat.NNMI)

    if(mfamily != 'gaussian') dat.NNMI[,k] <- factor(dat.NNMI[,k],labels=datalevels)


    # impute survival times for censored observations
    #### need to figure out whether needs a seperate covariates matrix, xc

    if (imputeCT) {

      X1.M<-dat.NNMI[,k]
      YY <- as.matrix(cbind(X1.M,xc))
      X1.M.B<-X1.M[ID.B]
      xc.B <- xc[ID.B,]
      YY.B <- as.matrix(cbind(X1.M.B,xc.B))
      #Cox working model for failure time
      cox.t <- survival::coxph(survival::Surv(To.B, Id.B) ~YY.B)
      coef.t <- cox.t$coef
      score.t <- YY.B %*% coef.t
      score.t.c <- (score.t - mean(score.t))/sqrt(as.numeric(stats::var(score.t)))
      #Cox working model for censoring time
      cox.c <- survival::coxph(survival::Surv(To.B, (1-Id.B)) ~YY.B)
      coef.c <- cox.c$coef
      score.c <- YY.B %*% coef.c
      score.c.c <- (score.c - mean(score.c))/sqrt(as.numeric(stats::var(score.c)))
      for(i in 1:N) {
        if(I.d[i] == 0) {
          score.ti <- YY[i,  ] %*% coef.t
          score.ti.c <- (score.ti - mean(score.t))/sqrt(as.numeric(stats::var(score.t)))
          score.ci <- YY[i,  ] %*% coef.c
          score.ci.c <- (score.ci - mean(score.c))/sqrt(as.numeric(stats::var(score.c)))
          score <-(w1*(score.t.c-as.numeric(score.ti.c))^2+w2*(score.c.c-as.numeric(score.ci.c))^2)^0.5
          To.B.s <- To.B[To.B > oT[i]]
          Id.B.s <- Id.B[To.B > oT[i]]
          score.s <- (score[To.B > oT[i]])
          To.B.s.r <- To.B.s[order(rank(score.s))]
          Id.B.s.r <- Id.B.s[order(rank(score.s))]
          nn.t <- min(NN.t,length(To.B.s.r))
          if(nn.t > 0) {
            dat.km <- cbind(To.B.s.r[1:nn.t], Id.B.s.r[1:nn.t])
            if(is.na(sum(dat.km[, 1])) == F) {
              km <- survival::survfit(survival::Surv(dat.km[, 1], dat.km[, 2])~1)
              F.km <- 1 - km$surv
              u1 <- stats::runif(1)
              if(u1 <= max(F.km) & sum(dat.km[, 2]) >0) {
                dat.T.NNMI.t[i,k] <- min(km$time[F.km == min(F.km[F.km >= u1])])
                dat.Id.NNMI.t[i,k] <- 1
              }
              if(u1 > max(F.km) & sum(dat.km[, 2]) > 0) {
                dat.T.NNMI.t[i,k] <- max(km$time)
                dat.Id.NNMI.t[i,k] <- 0
              }
              if(sum(dat.km[, 2]) == 0) {
                dat.T.NNMI.t[i,k] <- max(dat.km[, 1])
                dat.Id.NNMI.t[i,k] <- 0
              }
            }
          }
        }
      }
    } else {
      dat.T.NNMI.t <- dat.Id.NNMI.t <- NA
    }

  }

  ### need to figure out the output format
  return(list(dat.NNMI=dat.NNMI,
              dat.T.NNMI.t=dat.T.NNMI.t,
              dat.Id.NNMI.t=dat.Id.NNMI.t))

}
