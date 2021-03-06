InverseMillsRatio <- function(q, mean, sd) {
  x <- (q - mean) / sd
  pdf <- dnorm(x, log = TRUE)
  cdf <- pnorm(x, log = TRUE)
  d <- exp(pdf - cdf)
  d[is.na(d)] <- 0
  return(d)
}

Pi_Zj_Zcut_new <- function(q, mean, sd, wl0) {
  a <- wl0 * pnorm(q, mean, sd)
  if(sum(a) == 0) return(wl0)
  return(a / sum(a))
}



SeparateKRpkmNew2 <- function(x, n, q, err = 1e-10) {
  k <- 1
  q <- max(q, min(x))
  c <- sum(x < q)
  x <- x[which(x >= q)]
  if (length(x) <= k) {
    warning(sprintf("The length of x is %i. Sorry, too little conditions\n", length(x)))
    return(cbind(0, 0, 0))
  }
  mean <- c()
  for (i in 1:k) {
    mean <- c(mean, sort(x)[floor(i * length(x)/(k + 1))])
  }
  if(k>1)
  {
    mean[1] <- min(x) - 1
    mean[length(mean)] <- max(x) + 1
  }
  p <- rep(1/k, k)
  sd <- rep(sqrt(var(x)), k)
  pdf.x.portion <- matrix(0, length(x), k)
  for (i in 1:n) {
    p0 <- p
    mean0 <- mean
    sd0 <- sd
    pdf.x.all <- t(p0 * vapply(x, function(x) dnorm(x, mean0, sd0), rep(0, k)))
    pdf.x.portion <- pdf.x.all/rowSums(t(pdf.x.all))
    cdf.q <- pnorm(q, mean0, sd0)
    cdf.q.all <- p0 * cdf.q
    cdf.q.portion <- cdf.q.all/sum(cdf.q.all)
    cdf.q.portion.c <- cdf.q.portion * c
    denom <- colSums(t(pdf.x.portion)) + cdf.q.portion.c
    p <- denom/(nrow(t(pdf.x.portion)) + c)
    im <- dnorm(q, mean0, sd0)/cdf.q * sd0
    im[is.na(im)] <- 0
    mean <- colSums(crossprod(x, t(pdf.x.portion)) + (mean0 - sd0 * im) * cdf.q.portion.c)/denom
    sd <- sqrt((colSums((x - matrix(mean0, ncol = length(mean0), nrow = length(x), byrow = TRUE)) ^ 2
                        * t(pdf.x.portion)) + sd0 ^ 2 * (1 - (q - mean0) / sd0 * im) * cdf.q.portion.c) / denom)
    if (!is.na(match(NaN, sd))) {
      break
    }
    if ((mean(abs(p - p0)) <= err) &&
        (mean(abs(mean - mean0)) <= err) &&
        (mean(abs(sd - sd0)) <= err)) {
      break
    }
  }
  return(cbind(p, mean, sd))
}

#' Calcuate
#'
#' @param x data, example: x<-runif(100,0,1)
#' @param n rounds
#' @param q cutoff
#' @param k k=1..5
#' @param err
#'
#' @return a matrix contains pi, mean and sd
#' @export
#'
#' @examples
SeparateKRpkmNew <- function(x, n, q, k, err = 1e-10) {
  if (k == 1) return(SeparateKRpkmNew2(x, n, q, err = 1e-10))
  q <- max(q, min(x))
  c <- sum(x < q)
  x <- x[which(x >= q)]
  if (length(x) <= k) {
    warning(sprintf("The length of x is %i. Sorry, too little conditions\n", length(x)))
    return(cbind(0, 0, 0))
  }
  mean <- c()
  for (i in 1:k) {
    mean <- c(mean, sort(x)[floor(i * length(x) / (k + 1))])
  }
  mean[1] <- min(x) - 1  # What is those two lines for?
  mean[length(mean)] <- max(x) + 1  # Without them the result of mean[1] is slightly different.
  p <- rep(1 / k, k)
  sd <- rep(sqrt(var(x)), k)
  pdf.x.portion <- matrix(0, length(x), k)

  for (i in 1:n) {
    p0 <- p
    mean0 <- mean
    sd0 <- sd

    pdf.x.all <- t(p0 * vapply(x, function(x) dnorm(x, mean0, sd0), rep(0, k)))
    pdf.x.portion <- pdf.x.all / rowSums(pdf.x.all)
    cdf.q <- pnorm(q, mean0, sd0)
    cdf.q.all <- p0 * cdf.q
    cdf.q.portion <- cdf.q.all / sum(cdf.q.all)
    cdf.q.portion.c <- cdf.q.portion * c
    denom <- colSums(pdf.x.portion) + cdf.q.portion.c
    p <- denom / (nrow(pdf.x.portion) + c)
    im <- dnorm(q, mean0, sd0) / cdf.q * sd0
    im[is.na(im)] <- 0
    mean <- colSums(crossprod(x, pdf.x.portion) + (mean0 - sd0 * im) * cdf.q.portion.c) / denom
    sd <- sqrt((colSums((x - matrix(mean0, ncol = length(mean0), nrow = length(x),
                                    byrow = TRUE)) ^ 2 * pdf.x.portion) + sd0 ^ 2 * (1 - (q - mean0) / sd0 * im) *
                  cdf.q.portion.c) / denom)
    if (!is.na(match(NaN, sd))) {
      break
    }
    if ((mean(abs(p - p0)) <= err) && (mean(abs(mean - mean0)) <= err) &&
        (mean(abs(sd - sd0)) <= err)) {
      break
    }
  }
  return(cbind(p, mean, sd))
}


# x: data, x<-runif(100,0,1), n: 300, (cutoff for stop, diff of total up ABS<1e-6, q=0 for testing data, k=1...5)

SeparateKRpkmNewp <- function(x, n, q, k, err = 1e-10) {
  q <- max(q, min(x))
  c <- sum(x < q)
  x <- x[which(x >= q)]
  mean <- c()
  for (i in 1:k) {
    mean <- c(mean, sort(x)[floor(i * length(x) / (k + 1))])
  }
  mean[1] <- min(x) - 1  # What is those two lines for?
  mean[length(mean)] <- max(x) + 1  # Without them the result of mean[1] is slightly different.
  p <- rep(1 / k, k)
  sd <- rep(sqrt(var(x)), k)
  t <- matrix(0, length(x), k)

  for (i in 1:n) {
    p0 <- p
    mean0 <- mean
    sd0 <- sd

    for (row in 1:nrow(t)) {
      all <- p0 * dnorm(x[row], mean0, sd0)
      t[row, ] <- all / sum(all)
    }

    pZil <- Pi_Zj_Zcut_new(q, mean0, sd0, p0) ################################################################
    denom <- (colSums(t) + pZil * c)
    all <- denom / (nrow(t) + c)
    p <- all / sum(all)
    im <- InverseMillsRatio(q, mean0, sd0)
    mean <- colSums(crossprod(x, t) + (mean0 - sd0 * im) * pZil * c) / denom

    a <- rep(0, length(sd))

    for (col in 1:length(sd)) {
      a[col] <- sum((x - mean0[col])^2 * t[, col])
    }

    sd <- sqrt((a + (sd0)^2 * (1 - (q - mean0) / sd0 * im) * pZil * c) / denom)
    if ((mean(abs(p - p0)) <= err) && (mean(abs(mean - mean0)) <= err) && (mean(abs(sd - sd0)) <= err)) {
      break
    }
  }
  return (cbind(p, mean, sd))
}

LogSeparateKRpkmNew <- function(x, n, q, k, err = 1e-10) {
  return (SeparateKRpkmNew(log(x), n, log(q), k, err))
}
#' Calcuate the LTMG_2LR for some genes
#'
#' @param x data, a List of NumericVector
#' @param n rounds
#' @param q cutoff of the elements in x
#' @param r maximum value of the standard diversion
#' @param s minimum value of the standard diversionzz
#' @param k number of peaks, should be 2
#' @param err the upper bound on the absolute error
#'
#' @return a matrix contains pi, mean and sd
#' @export
#'
#' @examples
SeparateKRpkmNewLRPlus <- function(x, n, q, r, s = 0.05, k = 2, err = 1e-10, M = Inf, m = -Inf) {
  c <- sapply(x, function(x) sum(x < q), simplify = "array")

  x_r0 <- lapply(x, function(x) x[which(x < q)])
  x_r <- lapply(x, function(x) x[which(x >= q)])
  x_r0.length <- sapply(x_r0, length)
  x_r.length <- sapply(x_r, length)
  x_r.non.zero <- x_r.length > 0
  x_r.non.allpos<- x_r0.length > 0
  x_r.input<-(x_r.non.allpos*x_r.non.zero)>0

  x_r <- x_r[x_r.input]

  if (sum(x_r.non.zero*x_r.non.allpos) == 0) {
    warning("Completely all 0/all pos conditions\n")
    results_c<-list()
    for(i in 1:length(x))
    {
      ccc<-matrix(0,2,3)
      colnames(ccc)<-c("p","mean","sd")
      if(x_r.non.zero[i]==0)
      {
        ccc[1,1]<-1
        ccc[2,1]<-0
        ccc[1,2]<-m
        ccc[2,2]<-M
        ccc[1,3]<-s
        ccc[2,3]<-s
      }
      if(x_r.non.allpos[i]==0)
      {
        ccc[1,1]<-0
        ccc[2,1]<-1
        ccc[1,2]<-m
        ccc[2,2]<-mean(x[[i]])
        ccc[1,3]<-s
        ccc[2,3]<-sd(x[[i]])
      }
      results_c[[i]]<-ccc
    }
    return(results_c)
  }

  if(is.na(max(x_r.length[x_r.input])<5))
  {
    browser()
  }
  if(max(x_r.length[x_r.input])<5)
  {
    warning("Too little non-zero part, forced ZIG\n")
    results_c<-list()
    for(i in 1:length(x))
    {
      ccc<-matrix(0,2,3)
      colnames(ccc)<-c("p","mean","sd")
      if(x_r.non.zero[i]==0)
      {
        ccc[1,1]<-1
        ccc[2,1]<-0
        ccc[1,2]<-m
        ccc[2,2]<-M
        ccc[1,3]<-s
        ccc[2,3]<-s
      }
      if(x_r.non.allpos[i]==0)
      {
        ccc[1,1]<-0
        ccc[2,1]<-1
        ccc[1,2]<-m
        ccc[2,2]<-mean(x[[i]])
        ccc[1,3]<-s
        ccc[2,3]<-sd(x[[i]])
      }
      if((x_r.non.allpos[i]!=0)&(x_r.non.zero[i]!=0))
      {
        ccc[1,1]<-sum(x[[i]]<q)/length(x[[i]])
        ccc[2,1]<-sum(x[[i]]>=q)/length(x[[i]])
        ccc[1,2]<-m
        ccc[2,2]<-mean(x[[i]][which(x[[i]]>=q)])
        ccc[1,3]<-s
        ccc[2,3]<-sd(x[[i]][which(x[[i]]>=q)])
        if(is.na(ccc[2,3]))
        {
          ccc[2,3]<-s
        }
      }
      results_c[[i]]<-ccc
    }
    return(results_c)
  }

  c<- c[x_r.input]
  c_sum <- sum(c)

  ncol <- length(x_r)

  x_all <- c(x_r, recursive = TRUE)

  p <- matrix(1 / k, nrow = k, ncol = ncol)

  mean <- matrix(nrow = k, ncol = ncol)
  for (col in 1:ncol) {
    ni <- length(x_r[[col]])
    for (row in 1:k) {
      tg_ic <- floor(row * ni / (k + 1))
      cc <- sort(x_r[[col]])[tg_ic]
      if (tg_ic < 1) {
        cc <- sort(x_r[[col]])[1] - 0.5
      }
      if (tg_ic > ni) {
        cc <- sort(x_r[[col]])[ni] + 0.5
      }
      mean[row, col] <- cc
    }
  }

  sd <- matrix(sqrt(vapply(x_r, var, 0)), nrow = k, ncol = ncol, byrow = TRUE)

  if(anyNA(sd)) {
    sd[is.na(sd)] <- 1
  }

  sd[which(sd<s)]<-s

  p0 <- p
  mean0 <- mean
  sd0 <- sd

  t <- lapply(x_r, function(x) matrix(nrow = length(x), ncol = k))
  for (i in 1:n) {
    ccc <- matrix(nrow = k, ncol = ncol)
    wad <- rep(0, ncol + 1)
    mean_all <- rep(0, ncol + 1)
    sd_all <- rep(0, ncol + 1)
    sd_all_1_sum <- 0
    for (col in 1:ncol) {
      t_u <- t(p[, col] * vapply(x_r[[col]], function(x) dnorm(x, mean[, col], sd[, col]), rep(0, k)))
      t[[col]] <- t_u / rowSums(t_u)
      pZil2 <- Pi_Zj_Zcut_new(q, mean[, col], sd[, col], p[, col])
      denom2 <- (colSums(t[[col]]) + pZil2 * c[[col]])
      im <- InverseMillsRatio(q, mean[, col], sd[, col])
      mean0[, col] <- colSums(crossprod(x_r[[col]], t[[col]]) + (mean[, col] - sd[, col] * im) * pZil2 * c[[col]]) / denom2
      if(denom2[1]==0)
      {
        mean0[1, col]<-mean[1, col]
      }
      if(denom2[2]==0)
      {
        mean0[2, col]<-mean[2, col]
      }
      if (anyNA(denom2)) {
        warning("denom2 conttains NA\n")
        print(x_r)
        print(col)
        print(dnorm(x_r[[col]][1], mean[1, col], sd[1, col]))
        print(dnorm(x_r[[col]][2], mean[2, col], sd[2, col]))
        browser()
      }
      if (anyNA(mean0)) {
        warning("mean0 conttains NA\n")
      }
      if (mean0[1, col] > q) {#denom2[1] == 0 ||
        mean0[1, col] <- q
      }
      if (mean0[2, col] < q) {#denom2[2] == 0 ||
        mean0[2, col] <- q
      }
      sd0[, col] <- sqrt((colSums((x_r[[col]] - matrix(mean[, col], ncol = length(mean[, col]), nrow = length(x_r[[col]]),
                                                       byrow = TRUE)) ^ 2 * t[[col]]) + (sd[, col]) ^ 2 * (1 - (q - mean[, col]) / sd[, col] * im) * pZil2 * c[[col]]) / denom2)
      if (denom2[1] == 0 || sd0[1, col] > r) {
        sd0[1, col] <- r
      }
      if (denom2[2] == 0 || sd0[2, col] > r) {
        sd0[2, col] <- r
      }

      wl2 <- denom2 / (nrow(t[[col]]) + c[[col]])

      # Reorder by mean0
      tg_R <- order(mean0[, col])
      ccc[, col] <- wl2[tg_R] * (c[[col]] + length(x_r[[col]])) / (sum(c) + sum(lengths(x_r)))
      mean0[, col] <- mean0[, col][tg_R]
      sd0[, col] <- sd0[, col][tg_R]

      wad[[col + 1]] <- ccc[-1, col]
      wad[1] <- wad[1] + ccc[1, col]
      mean_all[[col + 1]] <- mean0[, col][-1]
      mean_all[1] <- mean_all[1] + mean0[1, col] * c[[col]] / c_sum
      sd_all[[col + 1]] <- sd0[, col][-1]
      sd_all_1_sum <- sd_all_1_sum + sd0[1, col] ^ 2 * c[[col]]
    }
    sd_all[1] <- sqrt(sd_all_1_sum / c_sum)

    t0 <- matrix(nrow = sum(vapply(x_r, length, 0)), ncol = k + ncol - 1)
    for (row in 1:nrow(t0)) {
      t0_u <- wad * dnorm(x_all[row], mean_all, sd_all)
      t0[row, ] <- t0_u / sum(t0_u)
    }

    pZil0 <- Pi_Zj_Zcut_new(q, mean_all, sd_all, wad)
    denom0 <- colSums(t0) + pZil0 * c_sum

    im1 <- InverseMillsRatio(q, mean_all[1], sd_all[1])
    mean0[1, ] <- (sum(x_all * t0[, 1]) + (mean_all[1] - sd_all[1] * im1) * pZil0[1] * c_sum) / denom0[1]
    sd0[1, ] <- sqrt((sum((x_all - mean_all[1]) ^ 2 * t0[, 1]) + sd_all[1] ^ 2 * (1 - (q - mean_all[1]) / sd_all[1] * im1) *
                        pZil0[1] * c_sum) / denom0[1])

    for (col in 1:ncol) {
      t_u <- t(p[, col] * vapply(x_r[[col]], function(x) dnorm(x, mean[, col], sd[, col]), c(0, 0)))
      t[[col]] <- t_u / rowSums(t_u)
      pZil <- Pi_Zj_Zcut_new(q, mean[, col], sd[, col], p[, col])
      p_u <- (colSums(t[[col]]) + pZil * c[[col]]) / (nrow(t[[col]]) + c[[col]])
      p0[, col] <- p_u / sum(p_u)
    }

    for (col in 1:ncol) {
      sd0[1, col] <- min(sd0[1, col], r)
      sd0[2, col] <- max(min(sd0[2, col], r), s)
      mean0[1, col] <- min(mean0[1, col], q)
      mean0[2, col] <- max(mean0[2, col], q)
    }

    #print(i)
    #print(p0)
    #print(mean0)
    #print(sd0)

    if (anyNA(sd0)) {
      warning("Found at least one NA in sd")
      break
    }

    if ((mean(abs(p - p0)) <= err) && (mean(abs(mean - mean0)) <= err) && (mean(abs(sd - sd0)) <= err)) {
      break
    }

    p <- p0
    mean <- mean0
    sd <- sd0
  }

  print(i)

  ret <- vector("list", length(x))

  mean_peak1 = mean[1, 1]
  sd_peak1 = sd[1, 1]

  col <- 1
  for (index in which((x_r.non.zero*x_r.non.allpos)==1)){
    ret[[index]] <- cbind(p = p[, col], mean = mean[, col], sd = sd[, col])
    col <- col + 1
  }
  for (index in 1:length(x)){
    if((x_r.non.zero[index]*x_r.non.allpos[index])==0)
    {
      ccc<-matrix(0,2,3)
      colnames(ccc)<-c("p","mean","sd")
      if(x_r.non.zero[index]==0)
      {
        ccc[1,1]<-1
        ccc[2,1]<-0
        ccc[1,2]<-mean_peak1
        ccc[2,2]<-M
        ccc[1,3]<-sd_peak1
        ccc[2,3]<-s
      }
      if(x_r.non.allpos[index]==0)
      {
        ccc[1,1]<-0
        ccc[2,1]<-1
        ccc[1,2]<-mean_peak1
        ccc[2,2]<-mean(x[[index]])
        ccc[1,3]<-sd_peak1
        ccc[2,3]<-sd(x[[index]])
      }
      ret[[index]] <- ccc
    }
  }
  return(list(ret, i))
}

SeparateKRpkmNewLR <- function(x, n, q, r, s = 0.05, k = 2, err = 1e-10, M = Inf, m = -Inf) {
  return (SeparateKRpkmNewLRPlus(x, n, q, r, s, k, err, M, m)[[1]])
}

LogSeparateKRpkmNewLR <- function(x, n, q, r, k = 2) {
  return(SeparateKRpkmNewLR(log(x), n, log(q), r, k))
}
