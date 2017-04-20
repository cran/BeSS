 bess.one = function(x, y, family = c("gaussian", "binomial", "cox"),
                  s = 1,
                  max.steps = 15,
                  glm.max = 1e6,
                  cox.max = 20,
                  normalize = FALSE)
{
  family <- match.arg(family)
  if(ncol(x)==1|is.vector(x)) stop("x should be two columns at least!")
  if(missing(family)) stop("Please input family!")
  if(family=="binomial")
  {
    if(length(unique(y))!=2)  stop("Please input binary variable!")else
      if(setequal(y_names<-unique(y),c(0,1))==FALSE)
      {
        y[which(y==unique(y)[1])]=0
        y[which(y==unique(y)[2])]=1
        y=as.numeric(y)
      }
  }
  if(family=="cox")
  {
    if(!is.matrix(y)) y=as.matrix(y)
    if(ncol(y)!=2) stop("Please input y with two columns!")
  }
  if(is.vector(y))
  {
    if(nrow(x)!=length(y)) stop("Rows of x must be the same as length of y!")
  }else{
    if(nrow(x)!=nrow(y)) stop("Rows of x must be the same as rows of y!")
  }

  beta0=rep(0,ncol(x))
  if(s>length(beta0))
  {stop("s is too large")}

  x=as.matrix(x)

  if(family=="gaussian")
  {
    out=bess.lm(x=x,y=y,beta0=beta0,s=s,max.steps=max.steps,normalize=normalize)
    class(out)="bess.one"
    return(out)
  }
  if(family=="binomial")
  {
    fit = bess.glm(x=x,y=y,beta0=beta0,
          intercept=0,s=s,
          max.steps=max.steps,
          glm.max=glm.max,
          normalize=normalize)
    if(!setequal(y_names,c(0,1)))
    {
     fit[[8]] = y_names
     names(fit[[8]])=y_names
     class(fit)="bess.one"
     return(fit)
    }else
    {
      class(fit)="bess.one"
      return(fit)
    }
  }
  if(family=="cox")
    {
     out=return(bess.cox(x=x,y=y,beta0=beta0,
                         s=s,
                         cox.max=cox.max,
                         max.steps=max.steps,
                         normalize=normalize))

     class(out)="bess.one"
     return(out)
    }
}



