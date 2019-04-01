rotate.data <- function(data, chans=NULL, theta=NULL)
{
  if (class(data)== "flowFrame" & !is.null(chans))
  {
    data.new <- exprs(data)[,chans]
    if (is.null(theta))
    {
      reg.slope <- atan(lm(data.new[,1] ~ data.new[,2])$coefficients[2])
      theta <- pi/2 - reg.slope
    }
    data.new <- data.new %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
    exprs(data)[,chans] <- data.new
  }else{
    data <- data %*% matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=T)
  }
  return(list(data=data,theta=theta))
}
