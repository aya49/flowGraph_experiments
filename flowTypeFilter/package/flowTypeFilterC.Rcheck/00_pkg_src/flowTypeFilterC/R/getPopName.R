##Given a phenotype and a list of marker names, extract the phenotype name as a combination of +/- markers.
getPopName <- function(digit, Markers){
  M<-length(Markers);
  a<-c(as.vector(matrix(0,1, M-length(digitsBase(digit,3)))),digitsBase(digit,3))
  a <- a-1
  str=''
  for (i in 1:length(a)){
    if (a[i]==0)
      ##str=paste(str, Markers[i], 'o', sep='');
      next;
    if (a[i]>0)
      str=paste(str, Markers[i], '+', sep='');
    if (a[i]<0)
      str=paste(str, Markers[i], '-', sep='');
  }
  return (str);
}
