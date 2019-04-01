coldel = c()

for (i in c(1:ncol(matrixPval)) ) {
  if (sum(matrixPval[,i]!=0) < 3)
  coldel = append(coldel, i)
}

rowdel = c()
for (i in c(1:nrow(matrixPval)) ) {
  if (sum(matrixPval[,i]!=0) < 3)
  rowdel = append(rowdel, i)
}