mkSingleThread<-function(nNodes,g)
# This function creates one thread for a fixed j-mer. 
{
  jj=1
  mat0=matrix(ncol=nNodes)
#    mat0=matrix(ncol=g$nZ)
  mat=mat0
#    print(mat0)
  mat0[1,]=g$mylets[1:nNodes]
#    print(mat0)
  if (nNodes>1) 
  {
    for (patchLen in 2:nNodes) {
      for (startThread in 1:(nNodes-patchLen+1)) {
#          cat("\nstartThread=",startThread," patchLen=",patchLen,"\n")
        mat[jj,]=mat0 # reset in case not properly reset
        mat[jj,startThread:(startThread+patchLen-1)]=mat0[startThread]  # fill whole band with start thread label
        mat=rbind(mat,mat0)
        jj=jj+1
      }
      if (patchLen == nNodes/2) {
        mat[jj,]=mat0 # reset in case not properly reset
        mat[jj,1:patchLen]=mat0[1]  # fill first half
        mat[jj,(patchLen+1):(2*patchLen)]=mat0[patchLen+1]  # fill second half
        mat=rbind(mat,mat0)
        jj=jj+1
      }
    }
  } else mat=mat0
  mat
}

