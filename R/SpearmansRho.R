SpearmansRho=function(InputDists,OutputDists){
  
  if(!is.matrix(InputDists)){
    warning('InputDists is not a matrix. Calling as.matrix()')
    InputDists=as.matrix(InputDists)
  }
  if(!is.matrix(OutputDists)){
    warning('OutputDists is not a matrix. Calling as.matrix()')
    OutputDists=as.matrix(OutputDists)
  }
  if(!mode(InputDists)=='numeric'){
    warning('InputDists is not a numeric matrix. Calling mode(InputDists)="numeric"')
    mode(InputDists)='numeric'
  }
  if(!mode(OutputDists)=='numeric'){
    warning('OutputDists is not a numeric matrix. Calling mode(OutputDists)="numeric"')
    mode(OutputDists)='numeric'
  }
  x=InputDists[lower.tri(InputDists, diag = FALSE)]
  y=OutputDists[lower.tri(OutputDists, diag = FALSE)]

  rho=cor(x,y,method = "spearman")
  
  return(rho)
}