SpearmanError = function(InputDists,OutputDists){
  # Calculates the error of a projection with spearman's rank correlation coefficient
  # 
  # INPUT
  # InputDists(1:n,1:n)            dissimilarities in Input Space between the n data points
  # OutputDists(1:n,1:n)           dissimilarities in Input Space between the n data points
  # OUTPUT
  # rho                                 rank correlation coefficient
  #
  # author: MT 01/2016
  # based on Diss
#EXAMPLE:
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
  
  VectorOfInputDists=InputDists[lower.tri(InputDists, diag = FALSE)]
  VectorOfOutputDists=OutputDists[lower.tri(OutputDists, diag = FALSE)]

  bInput = rank(VectorOfInputDists)
  bOutput = rank(VectorOfOutputDists)

  # Differenzen zwischen den Rangen
  diffRang = (bOutput - bInput)^2
  n=length(VectorOfInputDists)
  #kappa = n*(n-1)/2 #Wenn Distanzmatrix genommen werden wuerde
  kappa=n #im Falle der unteren Dreicksmatrix
  rho = 1 - ((1/ (kappa^3 - kappa)) * 6 * sum(diffRang))
  return(rho = rho)
}