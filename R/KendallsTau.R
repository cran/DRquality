KendallsTau=function(InputDists,OutputDists){
# tau=KendallsTau(InputDists,OutputDists)
# Berechnet den statistischen Zusammenhang nach Kendall
#
# INPUT
# InputDists             Matrize der Distanzen des Eingaberaumes
# OutputDists            Matrize der Distanzen des Ausgaberaumes 
# 
# OUTPUT
# tau										 numeric, Kendalls tau
# Author: MT 10/2015
   if (!requireNamespace('pcaPP', quietly = TRUE)) {
    message(
      'Subordinate package (pcaPP) is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return(
      list(
        Object = "Subordinate package (pcaPP) is missing.
                Please install the package which is defined in 'Suggests'."
      )
    )
  }
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
  
  #return(cor(x,y,method='kendall'))
  return(pcaPP::cor.fk(x,y))
  
}