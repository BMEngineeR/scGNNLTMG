CreateLTMGObject <- function(x = input_dir) {
  
  BRIC_Object<- new(Class = 'LTMG', InputData = as.matrix(x))
  return(BRIC_Object)
  
}


