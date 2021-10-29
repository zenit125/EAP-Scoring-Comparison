Code_to_SPSS_Label <- function(Data)
{
  return(
    as.factor(
    names(attributes(Data)$labels)[match(Data, attributes(Data)$labels)]
  ))
}