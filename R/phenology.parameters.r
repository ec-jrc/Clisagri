phenology.parameters = function(breeder=FALSE)
{
  if(!breeder)
    print(read.table(file="data/PARAM.DESC.csv", header = TRUE, sep=";"))
  else
    print(read.table(file="data/PARAM.DESC.breeder.csv", header = TRUE, sep=";"))
}
