library(dplyr)
library(binaryLogic)


toget_immune_column_combi = function(){
  ncol=10
  bitemp = as.binary(1:(2^(ncol)-1))
  
  toget_colcombi = vector()
  for (n in 1:length(bitemp)) {
    v = as.character(bitemp[[n]])
    v = paste0(v, collapse = "")
    v = sprintf("%010s", v)
    v = gsub(" ", "0", v)
    toget_colcombi[n] = v
  }
  
  return(toget_colcombi)
}

TIL_immune_column_combi = function(){
  ncol=6
  bitemp = as.binary(1:(2^(ncol)-1))
  
  TIL_colcombi = vector()
  for (n in 1:length(bitemp)) {
    v = as.character(bitemp[[n]])
    v = paste0(v, collapse = "")
    v = sprintf("%06s", v)          #only immune = 6, immune+clinical = 10
    v = gsub(" ", "0", v)
    TIL_colcombi[n] = v
  }
  
  return(TIL_colcombi)
}

PBMC_immune_column_combi = function(){
  ncol=6
  bitemp = as.binary(1:(2^(ncol)-1))
  
  PBMC_colcombi = vector()
  for (n in 1:length(bitemp)) {
    v = as.character(bitemp[[n]])
    v = paste0(v, collapse = "")
    v = sprintf("%06s", v)            #only immune = 6, immune+clinical = 10
    v = gsub(" ", "0", v)
    PBMC_colcombi[n] = v
  }
  
  return(PBMC_colcombi)
}


# dec2Ndigit = function(n, number, ncol) {
#   ## n: to_digit, if binary == 2
#   ## number: dec to be converted 
#   
#   ## to N digit
#   toNdigit = ""
#   while (number){
#     toNdigit = paste0(toNdigit, as.character(number %% n))
#     number = number %/% n
#   }
#   
#   ## void to 0
#   ncol = ncol - nchar(toNdigit)
#   add_void = paste0(rep("0", ncol), collapse = "")
#   
#   ## reverse to N digit
#   rev_toNdigit = vector()
#   for(i in nchar(toNdigit):1) { rev_toNdigit = append(rev_toNdigit, substr(toNdigit, i, i)) }
#   rev_toNdigit = paste0(rev_toNdigit, collapse = "")
#   
#   ## column_combination
#   column_combination = paste0(add_void, rev_toNdigit)
#   
#   return(column_combination)
# }


  



