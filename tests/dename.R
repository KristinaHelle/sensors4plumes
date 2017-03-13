dename = function(x, kind){
  if (missing(kind)){
    names(x) = NULL    
  } else {
    if (kind == "rownames"){
      rownames(x) = NULL
    }
    if (kind == "dimnames"){
      dimnames(x) = NULL
    }
  }
  return(x)
}
