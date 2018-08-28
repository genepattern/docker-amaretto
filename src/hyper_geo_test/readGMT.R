readGMT<-function (filename) 
{
    a = scan(filename, what = list("", ""), sep = "\t", quote = NULL, 
        fill = T, flush = T, multi.line = F,strip.white = TRUE)
    geneset.names = a[1][[1]]
    geneset.descriptions = a[2][[1]]
    dd = scan(filename, what = "", sep = "\t", quote = NULL,strip.white = TRUE)
    nn = length(geneset.names)
    n = length(dd)
    ox = rep(NA, nn)
    ii = 1
    for (i in 1:nn) {
        #print(i)
        while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] != 
            geneset.descriptions[i])) {
            ii = ii + 1
        }
        ox[i] = ii
        ii = ii + 1
    }
    genesets = vector("list", nn)
    
    for (i in 1:(nn - 1)) {
       # print(i, fill = T)
        if(nn>1){
        i1 = ox[i] + 2
        i2 = ox[i + 1] - 1
        geneset.descriptions[i] = dd[ox[i] + 1]
        current_geneset = dd[i1:i2]
        genesets[[i]] = toupper(current_geneset[current_geneset!=""])
      }else{
        current_geneset = dd[3:n]
        genesets[[1]] = toupper(current_geneset[current_geneset!=""])
      }
    }
    
    geneset.descriptions[nn] = dd[ox[nn] + 1]
    current_geneset = dd[(ox[nn] + 2):n]
    genesets[[nn]] = toupper(current_geneset[current_geneset!=""])
    out = list(genesets = genesets, geneset.names = geneset.names, 
        geneset.descriptions = geneset.descriptions)
    class(out) = "GSA.genesets"
    return(out)
}
