############ PLOTTING FUNCTION FOR HIPHOP PLOTS
p10 <- function(matI,x,sig = 1,ylab1 = "fitness defect",xlab1 = "gene",las1 = 2,font1 = 3,cex1 = 2,pch1 = 17,
                cex.main1 = 2,cex.axis1 = 1.5,cex.lab1=1.5,cex.text1 = 1.2,...)  {
  
  mycolors = 
    c(
      "darkorange",
      "dodgerblue",
      "limegreen"
      )
  
  col = ifelse (matI[, x] > sig, 1,      ifelse (matI[, x] < -sig, 3, 2))
  w <- which(matI[, x]  > sig | matI[, x]  < -sig, arr.ind = T)
  palette(mycolors)
  posw = ifelse (w > nrow(matI) - 0.1 * nrow(matI) ,2, ifelse (w <  nrow(matI) + 0.1 * nrow(matI), 4 ,  2))
  plot(
    matI[, x],
    col = col ,cex.lab = cex.lab1,
    main = paste(colnames(matI)[x], "fitness"),
    ylab = ylab1, xaxt = "n",
    xlab = xlab1,
    las = las1,
    cex = cex1,
    cex.axis = cex.axis1,
    cex.main = cex.main1,
    pch = pch1,
    ...
  )
  if (length(w != 0))
    text(
      w,
      matI[w, x],
      names(w),
      pos = posw,
      cex = cex.text1,
      font = font1,
      ...
    )
  abline(
    h = sig ,
    col = "red",
    lty = 2,
    lwd = 2
  )
  abline(
    h = -sig ,
    col = "red",
    lty = 2,
    lwd = 2
  )
}

######### PLOTTING FUNCTION
stringWINDOW = function(x, width = 25){
  
  strng = paste(strwrap(x,width = width), collapse="\n")
  
  strng
  
}






############ GRAPHING FUNCTION
meltDF = function(mat, row, df = phiphop) {
  mx = reshape2::melt(mat[row, ], value.name = "fitness_defect")
  mx = data.frame(
    screen = rownames(mx),
    gene = row,
    fitness_defect = mx$fitness_defect,
    stringsAsFactors = F
  )
  
  mx$gene = paste0("<a href=https://www.yeastgenome.org/locus/",mx$gene,">",mx$gene,"</a>")
  
  mx = mx[, c("gene", "screen", "fitness_defect")]
  m = match(mx$screen, df$name)
  table(is.na(m))
  
  
  
  mx$target = df$target_html[m]
  
  mx$signature = df$signature[m]
  mx$shape = 19
  
  
  g = grep(row,df$target)
  if(length(g) > 0) mx$shape[g] = 17
  
  mx$PCID = df$PCID[m]
  mx$FDA = factor(df$FDA[m])
  mx$mechanism = df$mechanism[m]
  n =   c("screen",
    "gene",
    "fitness_defect",
    "PCID",
    "mechanism",
    
    
    "target",
    "signature"   ,
     "FDA",
    "shape"
  )
  g = grep("fitness",colnames(mx))
 
  mx = mx[,n]
  mx
}


############ GENE ANNOTATION
# generates annotation for the HOP fitness plots for a single screen from a matrix
# mat - data matrix; each colname is a screen. values are defect scores, rownames are gene names
# fdat - gene annotation data.frame file with genes, ORFS, descriptors, etc.
# cmp = name of screen of interest; column name
# sgdlink - provides link to SGD for gene names if desired
# enrichtype = FD for FITNESS DEFECT; FA for FITNESS ADVANTAGE (resistant) 
# out put is shown on the HOP fitness tap underneath the plots; used as input for the ggfit plots below
####
geneAnno <-  function(mat,fdat = NULL,cmp,sgdlink = T, enrichtype = "FD"){
  library(dplyr)
  
  if(is.null(fdat)){
  fdat = read.delim("2020oct28_files/december4_fdat_updated.txt",stringsAsFactors = F,check.names = F)}
  mat = mat [order(rownames(mat)),]
  w1 = which(colnames(mat) %in% cmp)
  dmat = data.frame(GENE = rownames(mat),FD = mat[,w1],stringsAsFactors = F,check.names = F)
  
  dmat = dmat %>% arrange(desc(FD))
  m = match(dmat$GENE,fdat$sgd_gene)
  dmat[,c("ORF","descriptor")] = NA
  dmat[,c("ORF","descriptor")] = fdat[m,c("sgd_orf","descriptor")]
  
  if(sgdlink) dmat$GENE = paste0("<a href=https://www.yeastgenome.org/locus/",dmat$GENE,">",dmat$GENE,"</a>")
  
  dmat = dmat[,c("FD","ORF","GENE","descriptor")]
  
  if(enrichtype == "fa"){names(dmat)=c("FA","ORF","GENE","descriptor")
  dmat$FA = format(round(dmat$FA,2),nsmall = 1,scientific = F)} else {
    
    dmat$FD = format(round(dmat$FD,2),nsmall = 1,scientific = F)}
  
  dmat
  
}
###################### SIMILAR TO ABOVE FUNCTION BUT INPUT IS A DATAFRAME RATHER INSTEAD OF A MATRIX ##############################
dfgeneAnno <-  function(df,fdat = NULL,gene = "gene",sgdlink = T,keep = F){
  
  if(is.null(fdat)){
  fdat = read.delim("2020oct28_files/december4_fdat_updated.txt",stringsAsFactors = F,check.names = F)}
  gene1 = which(names(df) %in% gene)
  
  
  w1 = which(df[,gene1] %in% fdat$sgd_gene)             
  
  m = match(df[w1,gene1],fdat$sgd_gene)
  
  namegiven = names(df)
  
  nametouse = c("ORF","GENE","description")
  w = which(nametouse %in% namegiven)
  if(length(w) > 0) nametouse[w]=paste0(nametouse[w],":A")
  df[,nametouse]= NA
  
  df[w1,nametouse] = fdat[m,c("sgd_orf","sgd_gene","descriptor")]
  
  
  if(sgdlink) df$GENE[w1] = paste0("<a href=https://www.yeastgenome.org/locus/",fdat$sgd_gene[m],">",fdat$sgd_gene[m],"</a>")
  if(keep == F) { df = df[,-gene1] }
    
  
  df
}


