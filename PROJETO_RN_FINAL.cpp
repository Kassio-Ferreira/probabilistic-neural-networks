#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec PNN(mat patterns, uvec labels, mat testPattern, float sigma=0.5, int kernel=1){
  
  int j=0;
  int n = patterns.n_rows;
  int d = patterns.n_cols;
  uvec c1 = unique(labels);
  int c = c1.n_elem;
  mat W(n,d); 
  //W.zeros();
  mat A(n,c);
  A.zeros();
  vec classe(testPattern.n_rows);
  int f = 0;
  vec x_jk(d);
  
  do{
    W.row(j) = patterns.row(j)/conv_to<float>::from(sqrt(sum(pow(patterns.row(j), 2), 1)));
    //cout << W.row(j);
    f = labels.at(j);
    A(j,f-1) = 1;
    j++;
  } while (j < n);
  
  
  for(int i=0; i<testPattern.n_rows; i++){
    
    int k=0;
    int c2 = 0;
    rowvec gc(c);
    gc.zeros();
    float zk;
    
    do{
      
      zk = conv_to<float>::from(W.row(k) * testPattern.row(i).t());  
      //cout << W.row(k) * testPattern.row(i).t() << endl;
      
      c2 = index_max(A.row(k));
      
      if(kernel == 1){
        gc(c2) += exp((zk-1)/pow(sigma,2));  //gaussiana
      } else if(kernel == 2){
        gc(c2) += sqrt(1+pow(zk*sigma,2));   //multiquadrica
      } else if(kernel == 3){
        gc(c2) += pow(zk,2) * log(zk);       //spline poli-harmônica
      }
      
      //gc(c2) += sqrt(1+pow(zk*sigma,2));
      //gc(c2) += pow(zk,2) * log(zk);
      //cout << gc << endl;
      
      k++;
    } while (k < n);
    //cout << index_max(gc) << endl;
    classe(i) = index_max(gc);
  }
  
  return classe+1;
}

/*** R
library(caret) # para criar os folds
library(microbenchmark)

# PolSAR data -------------------------------------------------------------

# <1.> Abrir a imagem
# usar ImgMPAE
dados <- load(file.choose()) 
  
# >function-@.1<
imagematrix <- function(mat, type=NULL, ncol=dim(mat)[1], nrow=dim(mat)[2],
                        noclipping=FALSE) {
  if (is.null(dim(mat)) && is.null(type)) stop("Type should be specified.")
    if (length(dim(mat)) == 2 && is.null(type)) type <- "grey"
    if (length(dim(mat)) == 3 && is.null(type)) type <- "rgb"
    if (type != "rgb" && type != "grey") stop("Type is incorrect.")
      if (is.null(ncol) || is.null(nrow)) stop("Dimension is uncertain.")
        imgdim <- c(ncol, nrow, if (type == "rgb") 3 else NULL)
        if (length(imgdim) == 3 && type == "grey") {
  # force to convert grey image
          mat <- rgb2grey(mat)
        }
        if (noclipping == FALSE && ((min(mat) < 0) || (1 < max(mat)))) {
          warning("Pixel values were automatically clipped because of range over.") 
          mat <- clipping(mat)
        }
        mat <- array(mat, dim=imgdim)
          attr(mat, "type") <- type
          class(mat) <- c("imagematrix", class(mat))
          mat
}

# >function-@.2<
print.imagematrix <- function(x, ...) {
  x.dim <- dim(x)
  cat("size: ", x.dim[1], "x", x.dim[2], "\n")
  cat("type: ", attr(x, "type"), "\n")
}

# >function-@.3<
plot.imagematrix <- function(x, ...) {
  colvec <- switch(attr(x, "type"),
                   grey=grey(x),
                   rgb=rgb(x[,,1], x[,,2], x[,,3]))
  if (is.null(colvec)) stop("image matrix is broken.")
    colors <- unique(colvec)
    colmat <- array(match(colvec, colors), dim=dim(x)[1:2])
    image(x = 0:(dim(colmat)[2]), y=0:(dim(colmat)[1]),
          z = t(colmat[nrow(colmat):1, ]), col=colors,
          xlab="", ylab="", axes=FALSE, asp=1, ...)
}

# >function-@.4<
imageType <- function(x) {
  attr(x, "type")
}

# >function-@.5<
rgb2grey <- function(img, coefs=c(0.30, 0.59, 0.11)) {
  if (is.null(dim(img))) stop("image matrix isn't correct.")
    if (length(dim(img))<3) stop("image matrix isn't rgb image.")
      imagematrix(coefs[1] * img[,,1] + coefs[2] * img[,,2] + coefs[3] * img[,,3],
                  type="grey")
}

# >function-@.6<
clipping <- function(img, low=0, high=1) {
  img[img < low] <- low
  img[img > high] <- high
  img
}

# >function-@.7<
normalize <- function(img) {
  (img - min(img))/(max(img) - min(img))
}

# >function-@.8< Fun??o que equaliza uma imagem de uma banda
equalize <- function(imagem) {
  imagemeq <- ecdf(imagem)(imagem)
  dim(imagemeq) <- dim(imagem)
  return(imagemeq)
}

# Função showBand para plotar as imagens

showBand = function(img,Pmin=0.0,Pmax=.8){
# parameters: img: DataBase from PolSAR image | Pmin,Pmax: Percentile from data
# return: Plot of a PolSAR image
#library(rimage)
  Dimg<-dim(img)
  if(length(Dimg)!=3){
    imgP<-array(0,c(Dimg[1],Dimg[2],3))
    imgP[,,1]<-img
    imgP[,,2]<-imgP[,,3]<-imgP[,,1]
  }else{
    imgP<-img
}
  
MINimg<-quantile(imgP,Pmin)
MAXimg<-quantile(imgP,Pmax)
plot(imagematrix((imgP-MINimg)/(MAXimg-MINimg)))}


HH <-imagematrix( normalize( Re( ImgMPAE[,,1] ) ) )
HHeq <- equalize(HH)
  
HV <-imagematrix( normalize( Re( ImgMPAE[,,2] ) ) )
HVeq <- equalize(HV)
  
VV <-imagematrix( normalize( Re( ImgMPAE[,,3] ) ) )
VVeq <- equalize(VV)
  
  
b1 <- imagematrix( normalize(Re(ImgMPAE[,,1]+ImgMPAE[,,3])) )
b1e <- equalize(b1)
b2 <- imagematrix( normalize(Re(ImgMPAE[,,1]-ImgMPAE[,,3])) )
b2e <- equalize(b2)
b3 <- imagematrix( normalize(Re(ImgMPAE[,,2])) )
b3e <- equalize(b3)
  
CANAISequalizados <- imagematrix(
    c(HHeq, HVeq, VVeq),
    type="rgb",
    ncol=dim(b1e)[1], nrow=dim(b1e)[2]
)
  
Md = Re(ImgMPAE[,,1])
dim(Md)
DimY = dim(Md)
  
# amostra 
par(mar=c(0,0,0,0))
showBand(CANAISequalizados)
rect(217,171,266,220, border="darkgreen",lwd=5) # cidade
rect(197,101,246,150, border="darkblue",lwd=5) # floresta
rect(75,75,124,124, border="red",lwd=5) # mar
  
# cidade
cidade = CANAISequalizados[(DimY[1]-220):(DimY[1]-171),217:266,]
showBand(cidade)
dim(cidade)
  
amostraCidade = 
  cbind(
    rep(1, 2500),
    as.vector(cidade[,,1]),
    as.vector(cidade[,,2]),
    as.vector(cidade[,,3])
  )
  
# floresta
floresta = CANAISequalizados[(DimY[1]-150):(DimY[1]-101),197:246,]
showBand(floresta)
dim(floresta)
  
amostraFloresta = 
  cbind(
    rep(2, 2500),
    as.vector(floresta[,,1]),
    as.vector(floresta[,,2]),
    as.vector(floresta[,,3])
  )
  
  
# mar
mar = CANAISequalizados[(DimY[1]-124):(DimY[1]-75),75:124,]
showBand(mar)
dim(mar)
  
amostraMar = 
  cbind(
    rep(3, 2500),
    as.vector(mar[,,1]),
    as.vector(mar[,,2]),
    as.vector(mar[,,3])
  )
  
  
# amostra final
Exemplos = rbind(amostraCidade, amostraFloresta, amostraMar)
Exemplos = Exemplos[sample(1:nrow(Exemplos), nrow(Exemplos), replace = FALSE), ]
head(Exemplos)

# Função de seleção de exemplos -------------------------------------------

ExampleSelection = function(treino, k=20, classes=7, method="k-vizinhos"){
    
  grupos = list();
  for(i in 1:classes){
    grupos[[i]] = treino[treino[,1] == i, 2:ncol(treino)]
  }
    
# pega os k elementos com menor distância para o protótipo
  for(j in 1:classes){
    prototipo  = colMeans(grupos[[j]])
    distWithPrototype = dist(rbind(prototipo, grupos[[j]]))
    distWithPrototype = as.matrix(distWithPrototype); 
      
    if(method == "k-vizinhos"){
      bestElements = order(distWithPrototype[1,-1])[1:k]
    } else if(method == "sample"){
      bestElements = sample(1:length(distWithPrototype[1,-1]),k,
                            prob =distWithPrototype[1,-1]/sum(distWithPrototype[1,-1]))
    }
      
    grupos[[j]] = cbind(j,(grupos[[j]])[bestElements,])
        
  }
    
  newdata = grupos[[1]]
  for(i in 2:length(grupos)){
    newdata = rbind(newdata, grupos[[i]])
  }
    
  return(newdata)
}


# Filtro de voto majoritário ----------------------------------------------
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

windows_VM <- function(M, windowSize){
  
  #ambiente 1: cidade, 2:floresta, 3:mar
  
  windowSizeA = windowSize-1
  ss = seq(1,nrow(M),windowSize)
  tt = seq(1,ncol(M),windowSize)
  
  for(i in ss){
    for(j in tt){
      M[i:(i+windowSizeA), j:(j+windowSizeA)] <- getmode(M[i:(i+windowSizeA), j:(j+windowSizeA)])
    }
  }
  
  return(M)
}

# Validação cruzada para avaliar Kernels (também com novo método) ---------
  
mediaTotal_gauss = mediaTotal_spline = NULL
mediaTotal_gauss_NM = mediaTotal_spline_NM = NULL
  
pb <- txtProgressBar(min = 0, max = 30, style = 3)

for(j in 1:30){
      
    Sys.sleep(0.1)
    setTxtProgressBar(pb, j)
      
    foldsList = createFolds(Exemplos[,1], k=10, list=TRUE)  
    medias_gaussiana = medias_spline =  medias_multiquadrica = NULL
    medias_gaussiana_NM = medias_spline_NM =  medias_multiquadrica_NM = NULL
      
  # validação cruzada:
    for(i in 1:10){
        
      testPattern = Exemplos[foldsList[[i]],-1] # fold teste
      labelstest = Exemplos[foldsList[[i]],1] # rotulos de teste
      patterns = Exemplos[-foldsList[[i]],-1]  # fold treino
      labels = Exemplos[-foldsList[[i]],1] # rotulos de treino
        
    # PNN com RBF gaussiana
      resGauss = PNN(patterns, labels, testPattern, sigma=0.05, kernel=1)
      medias_gaussiana[i] = mean(resGauss == labelstest)
        
    # PNN com RBF spline poli-harmonica
      resSpline = PNN(patterns, labels, testPattern, kernel=3)
      medias_spline[i] = mean(resSpline == labelstest)
        
    # Reorganizar os dados para rodar o novo método com o melhor K para cada kernel
      patternsPOLSARinit = Exemplos[-foldsList[[i]], ]
      testPattern = Exemplos[-foldsList[[i]],-1]
      labelstest = Exemplos[-foldsList[[i]],1]
        
    # gaussiana
      patternsPOLSAR_gauss = ExampleSelection(patternsPOLSARinit, k=510, classes=3)
      patternsPOLSARLabels_gauss = patternsPOLSAR_gauss[ ,1]
      patternsPOLSAR_gauss = patternsPOLSAR_gauss[ ,-1] # conjunto de treinamento para novo método
    # usando a rede
      resGauss_NM = PNN(patternsPOLSAR_gauss, patternsPOLSARLabels_gauss, testPattern, sigma=0.05, kernel=1)
      medias_gaussiana_NM[i] = mean(resGauss_NM == labelstest)
        
    # spline 
      patternsPOLSAR_spline = ExampleSelection(patternsPOLSARinit, k=805, classes=3)
      patternsPOLSARLabels_spline = patternsPOLSAR_spline[ ,1]
      patternsPOLSAR_spline = patternsPOLSAR_spline[ ,-1] # conjunto de treinamento para novo método
        
      resSpline_NM = PNN(patternsPOLSAR_spline, patternsPOLSARLabels_spline, testPattern, kernel=3)
      medias_spline_NM[i] = mean(resSpline_NM == labelstest)
        
  }
      
  mediaTotal_gauss[j] = mean(medias_gaussiana)
  mediaTotal_spline[j] = mean(medias_spline)

        
  mediaTotal_gauss_NM[j] = mean(medias_gaussiana_NM)
  mediaTotal_spline_NM[j] = mean(medias_spline_NM)

}
    
    
mediaTotal_gauss
mediaTotal_spline

mediaTotal_gauss_NM
mediaTotal_spline_NM

# tabela
library(xtable)
  
validacao_cruzada = t(data.frame(Media_Gaussiano = mean(mediaTotal_gauss),
                                 Media_Spline = mean(mediaTotal_spline),
                                 media_Gaussiano_NM = mean(medias_gaussiana_NM),
                                 media_Spline_NM = mean(medias_spline_NM)))

colnames(validacao_cruzada) = "Acurácia"
tabelaVC = xtable(validacao_cruzada,sanitize.rownames.function = italic,
                  sanitize.colnames.function = large,
                  booktabs = TRUE) 
  
digits(tabelaVC) <- 4

# Obtenção do melhor K ----------------------------------------------------

# ir pegando os K mais próximos até que o p-valor do teste t seja menor que 0.05
# incrementar K até que p < 0.05. Pare e retorne K. (comparar com a rede default)
# com K exemplos eu consigo uma acurácia com diferença média de 5%.

# Kernel Spline -----------------------------------------------------------

medias_spline
  
foldsList = createFolds(Exemplos[,1], k=10, list=TRUE)  
medias_spline_NM = NULL

for(K in seq(5,2200, 50)){  
# validação cruzada:
  for(i in 1:10){
    
# Reorganizar os dados para rodar o novo método com o melhor K para cada kernel
    patternsPOLSARinit = Exemplos[-foldsList[[i]], ]
    testPattern = Exemplos[-kk,-1]
    labelstest = Exemplos[-kk,1]
    
    testPattern = Exemplos[foldsList[[i]],-1] # fold teste
    labelstest = Exemplos[foldsList[[i]],1] # rotulos de teste
    patterns = Exemplos[-foldsList[[i]],-1]  # fold treino
    labels = Exemplos[-foldsList[[i]],1] # rotulos de treino
    
# spline
    patternsPOLSAR_spline = ExampleSelection(patternsPOLSARinit, k=K, classes=3)
    patternsPOLSARLabels_spline = patternsPOLSAR_spline[ ,1]
    patternsPOLSAR_spline = patternsPOLSAR_spline[ ,-1] # conjunto de treinamento para New Method
# usando a rede
    resSpline_NM = PNN(patternsPOLSAR_spline, patternsPOLSARLabels_spline, testPattern, kernel=3)
    medias_spline_NM[i] = mean(resSpline_NM == labelstest)
  }
  
  if(t.test(medias_spline, medias_spline_NM, mu=0.05, alternative = "less")$p.value < 0.05)
    break;
  
  print(K)
}  

melhor_K = K-50

# classificar toda a imagem com novo método 
par(mfrow=c(1,2))
showBand(CANAISequalizados)
  
im = 
  cbind(
    as.vector(CANAISequalizados[,,1]),
    as.vector(CANAISequalizados[,,2]),
    as.vector(CANAISequalizados[,,3])
)

patternsPOLSAR = ExampleSelection(patternsPOLSARinit, k=melhor_K, classes=3)
patternsPOLSARLabels = patternsPOLSAR[ ,1]
patternsPOLSAR = patternsPOLSAR[ ,-1] # conjunto de treinamento para New Method; dim(patternsNM)

imagemClassificada = PNN(patternsPOLSAR, patternsPOLSARLabels, im, sigma=1, kernel=3)
par(mar=c(0,0,0,0))
showBand(matrix(imagemClassificada, nrow = nrow(CANAISequalizados), byrow = FALSE), Pmax = 1)
  
# voto majoritário com janela 3x3
img_total_VM = windows_VM(matrix(imagemClassificada, nrow = nrow(CANAISequalizados), byrow = FALSE)[1:228,1:435], 3)
showBand(img_total_VM, Pmax = 1)

# voto majoritário com janela 5x5
img_total_VM = windows_VM(matrix(imagemClassificada, nrow = nrow(CANAISequalizados), byrow = FALSE)[1:225,1:435], 5)
showBand(img_total_VM, Pmax = 1)

# tempo de processamento
gg = microbenchmark(
  tradicional = PNN(patterns, labels, im, kernel = 3),
  novo_metodo = PNN(patternsPOLSAR, patternsPOLSARLabels, im, kernel=3),
  times = 1L)
  
gg
  
# Kernel Gaussiano --------------------------------------------------------
medias_gaussiana
    
foldsList = createFolds(Exemplos[,1], k=10, list=TRUE)  
medias_gaussiana_NM = NULL
  
for(K in seq(10,2200, 50)){  
# validação cruzada:
  for(i in 1:10){
      
# Reorganizar os dados para rodar o novo método com o melhor K para cada kernel
    patternsPOLSARinit = Exemplos[-foldsList[[i]], ]
    testPattern = Exemplos[-kk,-1]
    labelstest = Exemplos[-kk,1]
      
    testPattern = Exemplos[foldsList[[i]],-1] # fold teste
    labelstest = Exemplos[foldsList[[i]],1] # rotulos de teste
    patterns = Exemplos[-foldsList[[i]],-1]  # fold treino
    labels = Exemplos[-foldsList[[i]],1] # rotulos de treino
      
# gaussiana
    patternsPOLSAR_gaussiana = ExampleSelection(patternsPOLSARinit, k=K, classes=3)
    patternsPOLSARLabels_gaussiana = patternsPOLSAR_gaussiana[ ,1]
    patternsPOLSAR_gaussiana = patternsPOLSAR_gaussiana[ ,-1] # conjunto de treinamento para New Method
# usando a rede
    resGaussiana_NM = PNN(patternsPOLSAR_gaussiana, patternsPOLSARLabels_gaussiana, testPattern, sigma=0.05 ,kernel=1)
    medias_gaussiana_NM[i] = mean(resGaussiana_NM == labelstest)
  }
    
  if(t.test(medias_gaussiana, medias_gaussiana_NM, mu=0.04, alternative = "less")$p.value < 0.05)
    break;
    
  print(K)
}  
  
melhor_K = K-50
  
# classificar toda a imagem com novo método 
par(mfrow=c(1,2))
showBand(CANAISequalizados)
  
patternsPOLSAR = ExampleSelection(patternsPOLSARinit, k=melhor_K, classes=3)
patternsPOLSARLabels = patternsPOLSAR[ ,1]
patternsPOLSAR = patternsPOLSAR[ ,-1] # conjunto de treinamento para New Method; dim(patternsNM)
  
imagemClassificada = PNN(patternsPOLSAR, patternsPOLSARLabels, im, sigma=0.05, kernel=1)
par(mar=c(0,0,0,0))
showBand(matrix(imagemClassificada, nrow = nrow(CANAISequalizados), byrow = FALSE), Pmax = 1)

# voto majoritário com janela 3x3
img_total_VM = windows_VM(matrix(imagemClassificada, nrow = nrow(CANAISequalizados), byrow = FALSE)[1:228,1:435], 3)
showBand(img_total_VM, Pmax = 1)

# voto majoritário com janela 5x5
img_total_VM = windows_VM(matrix(imagemClassificada, nrow = nrow(CANAISequalizados), byrow = FALSE)[1:225,1:435], 5)
showBand(img_total_VM, Pmax = 1)

# tempo de processamento
    
gg1 = microbenchmark(
  tradicional = PNN(patterns, labels, im, kernel = 3),
  novo_metodo = PNN(patternsPOLSAR, patternsPOLSARLabels, im, kernel=1),
  times = 1L)
    
gg1
  
*/
