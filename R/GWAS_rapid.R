#' Rapid GWAS with PCA
#'
#' @param pheno file with numeric phenotypic values
#' @param geno data.frame with genotype calls coded as 0,1,2.
#' @param Cov numeric data.frame with covariates values
#' @param GM genetic map of data with chr and position of each SNP
#' @param PCA.M number of principal components to use default is 3
#' @param cutoff  If cutoff is default, uses Bonferroni; 0.05/number of SNPs
#' @return GWAS P-values
GWASapply_rapid<- function(pheno=NULL, geno=NULL, Cov=NULL, GM=NULL, PCA.M=3,cutoff=NULL){
  GD=geno
  n=nrow(GD)
  m=ncol(GD)
  CV=Cov[,-1]
  y=pheno
  PCA=prcomp(GD)
  P=apply(GD,2, function(x)
    # CHECK FOR GENOTYPE DISTRIBUTION
    if(max(x)==min(x)){p=1
    P=p[length(p)]
    # IF MAX NOT EQUAL TO MIN
    }else{
      # CHECK FOR DEPENDENCE
      # IF NO CV INPUT
      if (is.null(CV)){X=cbind(1, PCA$x[,1:PCA.M],x)
      # IF THERE IS CV INPUT
      }else{fD <- fix_Dep(PCA$x[,1:PCA.M], as.matrix(CV))
      # WITH CV INPUT, CONDITION 1: NO DEPENDENCE
      if (is.null(fD)){X=cbind(1, PCA$x[,1:PCA.M],CV, x)

      # WITH CV INPUT, CONDITION 2: WITH DEPENDENCE
      }else {X=cbind(1, PCA$x[,1:PCA.M],x)}
      }# END FOR DEPENDECE
      # SOLVE THE LINEAR REGRESSION MATRIX:
      LHS=t(X)%*%X
      C=solve(LHS)
      RHS=t(X)%*%y
      b=C%*%RHS
      yb=X%*%b
      e=y-yb
      n=length(y)
      ve=sum(e^2)/(n-1)
      vt=C*ve
      t=b/sqrt(diag(vt))
      p=2*(1-pt(abs(t),n-2))
      P=p[length(p)]
    }# END FOR MAX NOT EQUAL TO MIN
    # COLLECT P-VALUES
  ) #end of appply function for markers
  P=t(matrix(P))
  P.value=P
  GWAS.Results=list(P.value.res=P.value)
  return(GWAS.Results)
}
