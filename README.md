# orthomr

Orthogonal polynomials for non-linear Mendelian randomization.

## Installation

You can install the development version of orthomr from GitHub:
``` r
# install.packages("devtools")
remotes::install_github("DrGBL/orthomr")
```

## Example

####### simulate #######
``` r
test<-sim_exp_out_non_linear(pleiotropy=TRUE,
                             n_s=10000,
                             n_ins=30,
                             func=as.function(polynomial(coef = c(1,10,3,0))),
                             func_conf_x=as.function(polynomial(coef = c(3,9,-6,1))),
                             func_conf_y=as.function(polynomial(coef = c(3,-2,4))))
Y<-test[["sims"]]$Y
X<-test[["sims"]]$X
Y_true<-test[["sims"]]$Y_true
X_true<-test[["sims"]]$X_true
G<-test[["geno"]]
af_ins<-test[["af_ins"]]


###### do the gwass for exposures and outcome #######
#orthogonal polynomials
n_deg<-5
polyX<-orthopol(degree=n_deg,X=X)

#exposures
gwas_exp<-list()
for(j in 1:n_deg){
  gwas_tmp<-c()
  for(i in 1:ncol(G)){
    mod<-lm(polyX[["values"]][,j]~G[,i])
    gwas_tmp<-data.frame(ins=paste0("X",i),
                         beta=coefficients(mod)[2],
                         se=coef(summary(mod))[, "Std. Error"][2],
                         pval=coef(summary(mod))[, "Pr(>|t|)"][2],
                         af=af_ins[i]) %>%
      bind_rows(gwas_tmp,.)
  }
  gwas_exp[[j]]<-gwas_tmp
}

#outcomes
gwas_out<-c()
for(i in 1:ncol(G)){
  mod<-lm(Y~G[,i])
  gwas_out<-data.frame(ins=paste0("X",i),
                           beta=coefficients(mod)[2],
                           se=coef(summary(mod))[, "Std. Error"][2],
                           pval=coef(summary(mod))[, "Pr(>|t|)"][2],
                           af=af_ins[i]) %>%
    bind_rows(gwas_out,.)
}

#putting them together
mvmr_exp<-merge_exp_inst(gwas_exp)

mvmr_out<-gwas_out %>%
  filter(ins %in% mvmr_exp$SNP) %>%
  rename(SNP=ins) %>%
  mutate(outcome="outcome") %>%
  mutate(id.outcome="outcome") %>%
  mutate(effect_allele.outcome="A") %>%
  mutate(other_allele.outcome="C") %>%
  rename(beta.outcome=beta) %>%
  rename(se.outcome=se) %>%
  rename(pval.outcome=pval) %>%
  rename(eaf.outcome=af)


#now do the actual MRs
sim_res<-mr_sim_res(mvmr_exp=mvmr_exp,
                    mvmr_out=mvmr_out,
                    test=test,
                    X=NULL,
                    inner_int=0.90,
                    intercept_choice=TRUE,
                    polyX_coef=polyX$coef)
sim_res$final_plots
```

## License

GPL (>= 3)
