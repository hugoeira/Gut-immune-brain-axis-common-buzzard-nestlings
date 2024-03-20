# Bayesian Structural Equation Modelling - modelling Haptoglobin immune assay



------

## Table of contents

[TOC]

------





## A) 16S rRNA (bacterial microbiota) SEM analysis



### 1. Define SEM for each diversity measurement  

```R
#Load Packages
library(brms)
library(rstan)
library(bayesplot)
library(bayestestR)
library(parallel)
library(svglite)
library(ggplot2)

#Load the data
metadata <- readRDS("16s_metadata_immune.rds")

#scale all predictors to range between 0-1 if they are not already naturally on that scale

#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("bci_two", "shannon_entropy", "faith_pd", "observed_features", "cort", "hapto", "age_days")


for(i in 1:ncol(metadata[,which(colnames(metadata)%in%scalecols)])){
  metadata[,which(colnames(metadata)%in%scalecols)][,i]<-range.use(metadata[,which(colnames(metadata)%in%scalecols)][,i],0,1)
}

# Define structural equation model paths

#Shannon
path1 <- bf(bci_two ~ shannon_entropy + cort + hapto+ age_days + (1|nest/ring_number))
path2 <- bf(hapto~  shannon_entropy + cort + age_days + (1|nest/ring_number))
path3 <- bf(shannon_entropy ~ cort + age_days + (1|nest/ring_number))
path4 <- bf(cort ~ age_days + (1|nest/ring_number)) + skew_normal()

sem_hapto_shannon <- path1 + path2 + path3 + path4

#Faith PD
path1 <- bf(bci_two ~ faith_pd + cort + hapto+ age_days + (1|nest/ring_number))
path2 <- bf(hapto~  faith_pd + cort + age_days + (1|nest/ring_number))
path3 <- bf(faith_pd ~ cort + age_days + (1|nest/ring_number))
path4 <- bf(cort ~ age_days + (1|nest/ring_number)) + skew_normal()

sem_hapto_faith <- path1 + path2 + path3 + path4

#N° of observed ASV's
path1 <- bf(bci_two ~ observed_features + cort + hapto+ age_days + (1|nest/ring_number))
path2 <- bf(hapto~  observed_features + cort + age_days + (1|nest/ring_number))
path3 <- bf(observed_features ~ cort + age_days + (1|nest/ring_number))
path4 <- bf(cort ~ age_days + (1|nest/ring_number)) + skew_normal()

sem_hapto_asv <- path1 + path2 + path3 + path4

```



### 2. Run brms

```R
ncores = detectCores()
options(mc.cores = parallel::detectCores())

#Shannon
model_hapto_shannon <-brm(sem_hapto_shannon + set_rescor(FALSE),
                         data = metadata,
                         warmup = 50000, iter = 100000,
                         control = list(adapt_delta = 0.99, max_treedepth = 15),
                         cores=ncores, chains=4, init=1000)
                         
#Faith PD                         
model_hapto_faith <-brm(sem_hapto_faith + set_rescor(FALSE),
                         data = metadata,
                         warmup = 50000, iter = 100000,
                         control = list(adapt_delta = 0.99, max_treedepth = 15),
                         cores=ncores, chains=4, init=1000)
                         
#N° of observed ASV's                     
model_hapto_asv <-brm(sem_hapto_asv + set_rescor(FALSE),
                         data = metadata,
                         warmup = 50000, iter = 100000,
                         control = list(adapt_delta = 0.99, max_treedepth = 15),
                         cores=ncores, chains=4, init=1000)
```



### 3. Model Diagnostics Shannon

#### 3.1 Model Summary

```R
#Model summary
summary_shannon<- summary(model_hapto_shannon)

#Bayes R2 
R2m_shannon <- bayes_R2(model_hapto_shannon,re_formula=NA)
R2c_shannon <- bayes_R2(model_hapto_shannon)
```



#### 3.2 Model diagnostics 

```R
# Model diagnostics 
diagnostic_shannon <- plot(model_hapto_shannon)

#Loop to save all diagnostic plots
diagnostic_plots <- list()
for (i in 1:length(diagnostic_shannon)) {
  diagnostic_plots[[i]] <- diagnostic_shannon[[i]]
  filename <- paste0("diagnostic", i, "_shannon")
  ggsave(filename = paste0(filename, ".png"), plot = diagnostic_plots[[i]], device = "png", dpi=300)
}
```

![16s_diagnostic_shannon](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/16s/Hapto/16s_diagnostic_shannon.png)



#### 3.3 Compare distribution of response variable to distributions of predicted response variable
```R
#Posterior predictive checks (one by one)
distribution_shannon_path1 <- pp_check(model_hapto_shannon, resp="bcitwo", ndraws=200)
distribution_shannon_path2 <- pp_check(model_hapto_shannon, resp="hapto", ndraws=200)
distribution_shannon_path3 <- pp_check(model_hapto_shannon, resp="shannonentropy", ndraws=200)
distribution_shannon_path4 <- pp_check(model_hapto_shannon, resp="cort", ndraws=200)

#Loop to save all distributions plot
responses <- c("bcitwo", "hapto", "shannonentropy", "cort")
response_names <- c("bci", "immune", "shannon", "cort")
for (i in seq_along(responses)) {
  pp_check_plot <- pp_check(model_hapto_shannon, resp = responses[i], ndraws = 200)
  filename <- paste0("distribution_shannon_", response_names[i])
  ggsave(filename = paste0(filename, ".svg"), plot = pp_check_plot, device = "svg", width = 8, height = 10)
}
```

![16s_posterior_shannon](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/16s/Hapto/16s_posterior_shannon.svg)



#### 3.4 Plot model posterior and credible intervals

```R
#Plot Model effects
plot1 <-mcmc_plot(model_hapto_shannon, type = "intervals",prob_outer=0.95, prob=0.95, 
                  variable =c("b_bcitwo_shannon_entropy", "b_bcitwo_cort", "b_bcitwo_hapto", "b_bcitwo_age_days",
							   "b_hapto_shannon_entropy", "b_hapto_cort", "b_hapto_age_days",
							   "b_shannonentropy_cort", "b_shannonentropy_age_days",
							   "b_cort_age_days")) 

plot1 <- plot1 + theme_classic() + geom_vline(xintercept = 0, linetype="dotted", color="blue")+
theme(axis.text.x = element_text(size = 16),  # Adjust the size as needed
axis.text.y = element_text(size = 16))+
theme(text = element_text(family = "Arial"))

ggsave(filename="16s_effect_sizes_shannon.svg", plot=plot1, device = "svg", width = 8, height = 10)
```



### 4. Model diagnostics - Faith PD

#### 4.1 Model summary

```R
#Model summary
summary_faith<- summary(model_hapto_faith)

#Bayes R2 
R2m_faith <- bayes_R2(model_hapto_faith,re_formula=NA)
R2c_faith <- bayes_R2(model_hapto_faith)
```



#### 4.2 Model diagnostics 

```R
# Model diagnostics
diagnostic_faith <- plot(model_hapto_faith)

#Loop to save all diagnostic plots
diagnostic_plots <- list()
for (i in 1:length(diagnostic_faith)) {
  diagnostic_plots[[i]] <- diagnostic_faith[[i]]
  filename <- paste0("diagnostic", i, "_faith")
  ggsave(filename = paste0(filename, ".png"), plot = diagnostic_plots[[i]], device = "png", dpi=200)
}
```

![16s_diagnostic_faith](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/16s/Hapto/16s_diagnostic_faith.png)



#### 4.3 Compare distribution of response variable to distributions of predicted response variable

```R
#Posterior predictive checks 
#Loop to save all distributions plot
responses <- c("bcitwo", "hapto", "faithpd", "cort")
response_names <- c("bci", "immune", "faith", "cort")
for (i in seq_along(responses)) {
  pp_check_plot <- pp_check(model_hapto_shannon, resp = responses[i], ndraws = 200)
  filename <- paste0("distribution_faith_", response_names[i])
  ggsave(filename = paste0(filename, ".svg"), plot = pp_check_plot, device = "svg", width = 8, height = 10)
}
```

![16s_posterior_faith](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/16s/Hapto/16s_posterior_faith.svg)



#### 4.4 Plot model posterior and credible intervals

```R
#Plot Model effects
plot2 <-mcmc_plot(model_hapto_shannon, type = "intervals",prob_outer=0.95, prob=0.95, 
                  variable =c("b_bcitwo_faith_pd", "b_bcitwo_cort", "b_bcitwo_hapto", "b_bcitwo_age_days",
							   "b_hapto_faith_pd", "b_hapto_cort", "b_hapto_age_days",
							   "b_faithpd_cort", "b_faithpd_age_days",
							   "b_cort_age_days")) 

plot2 <- plot2 + theme_classic() + geom_vline(xintercept = 0, linetype="dotted", color="blue")+
theme(axis.text.x = element_text(size = 16),  # Adjust the size as needed
axis.text.y = element_text(size = 16))+
theme(text = element_text(family = "Arial"))

ggsave(filename="16s_effect_sizes_faith.svg", plot=plot1, device = "svg", width = 8, height = 10)
```



### 5. Model diagnostics - N° of observed ASV's

#### 5.1 Model summary

```R
#Model summary
summary_asv<- summary(model_hapto_asv)

#Bayes R2 
R2m_asv <- bayes_R2(model_hapto_asv,re_formula=NA)
R2c_asv <- bayes_R2(model_hapto_asv)
```



#### 5.2 Model diagnostics 

```R
# Model diagnostics
diagnostic_asv <- plot(model_hapto_asv)

#Loop to save all diagnostic plots
diagnostic_plots <- list()
for (i in 1:length(diagnostic_asv)) {
  diagnostic_plots[[i]] <- diagnostic_asv[[i]]
  filename <- paste0("diagnostic", i, "_asv")
  ggsave(filename = paste0(filename, ".png"), plot = diagnostic_plots[[i]], device = "png", dpi= 300)
}
```

![16s_diagnostic_asv](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/16s/Hapto/16s_diagnostic_asv.png)



#### 5.3 Compare distribution of response variable to distributions of predicted response variable

```R
#Posterior predictive checks 
#Loop to save all distributions plot
responses <- c("bcitwo", "hapto", "observedfeatures", "cort")
response_names <- c("bci", "immune", "asv", "cort")
for (i in seq_along(responses)) {
  pp_check_plot <- pp_check(model_hapto_asv, resp = responses[i], ndraws = 200)
  filename <- paste0("distribution_asv_", response_names[i])
  ggsave(filename = paste0(filename, ".svg"), plot = pp_check_plot, device = "svg", width = 8, height = 10)
}
```

![16s_posterior_asv](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/16s/Hapto/16s_posterior_asv.svg)



#### 5.4 Plot model posterior and credible intervals

```R
#Plot Model effects
plot3 <-mcmc_plot(model_hapto_shannon, type = "intervals",prob_outer=0.95, prob=0.95, 
                  variable =c("b_bcitwo_observed_features", "b_bcitwo_cort", "b_bcitwo_hapto", "b_bcitwo_age_days",
							   "b_hapto_observed_features", "b_hapto_cort", "b_hapto_age_days",
							   "b_observedfeatures_cort", "b_observedfeatures_age_days",
							   "b_cort_age_days")) 

plot3 <- plot3 + theme_classic() + geom_vline(xintercept = 0, linetype="dotted", color="blue")+
theme(axis.text.x = element_text(size = 16),  # Adjust the size as needed
axis.text.y = element_text(size = 16))+
theme(text = element_text(family = "Arial"))

ggsave(filename="16s_effect_sizes_asv.svg", plot=plot1, device = "svg", width = 8, height = 10)
```





## B) 28S rRNA (eukaryotic microbiota) SEM analysis



### 1. Define SEM for each diversity measurement 

```R
#Load the data
metadata <- readRDS("28s_metadata_immune.rds")

#scale all predictors to range between 0-1 if they are not already naturally on that scale

#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("bci_two", "shannon_entropy", "faith_pd", "observed_features", "cort", "hapto", "age_days")


for(i in 1:ncol(metadata[,which(colnames(metadata)%in%scalecols)])){
  metadata[,which(colnames(metadata)%in%scalecols)][,i]<-range.use(metadata[,which(colnames(metadata)%in%scalecols)][,i],0,1)
}

# Define structural equation model paths

#Shannon
path1 <- bf(bci_two ~ shannon_entropy + cort + hapto+ age_days + (1|nest/ring_number))
path2 <- bf(hapto~  shannon_entropy + cort + age_days + (1|nest/ring_number))
path3 <- bf(shannon_entropy ~ cort + age_days + (1|nest/ring_number))
path4 <- bf(cort ~ age_days + (1|nest/ring_number)) + skew_normal()

sem_hapto_shannon <- path1 + path2 + path3 + path4

#Faith PD
path1 <- bf(bci_two ~ faith_pd + cort + hapto+ age_days + (1|nest/ring_number))
path2 <- bf(hapto~  faith_pd + cort + age_days + (1|nest/ring_number))
path3 <- bf(faith_pd ~ cort + age_days + (1|nest/ring_number))
path4 <- bf(cort ~ age_days + (1|nest/ring_number)) + skew_normal()

sem_hapto_faith <- path1 + path2 + path3 + path4

#N° of observed ASV's
path1 <- bf(bci_two ~ observed_features + cort + hapto+ age_days + (1|nest/ring_number))
path2 <- bf(hapto~  observed_features + cort + age_days + (1|nest/ring_number))
path3 <- bf(observed_features ~ cort + age_days + (1|nest/ring_number))
path4 <- bf(cort ~ age_days + (1|nest/ring_number)) + skew_normal()

sem_hapto_asv <- path1 + path2 + path3 + path4

```



### 2. Run brms

```R
ncores = detectCores()
options(mc.cores = parallel::detectCores())

#Shannon
model_hapto_shannon <-brm(sem_hapto_shannon + set_rescor(FALSE),
                         data = metadata,
                         warmup = 50000, iter = 100000,
                         control = list(adapt_delta = 0.99, max_treedepth = 15),
                         cores=ncores, chains=4, init=1000)
                         
#Faith PD                         
model_hapto_faith <-brm(sem_hapto_faith + set_rescor(FALSE),
                         data = metadata,
                         warmup = 50000, iter = 100000,
                         control = list(adapt_delta = 0.99, max_treedepth = 15),
                         cores=ncores, chains=4, init=1000)
                         
#N° of observed ASV's                     
model_hapto_asv <-brm(sem_hapto_asv + set_rescor(FALSE),
                         data = metadata,
                         warmup = 50000, iter = 100000,
                         control = list(adapt_delta = 0.99, max_treedepth = 15),
                         cores=ncores, chains=4, init=1000)
```





### 3. Model Diagnostics - Shannon

#### 3.1 Model Summary

```R
#Model summary
summary_shannon<- summary(model_hapto_shannon)

#Bayes R2 
R2m_shannon <- bayes_R2(model_hapto_shannon,re_formula=NA)
R2c_shannon <- bayes_R2(model_hapto_shannon)
```



#### 3.2 Model diagnostics 

```R
# Model diagnostics 
diagnostic_shannon <- plot(model_hapto_shannon)

#Loop to save all diagnostic plots
diagnostic_plots <- list()
for (i in 1:length(diagnostic_shannon)) {
  diagnostic_plots[[i]] <- diagnostic_shannon[[i]]
  filename <- paste0("diagnostic", i, "_shannon")
  ggsave(filename = paste0(filename, ".png"), plot = diagnostic_plots[[i]], device = "png", dpi=300)
}
```

![28s_diagnostic_shannon](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/28s/Hapto/28s_diagnostic_shannon.png)



#### 3.3 Compare distribution of response variable to distributions of predicted response variable

```R
#Posterior predictive checks (one by one)
distribution_shannon_path1 <- pp_check(model_hapto_shannon, resp="bcitwo", ndraws=200)
distribution_shannon_path2 <- pp_check(model_hapto_shannon, resp="hapto", ndraws=200)
distribution_shannon_path3 <- pp_check(model_hapto_shannon, resp="shannonentropy", ndraws=200)
distribution_shannon_path4 <- pp_check(model_hapto_shannon, resp="cort", ndraws=200)

#Loop to save all distributions plot
responses <- c("bcitwo", "hapto", "shannonentropy", "cort")
response_names <- c("bci", "immune", "shannon", "cort")
for (i in seq_along(responses)) {
  pp_check_plot <- pp_check(model_hapto_shannon, resp = responses[i], ndraws = 200)
  filename <- paste0("distribution_shannon_", response_names[i])
  ggsave(filename = paste0(filename, ".svg"), plot = pp_check_plot, device = "svg", width = 8, height = 10)
}
```

![28s_posterior_shannon](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/28s/Hapto/28s_posterior_shannon.svg)



#### 3.4 Plot model posterior and credible intervals

```R
#Plot Model effects
plot4 <-mcmc_plot(model_hapto_shannon, type = "intervals",prob_outer=0.95, prob=0.95, 
                  variable =c("b_bcitwo_shannon_entropy", "b_bcitwo_cort", "b_bcitwo_hapto", "b_bcitwo_age_days",
							   "b_hapto_shannon_entropy", "b_hapto_cort", "b_hapto_age_days",
							   "b_shannonentropy_cort", "b_shannonentropy_age_days",
							   "b_cort_age_days")) 

plot4 <- plot4 + theme_classic() + geom_vline(xintercept = 0, linetype="dotted", color="blue")+
theme(axis.text.x = element_text(size = 16),  # Adjust the size as needed
axis.text.y = element_text(size = 16))+
theme(text = element_text(family = "Arial"))

ggsave(filename="16s_effect_sizes_shannon.svg", plot=plot4, device = "svg", width = 8, height = 10)
```





### 4. Model diagnostics - Faith PD

#### 4.1 Model summary

```R
#Model summary
summary_faith<- summary(model_hapto_faith)

#Bayes R2 
R2m_faith <- bayes_R2(model_hapto_faith,re_formula=NA)
R2c_faith <- bayes_R2(model_hapto_faith)
```



#### 4.2 Model diagnostics 

```R
# Model diagnostics
diagnostic_faith <- plot(model_hapto_faith)

#Loop to save all diagnostic plots
diagnostic_plots <- list()
for (i in 1:length(diagnostic_faith)) {
  diagnostic_plots[[i]] <- diagnostic_faith[[i]]
  filename <- paste0("diagnostic", i, "_faith")
  ggsave(filename = paste0(filename, ".png"), plot = diagnostic_plots[[i]], device = "png", dpi=200)
}
```

![28s_diagnostic_faith](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/28s/Hapto/28s_diagnostic_faith.png)



#### 4.3 Compare distribution of response variable to distributions of predicted response variable

```R
#Posterior predictive checks 
#Loop to save all distributions plot
responses <- c("bcitwo", "hapto", "faithpd", "cort")
response_names <- c("bci", "immune", "faith", "cort")
for (i in seq_along(responses)) {
  pp_check_plot <- pp_check(model_hapto_faith, resp = responses[i], ndraws = 200)
  filename <- paste0("distribution_faith_", response_names[i])
  ggsave(filename = paste0(filename, ".svg"), plot = pp_check_plot, device = "svg", width = 8, height = 10)
}
```

![28s_posterior_faith](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/28s/Hapto/28s_posterior_faith.svg)



#### 4.4 Plot model posterior and credible intervals

```R
#Plot Model effects
plot5 <-mcmc_plot(model_hapto_shannon, type = "intervals",prob_outer=0.95, prob=0.95, 
                  variable =c("b_bcitwo_faith_pd", "b_bcitwo_cort", "b_bcitwo_hapto", "b_bcitwo_age_days",
							   "b_hapto_faith_pd", "b_hapto_cort", "b_hapto_age_days",
							   "b_faithpd_cort", "b_faithpd_age_days",
							   "b_cort_age_days")) 

plot5 <- plot5 + theme_classic() + geom_vline(xintercept = 0, linetype="dotted", color="blue")+
theme(axis.text.x = element_text(size = 16),  # Adjust the size as needed
axis.text.y = element_text(size = 16))+
theme(text = element_text(family = "Arial"))

ggsave(filename="16s_effect_sizes_faith.svg", plot=plot5, device = "svg", width = 8, height = 10)
```



### 5. Model diagnostics - N° of observed ASV's

#### 5.1 Model summary

```R
#Model summary
summary_asv<- summary(model_hapto_asv)

#Bayes R2 
R2m_asv <- bayes_R2(model_hapto_asv,re_formula=NA)
R2c_asv <- bayes_R2(model_hapto_asv)
```



#### 5.2 Model diagnostics 

```R
# Model diagnostics
diagnostic_asv <- plot(model_hapto_asv)

#Loop to save all diagnostic plots
diagnostic_plots <- list()
for (i in 1:length(diagnostic_asv)) {
  diagnostic_plots[[i]] <- diagnostic_asv[[i]]
  filename <- paste0("diagnostic", i, "_asv")
  ggsave(filename = paste0(filename, ".png"), plot = diagnostic_plots[[i]], device = "png", dpi= 300)
}
```

![28s_diagnostic_asv](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/28s/Hapto/28s_diagnostic_asv.png)



#### 5.3 Compare distribution of response variable to distributions of predicted response variable

```R
#Posterior predictive checks 
#Loop to save all distributions plot
responses <- c("bcitwo", "hapto", "observedfeatures", "cort")
response_names <- c("bci", "immune", "asv", "cort")
for (i in seq_along(responses)) {
  pp_check_plot <- pp_check(model_hapto_asv, resp = responses[i], ndraws = 200)
  filename <- paste0("distribution_asv_", response_names[i])
  ggsave(filename = paste0(filename, ".svg"), plot = pp_check_plot, device = "svg", width = 8, height = 10)
}
```

![28s_posterior_asv](/home/localadmin/microbiome-analysis/Path_analysis/alpha-diversity/28s/Hapto/28s_posterior_asv.svg)



#### 5.4 Plot model posterior and credible intervals

```R
#Plot Model effects
plot2 <-mcmc_plot(model_hapto_shannon, type = "intervals",prob_outer=0.95, prob=0.95, 
                  variable =c("b_bcitwo_faith_pd", "b_bcitwo_cort", "b_bcitwo_hapto", "b_bcitwo_age_days",
							   "b_hapto_faith_pd", "b_hapto_cort", "b_hapto_age_days",
							   "b_faithpd_cort", "b_faithpd_age_days",
							   "b_cort_age_days")) 

plot2 <- plot2 + theme_classic() + geom_vline(xintercept = 0, linetype="dotted", color="blue")+
theme(axis.text.x = element_text(size = 16),  # Adjust the size as needed
axis.text.y = element_text(size = 16))+
theme(text = element_text(family = "Arial"))

ggsave(filename="16s_effect_sizes_faith.svg", plot=plot1, device = "svg", width = 8, height = 10)
```



