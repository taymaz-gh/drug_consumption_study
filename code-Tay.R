

#======================================
#        Dimension Reduction
#======================================

rm(list=ls())
cat('\14')
dev.off()



drug <- read.csv('drug_consumption.data', header = F)
colnames(drug) <- c('id','age','gender','education','country','ethnicity',
                    'nscore','escore','oscore','ascore','cscore',
                    'impulsive','ss','alcohol','amphet',
                    'amyl','benzos','caff','cannabis','choc',
                    'coke','crack','ecstasy','heroin','ketamine','legalh','lsd',
                    'meth','mushrooms','nicotine','semer','vsa')
data <- drug


for(i in 2:6){data[,i] = as.factor(data[,i])}


data$age <- factor(data$age,
                   levels = c(-0.95197, -0.07854, 0.49788, 1.09449, 1.82213, 2.59171),
                   labels = c("18-24", "25-34", "35-44", "45-54", "55-64", "65+"))

unique(data$heroin)
unique(drug$age)
round(prop.table(table(data$age))*100,1)
round(prop.table(table(data$age, data$amphet))*100,1)
str(data)




#------------------------------------------------
# Step 1:
#-------------------------------------------

Addict_Stage <- data[,c(1,14:32)]

Addict_Stage[,-1] <- lapply(Addict_Stage[,-1], function(x) as.numeric(sub("CL", "", x)))





table(table(Addict_Stage$id))

addict_stage <- Addict_Stage[,-1]

library(nFactors)       
library(psych)           
library(EFA.dimensions) 
library(factoextra)     
library(cluster)        
library(tidyverse)      

#-------------------------------------------
# Step 2: 
# Performing PCA (with pre-scaled variables)
#-------------------------------------------


pca_result <- prcomp(addict_stage, scale = TRUE)   
pca_result$center
summary(pca_result)




#-------------------------------------------
# Step 3: 
# Explained Variance Summary
#-------------------------------------------

# Saving explained variance:
prop_var_expl <- summary(pca_result)$importance[2,]  # Proportion of variance explained
cumul_prop_var_expl <- summary(pca_result)$importance[3,]  # Cumulative proportion of variance explained



#-------------------------------------------
# Step 4: Scree Plot and Cumulative Plot
#-------------------------------------------


ls(pca_result)
ls(summary(pca_result))
rownames(summary(pca_result)$importance)
colnames(t(summary(pca_result)$importance))
pca_result_summary_imp = as.data.frame(t(summary(pca_result)$importance))
colnames(pca_result_summary_imp) = c('StdDev','PropVar','CumProp')
rownames(pca_result_summary_imp) = paste0("PC", 1:19)

barplot(pca_result_summary_imp$CumProp,
        names.arg = rownames(pca_result_summary_imp),
        col = "orange",
        main = "Scree Plot & Cummulative Plot",
        xlab = "Principle Component",
        ylab = "Explained Proportion of Variance",
        cex.names = 0.6,  # shrinks labels to 60% size
        las = 2)  ->  midpoints
# Adding another barplot for Proportion of Variance
barplot(pca_result_summary_imp$PropVar,
        col = rgb(0, 0, 1, 0.5),  # semi-transparent blue
        cex.names = 0.6,  # shrinks labels to 60% size
        las = 2,
        add = TRUE)
abline(h = 0.7 , lty = 3, col = 'orangered')
abline(h = 0.8 , lty = 2, col = 'red')
abline(h = 0.9 , lty = 3, col = 'orangered')
legend("topleft",
       legend = c("Cumulative Proportion", "Proportion of Variance"),
       fill = c("darkorange", "slateblue"),
       cex = 0.6,horiz = TRUE)
points(midpoints+0.1, pca_result_summary_imp$CumProp, pch = 19, col = "red")
points(midpoints-0.1, pca_result_summary_imp$PropVar, pch = 19, col = "blue")
text(x = par("usr")[2] + 1,         # a bit right of the plot edge
     y = mean(par("usr")[3:4] + .17),     # vertical center of the plot
     labels = "7 to 15 PCs",
     srt = 90,                      # rotate 90 degrees
     adj = 0,                       # bottom-aligned
     col = "tomato",
     font = 2,
     cex = .9,
     xpd = NA)                      # allow drawing outside plot

#-------------------------------------------
# Step 5: Elbow Method - Ratio Calculation
#-------------------------------------------

nElem <- length(cumul_prop_var_expl)  # Number of components

ratios <- (prop_var_expl[2:(nElem-1)]) / (prop_var_expl[3:nElem])

nComponents_selected <- which.max(ratios) + 1
nComponents_selected

plot(2:(nElem - 1), ratios,
     type = "b",
     xlab = "Principal Component",
     ylab = expression(Ratio == ~ lambda[~i] ~~ "/" ~~ lambda[~i + 1]),
     main = "Elbow Detection via Successive `Explained Variance` Ratios",
     pch = 19,
     col = "darkblue")
abline(v = nComponents_selected, col = "red", lty = 2)
abline(v = 7, col = "orangered", lty = 3)
axis(1, at = 7, labels = "7")
# The result shows that it is optimal to keep the first 10 PCs among 19.


#-------------------------------------------
# Step 6: Kaiser Rule - Eigenvalues > 1
#-------------------------------------------


ev <- eigen(cor(addict_stage))

ev$values

plot(ev$values, col='blue', pch = 19, type = "b", lty = 3,
     xlab = "Principle Component",
     ylab = "Eigenvalue :  λ",
     main = "Kaiser Rule: keep those with λ > 1",)
mtext("A Standardized PC with  λ < 1  explains less variance than a single variable.",
      side = 3, line = 0.5, cex = 0.75)
abline(h = 1 , lty = 2, col = 'red')
abline(h = 0.9 , lty = 3, col = 'orangered')
text(17, 1.25, expression(lambda==1), col = 'red', cex = 0.85)
text(18, 0.7, expression(lambda==.9), col = 'orangered', cex = 0.85)
axis(1, at = 6, labels = "6")
axis(1, at = 7, labels = "7")



#-------------------------------------------
# Step 7: Horn's Parallel Analysis
#-------------------------------------------

ap <- parallel(subject = nrow(addict_stage), var = ncol(addict_stage), rep = 1000)

nS <- nScree(x = ev$values, aparallel = ap$eigen$qevpea)

plotnScree(nS, main = "Parallel Analysis")



#-------------------------------------------
# Step 8: Velicer's MAP Test
#-------------------------------------------


map = MAP(addict_stage, corkind = 'pearson', verbose = TRUE)
ls(map)
map$avgsqrs

plot(map$avgsqrs[, "root"],
     map$avgsqrs[, "Avg.Corr.Sq."],
     type = "b",
     xlab = "Number of Components Removed",
     ylab = "Average Squared Partial Correlation",
     main = "Velicer's MAP Test",
     pch = 19,
     col = "darkblue")
mtext("PCA already removed Shared variance mostly - Unique variance or Noise remains",
      side = 3, line = 0.5, cex = 0.75, col = 'navy')
abline(v = 1 , lty = 2, col = 'red')
abline(v = 7 , lty = 3, col = 'orangered')
abline(h = 0.1 , lty = 4, col = 'steelblue')
axis(1, at = 1, labels = "1")
axis(1, at = 7, labels = "7")
text(x = par("usr")[2] - 6,             
     y = mean(par("usr")[3:4]),          
     labels = "Chance of Overfit",
     srt = 45,                            
     col = "deeppink",
     cex = 0.9,
     font = 1,
     xpd = NA)                            


#-------------------------------------------
# Step 9: Very Simple Structure (VSS) Method
#-------------------------------------------

fit_vss <- vss(addict_stage, n = 15, rotate = "oblimin",
               diagonal = FALSE, fm = "pc", complexity = 2)  # fm='pc' does not work with complexity>2
fit_vss

ls(fit_vss)

x <- 1:length(fit_vss$cfit.1)
vss1 <- fit_vss$cfit.1
vss2 <- fit_vss$cfit.2
# Plot with manual styling
plot(x, vss1, type = "b", col = "forestgreen", ylim = c(0, 1),
     pch = 19, lty = 1,
     xlab = "Number of Components",
     ylab = "VSS Fit",
     main = "Very Simple Structure Method (Complexity 1 to 2)")
lines(x, vss2, type = "b", col = "navy", pch = 17, lty = 2)
abline(v = 12 , lty = 2, col = 'red')
abline(v = 7 , lty = 3, col = 'orangered')
axis(1, at = 7, labels = "7")
text(x = par("usr")[2] - 2,             
     y = mean(par("usr")[3:4] - 0.06),           
     labels = "better fit - strongly \n correlated to one \n or a few factors",
     srt = 90,                           
     col = "deeppink",
     cex = 0.8,
     font = 1,
     xpd = NA)                            
mtext("Each variable may simply load only on up to 2 factors",
      side = 3, line = 0.5, cex = 0.75, col = 'navy')
mtext("Oblique rotation - allowing for non-orthogonal factors",
      side = 1, line = 4, cex = 0.75, col = 'navy')
# Add legend
legend("bottomright",
       legend = c("Complexity 1", "Complexity 2"),
       col = c("forestgreen", "navy"),
       pch = c(19, 17),
       lty = c(1, 2),
       cex = 0.8)




#-------------------------------------------
# Final decision about Number of PCs to keep
#-------------------------------------------

# Final Decision (according to steps 4 to 9) : PC1 to PC7

numPC = 7

#-------------------------------------------
# Step 10: Varimax Rotation (Optional but recommended)
#-------------------------------------------

eigvals <- pca_result$sdev[1:numPC]^2
eigvecs <- pca_result$rotation[,1:numPC]


unrotated_loadings <- sweep(pca_result$rotation[,1:numPC], 2, sqrt(eigvals), "*")

varimax_result <- varimax(unrotated_loadings)

varimax_result




rot.pca_mat <- as.matrix(varimax_result$loadings[, 1:numPC])

colnames(rot.pca_mat) <- paste("RC", 1:ncol(rot.pca_mat), sep="")
rot.pca_mat

rot_mat <- varimax_result$rotmat
rownames(rot_mat) <- paste("PC", 1:nrow(rot_mat), sep="")
colnames(rot_mat) <- paste("RC", 1:ncol(rot_mat), sep="")
rot_mat


for (i in 1:numPC) { print (rot.pca_mat[ (rot.pca_mat[,i])^2 > 0.11187 , i , drop=F ] ) } 



RC1 <- c(
  amphet = -0.4723091,
  cannabis = -0.7605636,
  coke = -0.3994801,
  ecstasy = -0.6966802,
  ketamine = -0.4437304,
  legalh = -0.7167084,
  lsd = -0.8052531,
  mushrooms = -0.8053280,
  nicotine = -0.4182750
)

RC2 <- c(
  amphet = 0.4059220,
  benzos = 0.6333221,
  coke = 0.4387975,
  crack = 0.7223068,
  heroin = 0.8094254,
  meth = 0.7286953,
  vsa = 0.3344768
)

RC3 <- c(
  caff = 0.7775604,
  nicotine = 0.5357334
)

RC4 <- c(
  amphet = -0.3760088,
  amyl = -0.8529194,
  coke = -0.5163405,
  ecstasy = -0.4318624,
  ketamine = -0.5421936
)

RC5 <- c(
  semer = 0.9707074
)

RC6 <- c(
  choc = 0.9207448
)

RC7 <- c(
  alcohol = 0.92874
)


# sum of squared loadings for each RC
Variance_Explained <- colSums(rot.pca_mat^2)
Variance_Explained
Prop_variance_explained <- colSums(rot.pca_mat^2) / nrow(rot.pca_mat)
Prop_variance_explained
sum(Prop_variance_explained)




weights <- c(
  alcohol     = (1.9 + 1.9) / 2,
  amphet      = (1.1 + 2.2) / 2,
  amyl        = (0.1 + 0.2) / 2,
  benzos      = (1.8 + 2.1) / 2,
  cannabis    = (1.0 + 1.7) / 2,
  choc        = (0.2 + 0.5) / 2,
  coke        = (2.3 + 2.8) / 2,
  caff        = (1.0 + 1.3) / 2,
  crack       = (3.0 + 3.0) / 2,
  ecstasy     = (0.2 + 1.1) / 2,
  heroin      = (3.0 + 3.0) / 2,
  ketamine    = (1.2 + 1.9) / 2,
  legalh      = (1.0 + 1.5) / 2,
  lsd         = (0.3 + 1.1) / 2,
  meth        = (2.8 + 2.2) / 2,
  mushrooms   = (0.2 + 0.7) / 2,
  nicotine    = (2.3 + 2.6) / 2,
  vsa         = (1.1 + 2.0) / 2,
  semer       = 0.0  # not in the table; set to 0 or NA
)


w <- weights[row.names(varimax_result$loadings)]

RC_weights <- colSums((varimax_result$loadings^2) * w) 

RC_weights <- RC_weights / sum(RC_weights)

RC_weights <- RC_weights[order(names(RC_weights))]  
RC_weights





##@@@@@@@@@@@@@@@@@@@@ VISUALIZATION @@@@@@@@@@@@@@@@@@@@@##

##      `graph1--.R`      

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##





#=====================================================
#     Determining the Optimal Number of Clusters
#=====================================================




rm(list=ls())
cat('\14')
dev.off()



drug <- read.csv('drug_consumption.data', header = F)
colnames(drug) <- c('id','age','gender','education','country','ethnicity',
                    'nscore','escore','oscore','ascore','cscore',
                    'impulsive','ss','alcohol','amphet',
                    'amyl','benzos','caff','cannabis','choc',
                    'coke','crack','ecstasy','heroin','ketamine','legalh','lsd',
                    'meth','mushrooms','nicotine','semer','vsa')
data <- drug


for(i in 2:6){data[,i] = as.factor(data[,i])}


data$age <- factor(data$age,
                   levels = c(-0.95197, -0.07854, 0.49788, 1.09449, 1.82213, 2.59171),
                   labels = c("18-24", "25-34", "35-44", "45-54", "55-64", "65+"))

Addict_Stage <- data[,c(1,14:32)]

Addict_Stage[,-1] <- lapply(Addict_Stage[,-1], function(x) as.numeric(sub("CL", "", x)))



addict_stage <- Addict_Stage[,-1]

library(nFactors)       
library(psych)           
library(EFA.dimensions)  
library(factoextra)     
library(cluster)        
library(tidyverse)      


numPC = 7


#----------------------------

 
qr(as.matrix(addict_stage))$rank  # this dataset only contains features of interest (not the response) 


pca.Results <- principal(scale(addict_stage),  
                         nfactors = numPC,     
                         rotate   = "varimax", 
                         scores   = TRUE)



ls(pca.Results)
str(pca.Results)
attributes(pca.Results$loadings)

summary_table <- rbind(
  `SS loadings`    = round(pca.Results$Vaccounted[1, , drop = FALSE], 3),
  `Proportion Var` = round(pca.Results$Vaccounted[2, , drop = FALSE], 3),
  `Cumulative Var` = round(pca.Results$Vaccounted[3, , drop = FALSE], 3) )

colnames(summary_table) <- colnames(pca.Results$loadings)
summary_table




pca_loads  <- as.data.frame(pca.Results$loadings[, 1:numPC]) 
pca_scores <- as.data.frame(pca.Results$scores  [, 1:numPC])

RC_scores = pca_scores[, sort(names(pca_scores))]   




MinVar =  as.matrix(pca_loads) |> abs() |> apply(1, max) |> which.min() |> names()

MaxLoad.of.MinVar = pca_loads[MinVar,] |> max()
threshold = MaxLoad.of.MinVar^2     

for (i in 1:numPC) { print (pca_loads[ (pca_loads[,i])^2 >= threshold , i , drop=F ] ) }  # LINE2: compare to LINE1




weights <- c(
  alcohol     = (1.9 + 1.9) / 2,
  amphet      = (1.1 + 2.2) / 2,
  amyl        = (0.1 + 0.2) / 2,
  benzos      = (1.8 + 2.1) / 2,
  cannabis    = (1.0 + 1.7) / 2,
  choc        = (0.2 + 0.5) / 2,
  coke        = (2.3 + 2.8) / 2,
  caff        = (1.0 + 1.3) / 2,
  crack       = (3.0 + 3.0) / 2,
  ecstasy     = (0.2 + 1.1) / 2,
  heroin      = (3.0 + 3.0) / 2,
  ketamine    = (1.2 + 1.9) / 2,
  legalh      = (1.0 + 1.5) / 2,
  lsd         = (0.3 + 1.1) / 2,
  meth        = (2.8 + 2.2) / 2,
  mushrooms   = (0.2 + 0.7) / 2,
  nicotine    = (2.3 + 2.6) / 2,
  vsa         = (1.1 + 2.0) / 2,
  semer       = 0.0  # not in the table; set to 0 or NA
)


w <- weights[row.names(pca.Results$loadings)]

RC_weights <- colSums((pca.Results$loadings^2) * w) 

RC_weights <- RC_weights / sum(RC_weights)

RC_weights <- RC_weights[order(names(RC_weights))]  


weighted_dist <- dist(scale(RC_scores) %*% diag(RC_weights), method = "manhattan")


as.matrix(weighted_dist)[1:5, 1:5]  
weighted_dist[1]  
attr(weighted_dist, "Size")  



#----------------------------------------------------
# Step 11:  
# Elbow Method to identify the optimal number of clusters
# for K-medoids (PAM) with Weighted Manhattan Distance
#----------------------------------------------------


set.seed(666) 

wcd = sapply(1:10, function(k){cluster::pam(weighted_dist, k = k)$objective[1]})

plot(1:10, wcd, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters K", ylab = "Total  WCD  (Within-Cluster Dissimilarity)",
     main = "Elbow Method for Determining Optimal K")
abline(v = 4, col = "red", lty = 2)
mtext("For K-medoids (PAM) with Weighted Manhattan Distance (L1 Norm - Absolute Value Metric)",
      side = 3, line = 0.5, cex = 0.7, col = 'navy')



#----------------------------------------------------
# Silhouette Method to determine the number of clusters for 
# K-medoids (PAM) and Hierarchical, with Weighted Manhattan Distance
#----------------------------------------------------


set.seed(666)


weighted_dist_matrix <- as.matrix(weighted_dist)

pam_custom <- function(x, k) {cluster::pam(x, k = k, diss = TRUE)}


fviz_nbclust(
  x = weighted_dist_matrix,
  FUNcluster = pam_custom,
  method = "silhouette",
  k.max = 10 ) +
  labs(title = "Silhouette Method for Optimal K in PAM Clustering \n (Weighted Manhattan Distance)")



library(cluster)
set.seed(666)


hc <- hclust(weighted_dist, method = "average")  # less sensetive to outliers than "complete" linkage

silh_scores <- sapply(2:10, function(k) { clus <- cutree(hc, k)
                              mean(silhouette(clus, dist = weighted_dist)[, 3]) })

plot(2:10, silh_scores, type = "b", pch = 19,
     xlab = "Number of Clusters K", ylab = "Average Silhouette Width",
     main = "Silhouette Method for Optimal K in Hierarchical Clustering \n (Weighted Manhattan Distance)")
abline(v = 2, col = "red", lty = 2)



#------------------------------------------------------
# Final Decision: K = 3 for PAM and K=2 for Hierarchical
#------------------------------------------------------


#----------------------------------------------------
# Step 12:  'K-Medoids (PAM)' Clustering 
#----------------------------------------------------


pam_result <- pam(
  x = weighted_dist,  
  k = 3,
  diss = TRUE)        


RC_scores$pam_cluster <- factor(pam_result$clustering)


library(ggplot2)

RC_scores$pam_cluster <- factor(pam_result$clustering)


ggplot(RC_scores, aes(x = RC1, y = RC2, color = pam_cluster)) +
  geom_point(alpha = 0.7) +
  stat_ellipse(type = "norm") +
  theme_minimal() +
  labs(title = "PAM Clustering (Weighted Manhattan, k = 3) - 7 Varimax-Rotated Components") +
  theme(plot.title = element_text(size = 10)) # control title size here
  



##@@@@@@@@@@@@@@@@@@@@ VISUALIZATION @@@@@@@@@@@@@@@@@@@@@##

##      `graph2.R`

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##




#----------------------------------------------------
# Step 13:  'Hierarchical' Clustering 
#----------------------------------------------------


hcA_result <- hclust(weighted_dist, method = "average")

hcC_result <- hclust(weighted_dist, method = "complete")

hcS_result <- hclust(weighted_dist, method = "single")



plot(hcA_result, labels = FALSE, hang = -1,
     main = "Hierarchical Clustering - Average Linkage")

rect.hclust(hcA_result, k = 2, border = "red")


plot(hcC_result, labels = FALSE, hang = -1,
     main = "Hierarchical Clustering - Complete Linkage")

rect.hclust(hcC_result, k = 2, border = "red")



plot(hcS_result, labels = FALSE, hang = -1,
     main = "Hierarchical Clustering - Single Linkage")

rect.hclust(hcS_result, k = 2, border = "red")


RC_scores$hcA_cluster <- factor(cutree(hcA_result, k = 2))


RC_scores$hcC_cluster <- factor(cutree(hcC_result, k = 2))


RC_scores$hcS_cluster <- factor(cutree(hcS_result, k = 2))


fviz_cluster(list(data = RC_scores[, c(1,3)], cluster = cutree(hcA_result, k = 2)),    
             geom = "point",
             ellipse.type = "norm",
             repel = TRUE,             
             ggtheme = theme_minimal()) +
  labs(title = "Hierarchical Clustering (Average, k = 2) on PCA Scores")


fviz_cluster(list(data = RC_scores[, c(1,3)], cluster = cutree(hcC_result, k = 2)),    
             geom = "point",
             ellipse.type = "norm",
             repel = TRUE,             
             ggtheme = theme_minimal()) +
  labs(title = "Hierarchical Clustering (Complete, k = 2) on PCA Scores")


fviz_cluster(list(data = RC_scores[, c(1,3)], cluster = cutree(hcS_result, k = 2)),    
             geom = "point",
             ellipse.type = "norm",
             repel = TRUE,             
             ggtheme = theme_minimal()) +
  labs(title = "Hierarchical Clustering (Single, k = 2) on PCA Scores")




#-----------------------------------------------
# Comparing results by Average Silhouette Score :
#-----------------------------------------------


pam3 <- pam(weighted_dist, k = 3, diss = TRUE)
silh_pam3 <- mean(silhouette(pam3$clustering, dist = weighted_dist)[, 3])

hcA <- hclust(weighted_dist, method = "average")
silh_avg2 <- mean(silhouette(cutree(hcA, k = 2), dist = weighted_dist)[, 3])

hcC <- hclust(weighted_dist, method = "complete")
silh_comp2 <- mean(silhouette(cutree(hcC, k = 2), dist = weighted_dist)[, 3])

hcS <- hclust(weighted_dist, method = "single")
silh_single2 <- mean(silhouette(cutree(hcS, k = 2), dist = weighted_dist)[, 3])

avg_silh <- c(
  PAM_k3       = silh_pam3,
  Hier_avg_k2  = silh_avg2,
  Hier_comp_k2 = silh_comp2,
  Hier_sing_k2 = silh_single2
)

print(round(avg_silh, 4))



library(gt)
library(tibble)

df_silh <- tibble(
  Method = names(avg_silh),
  Silhouette = round(as.numeric(avg_silh), 2) )

gt(df_silh) %>%
  tab_header(title = "Average Silhouette Widths by Clustering Method")



pam_tbl=table(RC_scores$pam_cluster)
hcA_tbl=table(RC_scores$hcA_cluster)
hcC_tbl=table(RC_scores$hcC_cluster)
hcS_tbl=table(RC_scores$hcS_cluster)



library(dplyr)
library(tibble)
library(gt)

pam_df <- as.data.frame(pam_tbl) %>% rename(Cluster = 1, PAM = 2)
hcA_df <- as.data.frame(hcA_tbl) %>% rename(Cluster = 1, `H - Average` = 2)
hcC_df <- as.data.frame(hcC_tbl) %>% rename(Cluster = 1, `H - Complete` = 2)
hcS_df <- as.data.frame(hcS_tbl) %>% rename(Cluster = 1, `H - Single` = 2)

df_counts <- pam_df %>%
  full_join(hcA_df, by = "Cluster") %>%
  full_join(hcC_df, by = "Cluster") %>%
  full_join(hcS_df, by = "Cluster") %>%
  arrange(Cluster)

gt(df_counts) %>%
  tab_header(title = "Cluster Membership Counts by Clustering Method")





#------------------------------------------------
# Investigating the meaning or interpretation of 
# each cluster based on drug usage patterns
#------------------------------------------------


Addict_Stage$pam_cluster <- RC_scores$pam_cluster ### THE RESPONCE VARIABLE


Addict_Stage[,-1] <- lapply(Addict_Stage[,-1], function(x) as.numeric(sub("CL", "", as.character(x))))
str(Addict_Stage)
#-------------------



library(dplyr)

RC_cluster_mad <- RC_scores %>%
  group_by(pam_cluster) %>%
  summarise(across(starts_with("RC"), ~ mean(abs(. - mean(.)))))

print(RC_cluster_mad)


library(reshape2)

melted <- melt(RC_cluster_mad, id.vars = "pam_cluster")

library(ggplot2)
# sum of absolute deviations (SAD)
ggplot(melted, aes(
       x = variable, y = pam_cluster,
       fill = cut(value,
                  breaks = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1,Inf),
                  labels = c("≤.1",".1–.2",".2–.3",".3–.4",".4–.5",
                             ".5–.6",".6–.7",".7–.8",".8–.9",".9–1",">1")))) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(value, 2)), size = 3) +
  scale_fill_brewer(palette = "YlOrRd", name = "SAD") +
  labs(title = "Within-Cluster  Sum of Absolute Deviation  per RC  -  PAM Clustering",
       x = "Rotated Component (RC)",
       y = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#------------------------

cluster1 = Addict_Stage[Addict_Stage$pam_cluster==1,]
cluster2 = Addict_Stage[Addict_Stage$pam_cluster==2,]
cluster3 = Addict_Stage[Addict_Stage$pam_cluster==3,]


drug_names <- colnames(cluster1)  

for (drug in drug_names) {
  cluster1[[drug]] <- factor(cluster1[[drug]], levels = 0:6)
  cluster2[[drug]] <- factor(cluster2[[drug]], levels = 0:6)
  cluster3[[drug]] <- factor(cluster3[[drug]], levels = 0:6)
}

drug_tables <- list()

for (drug in drug_names) {
  tbl <- rbind(
    round(prop.table(table(cluster1[[drug]])) * 100, 0),
    round(prop.table(table(cluster2[[drug]])) * 100, 0),
    round(prop.table(table(cluster3[[drug]])) * 100, 0)
  )
  rownames(tbl) <- c("Cluster 1", "Cluster 2", "Cluster 3")
  colnames(tbl) <- paste0("CL", 0:6)
  drug_tables[[drug]] <- tbl
}


print(drug_tables$heroin)


#-------------------------

barplot(weights[order(weights,decreasing=T)],cex.names=0.7,las=2)
abline(h=1); abline(h=2); abline(h=3)
range(weights)

DrugsA <- c("heroin", "crack", "coke", "ecstasy", "lsd", "meth", "amphet", "mushrooms")
DrugsB <- c("cannabis", "ketamine")
DrugsC <- c("benzos", "amyl")
DrugsU <- c("alcohol", "nicotine", "caff", "choc", "legalh", "vsa") # U for Unclassified

drug_class <- ifelse(names(weights) %in% DrugsA, "A",
                     ifelse(names(weights) %in% DrugsB, "B",
                            ifelse(names(weights) %in% DrugsC, "C", "U")))


df_weights <- data.frame( drug = names(weights),
                          weight = as.numeric(weights),
                          class = drug_class )

ggplot(df_weights, aes(x = reorder(drug, -weight), y = weight, fill = class)) +
  geom_bar(stat = "identity") +
  labs(x = "Drug", y = "Weight", title = "Drug Weights by Dependence Scores") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey2") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "grey2") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "grey2")


drugsA <- c("heroin","crack","coke","ecstasy","lsd","meth","amphet","mushrooms","nicotine")
drugsB <- c("cannabis","ketamine","alcohol","caff","legalh","vsa")
drugsC <- c("benzos","amyl","choc")


#----------------------------

# Step 1: 

wA <- weights[drugsA]
wB <- weights[drugsB]
wC <- weights[drugsC]


Weighted_Addict_Stage <- Addict_Stage



# Step 2:


for (drug in names(weights)) {
  if (drug %in% colnames(Weighted_Addict_Stage)) {
    Weighted_Addict_Stage[[drug]] <- Weighted_Addict_Stage[[drug]] * weights[drug] } }

w_cluster1 <- Weighted_Addict_Stage[Weighted_Addict_Stage$pam_cluster == 1, ]
w_cluster2 <- Weighted_Addict_Stage[Weighted_Addict_Stage$pam_cluster == 2, ]
w_cluster3 <- Weighted_Addict_Stage[Weighted_Addict_Stage$pam_cluster == 3, ]


wmatrix <- data.frame(
  cluster = 1:3,
  A = c(median(unlist(w_cluster1[, drugsA])),
        median(unlist(w_cluster2[, drugsA])),
        median(unlist(w_cluster3[, drugsA]))),
  B = c(median(unlist(w_cluster1[, drugsB])),
        median(unlist(w_cluster2[, drugsB])),
        median(unlist(w_cluster3[, drugsB]))),
  C = c(median(unlist(w_cluster1[, drugsC])),
        median(unlist(w_cluster2[, drugsC])),
        median(unlist(w_cluster3[, drugsC]))) )


pivot_longer(wmatrix, cols = -cluster, names_to = "drug_class", values_to = "w_median") %>%
ggplot(aes(x = drug_class, y = w_median, group = factor(cluster), color = factor(cluster))) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Weighted Median Drug Usage Score per Drug Class and Cluster",
       x = "Drug Class", y = "Weighted Median Score", color = "Cluster") +
  scale_color_manual(values = c('green','orange','red'))


