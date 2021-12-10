library(haven)
library(dplyr)
library(ggplot2)
library(stringr)
library(cluster)
library(tidyverse)

##### Import and create dataset nhanes

# Read data
demographics <- read_xpt("aging-biomarkers/DEMO_J.XPT")
biochemistry <- read_xpt("aging-biomarkers/BIOPRO_J.XPT")
cReactive <- read_xpt("aging-biomarkers/HSCRP_J.XPT")
glycohemoglobin <- read_xpt("aging-biomarkers/GHB_J.XPT")
bloodCount <- read_xpt("aging-biomarkers/CBC_J.XPT")
bloodPressure <- read_xpt("aging-biomarkers/BPX_J.XPT")

# Merge data
nhanes <- merge(demographics, biochemistry, by = "SEQN")
nhanes <- merge(nhanes, cReactive, by = "SEQN")
nhanes <- merge(nhanes, glycohemoglobin, by = "SEQN")
nhanes <- merge(nhanes, bloodCount, by = "SEQN")
nhanes <- merge(nhanes, bloodPressure, by = "SEQN")

# Get mean systolic blood pressure
nhanes$BPXSY_MEAN <-
  rowMeans(nhanes[, c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")], na.rm = T)

# Extract and rename variables of interest
nhanes <-
  nhanes %>% select(
    ID = "SEQN",
    Albumin = "LBXSAL",
    Alkaline_Phosphatase = "LBXSAPSI",
    Blood_Urea_Nitrogen = "LBXSBU",
    Creatinine = "LBXSCR",
    C_Reactive_Protein = "LBXHSCRP",
    Glucose = "LBXSGL",
    Glycohemoglobin = "LBXGH",
    Total_Cholesterol = "LBXSCH",
    Uric_Acid = "LBXSUA",
    White_Blood_Cell_Count = "LBXWBCSI",
    Lymphocyte_Percent = "LBXLYPCT",
    Mean_Cell_Volume = "LBXMCVSI",
    Red_Cell_Distribution_Width = "LBXRDW",
    Systolic_Blood_Pressure = "BPXSY_MEAN",
    Age = "RIDAGEYR",
    Gender = "RIAGENDR"
  )


##### Clean data

# Number of individuals before data cleaned
nrow(nhanes)

# Age missing is coded as ".".  Age is topcoded at 80.  Both these type of entries are removed.
# Entries containing NA are removed.
nhanes <- nhanes %>%
  filter(is.numeric(Age)) %>%
  na.omit() %>%
  filter(Age != 80)

# Gender is reclassified from 1 and 2 to Male and Female to make the data clearer.
nhanes$Gender <- nhanes$Gender %>% str_replace("1", "Male")
nhanes$Gender <- nhanes$Gender %>% str_replace("2", "Female")

# Number of individuals after data cleaned
nrow(nhanes)


##### Data distribution

# Age range of nhanes
range(nhanes$Age)

# How is Age distributed in nhanes?
nhanes %>% ggplot(aes(x = Age)) + geom_histogram(stat = "count")


##### Regression for each biomarker with age.

min_age = min(nhanes$Age)

# Function for scatterplot of variable vs. age
do_plot <- function(biomarker) {
  lower.cut = quantile(nhanes[, biomarker], 0.02)
  upper.cut = quantile(nhanes[, biomarker], 0.98)
  print(nhanes %>% ggplot(aes(x = Age, y = nhanes[, biomarker])) +
          geom_point(position = "jitter",
                     size = 0.7,
                     alpha = 0.5) +
          geom_smooth(method = "lm") + ylab(biomarker) +
          scale_x_continuous(breaks = seq(from = min_age, to = 79, by = 10)) +
          coord_cartesian(ylim = c(lower.cut * .9, upper.cut * 1.1)))
  print(paste(biomarker, "Correlation with Age:"))
  print(summary(lm(nhanes$Age ~ nhanes[, biomarker])))
}

# Generate scatterplots for all continuous variables
# Age and Gender are last 2 variables, don't plot those
# ID is the first variable, don't plot it
n <- ncol(nhanes) - 2
for (variable in colnames(nhanes[2:n])) {
  do_plot(variable)
}


##### 2nd cleaning and data prep

# This variable exhibits different behavior before and after age 20 so it's removed
nhanes_original <- nhanes
nhanes <- nhanes %>% select(-Alkaline_Phosphatase)

# Output cleaned nhanes data for use in Python
write.csv(nhanes, 'nhanes_data.csv')

# Make scaled version of nhanes
nhanes_scaled <- nhanes %>% mutate_at(vars(-Gender, -ID), scale)


##### Correlation Matrix

cormat <- nhanes %>% select(-Gender, -ID) %>% cor()

tidycor <- cormat %>%
  as.data.frame %>%
  rownames_to_column("var1") %>%
  pivot_longer(-1, names_to = "var2", values_to = "correlation")

tidycor %>% ggplot(aes(var1, var2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  geom_text(aes(label = round(correlation, 2)), color = "black", size = 2.5) + #overlays correlation values
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + #flips the x-axis labels
  coord_fixed()


##### Clustering: is clustering the best way to analyze the data?

# Do Blood Pressure and Blood Urea Nitrogen show obvious clusters on Age?
nhanes %>% ggplot(aes(Systolic_Blood_Pressure, Blood_Urea_Nitrogen, col = Age)) +
  geom_point()

# Do Blood Pressure and Blood Urea Nitrogen form true clusters at all?
clust_dat <- nhanes %>% select(Systolic_Blood_Pressure, Blood_Urea_Nitrogen)

sil_width <- vector() #empty vector to hold mean sil width

for (i in 2:10) {
  kms <- kmeans(clust_dat, centers = i) #compute k-means solution
  sil <- silhouette(kms$cluster, dist(clust_dat)) #get sil widths
  sil_width[i] <- mean(sil[, 3]) #take averages (higher is better)
}

# Silhouette plot: gives 2 as number of clusters
ggplot() +
  geom_line(aes(x = 1:10, y = sil_width)) +
  scale_x_continuous(name = "k", breaks = 1:10)

# Use result from silhouette plot.  How does 6 clusters look?
kmeans <- clust_dat %>% kmeans(2)
kmeansclust <-
  clust_dat %>% mutate(cluster = as.factor(kmeans$cluster))
kmeansclust %>% ggplot(aes(Age, Systolic_Blood_Pressure, color = cluster)) +
  geom_point()

##### PCA

# Systolic Blood Pressure is just one variable.  What if we look at all of them using PCA?
pca <- nhanes_scaled %>%
  select(-ID, -Age, -Gender) %>%
  prcomp(scale = TRUE)

# Loading score: how much original variable contributes to PC1
loading_scores <- pca$rotation[, 1]
loading_scores

# Ranks the contributors to PC 1, with the most significant at the top.
loading_score_ranked <-
  names(sort(abs(pca$rotation[, 1]), decreasing = TRUE))
loading_score_ranked

# Makes a scree plot for PC 1-10.  How much variance in the data does each PC explain
pca_var <- pca$sdev ^ 2
pca_var_per <- round(pca_var / sum(pca_var) * 100, 1)
x <- barplot(
  pca_var_per[1:10],
  main = "Scree Plot",
  xlab = "Principal Component",
  ylab = "Percent Variation",
  ylim = c(0, 25)
)
y <- pca_var_per[1:10]
text(x, y + 1, labels = as.character(y))

# Binds PCA data to original variables Age and Gender
pca_nhanes <- nhanes %>%
  select(ID, Age, Gender) %>%
  cbind(pca$x[, 1:3])

# Plots PC1 vs. PC2
# Does Age form obvious clusters in PC space?
ggplot(pca_nhanes, aes(PC1, PC2, col = Age)) +
  geom_point() +
  xlab(paste("PC1 - ", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca_var_per[2], "%", sep = "")) +
  ggtitle("PC1 vs PC2 Graph")


##### PCA Clustering

# Does the data form clusters in PC space?
clust_dat <- pca_nhanes %>% select(PC1, PC2)

sil_width <- vector() #empty vector to hold mean sil width

for (i in 2:10) {
  kms <- kmeans(clust_dat, centers = i) #compute k-means solution
  sil <- silhouette(kms$cluster, dist(clust_dat)) #get sil widths
  sil_width[i] <- mean(sil[, 3]) #take averages (higher is better)
}

ggplot() +
  geom_line(aes(x = 1:10, y = sil_width)) +
  scale_x_continuous(name = "k", breaks = 1:10)

kmeans <- clust_dat %>% kmeans(2)
kmeansclust <- clust_dat %>%
  mutate(cluster = as.factor(kmeans$cluster), Age = nhanes$Age)
kmeansclust %>% ggplot(aes(PC1, PC2, col = cluster)) +
  geom_point()

# Result: biomarkers don't cluster on age, use regression-based model
# Go to regression.ipynb, then randomforest.ipynb for models


##### Future work

# Experimented with removing outliers from dataset
# Not part of current model but may include in future models

# Rate of change for several biomarker levels different b/f and after 22 
nhanes_original_22plus <- nhanes_original %>% filter(Age >= 22)

# Alkaline Phosphatase: use dataset with young removed, remove outliers
nhanes_original_22plus %>%
  filter(Alkaline_Phosphatase < 200) %>%
  ggplot(aes(x = Age, y = Alkaline_Phosphatase)) +
  geom_point(position = "jitter", size = 0.7, alpha = 0.5) +
  geom_smooth(method = "lm") +
  ylab("Alkaline Phosphatase")

# Blood Urea Nitrogen: remove outliers
nhanes %>%
  filter(Blood_Urea_Nitrogen < 40) %>%
  ggplot(aes(x = Age, y =  Blood_Urea_Nitrogen)) +
  geom_point(position = "jitter", size = 0.7, alpha = 0.5) +
  geom_smooth(method = "lm") +
  ylab("Blood Urea Nitrogen")

# Creatinine: removed outliers
nhanes %>%
  filter(Creatinine < 3) %>%
  ggplot(aes(x = Age, y = Creatinine)) +
  geom_point(position = "jitter", size = 0.7, alpha = 0.5) +
  geom_smooth(method = "lm") +
  ylab("Creatinine")

# C Reactive Protein: removed outliers
nhanes %>%
  filter(C_Reactive_Protein < 20) %>%
  ggplot(aes(x = Age, y = C_Reactive_Protein)) +
  geom_point(position = "jitter", size = 0.7, alpha = 0.5) +
  geom_smooth(method = "lm") +
  ylab("C Reactive Protein")

# Glucose: removed outliers
nhanes %>%
  filter(Glucose < 200) %>%
  ggplot(aes(x = Age, y = Glucose)) +
  geom_point(position = "jitter", size = 0.7, alpha = 0.5) +
  geom_smooth(method = "lm") +
  ylab("Glucose")

# White Blood Cell Count: removed outliers
nhanes %>%
  filter(White_Blood_Cell_Count < 50) %>%
  ggplot(aes(x = Age, y = White_Blood_Cell_Count)) +
  geom_point(position = "jitter", size = 0.7, alpha = 0.5) +
  geom_smooth(method = "lm") +
  ylab("White Blood Cell Count")

# Lymphocyte Percent: use dataset with young removed
nhanes_original_22plus %>%
  ggplot(aes(x = Age, y = Lymphocyte_Percent)) +
  geom_point(position = "jitter", size = 0.7, alpha = 0.5) +
  geom_smooth(method = "lm") +
  ylab("Lymphocyte Percent")