library(haven)
library(dplyr)
library(ggplot2)
library(stringr)
library(cluster)
library(tidyverse)

##### Import and create data set nhanes

# Read data
demographics <- read_xpt("DEMO_J.XPT")
biochemistry <- read_xpt("BIOPRO_J.XPT")
cReactive <- read_xpt("HSCRP_J.XPT")
glycohemoglobin <- read_xpt("GHB_J.XPT")
bloodCount <- read_xpt("CBC_J.XPT")
bloodPressure <- read_xpt("BPX_J.XPT")
healthStatus <- read_xpt("HSQ_J.XPT")

# Merge data
nhanes <- merge(demographics, biochemistry, by = "SEQN")
nhanes <- merge(nhanes, cReactive, by = "SEQN")
nhanes <- merge(nhanes, glycohemoglobin, by = "SEQN")
nhanes <- merge(nhanes, bloodCount, by = "SEQN")
nhanes <- merge(nhanes, bloodPressure, by = "SEQN")
nhanes <- merge(nhanes, healthStatus, by = "SEQN")

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
#        Health_Status = "HSD010",
        Age = "RIDAGEYR",
        Gender = "RIAGENDR",
#        Income_Ratio = "INDFMPIR"   # Ratio of Income to Poverty, thought about adding but didn't
    )

##### Clean data

# Age missing is coded as ".".  Age is topcoded at 80.  Both these type of entries are removed.
# Entries containing NA's are removed.
nhanes <- nhanes %>% filter(is.numeric(Age)) %>% na.omit() %>% filter(Age != 80)

# Gender is reclassified from 1 and 2 to Male and Female to make the data clearer.
nhanes %>% mutate(Gender = as.factor(Gender))
nhanes$Gender <- nhanes$Gender %>% str_replace("1", "Male")
nhanes$Gender <- nhanes$Gender %>% str_replace("2", "Female")

#Exlude under 25 because ????????
nhanes_25plus <- nhanes %>% filter(Age > 25)

# Make scaled version of nhanes.  Keep raw data in case need to present.
nhanes_scaled <- nhanes %>% mutate_at(vars(-Gender, -ID), scale)

# Only include health status 1, 2, 3 (Excellent, Very Good, Good)
#write.csv(nhanes %>% filter(Health_Status < 4), "nhanes_data.csv", row.names = FALSE)

##### Correlation Matrix

cormat <- nhanes_scaled %>% select(-Gender, -ID) %>% cor()
cormat %>% as.data.frame %>% rownames_to_column("var1")

tidycor <- cormat %>%
    as.data.frame %>%
    rownames_to_column("var1") %>%
    pivot_longer(-1, names_to = "var2", values_to = "correlation")

tidycor %>% ggplot(aes(var1, var2, fill=correlation)) +
    geom_tile() +
    scale_fill_gradient2(low="red", mid="white", high = "blue") +
    geom_text(aes(label=round(correlation,2)),color = "black", size = 3.5)+ #overlays correlation values
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + #flips the x-axis labels
    coord_fixed()

##### Clustering: is clustering the best way to analyze the data?

# Does Blood Pressure cluster on Gender
nhanes_scaled %>% ggplot(aes(Age, Systolic_Blood_Pressure, col=Gender)) + geom_point()

# Does the data form clusters at all?
clust_dat<- nhanes_scaled %>% select(Age, Systolic_Blood_Pressure)

sil_width<-vector() #empty vector to hold mean sil width

for(i in 2:10){
    kms <- kmeans(clust_dat,centers=i) #compute k-means solution
    sil <- silhouette(kms$cluster, dist(clust_dat)) #get sil widths
    sil_width[i]<-mean(sil[,3]) #take averages (higher is better)
}

# Silhouette plot: gives 2 as number of clusters
ggplot() + 
geom_line(aes(x=1:10,y=sil_width)) +
scale_x_continuous(name="k",breaks=1:10)

# Use result from silhouette plot.  How does 2 clusters look?
kmeans1 <- clust_dat %>% kmeans(2)
kmeansclust <- clust_dat %>% mutate(cluster=as.factor(kmeans1$cluster))
kmeansclust %>% ggplot(aes(Age, Systolic_Blood_Pressure, color=cluster)) + geom_point()

##### PCA

# Systolic Blood Pressure is just one variable.  What if we look at all of them using PCA?
pca <- nhanes_scaled %>% select(-ID, -Age, -Gender) %>% prcomp(scale=TRUE) 

# Loading score: how much original variable contributes to PC1
loading_scores <- pca$rotation[,1]
loading_scores

# Ranks the contributors to PC 1, with the most significant at the top.
loading_score_ranked <- names(sort(abs(pca$rotation[,1]), decreasing=TRUE))
loading_score_ranked

# Makes a scree plot for PC 1-10.  How much variance in the data does each PC explain
pca_var <- pca$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)
x <- barplot(pca_var_per[1:10], 
             main="Scree Plot", xlab="Principal Component", ylab="Percent Variation",
             ylim=c(0,25))
y <- pca_var_per[1:10]
text(x,y+1,labels=as.character(y))

# Binds PCA data to original variables Age and Gender
pca_nhanes <- nhanes %>% select(ID, Age, Gender) %>% cbind(pca$x[,1:3])

# Plots PC1 vs. PC2
# Does Age form obvious clusters in PC space?
ggplot(pca_nhanes, aes(PC1, PC2, col=Age)) +
    geom_point() +
    xlab(paste("PC1 - ", pca_var_per[1], "%", sep="")) +
    ylab(paste("PC2 - ", pca_var_per[2], "%", sep="")) +
    ggtitle("PC1 vs PC2 Graph")

##### PCA Clustering

# Does the data form clusters in PC space?
clust_dat<- pca_nhanes %>% select(PC1,PC2)

sil_width<-vector() #empty vector to hold mean sil width

for(i in 2:10){
    kms <- kmeans(clust_dat,centers=i) #compute k-means solution
    sil <- silhouette(kms$cluster, dist(clust_dat)) #get sil widths
    sil_width[i]<-mean(sil[,3]) #take averages (higher is better)
}

ggplot() + 
    geom_line(aes(x=1:10,y=sil_width)) +
    scale_x_continuous(name="k",breaks=1:10)

kmeans1 <- clust_dat %>% kmeans(2)
kmeansclust <- clust_dat %>% mutate(cluster=as.factor(kmeans1$cluster), Age=nhanes$Age)
kmeansclust %>% ggplot(aes(PC1, PC2, col=cluster)) + geom_point()

##### Anthony

# How is Age distributed in nhanes?
nhanes %>% ggplot(aes(x=Age)) + geom_histogram(stat = "count")

# Age range of nhanes
range(nhanes$Age)

nhanes %>% ggplot(aes(x=Age, y=Systolic_Blood_Pressure)) + geom_point()
summary(lm(nhanes$Age~nhanes$Systolic_Blood_Pressure))

# hypotheses: older adults have greater BP and greater blood urea nitrogen

##### Regression for each variable

## Janice's work
# Create male and female subsets

# Function to do the correlations for a particular gender
do_correlations <- function(gender) {
    nhanes_male <- nhanes %>% filter(Gender == gender)
    
    # Albumin - slightly negative for male and female
    nhanes_male %>% ggplot(aes(x = Age, y = Albumin)) + geom_point() + geom_smooth(method =
                                                                                       "lm")
    
    ##... insert other correlations here...
}

do_correlations("Male")
do_correlations("Female")

#nhanes_male <- nhanes %>% filter(Gender == "Male")
#nhanes_female <- nhanes %>% filter(Gender == "Female")

do_plot <- function(biomarker) {
    nhanes %>% ggplot(aes(x=Age, y=nhanes[,biomarker])) + geom_point(position="jitter",size=0.7,alpha=0.5) + geom_smooth(method="lm") + ylab(biomarker)
    
}

do_plot("Albumin")


# Albumin - slightly negative for male and female
nhanes %>% ggplot(aes(x=Age, y=Albumin)) + geom_point() + geom_smooth(method="lm")
nhanes_male %>% ggplot(aes(x=Age, y=Albumin)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Albumin)) + geom_point() + geom_smooth(method="lm")

# Alkaline Phosphatase - negative for male and tiny bit negative for female
nhanes_male %>% ggplot(aes(x=Age, y=Alkaline_Phosphatase)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Alkaline_Phosphatase)) + geom_point() + geom_smooth(method="lm")

# Blood Urea Nitrogen - slightly positive for both
nhanes_male %>% ggplot(aes(x=Age, y=Blood_Urea_Nitrogen)) + geom_point() + geom_smooth(method="lm") + ggtitle("Blood Urea Nitrogen vs. Age - Male")
nhanes_female %>% ggplot(aes(x=Age, y=Blood_Urea_Nitrogen)) + geom_point() + geom_smooth(method="lm") + ggtitle("Blood Urea Nitrogen vs. Age - Female")

# Creatinine - slightly positive for male and female
nhanes_male %>% ggplot(aes(x=Age, y=Creatinine)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Creatinine)) + geom_point() + geom_smooth(method="lm")
nhanes_male %>%  filter(Creatinine < 5) %>% ggplot(aes(x=Age, y=Creatinine)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% filter(Creatinine < 3) %>% ggplot(aes(x=Age, y=Creatinine)) + geom_point() + geom_smooth(method="lm")

# C Reactive Protein - not really a correlation
nhanes_male %>% ggplot(aes(x=Age, y=C_Reactive_Protein)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=C_Reactive_Protein)) + geom_point() + geom_smooth(method="lm")

nhanes_male %>% filter(C_Reactive_Protein < 50) %>% ggplot(aes(x=Age, y=C_Reactive_Protein)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% filter(Creatinine < 50) %>% ggplot(aes(x=Age, y=C_Reactive_Protein)) + geom_point() + geom_smooth(method="lm")

# Glucose - slightly positive correlation for both
nhanes_male %>% ggplot(aes(x=Age, y=Glucose)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Glucose)) + geom_point() + geom_smooth(method="lm")

# Glycohemoglobin - positive for both
nhanes_male %>% ggplot(aes(x=Age, y=Glycohemoglobin)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Glycohemoglobin)) + geom_point() + geom_smooth(method="lm")

# Total Cholesterol - slightly positive for male and female
nhanes_male %>% ggplot(aes(x=Age, y=Total_Cholesterol)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Total_Cholesterol)) + geom_point() + geom_smooth(method="lm")

# Uric Acid - no correlation for male, positive for female
nhanes_male %>% ggplot(aes(x=Age, y=Uric_Acid)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Uric_Acid)) + geom_point() + geom_smooth(method="lm")

# White Blood Cell Count Visualization - no correlation or tiny bit negative for both
nhanes %>% arrange(desc(White_Blood_Cell_Count)) # to see highest/outliers for WBC
nhanes %>% ggplot(aes(x=Age)) + geom_point(aes(y=White_Blood_Cell_Count))
nhanes %>% filter(White_Blood_Cell_Count < 50) %>% ggplot(aes(x=Age)) + geom_point(aes(y=White_Blood_Cell_Count))
nhanes %>% filter(White_Blood_Cell_Count < 30) %>% ggplot(aes(x=Age, y=White_Blood_Cell_Count)) + geom_point() + geom_smooth(method="lm")
nhanes_male %>% filter(White_Blood_Cell_Count < 50) %>% ggplot(aes(x=Age, y=White_Blood_Cell_Count)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% filter(White_Blood_Cell_Count < 30) %>% ggplot(aes(x=Age, y=White_Blood_Cell_Count)) + geom_point() + geom_smooth(method="lm")

# Lymphocyte Percent - slightly negative for male and female
nhanes_male %>% ggplot(aes(x=Age, y=Lymphocyte_Percent)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Lymphocyte_Percent)) + geom_point() + geom_smooth(method="lm")

# Mean Cell Volume - positive for both
nhanes_male %>% ggplot(aes(x=Age, y=Mean_Cell_Volume)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Mean_Cell_Volume)) + geom_point() + geom_smooth(method="lm")

# Red Cell Distribution Width - slightly positive for male, no correlation for female
nhanes_male %>% ggplot(aes(x=Age, y=Red_Cell_Distribution_Width)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Red_Cell_Distribution_Width)) + geom_point() + geom_smooth(method="lm")

# Systolic Blood Pressure - positive for both
nhanes_male %>% ggplot(aes(x=Age, y=Systolic_Blood_Pressure)) + geom_point() + geom_smooth(method="lm") + ggtitle("Systolic Blood Pressure vs. Age - Male")
nhanes_female %>% ggplot(aes(x=Age, y=Systolic_Blood_Pressure)) + geom_point() + geom_smooth(method="lm") + ggtitle("Systolic Blood Pressure vs. Age - Female")
