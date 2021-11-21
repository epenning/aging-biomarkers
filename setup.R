library(haven)
library(dplyr)
library(ggplot2)
library(stringr)

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
        Health_Status = "HSD010",
        Age = "RIDAGEYR",
        Gender = "RIAGENDR"
    )

# James' Work (Plus, I added ID = "SEQN" to the code above)

# Age missing is coded as ".".  Age is topcoded at 80.  Both these type of entries are removed.
# Entries containing NA's are removed.
# Gender is reclassified from 1 and 2 to Male and Female to make the data clearer.
nhanes <- nhanes %>% filter(is.numeric(Age)) %>% filter(Age != 80) %>% filter(is.numeric(Gender)) %>% na.omit()
nhanes %>% mutate(Gender = as.factor(Gender))
nhanes$Gender <- nhanes$Gender %>% str_replace("1", "Male")
nhanes$Gender <- nhanes$Gender %>% str_replace("2", "Female")

# Only include health status 1, 2, 3 (Excellent, Very Good, Good)
write.csv(nhanes %>% filter(Health_Status < 4), "nhanes_data.csv", row.names = FALSE)

# HERE IS WHERE YOU ASSIGN GENDER.  Change below to "Female" to see the Female version.
nhanes <- nhanes %>% filter(Gender == "Male")

# PCA
pca <- nhanes %>% select(-ID, -Age, -Gender) %>% prcomp(scale=TRUE)

# Makes a scree plot for PC 1-10
pca_var <- pca$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)
x <- barplot(pca_var_per[1:10],
             main="Scree Plot", xlab="Principal Component", ylab="Percent Variation",
             ylim=c(0,25))
y <- pca_var_per[1:10]
text(x,y+1,labels=as.character(y))

# Cumulative sum of % variation accounted for by PC 1-10
# 3 PC's explain about 40% of variance, which is really good,
# but we remember we only entered 16 physiological measurements.
cumsum(pca_var_per[1:10])

#Creates Age Ranges
pca_nhanes <- nhanes %>% select(ID, Age) %>% cbind(pca$x[,1:3])
pca_nhanes <- pca_nhanes %>%
    mutate(Age_Range = case_when(
        Age <= 20 ~ "20 and under",
        Age > 20 & Age <= 40 ~ "21-40",
        Age > 40 & Age <= 60 ~ "41-60",
        Age > 60 & Age <= 80 ~ "61-79"
    )
    )

#### Clustering analysis


pca_nhanes

# Plots PC1 vs. PC2 and PC1 vs PC3.  Ellipse of each Age Range
ggplot(pca_nhanes, aes(PC1, PC2, fill = Age_Range, col = Age_Range)) +
    stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
    geom_point() +
    #geom_text(aes(label=ID, size=10)) +
    xlab(paste("PC1 - ", pca_var_per[1], "%", sep="")) +
    ylab(paste("PC2 - ", pca_var_per[2], "%", sep="")) +
    ggtitle("PC1 vs PC2 Graph")

ggplot(pca_nhanes, aes(PC1, PC3, fill = Age_Range, col = Age_Range)) +
    stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
    geom_point() +
    #geom_text(aes(label=ID, size=10)) +
    xlab(paste("PC1 - ", pca_var_per[1], "%", sep="")) +
    ylab(paste("PC3 - ", pca_var_per[3], "%", sep="")) +
    ggtitle("PC1 vs PC3 Graph")

# Ranks the contributors to PCA 1, with the most significant at the top.
loading_score_ranked <- names(sort(abs(pca$rotation[,1]), decreasing=TRUE))
loading_score_ranked

# Outputs correlation between top 10 physiological measurements and PCA 1-3.
pca_top_10_measurements <- nhanes %>% select(loading_score_ranked[1:10]) %>% cbind(pca$x[,1:3])
cor(pca_top_10_measurements[1:10], pca_top_10_measurements[, 11:13])

# Janice: Repeat PCA but with female
# Age missing is coded as ".".  Age is topcoded at 80.  Both these type of entries are removed.
# Entries containing NA's are removed.
# Gender is reclassified from 1 and 2 to Male and Female to make the data clearer.
nhanes <- nhanes %>% filter(is.numeric(Age)) %>% filter(Age != 80) %>% filter(is.numeric(Gender)) %>% na.omit()
nhanes %>% mutate(Gender = as.factor(Gender))
nhanes$Gender <- nhanes$Gender %>% str_replace("1", "Male")
nhanes$Gender <- nhanes$Gender %>% str_replace("2", "Female")

write.csv(nhanes, "nhanes_data.csv", row.names = FALSE)

# HERE IS WHERE YOU ASSIGN GENDER.  Change below to "Female" to see the Female version.
nhanes <- nhanes %>% filter(Gender == "Female")

# PCA
pca <- nhanes %>% select(-ID, -Age, -Gender) %>% prcomp(scale=TRUE)

# Makes a scree plot for PC 1-10
pca_var <- pca$sdev^2
pca_var_per <- round(pca_var/sum(pca_var)*100, 1)
x <- barplot(pca_var_per[1:10],
             main="Scree Plot", xlab="Principal Component", ylab="Percent Variation",
             ylim=c(0,25))
y <- pca_var_per[1:10]
text(x,y+1,labels=as.character(y))

# Cumulative sum of % variation accounted for by PC 1-10
# 3 PC's explain about 40% of variance, which is really good,
# but we remember we only entered 16 physiological measurements.
cumsum(pca_var_per[1:10])

#Creates Age Ranges
pca_nhanes <- nhanes %>% select(ID, Age) %>% cbind(pca$x[,1:3])
pca_nhanes <- pca_nhanes %>%
    mutate(Age_Range = case_when(
        Age <= 20 ~ "20 and under",
        Age > 20 & Age <= 40 ~ "21-40",
        Age > 40 & Age <= 60 ~ "41-60",
        Age > 60 & Age <= 80 ~ "61-79"
    )
    )

#### Clustering analysis


pca_nhanes

# Plots PC1 vs. PC2 and PC1 vs PC3.  Ellipse of each Age Range
ggplot(pca_nhanes, aes(PC1, PC2, fill = Age_Range, col = Age_Range)) +
    stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
    geom_point() +
    #geom_text(aes(label=ID, size=10)) +
    xlab(paste("PC1 - ", pca_var_per[1], "%", sep="")) +
    ylab(paste("PC2 - ", pca_var_per[2], "%", sep="")) +
    ggtitle("PC1 vs PC2 Graph")

ggplot(pca_nhanes, aes(PC1, PC3, fill = Age_Range, col = Age_Range)) +
    stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
    geom_point() +
    #geom_text(aes(label=ID, size=10)) +
    xlab(paste("PC1 - ", pca_var_per[1], "%", sep="")) +
    ylab(paste("PC3 - ", pca_var_per[3], "%", sep="")) +
    ggtitle("PC1 vs PC3 Graph")

# Ranks the contributors to PCA 1, with the most significant at the top.
loading_score_ranked <- names(sort(abs(pca$rotation[,1]), decreasing=TRUE))
loading_score_ranked

# Outputs correlation between top 10 physiological measurements and PCA 1-3.
pca_top_10_measurements <- nhanes %>% select(loading_score_ranked[1:10]) %>% cbind(pca$x[,1:3])
cor(pca_top_10_measurements[1:10], pca_top_10_measurements[, 11:13])




# Anthony work
dim(nhanes)
head(nhanes)


View(nhanes)
table(nhanes$RIAGENDR,nhanes$SDDSRVYR)
class(nhanes$RIDAGEYR)
nhanes %>% ggplot(aes(x="RIDAGEYR")) + geom_histogram(stat = "count")
range(nhanes$RIDAGEYR)
nhanes %>% ggplot(aes(x=Age, y=Systolic_Blood_Pressure)) + geom_point()
summary(lm(nhanes$Age~nhanes$Systolic_Blood_Pressure))

#11.21
library(tidyverse)
library(cluster)

nhanes <- nhanes %>% na.omit()

cormat <- nhanes %>% cor()

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

# hypotheses: older adults have greater BP and greater blood urea nitrogen


## Janice's work
# Create male and female subsets
nhanes_male <- nhanes %>% filter(Gender == "Male")
nhanes_female <- nhanes %>% filter(Gender == "Female")

# Albumin - slightly negative for male and female
nhanes_male %>% ggplot(aes(x=Age, y=Albumin)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Albumin)) + geom_point() + geom_smooth(method="lm")

# Alkaline Phosphatase - negative for male and tiny bit negative for female
nhanes_male %>% ggplot(aes(x=Age, y=Alkaline_Phosphatase)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Alkaline_Phosphatase)) + geom_point() + geom_smooth(method="lm")

# Blood Urea Nitrogen - slightly positive for both
nhanes_male %>% ggplot(aes(x=Age, y=Blood_Urea_Nitrogen)) + geom_point() + geom_smooth(method="lm")
nhanes_female %>% ggplot(aes(x=Age, y=Blood_Urea_Nitrogen)) + geom_point() + geom_smooth(method="lm")

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

