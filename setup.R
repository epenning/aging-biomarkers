library(haven)
library(dplyr)
library(ggplot2)

# Read data
demographics <- read_xpt("DEMO_J.XPT")
biochemistry <- read_xpt("BIOPRO_J.XPT")
cReactive <- read_xpt("HSCRP_J.XPT")
glycohemoglobin <- read_xpt("GHB_J.XPT")
bloodCount <- read_xpt("CBC_J.XPT")
bloodPressure <- read_xpt("BPX_J.XPT")

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
        Age = "RIDAGEYR"
    )

# Age is topcoded at 80 so 80 and above should be excluded
# Age missing is coded as "." which should be removed

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
