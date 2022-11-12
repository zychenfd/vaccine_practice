#==load package==
library(readr)
library(stringr)
library(dplyr)
library(reshape2)

#==read data==
gisaid <- read_delim("../data/gisaid_hcov-19.tsv", delim="\t")

#==clean data==
gisaid1 <- gisaid %>% filter(Host == "Human")
gisaid1$Lineage[gisaid1$Lineage %in% c("Q.1","Q.2","Q.3","Q.4","Q.5","Q.6","Q.7","Q.8")] <- "B.1.1.7"
gisaid1$Lineage[gisaid1$Lineage %in% c("B.1.351.1","B.1.351.2","B.1.351.3","B.1.351.4")] <- "B.1.351"
gisaid1$Lineage[str_detect(gisaid1$Lineage, "AY.")] <- "B.1.617.2"
gisaid1$Lineage[!gisaid1$Lineage %in% c("B.1.1.7","B.1.351","P.1","B.1.617.2")] <- "Non-VOC"
gisaid1$`Patient age`[gisaid1$`Patient age`== "12 days"] <- 12/365
gisaid1$`Patient age`[gisaid1$`Patient age`== "1 month"] <- 1/12
gisaid1$`Patient age`[gisaid1$`Patient age`== "4 months"] <- 4/12
gisaid1$`Patient age`[gisaid1$`Patient age`== "7 months"] <- 7/12
gisaid1$`Patient age`[gisaid1$`Patient age`== "9 months"] <- 9/12
gisaid1$`Patient age`[gisaid1$`Patient age`== "10 months"] <- 10/12
gisaid1$`Patient age` <- round(as.numeric(gisaid1$`Patient age`),2)

#==analysis data==
sex <- gisaid1 %>% group_by(Lineage, Gender) %>% summarise(case = n()) %>% filter(Gender != "unknown")
sex1 <- gisaid1  %>% filter(Gender != "unknown") %>% group_by(Lineage) %>% summarise(total = n())
sex2 <- left_join(sex, sex1)
sex2$Sex_prop <- paste(format(round(sex2$case * 100 /sex2$total,1), nsmall = 1),
                       "% (",sex2$case, "/", sex2$total,")",sep = "")
cut(gisaid1$`Patient age`, breaks = c(-Inf,17,65,Inf), labels = c("0-17","18-65","65+")) -> gisaid1$age_group
age <- gisaid1 %>% group_by(Lineage, age_group) %>% summarise(case = n()) %>% filter(!(is.na(age_group)))
age1 <- gisaid1  %>% filter(!(is.na(age_group))) %>% group_by(Lineage) %>% summarise(total = n())
age2 <- left_join(age, age1)
age2$Age_prop <- paste(format(round(age2$case * 100 /age2$total,1), nsmall = 1),
                       "% (",age2$case, "/", age2$total,")",sep = "")
sex3 <- dcast(sex2, Gender ~ Lineage, value.var = "Sex_prop")
age3 <- dcast(age2, age_group ~ Lineage, value.var = "Age_prop")
names(sex3)[1] <- "Vari"
names(age3)[1] <- "Vari"
total <- rbind(sex3, age3)

#==output==
write.csv(total,"../output/table_HK.csv", row.names = F)