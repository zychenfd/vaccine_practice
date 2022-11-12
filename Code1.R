# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("stringr")
# install.packages("ggsci")
# install.packages("ggpubr")
# install.packages("sf")
# install.packages("rgdal")
# install.packages("readr")
# install.packages("reshape2")
# install.packages("lubridate")

#==load packages==
library(dplyr)
library(ggplot2)
library(stringr)
library(ggsci)
library(ggpubr)
library(lubridate)

#==read data==
df <- read.csv("../data/sequences.csv") # read genome data

#====================
#==data clean==
#====================
#Step1. Remove sequences from non-human
df <- df %>% filter(Host == "Homo sapiens") 

#Step2. Remove sequences with incomplete date of collection (only year) or sampled before 1 December 2019
df$Collection_Date <- as.Date(df$Collection_Date)
df <- df %>% filter(!is.na(Collection_Date))  %>% filter(Collection_Date >= as.Date("2019-12-01"))

#Step3. Remove sequences with country information unavailable
df <- df %>% filter(Country %in% c("Kenya","Netherlands", "New Zealand", "Mexico", "India"))

#Step4. Variant Lineage
df <- df %>% filter(!Pangolin == "") %>% filter(!Pangolin == "unclassifiable")
df$Pangolin[df$Pangolin %in% c("Q.1","Q.2","Q.3","Q.4","Q.5","Q.6","Q.7","Q.8")] <- "B.1.1.7"
df$Pangolin[df$Pangolin %in% c("B.1.351.1","B.1.351.2","B.1.351.3","B.1.351.4")] <- "B.1.351"
df$Pangolin[df$Pangolin %in% c("P.1.1","P.1.2","P.1.3","P.1.4","P.1.5","P.1.6","P.1.7","P.1.8",
                                "P.1.9","P.1.10","P.1.10.2","P.1.11","P.1.12" )] <- "P.1"
df$Pangolin[str_detect(df$Pangolin, "AY.")] <- "B.1.617.2"
df$Pangolin[df$Pangolin %in% c("B","B.1","A.1","A")] <- "Ref"
df$Pangolin[!df$Pangolin %in% c("B.1.1.7","B.1.351","P.1","B.1.617.2","Ref")] <- "Others"
df1 <- df[,c("Accession","Pangolin","Country","Collection_Date","Release_Date")]
df1 <- df1 %>% filter(!(Pangolin == "B.1.1.7" & Collection_Date < as.Date("2020-09-20")))
df1 <- df1 %>% filter(!(Pangolin == "B.1.351" & Collection_Date < as.Date("2020-05-11")))
df1 <- df1 %>% filter(!(Pangolin == "P.1" & Collection_Date < as.Date("2020-11-03")))
df1 <- df1 %>% filter(!(Pangolin == "B.1.617.2" & Collection_Date < as.Date("2020-10-23")))

country <- "India" # choose which country to plot

#==============
#==data  plot==
#==============
#1. By week
df_week2 <- df1 %>% filter(Country == country) %>% group_by(Country, Collection_Date, Pangolin) %>% summarise(new = n())
df_week2 <- df_week2 %>% filter(Collection_Date <= as.Date("2021-09-26")) %>% filter(Collection_Date >= as.Date("2021-01-04"))
df_week2$week <- isoweek(df_week2$Collection_Date)
df_week2$year <- isoyear(df_week2$Collection_Date)
df_week2$week <- formatC(df_week2$week,flag = "0", width = 2)
df_week2$year_week <- paste(df_week2$year,"-", df_week2$week,sep = "")
df_week3 <- df_week2 %>% group_by(Country, Pangolin, year_week) %>% summarise(variant = sum(new))
total_seq <- df_week2 %>% group_by(Country, year_week) %>% summarise(num_sequenced = sum(new))
df_week3 <- left_join(df_week3, total_seq)

tran_data <- data.frame(date = seq(as.Date("2021-01-04"),as.Date("2021-09-26"), 7))
tran_data$year <- isoyear(tran_data$date)
tran_data$week <- isoweek(tran_data$date)
tran_data$week <- formatC(tran_data$week,flag = "0", width = 2)
tran_data$year_week <- paste(tran_data$year,"-", tran_data$week,sep = "")
tran_data$week <- as.numeric(tran_data$week)
tran_data <- tran_data[,c(1,4)]
df_week4 <- left_join(df_week3, tran_data)

#2. By month
df_month2 <- df1 %>% filter(Country == country) %>% group_by(Country, Collection_Date, Pangolin) %>% summarise(new = n())
df_month2 <- df_month2 %>% filter(Collection_Date <= as.Date("2021-09-26")) %>% filter(Collection_Date >= as.Date("2021-01-04"))
df_month2$month <- as.numeric(format(df_month2$Collection_Date, "%m"))
df_month2$year <- as.numeric(format(df_month2$Collection_Date, "%Y"))
df_month3 <- df_month2 %>% group_by(Country, Pangolin, year, month) %>% summarise(variant = sum(new))
total_seq <- df_month2 %>% group_by(Country, year, month) %>% summarise(num_sequenced = sum(new))
df_month4 <- left_join(df_month3, total_seq)
df_month4$month <- as.character(df_month4$month)

#==plot==
Sys.setlocale("LC_TIME", "US")
factor(df_week4$Pangolin, levels = c("Ref","Others","B.1.1.7","B.1.351","P.1","B.1.617.2")) -> df_week4$Pangolin
ggplot() + 
  geom_bar(data = df_week4, aes(x = date, y = variant, fill = Pangolin),
           position = "stack", stat = "identity",color = "black")+
  scale_y_continuous("No. of variants")+
  scale_x_date("Date of specimen collection",date_breaks = "1 month",  date_labels = "%b", expand = c(0.1,0))+
  coord_cartesian(clip = "off")+
  theme_bw()+
  labs(title = paste(country,"(by week)"))+
  theme(legend.position = "none")+
  scale_fill_manual("", values = c("grey80",pal_npg("nrc", alpha =1)(5)[5:1])) ->p1
p1

ggplot() + geom_bar(data = df_week4,aes(x = date, y = variant/num_sequenced, fill = Pangolin),
                    position = "stack", stat = "identity",color = "black")+
  scale_y_continuous("Proportion (%)", breaks = seq(0,1,0.25),labels = seq(0,100,25))+
  scale_x_date("Date of specimen collection",date_breaks = "1 month",  date_labels = "%b", expand = c(0.1,0))+
  coord_cartesian(ylim = c(0, 1),clip = "off")+
  theme_bw()+
  labs(title = paste(country,"(by week)"))+
  theme(legend.position = "none")+
  scale_fill_manual("", values = c("grey80",pal_npg("nrc", alpha =1)(5)[5:1])) ->p2
p2

factor(df_month4$Pangolin, levels = c("Ref","Others","B.1.1.7","B.1.351","P.1","B.1.617.2")) -> df_month4$Pangolin
ggplot() + 
  geom_bar(data = df_month4, aes(x = month, y = variant, fill = Pangolin),
           position = "stack", stat = "identity",color = "black",width = 0.6)+
  scale_y_continuous("No. of variants")+
  scale_x_discrete("Date of specimen collection", label = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep"),
                   expand = c(0.1,0))+
  coord_cartesian(clip = "off")+
  theme_bw()+
  labs(title = paste(country,"(by month)"))+
  # theme(axis.text.x = element_text(hjust = -0.3))+
  scale_fill_manual("", values = c("grey80",pal_npg("nrc", alpha =1)(5)[5:1])) ->p3
p3

ggplot() + geom_bar(data = df_month4,aes(x = month, y = variant/num_sequenced, fill = Pangolin),
                    position = "stack", stat = "identity",color = "black",width = 0.6)+
  scale_y_continuous("Proportion (%)", breaks = seq(0,1,0.25),labels = seq(0,100,25))+
  scale_x_discrete("Date of specimen collection", label = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep"),
                   expand = c(0.1,0))+
  coord_cartesian(ylim = c(0, 1),clip = "off")+
  theme_bw()+
  labs(title = paste(country,"(by month)"))+
  scale_fill_manual("", values = c("grey80",pal_npg("nrc", alpha =1)(5)[5:1])) ->p4
p4

#==output==
pdf("../output/Fig_1.pdf",width = 12, height = 10)
ggarrange(p1,p3,p2,p4,nrow  = 2, ncol = 2, widths = c(1,0.8))
dev.off()