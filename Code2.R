#==load packages==
library(dplyr)
library(ggplot2)
library(stringr)
library(ggsci)
library(ggpubr)
library(sf)
library(rgdal)

#==read data==
df <- read.csv("../data/sequences.csv") # read genome data
worldmap <- st_read("../data/World_Countries/World_Countries.shp") # read world map data

#====================
#==data clean==
#====================
#Step1. Remove sequences from non-human
df <- df %>% filter(Host == "Homo sapiens") 

#Step2. Remove sequences with incomplete date of collection (only year) or sampled before 1 December 2019
df$Collection_Date <- as.Date(df$Collection_Date)
df <- df %>%  filter(!is.na(Collection_Date))  %>% filter(Collection_Date >= as.Date("2019-12-01"))

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

#====================
#==Fig2. data clean==
#====================
st_crs(worldmap)
worldmap <- st_transform(worldmap, CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"))
tmp2 <- df1 %>% group_by(Country, Collection_Date, Pangolin) %>% summarise(new = n())
tmp2 <- tmp2 %>% filter(Collection_Date <= as.Date("2021-09-26")) %>% filter(Collection_Date >= as.Date("2021-01-04"))
tmp2$month <- as.numeric(format(tmp2$Collection_Date, "%m"))
tmp2$year <- as.numeric(format(tmp2$Collection_Date, "%Y"))
tmp3 <- tmp2 %>% group_by(Country, Pangolin, year, month) %>% summarise(variant = sum(new))
total_seq <- tmp2 %>% group_by(Country, year, month) %>% summarise(num_sequenced = sum(new))
tmp3 <- left_join(tmp3, total_seq)
tmp3$prevalence <- round(tmp3$variant * 100 / tmp3$num_sequenced,1)

worldmap1 <- left_join(worldmap, tmp3[tmp3$Pangolin == "B.1.617.2" & tmp3$month == 2,], by=c("COUNTRY" = "Country"))
worldmap2 <- left_join(worldmap, tmp3[tmp3$Pangolin == "B.1.617.2" & tmp3$month == 5,], by=c("COUNTRY" = "Country"))
worldmap3 <- left_join(worldmap, tmp3[tmp3$Pangolin == "B.1.617.2" & tmp3$month == 8,], by=c("COUNTRY" = "Country"))

theme1 <- theme(axis.ticks = element_blank(), 
                axis.line = element_blank(),
                axis.text = element_blank(),
                legend.justification=c(0,0), 
                legend.position = "bottom",
                legend.background = element_blank(),
                # legend.direction  = "horizontal",
                legend.key.size = unit(0.3,"cm"),
                legend.key.width = unit(3.5,"cm"),
                legend.title = element_text(size = 11), 
                legend.text = element_text(size = 10,lineheight=12), 
                legend.spacing = unit(5, "cm"), legend.spacing.x = NULL, # the spacing between legends (unit)
                legend.spacing.y = NULL,#the spacing between legends (unit)
                axis.title = element_blank(),
                plot.margin =  margin(0, 0, 0, 0, "cm"),
                panel.background = element_rect(fill = "white"),
                legend.box.background =element_blank(),
                legend.box.margin=  margin(0, 0, 0, 0, "cm"),
                panel.spacing = unit(0,"cm"),
                panel.border = element_rect(fill='transparent',colour="transparent"),
                legend.key = element_blank(),
                plot.title = element_text(hjust = 0.05, size = 15, vjust = 0))
#plot
ggplot() + 
  geom_sf(data = worldmap1, aes(fill = prevalence), size = 0.4,color = "black")+
  coord_sf(xlim = c(-170, 170), ylim = c(-57, 90))+
  theme(axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white"))+
  labs(title = "Prevalence of Delta variant in Feb")+
  scale_fill_gradientn("",colours=c("red","orange","yellow","green","blue"),
                       na.value="white",breaks = seq(0, 100, 25),
                       labels = c("0","25","50","75","100"),limits = c(0, 100))-> fig1
fig1

ggplot() + 
  geom_sf(data = worldmap2, aes(fill = prevalence), size = 0.4,color = "black")+
  coord_sf(xlim = c(-170, 170), ylim = c(-57, 90))+
  theme(axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white"))+
  labs(title = "Prevalence of Delta variant in May")+
  scale_fill_gradientn("",colours=c("red","orange","yellow","green","blue"),
                       na.value="white",breaks = seq(0, 100, 25),
                       labels = c("0","25","50","75","100"),limits = c(0, 100))-> fig2
fig2

ggplot() + 
  geom_sf(data = worldmap3, aes(fill = prevalence), size = 0.4,color = "black")+
  coord_sf(xlim = c(-170, 170), ylim = c(-57, 90))+
  theme(axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white"))+
  labs(title = "Prevalence of Delta variant in Aug")+
  scale_fill_gradientn("",colours=c("red","orange","yellow","green","blue"),
                       na.value="white",breaks = seq(0, 100, 25),
                       labels = c("0","25","50","75","100"),limits = c(0, 100))-> fig3
fig3

#==output==
pdf("../output/Fig_2.pdf",width = 10, height = 12.5)
ggarrange(fig1,fig2,fig3,nrow  = 3)
dev.off()
