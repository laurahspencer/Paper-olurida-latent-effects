# Adult size and mortality 

# read in data 
histo <- read.csv("data/histology.csv", header=T, stringsAsFactors = T, na.strings = c("NA", " NA ", "TBD"))

#Convert a few columns to factors & reorder
histo$Week <- as.factor(histo$Week)
histo$Bag.. <- as.factor(histo$Bag..)
histo$FEMSTAGE.COL <- as.factor(histo$FEMSTAGE.COL)
histo$MALSTAGE.COL <- as.factor(histo$MALSTAGE.COL)
histo$SEX <- factor(histo$SEX, levels=c("I", "M", "HPM", "H", "HPF", "F"))
histo$TREAT <- factor(histo$TREAT, levels=c("PRE", "A" , "B", "C", "D", "Wild"))
histo$Date <- factor(histo$Date, levels=c("11/30/17", "12/20/17", "1/23/18", "1/4/18", "2/27/18", "2/9/18", "3/13/18", "3/23/18"))

histo <- histo %>% 
  # add treatment information 
  mutate(TEMP=as.factor(str_replace(str_replace(TREAT, "A|C", "Cold"), "B|D", "Warm"))) %>% 
  mutate(FOOD=as.factor(str_replace(str_replace(TREAT, "A|B", "Low"), "C|D", "High"))) %>%
  filter(TREAT!="Wild") %>% droplevels() #remove those collected from the wild for this paper 

# Generate summary dataframe statistics 
size.adult <- histo %>%
  filter(!is.na(TREAT)) %>%
  group_by(TREAT, Week) %>%
  summarize(mean_length = mean(Length..cm., na.rm = TRUE), 
            sd_length = sd(Length..cm., na.rm = TRUE),
            mean_weight = mean(Est..Tissue.Weight..g., na.rm = TRUE),
            sd_weight = sd(Est..Tissue.Weight..g., na.rm = TRUE))

#add same "pre" values for all treatments  
size.adult[(nrow(size.adult)+1):(nrow(size.adult)+4),] <- size.adult[1,] 
size.adult[(nrow(size.adult)-3):(nrow(size.adult)),"TREAT"] <- as.factor(c("A", "B", "C", "D"))
# ====== WHY?

# How many samples per treatment & week? 
# Tissue wet weight 
histo %>% filter(Week!=0) %>% 
  group_by(Week, TEMP, TREAT) %>% 
  dplyr::select(Week, TEMP, FOOD, Est..Tissue.Weight..g.) %>%
  tally()

# Shell length  
histo %>% filter(Week!=0) %>% 
  group_by(Week, TEMP, TREAT) %>% 
  dplyr::select(Week, TEMP, FOOD, Length..cm.) %>% 
  tally()

# Examine weight and length over time 
# NOTE: balanced design for all weeks & treatments EXCEPT 
# tissue weight for week 16 - missing data for Cold+High Food 


# =============================================
# Does weight change over time, diff by treat? 
# ============================================

# plot weight over time by treatments 
ggplot(data=size.adult %>% filter(Week!=15 & Week !=16), aes(x=Week, y=mean_weight, group=TREAT, col=TREAT)) +
  #geom_errorbar(aes(ymin=mean_length-sd_length, ymax=mean_length+sd_length), width=.1) +
  geom_line()+
  geom_point()

# ANCOVA analysis using TEMP and FOOD as separate factors, Only use weeks during treatments 
hist(subset(histo, Week!=15 & Week!=16)$Est..Tissue.Weight..g.^(1/3)) #cube-root transform 
shapiro_test(subset(histo, Week!=15 & Week!=16)$Est..Tissue.Weight..g.^(1/3))

weight.ancova <- lm(Est..Tissue.Weight..g.^(1/3) ~ FOOD*TEMP*as.numeric(as.character(Week)), 
                    data=subset(histo, Week!=15 & Week!=16))
plot(weight.ancova, add.smooth = FALSE, which=1) #assess linearity - looks good 
plot(weight.ancova, which=2) #assess normality - looks pretty good 
plot(weight.ancova, add.smooth=FALSE, which=3) #assess variance - looks good 
anova(weight.ancova) # does weight change over time, and does treatment effect rate of change? 
# result: no effect of food or temp on weight, no sign. change over time (despite looking like it does)

anova(lm(Est..Tissue.Weight..g.^(1/3) ~ as.numeric(as.character(Week)), 
   data=subset(histo, Week!=15 & Week!=16))) #how about just with week? # no 

# Redo analysis as suggested by JEMBE editor - 3-way ANOVA with week as factor 
summary(aov(Est..Tissue.Weight..g.^(1/3) ~ FOOD*TEMP*as.factor(Week), data=subset(histo,Week!="0" & Week!=15 & Week!=16))) #no effect of 
# No effect of temp, food, or week on weight 

# plot mean weight over time 

pdf("results-figures/broodstock-weight.pdf", width = 6, height = 3.5)
ggplot(data=subset(size.adult, TREAT!="PRE" & mean_weight<3.5), 
       aes(x=Week, y=mean_weight/mean_length, group=TREAT, col=TREAT)) + 
  geom_line() + theme_bw(base_size = 11) + geom_point(size=3) + 
  ylab("Mean wet tissue weight / shell length (g/cm)") + 
  ggtitle(label = "Adult tissue weight (standardized by length) over time") + 
  theme(axis.title.x = element_blank(), legend.position = "none") +
  scale_color_manual(
    values=c("#0571b0","#92c5de","#ca0020","#f4a582","gray30"),
    name="Treatment",
    breaks=c("C", "A", "D", "B", "Wild"),
    labels=c("7°C+high-food", "7°C+low-food", "10°C+high-food", "10°C+low-food", "Wild")) +
  scale_x_discrete(labels= c("Nov 30", "Dec 20", "Jan 4", "Jan 23", "Feb 9", "Feb 27", "Mar 3", "Mar 23")) 
dev.off()


# =============================================
# Does length change over time, diff by treat? 
# ============================================

# plot mean length over time 
ggplot(data=size.adult %>% filter(Week!=15 & Week !=16), aes(x=Week, y=mean_length, group=TREAT, col=TREAT)) +
  #geom_errorbar(aes(ymin=mean_length-sd_length, ymax=mean_length+sd_length), width=.1) +
  geom_line()+
  geom_point()

hist(subset(histo, Week!=15 & Week!=16)$Length..cm.)  
shapiro_test(subset(histo, Week!=15 & Week!=16)$Length..cm.)
length.ancova <- lm(Length..cm. ~ FOOD*TEMP*as.numeric(as.character(Week)),
                    data=subset(histo,Week!=15 & Week!=16))

plot(length.ancova, add.smooth = FALSE, which=1) #assess linearity - looks good 
plot(length.ancova, which=2) #assess normality - looks good 
plot(length.ancova, add.smooth=FALSE, which=3) #assess variance - looks good 
anova(length.ancova) # does length change over time, and does treatment effect rate of change? 
# result: no effect of food or temp on length

anova(lm(Length..cm. ~ as.numeric(Week), data=subset(histo,Week!=15 & Week!=16)))
# Length does not change 

# Redo analysis as suggested by JEMBE editor - 3-way ANOVA with week as factor 
anova(lm(Length..cm. ~ Week*FOOD*TEMP, data=subset(histo,Week!="0" & Week!=15 & Week!=16)))
# No effect of temp, food, or week on shell length 


# =====================
# ADULT MORTALITY 
# =====================

brood.mortality <- read.csv("data/adult-mortality.csv", header=T, stringsAsFactors = F)
brood.mortality$Date <- as.Date(brood.mortality$Date, format = "%m/%d/%y")
brood.mortality <- melt(data = brood.mortality, id.vars = "Date", value.name = "Alive", variable.name = "TREAT.rep")
brood.mortality$TREAT <- brood.mortality$TREAT.rep
brood.mortality$TREAT <- gsub("1", "", brood.mortality$TREAT)
brood.mortality$TREAT <- as.factor(gsub("2", "", brood.mortality$TREAT))

# Survival Analysis 
survdiff(data = brood.mortality, formula = Surv(Alive) ~ TREAT)

# ======================
# ADULT MORTALITY FIGURE
# ======================

pdf(file="results-figures/broodstock-survival.pdf", width=7, height=3)
ggplot(data=brood.mortality, aes(x=Date, y=100*Alive, group=TREAT.rep, col=TREAT)) + 
  theme_bw(base_size = 13) + ggtitle("Adult survival over time") + 
  geom_step()+ geom_point() + xlab("Date") + ylab("Percent survival") + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), 
                     name=element_blank(), 
                     labels = c("Cold / High Food", "Cold / Low Food", "Warm / High Food", "Warm / Low Food")) + 
  scale_x_date(date_breaks = "3 week",date_labels ="%b %d", 
               limits = c(min(brood.mortality$Date), max(brood.mortality$Date))) +
  scale_y_continuous(breaks=c(50, 60, 70, 80, 90, 100)) +
  geom_vline(xintercept = as.numeric(as.Date("2018-02-28")), linetype="dashed", color = "gray50", size=.5) +
  theme(legend.position = "none", axis.title.x = element_blank())

dev.off()

