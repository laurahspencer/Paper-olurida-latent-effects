load(file = "data/larval-survival") #object = survival

Collection <- read.csv("data/Larvae-collection.csv", header=T, stringsAsFactors = F)
Collection <- Collection[c(1,2,3,4,5,6,7,8,9,10,11,12,13,22)]
Collection$Date.collected <- as.Date(Collection$Date.collected, format = "%m/%d/%y")
Collection$Treatment <- as.factor(Collection$Treatment)
Collection$Bag.. <- as.factor(Collection$Bag..)
Collection$Broodstock <- as.numeric(Collection$Broodstock)
Collection$Vol..for.count..mL. <- as.numeric(Collection$Vol..for.count..mL.)
Collection$Count.A.LIVE <- as.numeric(Collection$Count.A.LIVE)
Collection$Count.A.DEAD <- as.numeric(Collection$Count.A.DEAD)
Collection$Count.B.LIVE <- as.numeric(Collection$Count.B.LIVE)
Collection$Count.B.DEAD <- as.numeric(Collection$Count.B.DEAD)
Collection$Count.C.LIVE <- as.numeric(Collection$Count.C.LIVE)
Collection$Count.C.DEAD <- as.numeric(Collection$Count.C.DEAD)
Collection$Live.Larvae <- (rowMeans(subset(Collection, select = c("Count.A.LIVE","Count.C.LIVE","Count.C.LIVE")), na.rm = TRUE))/Collection$Vol..for.count..mL.*Collection$Total.Vol..mL.
Collection$Dead.Larvae <- (rowMeans(subset(Collection, select = c("Count.A.DEAD","Count.C.DEAD","Count.C.DEAD")), na.rm = TRUE))/Collection$Vol..for.count..mL.*Collection$Total.Vol..mL.
Collection$Live.Larvae.norm <- Collection$Live.Larvae/Collection$Broodstock
Collection <- merge(x=Collection, y=unique(survival[,c("TRT.REP", "TEMP", "FOOD")]), by.x = "Treatment", by.y = "TRT.REP", all.x=T, all.y=T)

Collection$TREAT <- as.factor(paste(Collection$TEMP, "-", Collection$FOOD))
colnames(Collection) <- c("Rep", "Date", "Bag", "Broodstock", "Group", "Vol.sampled", "Vol.total", "Live.A", "Dead.A", "Live.B", "Dead.B", "Live.C", "Dead.C", "Sample.number", "Live.Larvae", "Dead.Larvae", "Live.Larvae.norm" ,"TEMP", "FOOD", "TREAT")

Collection$Perc.live <- (Collection$Live.Larvae/(Collection$Dead.Larvae + Collection$Live.Larvae))*100
survival

#save collection data R object for later merge with larval size & survival 
save(Collection, file = "data/larvae-collection-data") 
saveRDS(Collection, file = "data/larvae-collection-data.rds") 

# Calculate summary stats by treatment including reps 
treat_total.rep <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~ Date + TREAT + TEMP + FOOD + Bag, data = Collection, sum, na.rm = TRUE) %>%
  group_by(TEMP, FOOD, TREAT, Bag) %>% 
  mutate(cum.total=cumsum(Live.Larvae),cum.percap = cumsum(Live.Larvae.norm),CalDay = format(Date,"%j")) %>% 
  group_by(TEMP, FOOD, TREAT, Bag) %>% dplyr::summarize(overall_Total = sum(Live.Larvae, na.rm = T), mean.Live.Larvae = mean(Live.Larvae, na.rm=T), mean.percap = mean(Live.Larvae.norm, na.rm=T), total.percap = sum(Live.Larvae.norm, na.rm=T), maxday = as.numeric(CalDay[which.max(Live.Larvae)]), max = max(Live.Larvae), max.percap = max(Live.Larvae.norm), first.big = as.numeric(CalDay[which(Live.Larvae > 10000)[1]])-58, release.days = as.numeric(length(CalDay[Live.Larvae > 10000])))

treat_total.rep <- merge(x=treat_total.rep, y=aggregate(Broodstock ~ Bag, data=Collection, median), by.x = "Bag", by.y = "Bag")[,-1]
treat_total.rep$trial <- c("13") 

# Estimate # and % broodstock that reproduced as females 
# treat_total.rep$noFemales <- treat_total.rep$overall_Total/215000 #estimate of # females 
# treat_total.rep$perc.spawn <- 100*(treat_total.rep$noFemales/treat_total.rep$Broodstock) #estimate of # females 
# ggplot(data=treat_total.rep, aes(x=TEMP:FOOD, y=perc.spawn)) + geom_boxplot() + geom_point()
# summary(treat_total.rep$perc.spawn)

# Treatment total (all reps combined), by date 
treat_total <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~ Date + TREAT + TEMP + FOOD, data = Collection, sum, na.rm = TRUE)


# ======================
# Statistical Analysis 
# ======================

# Larval release timing 

## ----------- First big day - RELEASE ONSET 

### visualize data 
ggplot(treat_total.rep, aes(x=FOOD:TEMP, y=first.big)) +    
  geom_violin() + geom_jitter(aes(color=FOOD:TEMP), size=3, alpha=0.8) +
  theme_minimal() +coord_flip() +
  scale_color_manual(values=c('#92c5de','#f4a582','#0571b0','#ca0020'), 
                     name=element_blank(), 
                     labels = c("Low-food+7°C","Low-food+10°C",
                                "High-food+7°C","High-food+10°C"),
                     guide = guide_legend(reverse = TRUE)) + 
  labs(title="Timing of larval release onset", 
       y=("No. days in spawning conditions")) + 
  theme(axis.title.y = element_blank(), legend.position = "none",
        plot.title = element_text(size=12)) +
  scale_x_discrete(labels=c("High:Warm"="High-food+10°C", 
                            "High:Cold"="High-food+7°C", 
                            "Low:Warm"="Low-food+10°C", 
                            "Low:Cold"="Low-food+7°C")) +
  scale_y_continuous(breaks=c(31,33,35,37,39), limits =c(30.5, 40))

### check for outliers 
treat_total.rep %>% group_by(FOOD, TEMP) %>% 
  select(FOOD, TEMP, first.big) %>% 
  identify_outliers(first.big) #no extreme outliers 

### check other ANOVA assumptions 
hist(treat_total.rep$first.big, breaks = 16, main="response variable distribution")
lm.big <- lm(first.big~FOOD*TEMP, data=treat_total.rep)
bc <- boxcox(first.big ~ FOOD*TEMP, data=treat_total.rep)
(lambda <- bc$x[which.max(bc$y)])
lm.big.bc <- lm(first.big^lambda~FOOD*TEMP, data=treat_total.rep)

hist(lm.big.bc$residuals, breaks=16, main="residual distribution")
hist(stdres(lm.big.bc), breaks=16, main="standardized residual distribution")
shapiro_test(stdres(lm.big.bc))
shapiro_test(lm.big.bc$residuals)

ggqqplot(stdres(lm.big.bc)) #residuals have some issues 
ggqqplot(residuals(lm.big.bc)) #residuals have some issues 

# run ANOVA, but doesn't meet assumptions (residuals are an issue); will do another "robust" anova in a sec
Anova(lm.big.bc, type="II") # run anova despite residuals' deviance from normality
treat_total.rep %>% levene_test(first.big ~ TEMP*FOOD) #looks good 

shapiro_test(residuals(lm.big)) # residuals not normal 
ggqqplot(treat_total.rep, "first.big",ggtheme=theme_bw()) + facet_grid(TEMP~FOOD) # not too bad, though  
treat_total.rep %>% levene_test(first.big ~ TEMP*FOOD) #looks good 

# Run 2-way ANOVA with robust estimation from WRS2 package 
pbad2way(first.big ~ FOOD*TEMP,
         data = treat_total.rep,
         est = "mom",    # modified M-estimator
         nboot = 5000)   # number of bootstrap samples
# as appears in the boxplot, food is a significant factor. Earlier release in low food groups. 

# try again with two-way ANOVA for trimmed means with interaction effects 
t2way(first.big ~ FOOD*TEMP,
      data = treat_total.rep, tr = 0.2) # again, no interaction 

# now use Kruskal Wallis test to examine chi-square for food as main effect 
kruskal.test(first.big ~ FOOD, data = treat_total.rep) 

# Check main effect of food, too, since the standard ANOVA indicated that it could be sign.  
kruskal.test(first.big ~ TEMP, data = treat_total.rep) 


## ----------- Max day - DATE OF MAXIMUM RELEASE 
ggplot(treat_total.rep, aes(x=FOOD:TEMP, y=maxday)) + 
  geom_violin() + geom_jitter(aes(color=FOOD:TEMP)) + theme_minimal()

treat_total.rep %>% group_by(FOOD, TEMP) %>% 
  select(FOOD, TEMP, maxday) %>% 
  identify_outliers(maxday) #no extreme outliers 

# check out ANOVA assumptions 
lm.max <- lm(maxday~FOOD*TEMP, data=treat_total.rep)
ggqqplot(residuals(lm.max)) #looks good 
shapiro_test(residuals(lm.max)) #good 
ggqqplot(treat_total.rep, "maxday",ggtheme=theme_bw()) + facet_grid(TEMP~FOOD) #looks good 
treat_total.rep %>% levene_test(maxday ~ TEMP*FOOD) #looks good 

# run ANOVA 
Anova(lm(maxday ~ TEMP*FOOD, treat_total.rep), type="II") #TEMP:FOOD interaction, but after rounding no sign. effect  


## -----------  Average daily no. of larvae released (~brood size)
ggplot(treat_total.rep, aes(x=FOOD:TEMP, y=mean.Live.Larvae)) + 
  geom_violin() + geom_jitter(aes(color=FOOD:TEMP)) + theme_minimal()

treat_total.rep %>% group_by(FOOD, TEMP) %>% 
  select(FOOD, TEMP, mean.Live.Larvae) %>% 
  identify_outliers(mean.Live.Larvae) #no extreme outliers 

# check out ANOVA assumptions 
lm.daily <- lm(mean.Live.Larvae~FOOD*TEMP, data=treat_total.rep)
ggqqplot(residuals(lm.daily)) #looks good 
shapiro_test(residuals(lm.daily)) #good 
ggqqplot(treat_total.rep, "mean.Live.Larvae",ggtheme=theme_bw()) + facet_grid(TEMP~FOOD) #looks good 
treat_total.rep %>% levene_test(mean.Live.Larvae ~ TEMP*FOOD) #looks good 

# run ANOVA 
anova(lm(mean.Live.Larvae ~ TEMP*FOOD, treat_total.rep))   


# -------- Total release not normalized
ggplot(treat_total.rep, aes(x=FOOD:TEMP, y=overall_Total)) + 
  geom_violin() + geom_jitter(aes(color=FOOD:TEMP)) + theme_minimal()
treat_total.rep %>% group_by(FOOD, TEMP) %>% 
  select(FOOD, TEMP, overall_Total) %>% 
  identify_outliers(overall_Total) #no  outliers 

# check out ANOVA assumptions 
lm.total <- lm(overall_Total~FOOD*TEMP, data=treat_total.rep)
ggqqplot(residuals(lm.total)) #looks good 
shapiro_test(residuals(lm.total)) #good 
ggqqplot(treat_total.rep, "overall_Total",ggtheme=theme_bw()) + facet_grid(TEMP~FOOD) #looks good 
treat_total.rep %>% levene_test(overall_Total ~ TEMP*FOOD) #looks good 

# run ANOVA 
anova(lm(overall_Total ~ TEMP*FOOD, treat_total.rep)) #no effects

# ---------- Total release normalized by broodstock - also check, but don't report 
ggplot(treat_total.rep, aes(x=FOOD:TEMP, y=total.percap)) + 
  geom_violin() + geom_jitter(aes(color=FOOD:TEMP)) + theme_minimal()
treat_total.rep %>% group_by(FOOD, TEMP) %>% 
  select(FOOD, TEMP, total.percap) %>% 
  identify_outliers(total.percap) #no  outliers 

# check out ANOVA assumptions 
lm.percap <- lm(total.percap~FOOD*TEMP, data=treat_total.rep)
ggqqplot(residuals(lm.percap)) #looks good 
shapiro_test(residuals(lm.percap)) #good 
ggqqplot(treat_total.rep, "total.percap",ggtheme=theme_bw()) + facet_grid(TEMP~FOOD) #looks good 
treat_total.rep %>% levene_test(total.percap ~ TEMP*FOOD) #looks good 

# run ANOVA 
anova(lm(total.percap ~ TEMP*FOOD, treat_total.rep)) #no effects 

# ===========================
# SUPPLEMENTAL MATERIALS 
# ===========================

# ===========================
# Larval Production Figures 
# ===========================

Collection.total <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~ TREAT + TEMP + FOOD + Rep + Bag, data = Collection, sum, na.rm = TRUE)

png(filename = "results-figures/total-larvae-released-points.png", width=400, height = 300)
ggplot(Collection.total, aes(x=TEMP:FOOD, y=Live.Larvae.norm)) + 
  geom_jitter(width=0.2, size=4, shape=21, aes(color=TEMP:FOOD), stroke=1.5) + 
  theme_bw() + labs(title="Total larvae released", y=("Total released (per adult)")) + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), 
                     name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", 
                                                      "Warm / High Food", "Warm / Low Food")) + 
  theme(text = element_text(size = 12)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  scale_y_continuous(labels=scales::comma,lim=c(30000,170000))
dev.off()


#Calculate cumulative larvae released through time for plots 
Collection.cumul <- aggregate(cbind(Live.Larvae, Live.Larvae.norm) ~ Date + TREAT + TEMP + FOOD, data = Collection, sum, na.rm = TRUE) %>%
  group_by(TEMP, FOOD, TREAT) %>% 
  mutate(cum.total=cumsum(Live.Larvae),cum.percap = cumsum(Live.Larvae.norm),CalDay = format(Date,"%j")) %>% 
  arrange(Date) %>% dplyr::select(Date,CalDay,TREAT,TEMP,FOOD,Live.Larvae,Live.Larvae.norm, cum.total,cum.percap)

# Relevel food factors such that Low comes before High
#Collection.cumul$FOOD <- factor(Collection.cumul$FOOD, levels = rev(levels(Collection.cumul$FOOD)))

TREATS <- levels(Collection$TREAT)
TREATS.col <- c("#0571b0","#92c5de","#ca0020","#f4a582")

p <- list()
p[[1]] <- ggplot(data=Collection.cumul, aes(x=Date, y=cum.percap, group=TEMP:FOOD, color=TEMP:FOOD)) + 
  geom_line(size=.75) + 
  scale_color_manual(values=TREATS.col, name=element_blank(),
                     labels = c("Cold / High Food", "Cold / Low Food",
                                "Warm / High Food", "Warm / Low Food")) +
  theme_classic(base_size = 14) + 
  scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", 
               limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + 
  theme(legend.position = "top", axis.title.x = element_blank(), axis.title.y = element_blank(), 
        plot.title = element_text(size = 14, hjust = 0), title = element_blank(), 
        axis.line = element_line(size = .5, colour = "gray50")) + 
  scale_y_continuous(limits = c(0, 400000), position = "right") +
  guides(color=guide_legend(nrow=2,byrow=TRUE))

for (i in 1:length(TREATS)) {
  p[[i+1]] <- ggplot(data=subset(Collection.cumul, TREAT==TREATS[i]), aes(x=Date, y=Live.Larvae, fill=TREAT)) + 
    geom_bar(fill=TREATS.col[i], stat="identity",width=.5, position = position_dodge(width=2)) + 
    theme_classic(base_size = 14) +  
    theme(plot.title = element_text(size = 14, hjust = 0), axis.title.y = element_blank(), 
          axis.title.x = element_blank(),  axis.title = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none",
          axis.line = element_line(size = .5, colour = "gray50")) + 
    scale_x_date(date_breaks = "1 week",date_labels ="%b-%d", 
                 limits=c(as.Date("2018-03-30"), as.Date("2018-04-27"))) + 
    scale_y_continuous(labels = scales::scientific, 
                       limits = c(0,3694094), breaks=c(2500000, 5000000), position = "right") +
    geom_smooth(method="loess", color="gray60", size=0.6, 
                linetype="dashed", span = 0.25, method.args = list(degree=1), se = FALSE)
}
pdf(file = "results-figures/larval-production-timeline.pdf", width = 5, height = 8)
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], align = "v", nrow = 5, rel_heights = c(2/6, 1/6, 1/6, 1/6, 1/6))
dev.off()


# ---------------- Estimated % spawn as females   
shapiro.test(treat_total.rep$perc.spawn) # normal 
hist(treat_total.rep$perc.spawn) # normal 
summary(aov(perc.spawn ~ TEMP*FOOD, data = treat_total.rep)) #No diff 
aggregate(perc.spawn ~ TEMP+FOOD, data=treat_total.rep, mean) #mean percent by treatments  
mean(treat_total.rep$perc.spawn)
sd(treat_total.rep$perc.spawn)

# --------------- Summary stats on larval collection 
mean(Collection$Live.Larvae)

aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, mean, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, median, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, max, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, sum, na.rm = TRUE)
aggregate(Broodstock ~ TEMP + FOOD, data = Collection, median, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD, data = Collection, sd, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD + Bag, data = Collection, sum, na.rm = TRUE)
aggregate(Live.Larvae ~ TEMP + FOOD + Bag, data = Collection, sum, na.rm = TRUE)
