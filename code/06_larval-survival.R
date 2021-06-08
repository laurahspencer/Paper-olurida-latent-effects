survival <- read.csv("data/Survival.csv", header=T, na.strings = "NA", stringsAsFactors = F, colClasses=
                       c(rep("character", times=2), "numeric", rep("factor", times=4), rep("numeric", times=4), "character", rep("numeric", times=5), rep("character", times=6)))
survival <- survival[c(-18:-23)]
survival$Date.stocked <- as.Date(survival$Date.stocked, format = "%m/%d/%y")
survival$Date.initial.count <- as.Date(survival$Date.initial.count, format = "%m/%d/%y")
survival$Date.imaged <- as.Date(survival$Date.imaged, format = "%m/%d/%y")
survival$TEMP <- survival$TRT
survival$TEMP <- gsub("\\<A\\>|\\<C\\>", "Cold", survival$TEMP)
survival$TEMP <- gsub("\\<B\\>|\\<D\\>", "Warm", survival$TEMP)
survival$FOOD <- survival$TRT
survival$FOOD <- gsub("\\<A\\>|\\<B\\>", "Low", survival$FOOD)
survival$FOOD <- gsub("\\<C\\>|\\<D\\>", "High", survival$FOOD)
survival$FOOD <- as.factor(survival$FOOD)
survival$TEMP <- as.factor(survival$TEMP)
survival$Dead.50.days <- (3*800)-survival$Live.50.days
survival$Dead.35.days <- 800-survival$Live.35.days
survival$TREAT <- as.factor(paste(survival$TEMP, "-", survival$FOOD))

# Save larval survival R object for use in collection analysis (just for treatment reps)
save(survival, file = "data/larval-survival") #object = survival

# General summary stats 
summary(100*(na.omit(survival$Live.35.days)/800)) #survival to 35 days 

# How long was larval collection period? 
max(survival$Date.stocked)-min(survival$Date.stocked)


# Create new dataframe that averages survival for each larval group (both percentage, and count)
# NOTE: survival was assessed at day 35 for each triplicate silo per biol. replicate, then again at day 50 for each biol. rep. 

Survival.family.35 <- do.call(data.frame, aggregate(Live.35.days ~ Family+TEMP+FOOD+TRT.REP+TREAT, data = survival, FUN = function(x) c(mean = mean(x)/(800)*100, sd = sd(x)/(800)*100, mean.live = mean(x), mean.dead = (800*3)-mean(x))))

names(Survival.family.35) <- c("Family", "TEMP", "FOOD", "TRT.REP", "TREAT", "Mean.Live.35", "SD.Live.35", "mean.live", "mean.dead")

# Plot survival again, this time  one point showing average survival of each group 
ggplot(Survival.family.35, aes(x=FOOD:TEMP, y=Mean.Live.35)) + 
  geom_boxplot(aes(color=FOOD:TEMP), outlier.shape = NA) + 
  geom_jitter(width=0.35, size=2, aes(color=FOOD:TEMP)) + 
  theme_bw(base_size = 12) + 
  labs(title="Mean % survival, by parental treatment", y=("% survival"), x="Treatment") + 
  scale_color_manual(values=c('#92c5de','#f4a582','#0571b0','#ca0020'), name=element_blank(), 
                      labels = c("Cold+\nLow Food", "Warm+\nLow Food", "Cold+\nHigh Food", "Warm+\nHigh Food")) + theme(text = element_text(size = 12), legend.position = "none") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())

# ======================
# Statistical analysis 
# ======================

# --------- Test for effects of treatment on larval survival (using average survival for each biological replicate)
Anova(glm(cbind(mean.live, mean.dead) ~ FOOD*TEMP, data=Survival.family.35, quasibinomial)) 

# ---------- Test for other factors 

# Read in larval production & size data to merge with survival data to assess various factors 
load(file = "data/larvae-collection-data") #save for later use 
load(file = "data/larval.size") #object = average.all

# Merge survival data with collection data 
survival.collect <- merge(x=survival, y=Collection[c(5,14,15,16,20)], by.x="Family", by.y="Group", all.x =TRUE, all.y=FALSE)

# Create a master dataframe of larval survival, larval size, and collection data to test for various effects on survival 
survival.collect$Sample.number <- as.factor(survival.collect$Sample.number) #first convert sample # to factor 
master <- merge(x=average.all[,c("Sample.number", "Length", "Width")], 
                y=survival.collect, by="Sample.number", all.x=T, all.y=T) %>% 
  mutate(Sample.number=as.numeric(Sample.number)) %>%
  mutate(ethanol=ifelse(Sample.number < 50, 'YES', 'NO')) %>%
mutate(length.cor = ifelse(Sample.number<50, Length-5.620869, Length)) 

# Run GLMS to test for effects of other factors on survival  
# Use quasibinomial model; design is balanced 

# ----------- Test ALL possible factors in one model - same result 
Anova(glm.all <- glm(cbind(Live.35.days, Dead.35.days) ~ 
                       FOOD + TEMP + FOOD:TEMP + Live.Larvae + length.cor + Date.stocked,
                     data=master, quasibinomial))

Anova(glm.all <- glm(cbind(Live.50.days, Dead.50.days) ~ 
                       FOOD + TEMP + FOOD:TEMP + Live.Larvae + length.cor + Date.stocked,
                     data=master, quasibinomial))

# --------- Date stocked sign. factor
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=master, quasibinomial)) 
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Date.stocked, data=master, quasibinomial)) #positive estimate (surv. increased over time)

# Plot survival by date larvae were released factor 
ggplot(master, aes(x=Date.stocked, y=100*Live.35.days/(800))) + 
  geom_point(size=2.5) + theme_bw(base_size = 14) + 
  labs(title="% larval survival ~ Date stocked", 
       y=("% survival"), 
       x="# live larvae collected") + 
  geom_smooth(method = "auto", color="gray50")

# ---------- No. of live larvae collected sign. factor
Anova(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae, data=master, quasibinomial)) 
summary(glm(cbind(Live.35.days, Dead.35.days) ~ Live.Larvae, data=master, quasibinomial)) #positive estimate (surv. higher from large batches)
ggplot(master, aes(x=Live.Larvae, y=100*Live.35.days/(800))) + 
  geom_point(size=2.5) + theme_bw(base_size = 14) + 
  labs(title="% larval survival ~ # larvae collected", 
       y=("% survival"), 
       x="# live larvae collected") + 
  geom_smooth(method = "auto", color="gray50")

# ----------- Larval length (aka width)

#test again using random effects model to include use of ethanol as preservation - definitely not sign.
Anova(lme(cbind(Live.35.days, Dead.35.days) ~ Length, random=~1|ethanol, data=na.omit(master))) 
Anova(lme(cbind(Live.50.days, Dead.50.days) ~ Length, random=~1|ethanol, data=na.omit(master))) 

# =======================
# LARVAL SURVIVAL FIGURE
# ======================= 

# first reverse order of food factor 
survival$FOOD <- factor(survival$FOOD, levels = rev(levels(survival$FOOD)))

# Calculate average larval survival by parental winter treatment - to add to boxplots 
survival.mean <- setNames(aggregate(100*Live.35.days/(800) ~ FOOD:TEMP, data=survival, mean), c("FOOD", "TEMP", "surv.mean"))

pdf(file="results-figures/larval-survival.pdf", width =4, height = 7.5)
ggplot(survival, aes(x=TEMP:FOOD, y=100*Live.35.days/(800))) + 
  geom_boxplot(aes(color=TEMP:FOOD), outlier.shape = NA) + 
  geom_jitter(width=0.35, size=2, aes(color=TEMP:FOOD)) + 
  theme_bw(base_size = 12) + 
  labs(title="% Survival, by parental treatment", y=("% Survival"), x="Treatment replicate") + 
  scale_color_manual(values=c('#92c5de', '#0571b0','#f4a582', '#ca0020'), name=element_blank(),
                     labels = c( "7°C+Low-food", "7°C+High-food","10°C+Low-food","10°C+High-food")) +
  theme(text = element_text(size = 12), legend.position = "bottom") + guides(color=guide_legend(nrow=2,byrow=F)) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  geom_point(data=survival.mean, aes(y=surv.mean, fill=TEMP:FOOD), col="gray40", shape=22, size=4, show.legend=F) +
  scale_fill_manual(values=c('#92c5de','#0571b0','#f4a582','#ca0020'))
dev.off()


# ==============================
# SUPPLEMENTAL MATERIALS FIGURES 
# ==============================

# Relevel food factors such that Low comes before High
master$FOOD <- factor(master$FOOD, levels = rev(levels(master$FOOD)))

# ------------- Plot larval survival ~ length (corrected), with fitted lm  
png(file="results-figures/larval-surv-size.png", width = 600, height = 325)
subset(master, ethanol=="NO") %>% dplyr::select(Length, Live.35.days, TREAT.x) %>% 
  filter(!is.na(TREAT.x), !is.na(Length)) %>% 
  ggplot(aes(x=Length, y=100*Live.35.days/(800))) + 
  geom_point(size=2, color="gray40") + 
  labs(title="Larval survival ~ shell width upon release", 
       y=("% Survival"), 
       x="Shell Length (µm)") + 
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  geom_smooth(method='lm', color="gray50") +
  annotate("text", x=175, y=40, color="gray30", 
            label="Adjusted R-squared=0.012\np-value=0.162")
dev.off()

# ------------plot larval survival by date collecte/releasd 
png(file="results-figures/larval-surv-time-col.png", width = 600, height = 325)
ggplot(master, aes(x=Date.stocked, y=100*(Live.35.days/800))) + theme_bw(base_size = 14) + 
  geom_point(size=2.5, aes(color=TEMP:FOOD)) + 
  labs(title="% Survival by treatment and date larvae released", y=("Percent Survival"), x=("Date Released")) + 
  ylim(0,60) + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), 
                     name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", 
                                                      "Warm / High Food", "Warm / Low Food"))
dev.off()

# ------------ plot larval survival by no. of larvae collected 
png(file="results-figures/larval-surv-number-col.png", width = 600, height = 325)
ggplot(master, aes(x=Live.Larvae, y=100*(Live.35.days/800))) + theme_bw(base_size = 14) + 
  geom_point(size=2.5, aes(color=TEMP:FOOD)) + ylim(0,60) + 
  labs(title="% Survival by treatment and no. larvae released in group", 
                                                    y=("Percent Survival"), x=("No. larvae released")) + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), 
                     name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", 
                                                      "Warm / High Food", "Warm / Low Food"))
dev.off()


# ---------- plot larval survival by variance in larval length within a sample 
png(file="results-figures/larval-surv-time-col.png", width = 600, height = 325)
ggplot(master, aes(x=Date.stocked, y=100*(Live.35.days/800))) + theme_bw(base_size = 14) + 
  geom_point(size=2.5, aes(color=TEMP:FOOD)) + 
  labs(title="% Survival by treatment and date larvae released", y=("Percent Survival"), x=("Date Released")) + 
  ylim(0,60) + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), 
                     name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", 
                                                      "Warm / High Food", "Warm / Low Food"))
dev.off()

# ------------ plot larval survival by no. of larvae collected 
png(file="results-figures/larval-surv-number-col.png", width = 600, height = 325)
ggplot(master, aes(x=Live.Larvae, y=100*(Live.35.days/800))) + theme_bw(base_size = 14) + 
  geom_point(size=2.5, aes(color=TEMP:FOOD)) + ylim(0,60) + 
  labs(title="% Survival by treatment and no. larvae released in group", 
       y=("Percent Survival"), x=("No. larvae released")) + 
  scale_color_manual(values=c('#0571b0','#92c5de','#ca0020','#f4a582'), 
                     name=element_blank(), labels = c("Cold / High Food", "Cold / Low Food", 
                                                      "Warm / High Food", "Warm / Low Food"))
dev.off()


