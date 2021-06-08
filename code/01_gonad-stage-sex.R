
# HISTOLOGY DATA
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


# How many samples total? How many were each sex?
histo %>% count() #total number sampled 
histo %>% count(SEX) %>% mutate(percent=n/(sum(n)))

# This is how many were some flavor of hermaph.: 
histo %>% count(SEX) %>% mutate(percent=n/(sum(n))) %>%
  filter(SEX %in% c("HPM", "H", "HPF")) %>% dplyr::select(percent) %>% sum()


# ========================================
# GONAD DEVELOPMENT (STAGE AND SEX RATIO)
# ========================================

# ============= All treatment weeks prior to reproductive conditioning combined (3-12):

print(sex.all <- table(subset(histo, Week!="0" & Week!="15" & Week!="16")$TREAT, subset(histo, Week!="0" & Week!="15" & Week!="16")$SEX))
chisq.test(sex.all[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.all.col <- table(subset(histo, Week!="0" & Week!="15" & Week!="16")$TREAT, subset(histo, Week!="0" & Week!="15" & Week!="16")$SEX.COL))
chisq.test(sex.all.col[2:5,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.all <- table(subset(histo, Week!="0" & Week!="15" & Week!="16")$TREAT, subset(histo, Week!="0" & Week!="15" & Week!="16")$MALSTAGE.COL))
chisq.test(malstage.all[2:5,], simulate.p.value = T, B = 10000)  #different, p=0.009899 How?
pairwiseNominalIndependence(malstage.all[2:5,],fisher = TRUE,gtest  = FALSE, chisq  = FALSE, digits = 3)

# Female gonad stage 
print(femalstage.all <- table(subset(histo, Week!="0" & Week!="15" & Week!="16")$TREAT, subset(histo, Week!="0" & Week!="15" & Week!="16")$FEMSTAGE.COL))
chisq.test(femalstage.all[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff


# ============= Week 12/13 (end of treatment)

print(sex.wk13 <- table(subset(histo, Week=="13")$TREAT, subset(histo, Week=="13")$SEX))
chisq.test(sex.wk13[2:5,], simulate.p.value = T, B = 10000)  #can't use chi-sq
fisher.test(sex.wk13[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk13.col <- table(subset(histo, Week=="13")$TREAT, subset(histo, Week=="13")$SEX.COL))
chisq.test(sex.wk13.col[2:5,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk13 <- table(subset(histo, Week=="13")$TREAT, subset(histo, Week=="13")$MALSTAGE.COL))
chisq.test(malstage.wk13[2:5,], simulate.p.value = T, B = 10000)  # diff! 
pairwiseNominalIndependence(malstage.wk13[2:5,],fisher = TRUE,gtest  = FALSE, chisq  = FALSE, digits = 3)
malstage.wk13[2:5,] %>% sum() #how many sampled? 

# Female gonad stage 
print(femalstage.wk13 <- table(subset(histo, Week=="13")$TREAT, subset(histo, Week=="13")$FEMSTAGE.COL))
chisq.test(femalstage.wk13[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff



# ============== Week 15 (first reproductive conditioning sampling, 2 weeks in conditioning)

print(sex.wk15 <- table(subset(histo, Week=="15")$TREAT, subset(histo, Week=="15")$SEX))
chisq.test(sex.wk15[2:5,], simulate.p.value = T, B = 10000)  #cant use

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk15.col <- table(subset(histo, Week=="15")$TREAT, subset(histo, Week=="15")$SEX.COL))
chisq.test(sex.wk15.col[2:5,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk15 <- table(subset(histo, Week=="15")$TREAT, subset(histo, Week=="15")$MALSTAGE.COL))
chisq.test(malstage.wk15[2:5,], simulate.p.value = T, B = 10000)  #Yes diff. p-value = 0.0189
malstage.wk15[2:5,] %>% sum() #how many sampled?

# Female gonad stage 
print(femalstage.wk15 <- table(subset(histo, Week=="15")$TREAT, subset(histo, Week=="15")$FEMSTAGE.COL))
chisq.test(femalstage.wk15[2:5,-5], simulate.p.value = T, B = 10000)  #no sign. diff



# ============= Week 16 (last sampling, 3rd week of reproductive conditioning)

print(sex.wk16 <- table(subset(histo, Week=="16")$TREAT, subset(histo, Week=="16")$SEX))
chisq.test(sex.wk16[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff

# Collapsed sex designations (HPF = F, etc.)
print(sex.wk16.col <- table(subset(histo, Week=="16")$TREAT, subset(histo, Week=="16")$SEX.COL))
chisq.test(sex.wk16.col[2:5,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wk16 <- table(subset(histo, Week=="16")$TREAT, subset(histo, Week=="16")$MALSTAGE.COL))
chisq.test(malstage.wk16[2:5,], simulate.p.value = T, B = 10000)  #Yes diff. p-value = 0.008699
pairwiseNominalIndependence(malstage.wk16[2:5,],fisher = TRUE,gtest  = FALSE, chisq  = FALSE, digits = 3)
malstage.wk16[2:5,] %>% sum() # how many sampled?

# Female gonad stage 
print(femalstage.wk16 <- table(subset(histo, Week=="16")$TREAT, subset(histo, Week=="16")$FEMSTAGE.COL))
chisq.test(femalstage.wk16[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff



# ============= Weeks 15 & 16 (2 conditioning weeks) 

print(sex.wkcond <- table(subset(histo, Week=="15" | Week=="16")$TREAT, subset(histo, Week=="15" | Week=="16")$SEX))
chisq.test(sex.wkcond[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff
 
# Collapsed sex designations (HPF = F, etc.)
print(sex.wkcond.col <- table(subset(histo, Week=="15" | Week=="16")$TREAT, subset(histo, Week=="15" | Week=="16")$SEX.COL))
chisq.test(sex.wkcond.col[2:5,-1], simulate.p.value = T, B = 10000)  #no sign. diff

# Male gonad stage 
print(malstage.wkcond <- table(subset(histo, Week=="15" | Week=="16")$TREAT, subset(histo, Week=="15" | Week=="16")$MALSTAGE.COL))
chisq.test(malstage.wkcond[2:5,], simulate.p.value = T, B = 10000)  #Yes diff p-value = 0.0014
pairwiseNominalIndependence(malstage.wkcond[2:5,],fisher = TRUE,gtest  = FALSE, chisq  = FALSE, digits = 3)
malstage.wkcond[2:5,] %>% sum() #how many sampled? 
(6+8)/(2+4+6+8+4) # % of D treatment (Warm/High) with late-stage sperm (2 or 3)
(5+0)/(15+3+5+1) # % of C treatment (Cold/High) with late-stage sperm (2 or 3)

2/(2+4+6+8+4) # % of D treatment (Warm/High) that lacked sperm 
15/(15+3+5+1) # % of C treatment (Cold/High) that lacked sperm 

# Female gonad stage 
print(femalstage.wkcond <- table(subset(histo, Week=="15" | Week=="16")$TREAT, subset(histo, Week=="15" | Week=="16")$FEMSTAGE.COL))
chisq.test(femalstage.wkcond[2:5,], simulate.p.value = T, B = 10000)  #no sign. diff



# ===========================
# SPERMATOCYTE STAGE FIGURE 
# ===========================

treats <- as.character(unique(subset(histo, TREAT!="PRE")$TREAT))

# For both male and female gamete stage, a column that combines stages 2 and 3
histo <- histo %>% 
  mutate(MALSTAGE.COL2 = droplevels(replace(MALSTAGE.COL, MALSTAGE.COL=="2", "3")), 
         FEMSTAGE.COL2 = droplevels(replace(FEMSTAGE.COL, FEMSTAGE.COL=="2", "3"))) #%>%

# Calculate # of samples in each category, and plot 
malestage.treats <- list()
for (i in 1:length(treats)) {
  malestage.treats[[i]] <- as.data.frame(table(subset(histo, c(TREAT==treats[i] | TREAT=="PRE"))$Week, subset(histo, c(TREAT==treats[i] | TREAT=="PRE"))$MALSTAGE.COL2), margin=1)
  malestage.treats[[i]]["Var2"] <- droplevels(malestage.treats[[i]]["Var2"])
}
mains=c("Cold/Low", "Warm/Low", "Cold/High", "Warm/High")

# Generate figures 
p <- list()
for (i in 1:3){
  p[[i]] <- ggplot(data=subset(malestage.treats[[i]], Var1!="0" & Var2!="0"), 
                   aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + 
    theme_bw(base_size = 12) + geom_point(size=2) + ylim(c(0,10)) + 
    geom_vline(xintercept = 5.1, linetype="dashed", color = "gray50", size=.5) + 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
          legend.position = "none", 
          axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
    scale_color_manual(values=c("gray45", "gray10", "gray75")) +
    ggtitle(mains[i])
}
p[[4]] <- ggplot(data=subset(malestage.treats[[4]], Var1!="0" & Var2!="0"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + 
  theme_bw(base_size = 12) + geom_point(size=2) + ylim(c(0,10)) + 
  geom_vline(xintercept = 5.1, linetype="dashed", color = "gray50", size=.5) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), 
        legend.position = "bottom", 
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank(), 
        legend.title=element_text(size=14), legend.text = element_text(size=14))+ 
  guides(color = guide_legend(override.aes = list(size=6))) + 
  scale_x_discrete(labels= c("Dec 20", "Jan 4", "Jan 23", "Feb 9", "Feb 27", "Mar 13", "Mar 23")) + 
  scale_color_manual(values=c("gray45", "gray10", "gray75"),
                     name=element_blank(), labels = c("Early", "Advanced/Ripe", "Spawned")) +
  ggtitle(mains[4])

# Order of plots is: low/cold (A), low/warm (B), high/cold (C), high/warm (D)
pdf(file = "results-figures/malestage-line-plots.pdf", width = 4.5, height = 8)
do.call(grid.arrange, c(c(p[1], p[2], p[3], p[4]), list(ncol=1, nrow=4))) #would need to resize for pub. 
dev.off()



# ============================
# EXTRA PLOTS FOR SUPPLEMENTAL 
# ============================


# ===========================
# SEX RATIO FIGURE 
# ===========================

sexcol.treats <- list()
for (i in 1:length(treats)) {
  sexcol.treats[[i]] <- as.data.frame(table(subset(histo, TREAT==treats[i] | TREAT=="PRE" )$Week, subset(histo, TREAT==treats[i] | TREAT=="PRE" )$SEX.COL), margin=1)
}

p <- list()
for (i in 1:3){
  p[[i]] <- ggplot(data=subset(sexcol.treats[[i]], Var1!=0 & (Var2=="F" | Var2=="M" | Var2=="H")), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + 
    scale_y_continuous(position = "right", limits = c(0,10), breaks=c(0,5,10)) + 
    geom_line() + theme_bw(base_size = 12) + geom_point(size=2.5) + 
    geom_vline(xintercept = 5.1, linetype="dashed", color = "gray50", size=.5) + 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", 
          axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
    scale_color_manual(values=c("gray10", "gray45", "gray75")) +
    ggtitle(mains[i])
}

p[[4]] <- ggplot(data=subset(sexcol.treats[[4]], Var1!=0 & (Var2=="F" | Var2=="M" | Var2=="H")), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + 
  scale_y_continuous(position = "right", limits = c(0,10), breaks=c(0,5,10)) + 
  geom_line() + theme_bw(base_size = 12) + geom_point(size=2.5) + 
  geom_vline(xintercept = 5.1, linetype="dashed", color = "gray50", size=.5) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "bottom", 
        axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        legend.title=element_text(size=13), legend.text = element_text(size=13)) + 
  guides(color = guide_legend(override.aes = list(size=6))) + 
  scale_x_discrete(labels= c("Dec 20", "Jan 4", "Jan 23", "Feb 9", "Feb 27", "Mar 13", "Mar 23")) + 
  scale_color_manual(values=c("gray10", "gray45", "gray75"), 
                     name=element_blank(), labels = c("Female / HPF", "Male / HPM", "Hermaph."))+
  ggtitle(mains[4])

# Order of plots is: low/cold (A), low/warm (B), high/cold (C), high/warm (D)
pdf(file = "results-figures/gonad-sex-line-plots.pdf", width = 4.8, height = 8)
do.call(grid.arrange, c(c(p[1], p[2], p[3], p[4]), list(ncol=1, nrow=4))) # would need to adjust ratios for pub.
dev.off() 


# ===========================
# OOCYTE STAGE FIGURE 
# ===========================

# plot female stage over time 

femalestage.treats <- list()
for (i in 1:length(treats)) {
  femalestage.treats[[i]] <- as.data.frame(table(subset(histo, TREAT==treats[i] | TREAT=="PRE" )$Week, subset(histo, TREAT==treats[i] | TREAT=="PRE" )$FEMSTAGE.COL2), margin=1)
}

p <- list()

for (i in 1:3){
  p[[i]] <- ggplot(data=subset(femalestage.treats[[i]], Var1!="0" & Var2!="0"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + theme_bw(base_size = 12) + geom_point(size=2) + ylim(c(0,10)) + 
    geom_vline(xintercept = 5.1, linetype="dashed", color = "gray50", size=.5) + 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "none", 
          axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    scale_color_manual(values=c("gray45", "gray10", "gray75")) +
    ggtitle(mains[i]) 
}

p[[4]] <- ggplot(data=subset(femalestage.treats[[4]], Var1!="0" & Var2!="0"), aes(x=Var1, y=Freq, group=Var2, col=Var2)) + geom_line() + 
  theme_bw(base_size = 12) + geom_point(size=2) + ylim(c(0,10)) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position = "bottom", 
        axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),legend.title=element_text(size=14), legend.text = element_text(size=14)) + 
  guides(color = guide_legend(override.aes = list(size=6))) + 
  scale_x_discrete(labels= c("Dec 20", "Jan 4", "Jan 23", "Feb 9", "Feb 27", "Mar 13", "Mar 23")) + 
  scale_color_manual(values=c("gray45", "gray10", "gray75"),
                     name=element_blank(), labels = c("Early", "Advanced/Ripe", "Spawned")) +
  ggtitle(mains[4])

pdf(file = "results-figures/femalestage-line-plots.pdf", width = 4.5, height = 8)
do.call(grid.arrange, c(c(p[1], p[2], p[3], p[4]), list(ncol=1, nrow=4)))
dev.off()


# ====================================
# STACKED BARPLOTS FOR EACH TREATMENT 
# ====================================

# Prepare contingency tables 

# SEX
print(CT.sex.all <- table(histo$TREAT, histo$SEX))
print(CT.sex.Week <- table(histo$TREAT, histo$SEX, histo$Week))
print(CT.sex.FL <- table(subset(histo, Week=="0" | Week=="13")$TREAT, subset(histo, Week=="0" | Week=="13")$SEX))

# MALE STAGE
print(CT.MALSTAGE.all <- table(histo$TREAT, histo$MALSTAGE.COL))
print(CT.MALSTAGE.Week <- table(histo$TREAT, histo$MALSTAGE.COL, histo$Week))
print(CT.MALSTAGE.FL <- table(subset(histo, Week=="0" | Week=="13")$TREAT, subset(histo, Week=="0" | Week=="13")$MALSTAGE.COL))

# FEMALE STAGE
print(CT.fem.stage.all <- table(histo$TREAT, histo$FEMSTAGE.COL))
print(CT.fem.stage.Week <- table(histo$TREAT, histo$FEMSTAGE.COL, histo$Week))
print(CT.fem.stage.FL <- table(subset(histo, Week=="0" | Week=="13")$TREAT, subset(histo, Week=="0" | Week=="13")$FEMSTAGE.COL))

# CT for Treatment A only 
CT.A.sex <- table(subset(histo, TREAT=="A" | TREAT=="PRE")$Week, subset(histo, TREAT=="A" | TREAT=="PRE")$SEX)
CT.A.MALSTAGE <- table(subset(histo, TREAT=="A" | TREAT=="PRE")$Week, subset(histo, TREAT=="A" | TREAT=="PRE")$MALSTAGE.COL)
CT.A.fem.stage <- table(subset(histo, TREAT=="A" | TREAT=="PRE")$Week, subset(histo, TREAT=="A" | TREAT=="PRE")$FEMSTAGE.COL)

# CT for Treatment B only 
CT.B.sex <- table(subset(histo, TREAT=="B" | TREAT=="PRE")$Week, subset(histo, TREAT=="B" | TREAT=="PRE")$SEX)
CT.B.MALSTAGE <- table(subset(histo, TREAT=="B" | TREAT=="PRE")$Week, subset(histo, TREAT=="B" | TREAT=="PRE")$MALSTAGE.COL)
CT.B.fem.stage <- table(subset(histo, TREAT=="B" | TREAT=="PRE")$Week, subset(histo, TREAT=="B" | TREAT=="PRE")$FEMSTAGE.COL)

# CT for Treatment C only 
CT.C.sex <- table(subset(histo, TREAT=="C" | TREAT=="PRE")$Week, subset(histo, TREAT=="C" | TREAT=="PRE")$SEX)
CT.C.MALSTAGE <- table(subset(histo, TREAT=="C" | TREAT=="PRE")$Week, subset(histo, TREAT=="C" | TREAT=="PRE")$MALSTAGE.COL)
CT.C.fem.stage <- table(subset(histo, TREAT=="C" | TREAT=="PRE")$Week, subset(histo, TREAT=="C" | TREAT=="PRE")$FEMSTAGE.COL)

# CT for Treatment D only 
CT.D.sex <- table(subset(histo, TREAT=="D" | TREAT=="PRE")$Week, subset(histo, TREAT=="D" | TREAT=="PRE")$SEX)
CT.D.MALSTAGE <- table(subset(histo, TREAT=="D" | TREAT=="PRE")$Week, subset(histo, TREAT=="D" | TREAT=="PRE")$MALSTAGE.COL)
CT.D.fem.stage <- table(subset(histo, TREAT=="D" | TREAT=="PRE")$Week, subset(histo, TREAT=="D" | TREAT=="PRE")$FEMSTAGE.COL)


# Plot sex and male/female stage over time for each treatment

# ----------------------
# A = 7C low nutrition
# ----------------------

jpeg(file = "results-figures/gonad-barplots-7C-low.jpeg", width = 1015, height = 550)

par(mfrow = c(1, 3), mar=c(1, 2, 3, 0), oma=c(4,4,4,0), col="gray30")

print(barplot(t(prop.table(CT.A.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))

colnames(CT.A.MALSTAGE) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.A.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))

colnames(CT.A.fem.stage) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.A.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))

mtext("7°C - Low Nutrition", side=3, line = .5, outer = T, cex = 2, col = "gray30")
mtext("% Sampled", side=2, line = 1, outer = T, cex = 1.5, col = "gray30")
mtext("Week Sampled", side=1, line = 2.5, outer = T, cex = 1.5, col = "gray30")

dev.off()

# ----------------------
# C = 7C high nutrition
# ----------------------

jpeg(file = "results-figures/gonad-barplots-7C-high.jpeg", width = 1015, height = 550)

par(mfrow = c(1, 3), mar=c(1, 2, 3, 0), oma=c(4,4,4,0), col="gray30")

print(barplot(t(prop.table(CT.C.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
chisq.test(CT.C.sex, simulate.p.value = T, B = 10000) 

colnames(CT.C.MALSTAGE) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.C.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.C.MALSTAGE, simulate.p.value = T, B = 10000) 

colnames(CT.C.fem.stage) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.C.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.C.fem.stage, simulate.p.value = T, B = 10000) 

mtext("7°C - High Nutrition", side=3, line = .5, outer = T, cex = 2, col = "gray30")
mtext("% Sampled", side=2, line = 1, outer = T, cex = 1.5, col = "gray30")
mtext("Week Sampled", side=1, line = 2.5, outer = T, cex = 1.5, col = "gray30")

dev.off()

# ----------------------
# B = 10C low nutrition
# ----------------------

jpeg(file = "results-figures/gonad-barplots-10C-low.jpeg", width = 1015, height = 550)

par(mfrow = c(1, 3), mar=c(1, 2, 3, 0), oma=c(4,4,4,0), col="gray30")

print(barplot(t(prop.table(CT.B.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
chisq.test(CT.B.sex, simulate.p.value = T, B = 10000) 

colnames(CT.B.MALSTAGE) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.B.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.B.MALSTAGE, simulate.p.value = T, B = 10000) 

colnames(CT.B.fem.stage) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.B.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.B.fem.stage, simulate.p.value = T, B = 10000) 

mtext("10°C - Low Nutrition", side=3, line = .5, outer = T, cex = 2, col = "gray30")
mtext("% Sampled", side=2, line = 1, outer = T, cex = 1.5, col = "gray30")
mtext("Week Sampled", side=1, line = 2.5, outer = T, cex = 1.5, col = "gray30")

dev.off()

# ----------------------
# D = 10C high nutrition
# ----------------------

jpeg(file = "results-figures/gonad-barplots-10C-high.jpeg", width = 1015, height = 550)

par(mfrow = c(1, 3), mar=c(1, 2, 3, 0), oma=c(4,4,4,0), col="gray30")

print(barplot(t(prop.table(CT.D.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
chisq.test(CT.D.sex, simulate.p.value = T, B = 10000)  #sign.

colnames(CT.D.MALSTAGE) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.D.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.D.MALSTAGE, simulate.p.value = T, B = 10000) 

colnames(CT.D.fem.stage) <- c("None present (0)", "Early (1)", "Advanced (2)", "Ripe (3)", "Spawned (4)")
print(barplot(t(prop.table(CT.D.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.5, cex.names = 1.7, legend.text = F, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-.5, 0), title="Gonad Stage", cex=1.5)))
chisq.test(CT.D.fem.stage, simulate.p.value = T, B = 10000) 

mtext("10°C - High Nutrition", side=3, line = .5, outer = T, cex = 2, col = "gray30")
mtext("% Sampled", side=2, line = 1, outer = T, cex = 1.5, col = "gray30")
mtext("Week Sampled", side=1, line = 2.5, outer = T, cex = 1.5, col = "gray30")

dev.off()


# GENERATE PLOTS JUST TO USE THE LEGENDS 

pdf(file = "results-figures/gonad-barplots-sex-4legend.pdf", width = 8.5, height = 6.5)
par(mar=c(1, 2, 3, 20), col="gray30")
colnames(CT.D.sex) <- c("Undifferentiated", "Male", "Male dominant", "Hermaphroditic", "Female dominant", "Female")
print(barplot(t(prop.table(CT.D.sex, 1)), main="Dominant gonad sex", xlab="Week Sampled", ylab="% Sampled", las=1, col=c("gray75", "#08519c", "#3182bd", "purple3","mediumorchid3", "#df65b0"),  cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.3, cex.names = 1.3, legend.text = T, args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Gonad Sex", cex=1.5)))
dev.off()

pdf(file = "results-figures/gonad-barplots-male-stage-4legend.pdf", width = 8.5, height = 6.5)
par(mar=c(1, 2, 3, 20), col="gray30")
print(barplot(t(prop.table(CT.D.MALSTAGE, 1)), main="Male", xlab=NA, ylab=NA, las=1, col=c("#eff3ff","#6baed6","#3182bd","#08519c","#bdd7e7"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.3, cex.names = 1.3, legend.text = T, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Male Gonad Stage", cex=1.5)))
dev.off()

pdf(file = "results-figures/gonad-barplots-female-stage-4legend.pdf", width = 8.5, height = 6.5)
par(mar=c(1, 2, 3, 20), col="gray30")
print(barplot(t(prop.table(CT.D.fem.stage, 1)), main="Female", xlab=NA, ylab=NA, las=1, col=c("#f1eef6","#df65b0","#dd1c77","#980043","#d7b5d8"), cex.main=1.7, col.axis = 'gray30', col.main = 'gray30', cex.lab=1.5, cex.axis = 1.3, cex.names = 1.3, legend.text = T, yaxt='n', args.legend = list(x = "topright", bty = "n", inset=c(-1, 0), title="Female Gonad Stage", cex=1.5)))
dev.off()


