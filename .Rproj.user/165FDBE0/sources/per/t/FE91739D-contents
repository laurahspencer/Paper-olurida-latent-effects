
# ==============================
# RIPE OOCYTE SIZE 
# ==============================

oocyte <- read.csv("data/Oocyte-length.csv", header=T, stringsAsFactors = T)
#oocyte.meanlength <- aggregate(Length ~ Sample+TEMP+FOOD+Week+TREAT+TREAT.NAME, oocyte, mean)

oocyte.size <- oocyte %>%
    group_by(Sample, Week, TEMP, FOOD, TREAT.NAME) %>%
    summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE)) %>%
    filter(TEMP!="Wild") %>% droplevels()
oocyte.size$cv <- oocyte.size$sd_length/oocyte.size$mean_length 

# check out ANOVA assumptions 
hist(oocyte.size$mean_length)
lm.oocyte <- lm(mean_length ~ TEMP*FOOD, data=oocyte.size)
ggqqplot(residuals(lm.oocyte)) #looks okay 
ggqqplot(oocyte.size, "mean_length",ggtheme=theme_bw()) + facet_grid(TEMP~FOOD) #looks good 
shapiro_test(residuals(lm.oocyte)) # not normal. try transforming. 

hist(oocyte.size$mean_length^2)
lm.oocyte2 <- lm(mean_length^2 ~ TEMP*FOOD, data=oocyte.size)
ggqqplot(residuals(lm.oocyte2)) #looks okay 
shapiro_test(residuals(lm.oocyte2)) # good  
leveneTest(lm.oocyte2) #assess variance 
plot(lm.oocyte2, which=1) 

# Compare oocyte length by treatment, using mean size per individual 
Anova(lm(mean_length^2 ~ TEMP*FOOD, data=oocyte.size), type = "II") #Use Car package with type II since unbalanced design 
summary(lm(mean_length^2 ~ TEMP*FOOD, data=oocyte.size)) #check intercept and coefficients 



# ===============
# OOCYTE FIGURE 
# ===============

# Relevel food factors such that Low comes before High
oocyte.size$FOOD <- factor(oocyte.size$FOOD, levels = rev(levels(oocyte.size$FOOD)))

# Calculate average oocyte size by parental winter treatment - to add to boxplots 
oocyte.length.mean <- aggregate(mean_length ~ FOOD:TEMP, data=oocyte.size, mean)

# FIGURE FOR PAPER! As of March 1, 2021 
pdf(file="results-figures/Stage3-oocyte-size.pdf", width=5.5, height=6)
ggplot(data=oocyte.size, 
       aes(x=TEMP:FOOD, y=mean_length, col=TEMP:FOOD)) + ylim(c(66, 109)) +
    geom_boxplot(lwd=0.5) + theme_bw(base_size = 11) + 
    ggtitle(label = "Ripe oocyte size") + ylab("Maximum oocyte length (µm)") + 
    scale_color_manual(values=c('#92c5de','#0571b0','#f4a582','#ca0020'), name=element_blank(),
                       labels = c( "7°C+Low-food", "7°C+High-food","10°C+Low-food","10°C+High-food")) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), legend.position = "left") + 
    geom_jitter(width = .2,aes(shape=as.factor(Week)),size=2) + 
    scale_shape_manual(values=c(15,16,17), name="Date oocyte sampled", labels=c("February 27", "March 13", "March 23")) +
    geom_point(data=oocyte.length.mean, aes(y=mean_length, fill=TEMP:FOOD), col="gray40", shape=22, size=4, show.legend=F) +
    scale_fill_manual(values=c('#92c5de','#0571b0','#f4a582','#ca0020'))
dev.off()

