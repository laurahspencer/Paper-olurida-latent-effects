larvalsize.path <- here::here("data", "larval-size.xlsx")
larvalsize <- larvalsize.path %>%
    excel_sheets() %>%
    set_names() %>%
    purrr::map_df(~ read_excel(path = larvalsize.path, sheet = .x), .id = "sheet") %>%
    dplyr::mutate(Sample.number = as.numeric(sheet)) %>% 
    inner_join(
        readRDS(file=here::here("data/larvae-collection-data.rds")) %>% 
            as_tibble() %>%
            dplyr::select(Rep, Date, Bag, Group, Sample.number, Live.Larvae, TEMP, FOOD, TREAT) %>%
            dplyr::mutate(Sample.number = as.numeric(Sample.number)))

# How many larvae were measured per sample? 
ngroups <- larvalsize %>% group_by(sheet) %>% tally()

# Create a dataframe with average shell measurements by sample number (aka larval batch)
average.all <- larvalsize %>% 
    dplyr::select(sheet, MaxFeret, MinFeret, Length, Width, Circularity, Sample.number, Date, TEMP, FOOD, Group, Live.Larvae) %>%
    mutate(sheet=as.numeric(larvalsize$sheet)) %>%
    group_by(sheet) %>%
    mutate_at(c("MaxFeret", "MinFeret", "Circularity", "Length", "Width", "Live.Larvae"), mean, na.rm=TRUE) %>%
    distinct() %>%
    mutate(ethanol=ifelse(sheet < 50, 'YES', 'NO')) %>%
    drop_na("sheet")

# Save mean shell length for each larval group to later merge with larval survival data 
save(average.all, file = "data/larval.size")

#How many larval groups counted per treatment? 
average.all %>% group_by(TEMP, FOOD) %>% tally()

# Plot shell width (aka length in data) by treatment, and by ethanol preservation 
# Note: I killed larvae using ethanol through sample 50, then stopped. This significantly altered the integrity of the larval samples (ethanol = well preserved), which altered the general size measurement. I therefore will include "ethanol" as a random variable in my models when I test effect of parental treatments 

average.all %>% 
    group_by(FOOD, TEMP) %>% 
    ggplot(aes(x=FOOD:TEMP, y=Length)) + 
    geom_boxplot(aes(col=FOOD:TEMP)) + 
    ggtitle("Length by treatment (and preservation method)") + 
    theme_minimal() +
    geom_jitter(width=0.2, aes(col=FOOD:TEMP)) +
    facet_wrap(~ethanol)

hist(average.all$Length)
shapiro.test(average.all$Length) #Data is normal 

bartlett.test(average.all$Length ~ average.all$FOOD)
bartlett.test(average.all$Length ~ average.all$TEMP) # variance not homogenous between temperatures
leveneTest(Length~FOOD*TEMP, data=average.all)  # Use Levene test to compare across all 4 treatments. Variance doesn't differ. 

# Balanced design? No. 
table(average.all$TEMP:average.all$FOOD) #how many larval groups per treatment?
sum(as.data.frame(table(average.all$TEMP:average.all$FOOD))$Freq) #how many larval groups total? 

# Include the "ethanol" random factor, since that presrevation method resulted in varying size
# Do type II ANOVA, since unbalanced design. and unequal variances. No diff. 
Anova(lme(Length ~  FOOD*TEMP, random=~1|ethanol,data=average.all), type="II", white.adjust=TRUE) # <-------------

# Posthoc on mixed model
emmeans(lme(Length ~  FOOD*TEMP, random=~1|ethanol,data=average.all), pairwise ~ FOOD+TEMP)

# Look at effect directly 
summary(lme(Length ~  FOOD*TEMP, random=~1|ethanol,data=average.all))

# Effect of ethanol preservation: -5.620869 
# Create new column that corrects data for differences in ethanol preservation 
average.all <- mutate(average.all, Length.corr=if (ethanol=="NO"){Length-5.620869}
                      else if (ethanol=="YES"){Length})

# Calculate average size by parental winter treatment - to add to boxplots 
length.mean <- aggregate(Length ~ FOOD:TEMP, data=average.all, mean)

# Relevel food factors such that Low comes before High
average.all$FOOD <- factor(average.all$FOOD, levels = rev(levels(average.all$FOOD)))

# ===================
# Larval Size Figure
# ===================

pdf(file=here::here("results-figures","larval-shell-length.pdf"), width = 4, height = 7)
average.all %>%
    ggplot(aes(x=TEMP:FOOD, y=Length)) + 
    geom_boxplot(lwd=0.5, aes(col=TEMP:FOOD)) + 
    theme_bw(base_size = 12) + 
    labs(title="Larval shell width upon release", y=("shell width (µm)"), x="treatment") + 
    scale_color_manual(values=c('#92c5de','#0571b0','#f4a582','#ca0020'), name=element_blank(), 
                       labels = c("7°C+Low-food","7°C+High-food","10°C+Low-food","10°C+High-food")) + 
    theme(text = element_text(size = 12.5)) + 
    theme(axis.title.x=element_blank(), 
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none") + 
    ylim(c(164.5, 207)) + #don't make ylim higher than 164.5, it omits a data point 
    geom_jitter(width=0.2, aes(col=TEMP:FOOD)) +
    geom_point(data=length.mean, aes(y=Length, fill=TEMP:FOOD), col="gray40", shape=22, size=4) +
    scale_fill_manual(values=c('#92c5de','#0571b0','#f4a582','#ca0020'))
dev.off()