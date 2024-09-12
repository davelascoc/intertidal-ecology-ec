                    ###  <<< R FOR MARINE BIOLOGY >>>  ###

         ###  <<< INTERTIDAL COMMUNITIES IN THE ECUADORIAN COAST>>>  ###


## Fundamentos de R - Daniel Velasco
## 12/09/2024


#data from the class of Techniques of Marine Research I

#load packages
library(readxl)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(ggtext)
library(vegan)
library(indicspecies)
library(RVAideMemoire)
library(rstatix)



# << COMPOSITION >> ##

# DATA MANAGEMENT

##import dataset
i           <-  read_excel("Datos/Datos Analizados/Master Database VFinal.xlsx", 
                               sheet = "Sessile")

##factors
i           <-  i   %>% mutate_at(c("Sector", "Code1", "Code2",
                                    "Initials"), factor)          %>%
                        mutate(Season = factor(Season,
                                            levels = c("Warm", "Cold")),
                        Zone = factor(Zone, levels = c("H", "M", "L")))

i_2023      <- i    %>% filter(Year == 2023, Season == "Warm", Zone == "L") %>%
                        mutate(Site = factor(Site,
                                      levels = c("Ballenita", "Machalilla",
                                                 "La Playita", "La Rinconada")),
                               Code2 = factor(Code2, levels =c(
                                 "BAN", "BAS", "MAN", "MAS",
                                 "LPN", "LPS", "LRN", "LRS")))

##separate data in sample information and community
saminfo     <- i_2023 %>%  select( 1:15) 
cominfo     <- i_2023 %>%  select(!1:15)

##NAs -> 0
cominfo[is.na(cominfo)]  <-  0

##summary
summary(saminfo)


##Hellinger: standardization method for community ecologists
hell          <-  decostand(cominfo, method = "hellinger")

##dissimilarity indices
diss          <-  vegdist(  hell, method = "bray", binary = F) # Bray-Curtis


# NON-METRIC MULTIDIMENSIONAL SCALING (NMDS)

##run mds
set.seed(200)
NMDS          <-  metaMDS(diss, k = 2, try = 200, trymax = 1000,
                            autotransform = F)
NMDS
NMDS$stress
saminfo[,16:17]   <- as.data.frame(scores(NMDS, display = "sites"))

stress  <-  data.frame(stress = NMDS$stress)
stress

saminfo[,16:17]  <-  data.frame(NMDS$points)

# PERMANOVA (Permutational multivariate analysis of variance)

##run permanova
permanova     <-  adonis2(diss ~ saminfo$Site)
permanova

pw            <-  pairwise.perm.manova(diss, saminfo$Code1)
pw

pv            <-  permanova$`Pr(>F)`[1]
r2            <-  permanova$R2[1]
Fv            <-  permanova$F[1]

##NMDS plotting
NMDSplot     <-   ggplot(saminfo, aes(x = -NMDS1, y = NMDS2, col = Site)) +
                  geom_point(aes(shape = Site), size = 2.7) +
                  annotate(geom = "text", family = "serif", x = -Inf, y = Inf,
                  hjust = -.1, vjust = 1.5, label = paste("2D Stress =",
                          format(round(NMDS$stress, 2), nsmall = 2))) +
                  stat_ellipse(type = "t", linetype = 1, level = 0.95,
                               show.legend = F) +
                  stat_ellipse(type = "norm", linetype = 3, level = 0.95,
                               show.legend = F) +
                  labs(x = "MDS1", y = "MDS2", col = "Site", shape = "Site",
                  subtitle = paste("PERMANOVA, F =", round(Fv, 3),
                                   ", R² =", round(r2, 3),
                                   ", p-value =", round(pv, 3)),
                  caption = "Low Intertidal Zone Composition - 2023") +
                  theme_test() +
                  scale_shape_manual(values = c(17, 15, 16, 8)) +
                  theme(axis.text = element_blank(),
                      axis.ticks = element_blank(),
                      text = element_text(family = "serif", size = 12)) +
                  scale_color_hue(l = 50)
NMDSplot

#ggsave("R_outputs/NMDS_2023_L.png", width = 7, height = 4, dpi = 1000)


##playing with ggplot
ggplot(saminfo, aes(x = -NMDS1, y = NMDS2, color = Site)) +
                    geom_point(aes(shape = Sector), size = 2.7) +
                    annotate(geom = "text", family = "serif", x = -Inf, y = Inf,
                             hjust = -.1, vjust = 1.5,
                             label = paste("2D Stress =",
                         format(round(NMDS$stress, 2), nsmall = 2))) +
                    stat_ellipse(type = "t", linetype = 1, level = 0.95,
                                 show.legend = F) +
                    stat_ellipse(type = "norm", linetype = 3, level = 0.95,
                                 show.legend = F) +
                    labs(x = "MDS1", y = "MDS2", col = "Site", shape = "Site",
                    subtitle = paste("PERMANOVA, F =", round(Fv, 3),
                                     ", R² =", round(r2, 3),
                                     ", p-value =", round(pv, 3)),
                    caption = "Low Intertidal Zone Composition - 2023") +
                    scale_shape_manual(values = c(17, 6)) +
                    theme_test() + theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    text = element_text(family = "serif", size = 12)) +
                    scale_color_hue(l = 50)

ggplot(saminfo, aes(x = -NMDS1, y = NMDS2, group = Code2, color = Code2)) +
                    geom_point(aes(shape = Code2), size = 2.7) +
                    annotate(geom = "text", family = "serif", x = -Inf, y = Inf,
                             hjust = -.1, vjust = 1.5,
                             label = paste("2D Stress =",
                         format(round(NMDS$stress, 2), nsmall = 2))) +
                    stat_ellipse(type = "t", linetype = 1, level = 0.95,
                                 show.legend = F) +
                    labs(x = "MDS1", y = "MDS2", col = "Site", shape = "Site",
                    subtitle = paste("PERMANOVA, F =", round(Fv, 3),
                                     ", R² =", round(r2, 3),
                                     ", p-value =", round(pv, 3)),
                    caption = "Low Zone Composition - 2023") +            
                    scale_shape_manual(values = c(17, 2, 15, 0, 16, 1, 18, 5)) +
                    theme_test() + theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    text = element_text(family = "serif", size = 12)) +
                    scale_color_hue(l = 50)


# HOMOGENEITY OF MULTIVARIATE DISPERSION
#------------------------------------------------------------------------------
dishom  <- betadisper(vegdist(cominfo, "bray", binary = T), saminfo$Site)
                     # NO homogeneity if p < 0.01
dishom

anova(dishom)
TukeyHSD(dishom)

saminfo$Distance <- dishom[["distances"]]
 
displot          <- ggplot(saminfo, aes(x = Site, y = Distance)) +
                    stat_boxplot(geom = "errorbar", width = 0.1, color = 1) +
                    geom_boxplot(fill = "white", col = 1, outlier.shape = NA) +
                    geom_jitter(aes(shape = Site, col = Site),
                    width = 0.25, size = 2) + xlab("") +
                    ylab("Distance to Centroid") +
                    stat_summary(fun = 'mean', geom = 'point', shape = 4) +
                  #  scale_shape_manual(values = c(17, 6)) +
                    scale_shape_manual(values = c(17, 15, 16, 8)) +
                    theme_test() +
                    theme(legend.position = "none",
                          text = element_text(family = "serif", size = 12)) +
                    scale_color_hue(l = 50)
displot

ggplot(saminfo, aes(x = Sector, y = Distance, color = Site)) +
                    stat_boxplot(geom = "errorbar", width = 0.15, color = 1) +
                    geom_boxplot(fill = "white", col = 1, outlier.shape = NA) +
                    geom_jitter(aes(shape = Sector), width = .25, size = 2) +
                    xlab("") + ylab("Distance to Centroid") +
                    stat_summary(fun = 'mean', geom = 'point', shape = 4) +
                    scale_shape_manual(values = c(17, 6)) +
                    facet_wrap(~Site, ncol = 4) + theme_test() +
                    theme(legend.position = "none",
                          text = element_text(family = "serif", size = 12)) +
                    scale_color_hue(l = 50)
#------------------------------------------------------------------------------

# OTHER STATS #
#------------------------------------------------------------------------------
#ANOSIM - Analysis of Similarities

        ## Test estadístico no paramétrico muy empleado para el análisis de
        ## comunidades. No tiene el requisito de la normalidad
        ## (podemos tener muchos ceros en nuestra matriz).

anosim(cominfo, saminfo$Site, distance = "bray", permutations = 999)

        ## https://stat.ethz.ch/pipermail/r-sig-ecology/2013-June/003864.html
        ## Pairwise tests are not possible in vegan.

#Indicator species analysis

        ## que especies difieren entre los tramos

summary(multipatt(cominfo, saminfo$Site, func = "r.g",
                  control = how(nperm = 9999)), alpha = 0.05)

#SIMPER - Similarity Percentages

simper(comm = cominfo, group = saminfo$Site, permutations = 999, ordered = T)
#------------------------------------------------------------------------------

rm(cominfo, dishom, displot, hell, NMDS, NMDSplot, permanova, pw, saminfo,
   stress, diss, Fv, pv, r2)



## << RICHNESS >> ##

# DATA MANAGEMENT

##ordenar por media de riqueza
rich_mean        <- i                      %>%
                    group_by(Site, Code1)  %>%
                    summarise(Mean = mean(Richness))  %>%  arrange(Mean)
rich_mean

i                <- i               %>%
                    select(1:15)    %>%
                    mutate(Site  = factor(Site,  levels = rich_mean$Site),
                           Code1 = factor(Code1, levels = rich_mean$Code1))

summary(i)
summary(i$Site)
by(i$Richness, i$Site, mean)

# IS OUR DATA PARAMETRIC?
#------------------------------------------------------------------------------
#identify outliers
i %>%       group_by(Site, Zone)               %>%  select(Richness)    %>%
            identify_outliers(Richness)    #    %>%  filter(is.extreme == T)

#residuals
normal    <- lm(Richness  ~ Site, i)

#qq-plot
ggqqplot(residuals(normal)) + theme_test()

#by sites
ggqqplot(i, "Richness", facet.by = "Site",  size = 1) + theme_test()
ggqqplot(i, "Richness", facet.by = "Code2", size = 1) + theme_test()

##normality Test: Shapiro-Wilk
shapiro_test(residuals(normal))     # no hay normalidad

#by sites
i  %>%  group_by(Site)  %>% shapiro_test(Richness)
i  %>%  group_by(Code2) %>% shapiro_test(Richness)


# Data distribution
ggplot(i, aes(x = Richness, fill = Zone)) + geom_histogram(colour = 1,
                    binwidth = 1) + facet_grid(Zone ~ Site) + theme_bw() +
                    theme(legend.position = "none") + ylab("Frequency")

ggplot(i, aes(x = Richness, fill = Zone)) + ylab("Frequency") +
                    geom_histogram(aes(y = ..density..), colour = 1,
                    binwidth = 1) + facet_grid(Zone ~ Site) + theme_bw() + 
                    geom_density(lwd = 1, colour = "black", alpha = 0.25)  +
                    theme(legend.position = "none")

# Transform Data

i$logRichness <- log(i$Richness + 1)

shapiro_test(residuals(lm(Richness     ~ Site, i)))
shapiro_test(residuals(lm(logRichness  ~ Site, i)))

ggplot(i, aes(x = logRichness, fill = Zone)) + ylab("Frequency") +
                    geom_histogram(aes(y = ..density..), colour = 1,
                    binwidth = .25) + facet_grid(Zone ~ Site) + theme_bw() +
                    geom_density(lwd = 1, colour = "black", alpha = 0.25) +
                    theme(legend.position = "none")

#homogeneity of variance
plot(normal, 1)

#Bartlett's test (more sensitive to deviations from normality)
bartlett.test(Richness ~ Site, i)

#Levene's test (less sensitive)
i  %>%  levene_test(Richness  ~ Site)
#------------------------------------------------------------------------------

# BOXPLOTS

#ggplot
ggplot(i, aes(x = Site, y = Richness, color = Site)) +
                    geom_jitter(width = .25, size = .6) +
                    geom_boxplot(fill = NA, col = 1, outlier.shape = NA) +
                    xlab("") + facet_wrap(~ Zone) + theme_test() +
                    ggtitle("Intertidal Zone Richness") +
                    stat_summary(fun = 'mean', geom = 'point',
                                 shape = 4, col = 1) +
                    stat_compare_means(method = "kruskal.test",
                    label = "p.signif", label.x = Inf,
                    label.y = Inf, hide.ns = F, hjust = 1.5, vjust = 1.5) +
                    theme(legend.position = "none",
                    text = element_text(family = "serif", size = 12),
                    axis.text.x = element_text(angle = 90,
                    vjust = 0.5, hjust=1),
                    plot.title = element_text(hjust = 0.5)) +
                    scale_color_frontiers()

#ggsave("R_outputs/richness_coast.png", width = 7, height = 4, dpi = 1200)

ggplot(i, aes(x = Zone, y = Richness, color = Zone)) +
                    geom_jitter(width = 0.25, size = 1) + xlab("") +
                    geom_boxplot(fill = NA, color = "black",
                    outlier.shape = NA) + stat_boxplot(geom = "errorbar",
                    width = 0.15, color = "black") + facet_wrap(~ Code1,
                    nrow = 1) + theme_test() + theme(legend.position = "none") 


#Kruskal-Wallis test
kruskal_test(Richness ~ Site, data = i)  #  entre sitios
kruskal_test(Richness ~ Zone, data = i)  #  entre zonas

##Dunn’s test (post hoc)
dunn_test(Richness ~ Site, data = i)
dunn_test(Richness ~ Zone, data = i)


#https://www.youtube.com/watch?v=pWs-mlCmOi8
#https://r-charts.com/es/r-base/titulo/

