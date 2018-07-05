################################################
# Author: Ni Jiang (njiang@danforthcenter.org)
# Update on July 5, 2018
################################################



library("ggplot2")
library("ggpubr")
library("gridExtra")
library("LaplacesDemon")

####################################################################
########### Load traits from DynamicRoots outputs ##################
####################################################################
#read the list of files from a local directory

TraitFilesPath <- "E:\\time series\\traits_data\\data"
TraitFileNames <- list.files(path = TraitFilesPath, pattern = ".txt", full.names = TRUE)
n_plants <- length(TraitFileNames)
total_volume <- matrix(data = NA, nrow = n_plants, ncol = 42)
total_length <- matrix(data = NA, nrow = n_plants, ncol = 42)
total_depth <- matrix(data = NA, nrow = n_plants, ncol = 42)
branch_number <- matrix(data = NA, nrow = n_plants, ncol = 42)
volume_growth_rate <- matrix(data = NA, nrow = n_plants, ncol = 41)
length_growth_rate <- matrix(data = NA, nrow = n_plants, ncol = 41)
depth_growth_rate <- matrix(data = NA, nrow = n_plants, ncol = 41)
number_growth_rate <- matrix(data = NA, nrow = n_plants, ncol = 41)
mean_length_growth_rate <- matrix(data = NA, nrow = n_plants, ncol = 41)
mean_volume_growth_rate <- matrix(data = NA, nrow = n_plants, ncol = 41)
mean_depth_growth_rate <- matrix(data = NA, nrow = n_plants, ncol = 41)
mean_tortuosity <- matrix(data = NA, nrow = n_plants, ncol = 42)
mean_radius <- matrix(data = NA, nrow = n_plants, ncol = 42)
mean_volume <- matrix(data = NA, nrow = n_plants, ncol = 42)
mean_length <- matrix(data = NA, nrow = n_plants, ncol = 42)
tS <- timeSequence(from = "4 06:00", format = "%d %H:%M", by = "4 hours", length.out = 42)
tS <- paste("day",format(tS, "%d %H:%M"), sep="")
FileNames <- c()

for ( i in 1:n_plants )
{
  data <- read.csv(TraitFileNames[i], header = FALSE, sep = "\t")
  data <- data[-ncol(data)]
  
  fnamex <- strsplit(TraitFileNames[i], "\\.")[[1]][1];
  FileNames[i] <- substring(sub(".*TIM", "", fnamex), 1, 5)
  
  traits <- data[, seq(2, ncol(data), 2)]
  TraitNames <- data[1, seq(1, ncol(data), 2)]
  colnames(traits) <- as.matrix(TraitNames)
  
  n_timepoints <- (ncol(traits)-5) / 11
  
  volume <- traits[, seq(6, ncol(traits), 11)]
  length <- traits[, seq(7, ncol(traits), 11)]
  depth <- traits[, seq(9, ncol(traits), 11)]
  tortuosity <- traits[, seq(11, ncol(traits), 11)]
  radius <- traits[, seq(12, ncol(traits), 11)]
  ### correction for Av_radius 1.#J
  for (j in 1:ncol(radius))
  {
    if(is.factor(radius[, j]))
      radius[, j] <- as.numeric(levels(radius[, j])[radius[, j]])
  }
  
  # remove short branches
  volume[length<=5] <- NA
  depth[length<=5] <- NA
  length[length<=5] <- NA
  tortuosity[length<=5] <- NA
  radius[length<=5] <- NA
  
  volume1 <- volume
  length1 <- length
  depth1 <- depth
  
  growth_rate <- volume1[, -1] - volume1[, -n_timepoints]
  growth_rate <- data.frame(rep(NA, nrow(growth_rate)), growth_rate)  
  length1[growth_rate<=0] <- NA
  volume1[growth_rate<=0] <- NA
  depth1[growth_rate<=0] <- NA
  growth_rate <- volume1[, -1] - volume1[, -ncol(volume)]
  check_growth <- (!apply(growth_rate, 1, function(y) all(is.na(y))))
  indices <- which(check_growth==TRUE|is.na(check_growth))
  filterd_length <- length1[indices, ]
  filterd_volume <- volume1[indices, ]
  filterd_depth <- depth1[indices, ]
  volume <- volume[indices, ]
  length <- length[indices, ]
  depth <- depth[indices, ]
  tortuosity <- tortuosity[indices, ]
  radius <- radius[indices, ]
  
  #compute total root traits
  total_volume[i, 1:n_timepoints] <- colSums(volume, na.rm = TRUE)
  volume_growth_rate[i, ] <- total_volume[i, -1] - total_volume[i, -n_timepoints]
  total_length[i, 1:n_timepoints] <- colSums(length, na.rm = TRUE)
  length_growth_rate[i, ] <- total_length[i, -1] - total_length[i, -n_timepoints]
  total_depth[i, 1:n_timepoints] <- colSums(depth, na.rm = TRUE)
  depth_growth_rate[i, ] <- total_depth[i, -1] - total_depth[i, -n_timepoints]
  branch_number[i, 1:n_timepoints] <- colSums(!is.na(volume))
  number_growth_rate[i, ] <- branch_number[i, -1] - branch_number[i, -n_timepoints]
  
  mean_length_growth_rate[i, 1:(n_timepoints-1)] <- colMeans(filterd_length[, -1]-filterd_length[, -ncol(filterd_length)], na.rm = TRUE)
  mean_volume_growth_rate[i, 1:(n_timepoints-1)] <- colMeans(filterd_volume[, -1]-filterd_volume[, -ncol(filterd_volume)], na.rm = TRUE)
  mean_depth_growth_rate[i, 1:(n_timepoints-1)] <- colMeans(filterd_depth[, -1]-filterd_depth[, -ncol(filterd_depth)], na.rm = TRUE)
  
  mean_tortuosity[i, 1:n_timepoints] <- colMeans(tortuosity, na.rm = TRUE)
  mean_radius[i, 1:n_timepoints] <- colMeans(radius, na.rm = TRUE)
  mean_length[i, 1:n_timepoints] <- colMeans(length, na.rm = TRUE)
  mean_volume[i, 1:n_timepoints] <- colMeans(volume, na.rm = TRUE)
  
}

#######################################
############   modeling  ##############
#######################################
##e.g. total root volume
scaledVolume <- total_volume*0.5384*0.5384*0.5384
coefVal <- c()
v_predit <- matrix(data = NA, nrow = 38, ncol = 42)
for (i in 1:38)
{
  v <- data.frame(t = c(0:41)*4, v = scaledVolume[i, ])
  v <- v[!is.na(v[,2]),]
  maxVal <- max(scaledVolume[i, ], na.rm = TRUE)
  minVal <- min(scaledVolume[i, ], na.rm = TRUE)
  coef_initial <- coef(lm(logit(v/(maxVal+minVal)) ~ t,data = v))
  fit<-nls(v ~ phi1/(1+exp(-(phi2+phi3*t))),
           start = list(phi1 = maxVal+minVal, phi2 = coef_initial[1], phi3 = coef_initial[2]),
           data = v)#, trace = TRUE)
  coefVal <- rbind(coefVal, coef(fit))
  v_predit[i, ] <- coef(fit)[1]/(1+exp(-(coef(fit)[2]+coef(fit)[3]*c(0:41)*4)))
}

################ boxplot for t-tests comparing model parameters ################
forBoxplot <- data.frame(coefVal, ip = -coefVal[,2]/coefVal[,3], genotype = c(rep("B73", 12), rep("Mo17", 14),
                                                                              rep("hybrid", 12)))

forBoxplot$genotype <- factor(forBoxplot$genotype, levels = c("B73", "Mo17", "hybrid"))

my_comparisons <- list( c("B73", "Mo17"), c("Mo17", "hybrid"), c("B73", "hybrid") )

p1 <- ggviolin(forBoxplot, x = "genotype", y = "phi1",
                fill = "genotype", palette = c("#E69F00", "#56B4E9", "#CC79A7"),
                xlab = FALSE, ylab = expression(alpha[1]), add = "boxplot",
                add.params = list(fill = "white", size = 1), size = 1) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", 
                     label = "p.signif", size = 1, label.y = c(max(forBoxplot$phi1))) + 
  theme_pubr(base_size = 24)
p1$layers[[3]]$aes_params$textsize <- 8 


p2 <- ggviolin(forBoxplot, x = "genotype", y = "ip",
                fill = "genotype", palette = c("#E69F00", "#56B4E9", "#CC79A7"),
                xlab = FALSE, ylab = expression(alpha[2]), add = "boxplot",
                add.params = list(fill = "white", size = 1), size = 1) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif", size = 1) + 
  theme_pubr(base_size = 24)
p2$layers[[3]]$aes_params$textsize <- 8 

p3 <- ggviolin(forBoxplot, x = "genotype", y = "phi3",
                fill = "genotype", palette = c("#E69F00", "#56B4E9", "#CC79A7"),
                xlab = FALSE, ylab = expression(alpha[3]), add = "boxplot",
                add.params = list(fill = "white", size = 1), size = 1) + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif", size = 1) + 
  theme_pubr(base_size = 24)
p3$layers[[3]]$aes_params$textsize <- 8 

p <- grid.arrange(p1, p2, p3, ncol = 3)

##############   growth curves   #######################
rects <- data.frame(xstart = seq(1, 41, 6)*4, xend = seq(3, 41, 6)*4)
df <- data.frame(TRL = as.vector(t(scaledVolume)), 
                 timepoint = 4*rep(c(0:41), 38), genotype = c(rep("B73", 12*42), rep("Mo17", 14*42), rep("hybrid", 12*42)))
df$genotype <- factor(df$genotype, levels = c("B73", "Mo17", "hybrid"))

ggplot() + geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), color = "gray", alpha = 0.2) + 
  geom_point(data = df, aes(x = timepoint, y = TRL, color = genotype), alpha = 0.5, size = 2) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7")) +
  geom_line(aes(x = 4*c(0:41), y = colMeans(v_predit[1:12,], na.rm = TRUE)), color = "#E69F00", size = 2) + 
  geom_line(aes(x = 4*c(0:41), y = colMeans(v_predit[13:26,], na.rm = TRUE)), color = "#56B4E9", size = 2) + 
  geom_line(aes(x = 4*c(0:41), y = colMeans(v_predit[27:38,], na.rm = TRUE)), color = "#CC79A7", size = 2) + 
  scale_x_continuous(breaks=seq(0,164,24)) +
  theme_classic(base_size = 32) + labs(x = "hours", y = expression(total~root~volume(mm^{3})))

################# growth rates and accelerations ##############
growth_rate <- volume_growth_rate*0.5384*0.5384*0.5384/4
growth_rate[growth_rate>100] <- NA

df1 <- data.frame(TRL = as.vector(t(growth_rate)), name = rep(FileNames, each = 41), 
                  timepoint = 4*rep(c(1:41), 38), genotype = c(rep("B73", 12*41), 
                                                               rep("Mo17", 14*41),
                                                               rep("hybrid", 12*41)))
df1$genotype <- factor(df1$genotype, levels = c("B73", "Mo17", "hybrid"))

df2 <- summarySE(df1, measurevar="TRL", groupvars=c("timepoint", "genotype"), na.rm=TRUE)
df2$genotype <- factor(df2$genotype, levels = c("B73", "Mo17", "hybrid"))

acc <- (growth_rate[, -1] - growth_rate[, -41])/4

df3 <- data.frame(TRL = as.vector(t(acc)), name = rep(FileNames, each = 40), 
                  timepoint = 4*rep(c(2:41), 38), genotype = c(rep("B73", 12*40), rep("Mo17", 14*40),
                                                               rep("hybrid", 12*40)))

df4 <- summarySE(df3, measurevar="TRL", groupvars=c("timepoint", "genotype"), na.rm=TRUE)
df4$genotype <- factor(df4$genotype, levels = c("B73", "Mo17", "hybrid"))


p <- ggplot() + geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), color = "gray", alpha = 0.2) +     
  geom_line(data = df2, aes(x=timepoint, y=TRL, colour=genotype), size = 1) +
  geom_point(data = df2, aes(x=timepoint, y=TRL, colour=genotype), size = 2) + 
  geom_errorbar(data = df2, aes(x=timepoint, y=TRL, colour=genotype, ymin=(TRL-se), ymax=(TRL+se)), size = 1, width = 0) +
  theme_classic(base_size = 32) + labs(x = "hours", y = expression(total~root~volume~growth~rate(mm^{3}/h))) +  
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#CC79A7"))

p <- p +  scale_y_continuous(sec.axis = sec_axis(~.*1, name = expression(acceleration(mm^{3}/h^{2}))))

p <- p + geom_bar(data = df4, aes(x=timepoint-2, y=TRL, colour=genotype, fill = genotype), alpha = 0.5, size = 1, stat = "identity", position = "identity") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#CC79A7")) + scale_x_continuous(breaks=seq(0,164,24))
















###################################################
################ Summarized data ##################
###################################################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}


  

