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

######################################
########## density plot ############## 
######################################
FileNames <- c()
B73_length_com <- c()
B73_soilangle_com <- c()
B73_branchingangle_com <- c()
Mo17_length_com <- c()
Mo17_soilangle_com <- c()
Mo17_branchingangle_com <- c()
hybrid_length_com <- c()
hybrid_soilangle_com <- c()
hybrid_branchingangle_com <- c()

for ( i in 1:12 )
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
  soilAngle <- traits[, seq(13, ncol(traits), 11)] - 90
  branchingAngle <- traits[, seq(16, ncol(traits), 11)]
  
  volume[length<=5] <- NA
  soilAngle[length<=5] <- NA 
  branchingAngle[length<=5] <- NA
  length[length<=5] <- NA
  
  volume1 <- volume
  length1 <- length
  
  growth_rate <- volume1[, -1] - volume1[, -n_timepoints]
  growth_rate <- data.frame(rep(NA, nrow(growth_rate)), growth_rate)  
  length1[growth_rate<=0] <- NA
  volume1[growth_rate<=0] <- NA
  soilAngle[growth_rate<=0] <- NA
  growth_rate <- volume1[, -1] - volume1[, -ncol(volume)]
  check_growth <- (!apply(growth_rate, 1, function(y) all(is.na(y))))
  indices <- which(check_growth==TRUE|is.na(check_growth))
  soilAngle <- soilAngle[indices, ]
  length <- length[indices, ]
  volume <- volume[indices, ]
  branchingAngle <- branchingAngle[indices, ]
  
  temp <- data.frame(TRL = as.vector(t(length))%/%(0.1*max(length, na.rm = TRUE)), timepoint = rep(c(0:(n_timepoints-1)), nrow(length)))
  B73_length_com <- rbind(temp[!is.na(temp[, 1]), ], B73_length_com)
  temp <- data.frame(TRL = as.vector(t(volume))%/%(0.1*max(volume, na.rm = TRUE)), timepoint = rep(c(0:(n_timepoints-1)), nrow(volume)))
  temp <- data.frame(TRL = as.vector(t(soilAngle))%/%10, timepoint = rep(c(0:(n_timepoints-1)), nrow(soilAngle)))
  temp <- temp[temp[, 1]>0, ]
  B73_soilangle_com <- rbind(temp[!is.na(temp[, 1]), ], B73_soilangle_com)
  temp <- data.frame(TRL = as.vector(t(branchingAngle))%/%10, timepoint = rep(c(0:(n_timepoints-1)), nrow(branchingAngle)))
  B73_branchingangle_com <- rbind(temp[!is.na(temp[, 1]), ], B73_branchingangle_com)
 
}

for ( i in 13:26 )
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
  soilAngle <- traits[, seq(13, ncol(traits), 11)] - 90
  branchingAngle <- traits[, seq(16, ncol(traits), 11)]
  
  volume[length<=5] <- NA
  soilAngle[length<=5] <- NA 
  branchingAngle[length<=5] <- NA
  length[length<=5] <- NA
  
  volume1 <- volume
  length1 <- length
  
  growth_rate <- volume1[, -1] - volume1[, -n_timepoints]
  growth_rate <- data.frame(rep(NA, nrow(growth_rate)), growth_rate)  
  length1[growth_rate<=0] <- NA
  volume1[growth_rate<=0] <- NA
  soilAngle[growth_rate<=0] <- NA
  growth_rate <- volume1[, -1] - volume1[, -ncol(volume)]
  check_growth <- (!apply(growth_rate, 1, function(y) all(is.na(y))))
  indices <- which(check_growth==TRUE|is.na(check_growth))
  soilAngle <- soilAngle[indices, ]
  length <- length[indices, ]
  volume <- volume[indices, ]
  branchingAngle <- branchingAngle[indices, ]
  
  temp <- data.frame(TRL = as.vector(t(length))%/%(0.1*max(length, na.rm = TRUE)), timepoint = rep(c(0:(n_timepoints-1)), nrow(length)))
  Mo17_length_com <- rbind(temp[!is.na(temp[, 1]), ], Mo17_length_com)
  temp <- data.frame(TRL = as.vector(t(volume))%/%(0.1*max(volume, na.rm = TRUE)), timepoint = rep(c(0:(n_timepoints-1)), nrow(volume)))
  temp <- data.frame(TRL = as.vector(t(soilAngle))%/%10, timepoint = rep(c(0:(n_timepoints-1)), nrow(soilAngle)))
  temp <- temp[temp[, 1]>0, ]
  Mo17_soilangle_com <- rbind(temp[!is.na(temp[, 1]), ], Mo17_soilangle_com)
  temp <- data.frame(TRL = as.vector(t(branchingAngle))%/%10, timepoint = rep(c(0:(n_timepoints-1)), nrow(branchingAngle)))
  Mo17_branchingangle_com <- rbind(temp[!is.na(temp[, 1]), ], Mo17_branchingangle_com)
  
}

for ( i in 27:38 )
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
  soilAngle <- traits[, seq(13, ncol(traits), 11)] - 90
  branchingAngle <- traits[, seq(16, ncol(traits), 11)]
  
  volume[length<=5] <- NA
  soilAngle[length<=5] <- NA 
  branchingAngle[length<=5] <- NA
  length[length<=5] <- NA
  
  volume1 <- volume
  length1 <- length
  
  growth_rate <- volume1[, -1] - volume1[, -n_timepoints]
  growth_rate <- data.frame(rep(NA, nrow(growth_rate)), growth_rate)  
  length1[growth_rate<=0] <- NA
  volume1[growth_rate<=0] <- NA
  soilAngle[growth_rate<=0] <- NA
  growth_rate <- volume1[, -1] - volume1[, -ncol(volume)]
  check_growth <- (!apply(growth_rate, 1, function(y) all(is.na(y))))
  indices <- which(check_growth==TRUE|is.na(check_growth))
  soilAngle <- soilAngle[indices, ]
  length <- length[indices, ]
  volume <- volume[indices, ]
  branchingAngle <- branchingAngle[indices, ]
  
  temp <- data.frame(TRL = as.vector(t(length))%/%(0.1*max(length, na.rm = TRUE)), timepoint = rep(c(0:(n_timepoints-1)), nrow(length)))
  hybrid_length_com <- rbind(temp[!is.na(temp[, 1]), ], hybrid_length_com)
  temp <- data.frame(TRL = as.vector(t(volume))%/%(0.1*max(volume, na.rm = TRUE)), timepoint = rep(c(0:(n_timepoints-1)), nrow(volume)))
  temp <- data.frame(TRL = as.vector(t(soilAngle))%/%10, timepoint = rep(c(0:(n_timepoints-1)), nrow(soilAngle)))
  temp <- temp[temp[, 1]>0, ]
  hybrid_soilangle_com <- rbind(temp[!is.na(temp[, 1]), ], hybrid_soilangle_com)
  temp <- data.frame(TRL = as.vector(t(branchingAngle))%/%10, timepoint = rep(c(0:(n_timepoints-1)), nrow(branchingAngle)))
  hybrid_branchingangle_com <- rbind(temp[!is.na(temp[, 1]), ], hybrid_branchingangle_com)
  
}

############ density plot ##############
my.colors <- colorRampPalette (c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'))


df <- data.frame(rbind(B73_length_com, 
                       Mo17_length_com,
                       hybrid_length_com), 
                 genotype = c(rep("B73", nrow(B73_length_com)), 
                              rep("Mo17", nrow(Mo17_length_com)),
                              rep("hybrid", nrow(hybrid_length_com))))
df$TRL[df$TRL==10] <- 9
df$genotype <- factor(df$genotype, levels = c("B73", "Mo17", "hybrid"))

label_text <- seq(0, 0.1*max(df[, 1]), 0.1)
label_text <- paste(label_text, "~", label_text+0.1)

df$TRL <- factor(df$TRL, levels = rev(levels(factor(df$TRL))))
ggplot(df, aes(x=factor(timepoint*4), fill=TRL)) + 
  geom_bar(position="fill") + theme_grey(base_size = 24) + labs(x = "hours", y = "root proportion by length") + 
  scale_x_discrete(breaks = 4*seq(0, 41, 6)) +  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~genotype, scales = "free_x") +
  theme_classic(base_size = 32) + scale_fill_manual(name = "length", labels = label_text[10:1], values = my.colors(10)[1:10]) +
  theme(strip.background = element_rect(colour="white"), legend.position = "none")

#############
df <- data.frame(rbind(B73_soilangle_com, Mo17_soilangle_com,
                       hybrid_soilangle_com), genotype = c(rep("B73", nrow(B73_soilangle_com)), 
                                                           rep("Mo17", nrow(Mo17_soilangle_com)),
                                                           rep("hybrid", nrow(hybrid_soilangle_com))))	 
df$genotype <- factor(df$genotype, levels = c("B73", "Mo17", "hybrid"))

label_text <- seq(0, 10*max(df[, 1]), 10)
label_text <- paste(label_text, "~", label_text+10)
df$TRL <- factor(df$TRL, levels = rev(levels(factor(df$TRL))))
ggplot(df, aes(x=factor(timepoint*4), fill=TRL)) + 
  geom_bar(position="fill") + theme_grey(base_size = 24) + labs(x = "hours", y = "root proportion by soil angle") + 
  scale_x_discrete(breaks = 4*seq(0, 41, 6)) +  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~genotype, scales = "free_x") +
  theme_classic(base_size = 32) + scale_fill_manual(name = "soil angle", labels = label_text[8:1], values = my.colors(8)[1:8]) +
  theme(strip.background = element_rect(colour="white"), legend.position = "none")

###############
df <- data.frame(rbind(B73_branchingangle_com, 
                       Mo17_branchingangle_com,
                       hybrid_branchingangle_com), genotype = c(rep("B73", nrow(B73_branchingangle_com)),
                                                                rep("Mo17", nrow(Mo17_branchingangle_com)),
                                                                rep("hybrid", nrow(hybrid_branchingangle_com))))	 
df$genotype <- factor(df$genotype, levels = c("B73", "Mo17", "hybrid"))

label_text <- seq(0, 10*max(as.numeric(df[, 1])), 10)
label_text <- paste(label_text, "~", label_text+10)
df$TRL <- factor(df$TRL, levels = rev(levels(factor(df$TRL))))
ggplot(df, aes(x=factor(timepoint*4), fill=TRL)) + 
  geom_bar(position="fill") + labs(x = "hours", y = "root proportion by branching angle") + 
  scale_x_discrete(breaks = 4*seq(0, 41, 6)) +  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~genotype, scales = "free_x") +
  theme_classic(base_size = 32) + scale_fill_manual(name = "branching angle", labels = label_text[13:1], values = my.colors(13)[1:13]) +
  theme(strip.background = element_rect(colour="white"), legend.position = "right")

####################################################################                          
############# Developmental analysis along primary root ############
####################################################################
branching_density <- matrix(data = NA, nrow = n_plants, ncol = 42)
branching_depth_com <- c()
N <- c()
FileNames <- c()
PL <- c()
direction_angle_hist_com <- c()

for (i in 1:n_plants)
{
  data <- read.csv(TraitFileNames[i], header = FALSE, sep = "\t")
  data <- data[-ncol(data)]
  
  fnamex <- strsplit(TraitFileNames[i], "\\.")[[1]][1];
  FileNames[i] <- substring(sub(".*TIM", "", fnamex), 1, 5)
  
  traits <- data[, seq(2, ncol(data), 2)]
  TraitNames <- data[1, seq(1, ncol(data), 2)]
  colnames(traits) <- as.matrix(TraitNames)
  
  n_timepoints <- (ncol(traits)-5) / 11
  
  #   length_endpoint <- traits[ncol(traits)-9]
  #   volumes_endpoint <- traits[ncol(traits)-10]
  #   traits <- traits[length_endpoint > 5 & volumes_endpoint > 50, ]
  
  volume <- traits[, seq(6, ncol(traits), 11)]
  length <- traits[, seq(7, ncol(traits), 11)]
  soilAngle <- traits[, seq(13, ncol(traits), 11)] - 90
  branchingAngle <- traits[, seq(16, ncol(traits), 11)]
  
  volume[length<=5] <- NA
  soilAngle[length<=5] <- NA
  length[length<=5] <- NA
  branchingAngle[length<=5] <- NA
  
  volume1 <- volume
  length1 <- length
  
  growth_rate <- volume1[, -1] - volume1[, -n_timepoints]
  growth_rate <- data.frame(rep(NA, nrow(growth_rate)), growth_rate)  
  length1[growth_rate<=0] <- NA
  volume1[growth_rate<=0] <- NA
  soilAngle[growth_rate<=0] <- NA
  growth_rate <- volume1[, -1] - volume1[, -ncol(volume)]
  check_growth <- (!apply(growth_rate, 1, function(y) all(is.na(y))))
  indices <- which(check_growth==TRUE|is.na(check_growth))
  soilAngle <- soilAngle[indices, ]
  length <- length[indices, ]
  traits <- traits[indices, ]
  children <- aggregate(traits[, 2], list(traits[, 2]), length)
  id <- children[which.max(children[, 2]), 1]
  traits_children <- traits[traits[, 2]==id & traits[, 1]!=id, ]
  traits_parent <- traits[traits[, 1]==id, ]
  volume_children <- traits_children[, seq(6, ncol(traits_children), 11)]
  N_children <- colSums(!is.na(volume_children))
  length_parent <- traits_parent[, seq(7, ncol(traits_parent), 11)]
  density <- N_children/length_parent
  branching_density[i, 1:n_timepoints] <- as.vector(t(density))
  fork_g_depth <- traits_children[, 4]/length_parent[1, n_timepoints]
  fork_g_depth[fork_g_depth>1] <- 1
  hv <- hist(fork_g_depth, breaks = seq(0, 1, 0.1), plot = FALSE)
  branching_depth <- hv$counts
  branching_depth_com <- rbind(branching_depth_com, branching_depth)
  N[i] <- length(fork_g_depth)
  PL[i] <- length_parent[1, n_timepoints]
  
  ###########get radial angle at the emergence time###########
  tipXYZ <- traits_children[, seq(10, ncol(traits_children), 11)]
  #find first branching angle
  tipXYZx1 <- matrix(data = NA, nrow = nrow(tipXYZ), ncol = 4)
  tipXYZx2 <- matrix(data = NA, nrow = nrow(tipXYZ), ncol = 4)
  for (j in 1:nrow(tipXYZ))
  {
    nonNA_index <- which(!is.na(volume_children[j, ]))
    if (length(nonNA_index) > 1)
    {
      tipXYZx1[j, ] <- unlist(strsplit(as.character(tipXYZ[j, nonNA_index[1]]), "[^-?0-9.]+"))
      tipXYZx2[j, ] <- unlist(strsplit(as.character(tipXYZ[j, nonNA_index[2]]), "[^-?0-9.]+"))
    }
  }
  
  x <- as.numeric(tipXYZx2[, 2]) - as.numeric(tipXYZx1[, 2])
  y <- as.numeric(tipXYZx2[, 3]) - as.numeric(tipXYZx1[, 3])
  len <- x*x + y*y
  
  diff_angle <- data.frame(x = x, y = y, depth = traits_children[, 4])
  diff_angle <- diff_angle[len > 0, ]
  diff_angle <- diff_angle[order(diff_angle$depth), ]
  angle <- atan2(diff_angle$y, diff_angle$x)*180/pi
  diff_value <- angle[-1] - angle[-length(angle)]
  diff_value <- (diff_value + 360) %% 360
  direction_angle_hist <- hist(diff_value, breaks = seq(0, 360, 36), plot = FALSE)  
  direction_angle_hist_com <- rbind(direction_angle_hist_com, direction_angle_hist$counts)
                          
}    
 
###################
branching_depth_com <- branching_depth_com[-c(14, 17), ]
Mo17_ave <- colMeans(branching_depth_com[13:24, ])
B73_ave <- colMeans(branching_depth_com[1:12, ])
hybrid_ave <- colMeans(branching_depth_com[25:36, ])

plant_id <- c(paste("B", c(sprintf("%02d", 1:12)), sep = "-"), 
              paste("M", c(sprintf("%02d", 1:12)), sep = "-"),
              paste("h", c(sprintf("%02d", 1:12)), sep = "-"),
              paste("a", c(sprintf("%02d", 1:3)), sep = "-"))

freq <- rbind(branching_depth_com, B73_ave, Mo17_ave, hybrid_ave)
df <- data.frame(x = rep(plant_id, each = 10), y = rep(seq(0, 1, 0.1)[-1], 39), 
                 N_roots = as.vector(t(freq)), genotype = c(rep("B73", 10*12), 
                                                            rep("Mo17", 10*12),
                                                            rep("hybrid", 10*12),
                                                            rep("mean", 10*3)))
df$genotype <- factor(df$genotype, levels = c("B73", "Mo17", "hybrid", "mean"))

p <- ggplot(data =  df, aes(x = x, y = y - 0.05)) + 
  geom_tile(aes(fill = N_roots), colour = "white") +
  #scale_fill_gradient(low = "white", high = "steelblue") +
  scale_fill_distiller(direction = 1, palette = "Blues") + 
  scale_y_continuous(expand = c(0.01, 0), breaks = seq(0, 1, 0.1), 
                     labels = paste(seq(0, 100, 10), "%", sep = ""), trans = 'reverse') + 
  facet_grid(~genotype, scales = "free_x", space = "free_x") +
  theme_classic(base_size = 24) + labs(y = "percent distance", fill = "#root") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_rect(colour="white"))                          
  
##########################################################
direction_angle_hist_com <- direction_angle_hist_com[-c(14, 17), ]
B73_ave <- colMeans(direction_angle_hist_com[1:12, ])
Mo17_ave <- colMeans(direction_angle_hist_com[13:24, ])
hybrid_ave <- colMeans(direction_angle_hist_com[25:36, ])

plant_id <- c(paste("B", c(sprintf("%02d", 1:12)), sep = "-"), 
              paste("M", c(sprintf("%02d", 1:12)), sep = "-"),
              paste("h", c(sprintf("%02d", 1:12)), sep = "-"),
              paste("a", c(sprintf("%02d", 1:3)), sep = "-"))
freq <- rbind(direction_angle_hist_com, B73_ave, Mo17_ave, hybrid_ave)
df <- data.frame(x = rep(plant_id, each = 10), y = rep(seq(0, 360, 36)[-1], 39), 
                 N_roots = as.vector(t(freq)), genotype = c(rep("B73", 10*12), 
                                                            rep("Mo17", 10*12),
                                                            rep("hybrid", 10*12),
                                                            rep("mean", 10*3)))
df$genotype <- factor(df$genotype, levels = c("B73", "Mo17", "hybrid", "mean"))

saveName <- paste(savePath, "radial_anlge_gap_new.pdf", sep = "")


p <- ggplot(data =  df, aes(x = x, y = y - 18)) + 
  geom_tile(aes(fill = N_roots), colour = "white") +
  scale_fill_distiller(direction = 1, palette = "Blues") + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 360, 36), trans = 'reverse') +
  facet_grid(~genotype, scales = "free_x", space = "free_x") +
  theme_classic(base_size = 24) + labs(y = "radial angle", fill = "#root") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_rect(colour="white"))
