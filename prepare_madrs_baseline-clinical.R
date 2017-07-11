#Takes prepare_madrs.csv file (output from prepare_clinical.py), runs stats and extracts useful info for write-up.

library('car') #VIF
library('psych') #Common stats functions
library('boot') #Bootstrapping
library('caret') #Model generation / evaluation
library('plyr') #Rename function
library('tidyr')
library('doMC') #Parallel processing
library('reshape2')
library('rfPermute')
library('randomForest')


#load data
#in_data <- read.csv('/home/orcasha/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs_vox_falff.csv',header=TRUE, sep=',') #Load data
#in_data <- read.csv('/home/orcasha/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs_vlsm_falff.csv', header = TRUE, sep = ',') #Load data
#in_data <- read.csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs_vlsm_sum_falff.csv', header = TRUE, sep = ',') #Load data
#in_data <- read.csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs_vlsm_falff.csv', header = TRUE, sep = ',') #Load data
#in_data <- read.csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs_vlsm_sum_falff.csv',header=TRUE, sep=',') #Load data
in_data <- read.csv('/home/peter/Dropbox/post_doc/florey_leeanne/study_scripts/prepare/prepare_madrs_falff.csv',header=TRUE, sep=',') #Load data

####Housekeeping - Define factors, rename, drop columns etc.####
in_data <- subset(in_data, select = -c(X, 
                                       site_id, 
                                       p_participant_num, 
                                       initials, Chol_3mth, 
                                       Trig_3mth, 
                                       HDLC_3mth, 
                                       LDLC_3mth, 
                                       sFol_3mth, 
                                       B12_3mth, 
                                       CRP_3mth, 
                                       madrs_score_stand, 
                                       Unnamed..0))

#Clean rec_treat_medic
in_data$rec_treat_medic[is.na(in_data$rec_treat_medic) == 1] = 0

in_data$gender <- factor(in_data$gender)
#in_data$mask_vlsm_overlap <- factor(in_data$mask_vlsm_overlap)
#in_data$madrs_group <- factor(in_data$madrs_group)
in_data <- rename(in_data, replace = c('MADRS_score_3mth' = 'madrs_score')) #Rename variable.
in_data <- rename(in_data, replace = c('Depress_history_3mths' = 'hist_dep')) #Rename variable.
in_data <- rename(in_data, replace = c('NIHSS_score_3mth' = 'nihss_score')) #Rename variable.
in_data$hist_dep[is.na(in_data$hist_dep)] <- 0 #Replace NAs with 0 (no history)
in_data$hist_dep <- factor(in_data$hist_dep)


####STANDARD STATS####
summary_d <- describeBy(in_data, 
                        group = in_data$madrs_group) #Split data by MADRS score.
summary_m <- describeBy(in_data, 
                        group = in_data$madrs_group, mat=TRUE) #Describe data by group, output as matrix for writing to spreadsheet.

write.csv(summary_m, file = '/home/peter/Dropbox/post_doc/florey_leeanne/papers/prepare/madrs_baseline/group_summaries.csv', sep=",")

#Write text files with low and high participant IDs.
#Note: Don't forget to go through output and add zeros for subs < 10
write(sprintf("%i-%i-%s", 
              in_data$site_id[in_data$madrs_group == 0], 
              in_data$p_participant_num[in_data$madrs_group == 0], 
              in_data$initials[in_data$madrs_group == 0]), 
              "/home/peter/Desktop/prepare/rest/output/prepare_low_madrs.txt", sep = " ")

write(sprintf("%i-%i-%s", 
              in_data$site_id[in_data$madrs_group == 1], 
              in_data$p_participant_num[in_data$madrs_group == 1], 
              in_data$initials[in_data$madrs_group == 1]), 
              "/home/peter/Desktop/prepare/rest/output/prepare_low_madrs.txt", sep = " ")


tbl <- table(in_data$madrs_group == 0, in_data$gender) #Create table of group x sex
print(tbl)
chi_tests.sex <- chisq.test(tbl,simulate.p.value = TRUE, B = 10000) #Chi test @ 10,0000 reps
chi_tests.sex$v <- sqrt(chi_tests.sex$statistic / sum(tbl)) #Effect size (R)

tbl=table(in_data$madrs_group == 0, in_data$hist_dep) #Create table of group x history of depression
print(tbl)
chi_tests.dep_hist <- chisq.test(tbl,simulate.p.value = TRUE,B = 10000) #Chi test @ 10,0000 reps
chi_tests.dep_hist$v <- sqrt(chi_tests.dep_hist$statistic / sum(tbl))
print(chi_tests.dep_hist)

tbl <- table(in_data$madrs_group == 1, in_data$depress) #Create table of group x reported sadness
print(tbl)
chi_tests.dep <- chisq.test(tbl,simulate.p.value = TRUE,B = 10000) #Chi test @ 10,0000 reps
chi_tests.dep$v <- sqrt(chi_tests.dep$statistic / sum(tbl))
print(chi_tests.dep)

tbl <- table(in_data$madrs_group == 1, in_data$discourage) #Create table of group x reported discouragement
print(tbl)
chi_tests.discourage <- chisq.test(tbl,simulate.p.value = TRUE,B = 10000) #Chi test @ 10,0000 reps
chi_tests.discourage$v <- sqrt(chi_tests.discourage$statistic / sum(tbl))
print(chi_tests.discourage)


tbl <- table(in_data$madrs_group == 1, in_data$lost_int) #Create table of group x reported loss of interest
print(tbl)
chi_tests.lost_int <- chisq.test(tbl,simulate.p.value = TRUE,B = 10000) #Chi test @ 10,0000 reps
chi_tests.lost_int$v <- sqrt(chi_tests.lost_int$statistic / sum(tbl))
print(chi_tests.lost_int)


tbl <- table(in_data$madrs_group, in_data$rec_treat_medic)
print(tbl)
chi_tests.meds <- chisq.test(tbl,simulate.p.value = TRUE,B = 10000) #Chi test @ 10,0000 reps
print(chi_tests.meds)

#Subsets for dep history

dep <- in_data[in_data$madrs_group == 1,]
dep_pos <- dep[dep$hist_dep == 1,]
dep_neg <- dep[dep$hist_dep == 0,]
t_tests.dep_madrs <- t.test(dep_pos$madrs_score, 
                            dep_neg$madrs_score, 
                            alternative = 'two.sided', 
                            paired = FALSE,
                            var.equal = FALSE)




#Run t-tests
t_tests.nihss <- t.test(in_data$nihss_score[in_data$madrs_group == 0], 
                        in_data$nihss_score[in_data$madrs_group == 1], 
                        alternative = c('two.sided'), 
                        paired = FALSE, 
                        var.equal = FALSE)

t_tests.age <- t.test(in_data$age_assessment[in_data$madrs_group == 0], 
                      in_data$age_assessment[in_data$madrs_group == 1], 
                      alternative = c('two.sided'), 
                      paired = FALSE, 
                      var.equal = FALSE)

t_tests.madrs <- t.test(in_data$madrs_score[in_data$madrs_group == 0], 
                        in_data$madrs_score[in_data$madrs_group == 1], 
                        alternative = c('two.sided'), 
                        paired = FALSE, 
                        var.equal = FALSE)



###RANDOM FOREST###


# #Test assumptions (multicollinearity)
# corrmat <- subset(in_data, select = c(falff_analysis_clustval_1, 
#                                       falff_analysis_clustval_2,
#                                       falff_analysis_clustval_3, 
#                                       falff_slow_4_analysis_clustval_1, 
#                                       falff_slow_4_analysis_clustval_2))
# 
# r <- cor(corrmat)
# cor.plot(r, upper = FALSE, numbers = TRUE, cex = 1)
# 
# #Due to highly correlated variables, PCA is employed.
# falff_main_all <- prcomp(subset(in_data, 
#                                 select = c(falff_analysis_clustval_1, 
#                                            falff_analysis_clustval_2, 
#                                            falff_analysis_clustval_3)))
# 
# falff_sub_all <- prcomp(subset(in_data, 
#                                select = c(falff_slow_4_analysis_clustval_1, 
#                                           falff_slow_4_analysis_clustval_2)))
# 
# #Cumulative % variance explained.
# falff_main_all_var <- cumsum((falff_main_all$sdev)^2) / sum(falff_main_all$sdev^2)
# falff_sub_all_var <- cumsum((falff_sub_all$sdev)^2) / sum(falff_sub_all$sdev^2)
# 
# biplot(falff_main_all)
# biplot(falff_sub_all)
# 
# in_data$falff_main_pca <- falff_main_all$x[,1] #Take 1st eigenvector
# in_data$falff_sub_pca <- falff_sub_all$x[,1] #Take 1st eigenvector
# 
# subdata <- subset(in_data, select = -c(falff_analysis_clustval_1, 
#                                        falff_analysis_clustval_2, 
#                                        falff_analysis_clustval_3, 
#                                        falff_slow_4_analysis_clustval_1, 
#                                        falff_slow_4_analysis_clustval_2))
# 
# ####Prepare data for random forest analysis####
# #Training data
# set.seed(411)
# train_idx <- createDataPartition(in_data$madrs_score, 
#                                  p = 0.8, 
#                                  list = FALSE) #Split 80 / 20
# 
# train_data <- in_data[train_idx,]
# 
# #Testing data with madrs score dropped
# test_data <- in_data[-train_idx,]
# madrs_obs <- test_data$madrs_score
# 
# ####ML####
# #Establish mtry parameter for tuning
# #k = length(train_data)
# k = 6
# mtries = k ^ 2 - k - 1
# tGrid <- expand.grid(mtry = 1:127)
# 
# #Establish training parameters
# params <- trainControl(method = 'repeatedcv', number = 10, repeats = 10)
# 
# #Set number of cores
# ncores <- detectCores()
# registerDoMC(cores = ncores - 1)
# 
# #Use caret to get best tune
# set.seed(411)
# fit_cv_rf <- train(madrs_score ~ gender + mask_vlsm_overlap + falff_analysis_clustval_1 + falff_analysis_clustval_2 + falff_analysis_clustval_3 + falff_slow_4_analysis_clustval_1 + falff_slow_4_analysis_clustval_2, 
#                      data = in_data, 
#                      method = 'rf',
#                      trControl = params,
#                      metric = 'RMSE',
#                      tuneGrid = tGrid,
#                      importance = TRUE,
#                      ntrees = 2000)
# 
# #Run RF + permutation of coefficents
# set.seed(411)
# fit_perms <- rfPermute(madrs_score ~ gender + mask_vlsm_overlap + falff_analysis_clustval_1 + falff_analysis_clustval_2 + falff_analysis_clustval_3 + falff_slow_4_analysis_clustval_1 + falff_slow_4_analysis_clustval_2, 
#                        data = in_data, 
#                        ntree = 2000,
#                        mtry = fit_cv_rf$bestTune[1,],
#                        nrep = 1000,
#                        num.cores = ncores)
# 
# 
# set.seed(411)
# fit_perms <- rfPermute(madrs_score ~ gender + mask_vlsm_overlap + falff_main_pca + falff_sub_pca, 
#                        data = in_data, 
#                        ntree = 2000,
#                        mtry = fit_cv_rf$bestTune[1,],
#                        nrep = 1000,
#                        num.cores = ncores)
# 
# print(fit_perms) #Stability of solution over number of trees
# print(rp.importance(fit_perms))
# plot(rp.importance(fit_perms), scale = TRUE) #Plot %incMSE (percent increase in MSE - higher = better) + IncNodePurity (more useful variables increase). Sig in red.



