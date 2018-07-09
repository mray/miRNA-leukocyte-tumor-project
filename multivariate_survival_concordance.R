library(survival)
library(Matrix)

# import survival data
survival_all = read.delim("xena_survival.txt", header=T)

# import mRNA cibersort results
mrna_cibersort_result = read.delim("mrna_cibersort_result.txt", header=T)

# import miRNA cibersort results
mirna_cibersort_result = read.delim("mirna_cibersort_result.txt", header=T)

# subset all files by matched samples
common_patients = Reduce(
  intersect,
  list(
    row.names(mrna_cibersort_result),
    row.names(mirna_cibersort_result),
    row.names(survival_all)
  )
)

mrna_cibersort_matrix = data.matrix(mrna_cibersort_result[common_patients,])
mirna_cibersort_matrix = data.matrix(mirna_cibersort_result[common_patients,])
survival_df = survival_all[common_patients,]

# randomize sample order
shuffled_order = sample(dim(survival_df)[1])
shuffled_survival_df = survival_df[shuffled_order,]
row.names(shuffled_survival_df) = row.names(survival_df)

surv = Surv(survival_df$survival_time, survival_df$dead)

k = 5
kfl = kfold_indexes(length(common_patients), k)


# MULTIVARIATE SURVIVAL ANALYSIS #

multivariate_cibersort_surv_df = cbind(surv, data.frame(cibersort_result))

multivariate_surv_cph = coxph(surv ~ ., data=multivariate_cibersort_surv_df)
multivariate_surv_cph_summary = summary.coxph.custom(multivariate_surv_cph)

multivariate_pvalues = -log10(multivariate_surv_cph_summary$coefficients[,5])

save.fig(append.extension(output.dir('multivariate_cox_regr_coef_pvalues')), width=10, height=7)
par(mar=c(13, 4, 4, 2) + 0.1)
barplot(
  multivariate_pvalues,
  beside=TRUE,
  col=matplotlib_colors[1],
  las=3,
  main=paste(c('Cox regression coefficient P-values: ', toupper(disease)), collapse=''),
  ylab='-log10(Cox regression coef. P-value)'
)
dev.off()

multivariate_coefs = multivariate_surv_cph_summary$coefficients[,1]

save.fig(append.extension(output.dir('multivariate_cox_regr_coefs')), width=10, height=7)
par(mar=c(13, 4, 4, 2) + 0.1)
barplot(
  multivariate_coefs,
  beside=TRUE,
  col=matplotlib_colors[1],
  las=3,
  main=paste(c('Cox regression coefficients: ', toupper(disease)), collapse=''),
  ylab='Cox regression coef.'
)
dev.off()

# /P-values from multivariate fit

descs = c('mRNA', 'miRNA', 'combined')
shuffled_descs = paste('Shuffled', descs)

sv_eval_results = list()
shuffled_sv_eval_results = list()
for (desc in descs) {
  sv_eval_results[[desc]] = list()
}
for (desc in shuffled_descs) {
  shuffled_sv_eval_results[[desc]] = list()
}

concordance_values = data.frame(
  matrix(0, ncol=length(descs), nrow=k)
)
names(concordance_values) = descs

shuffled_concordance_values = data.frame(
  matrix(0, ncol=length(shuffled_descs), nrow=k)
)
names(shuffled_concordance_values) = shuffled_descs

combined_concordance = data.frame(
  matrix(0, ncol=length(sm_choices), nrow=k)
)
names(combined_concordance) = sm_choices

surv_eval = function(surv_train, surv_test, data_train, data_test) {
  sv_train_df = data.frame(sv=surv_train, data_train)
  sv_test_df = data.frame(sv=surv_test, data_test)
  cc = coxph.control(iter.max=100)
  cph_train = coxph(sv ~ ., data=sv_train_df, control=cc)
  risk_test = predict(cph_train, newdata=sv_test_df, type='risk')
  concordance = survConcordance(surv_test ~ risk_test)$concordance
  
  list(
    cph_train=cph_train,
    risk_test=risk_test,
    concordance=concordance
  )
}

for (i in 1:k) {
  train_test = build.sig.env$kfl[[i]]
  train_patients = row.names(mirna_cibersort_result)[train_test$train]
  test_patients = row.names(mirna_cibersort_result)[train_test$test]
  surv_train = Surv(survival_df[train_patients, 'survival_time'], survival_df[train_patients, 'dead'])
  surv_test = Surv(survival_df[test_patients, 'survival_time'], survival_df[test_patients, 'dead'])
  shuffled_surv_train = Surv(shuffled_survival_df[train_patients, 'survival_time'], shuffled_survival_df[train_patients, 'dead'])
  shuffled_surv_test = Surv(shuffled_survival_df[test_patients, 'survival_time'], shuffled_survival_df[test_patients, 'dead'])
  
  # mRNA
  mrna_prop_train = data.matrix(mrna_cibersort_result[train_patients,])
  mrna_prop_test = mrna_cibersort_result[test_patients,]
  mrna_sv_results = surv_eval(surv_train, surv_test, mrna_prop_train, mrna_prop_test)
  sv_eval_results[['mRNA']][[i]] = mrna_sv_results
  concordance_values[i, 'mRNA'] = mrna_sv_results$concordance
  shuffled_mrna_sv_results = surv_eval(shuffled_surv_train, shuffled_surv_test, mrna_prop_train, mrna_prop_test)
  shuffled_sv_eval_results[['Shuffled mRNA']][[i]] = shuffled_mrna_sv_results
  shuffled_concordance_values[i, 'Shuffled mRNA'] = shuffled_mrna_sv_results$concordance
  # /mRNA
  
  # miRNA
  mirna_prop_train = mirna_prop_full[train_patients,]
  mirna_prop_test = mirna_prop_full[test_patients,]
  mirna_sv_results = surv_eval(surv_train, surv_test, mirna_prop_train, mirna_prop_test)
  sv_eval_results[['miRNA']][[i]] = mirna_sv_results
  concordance_values[i, 'miRNA'] = mirna_sv_results$concordance
  shuffled_mirna_sv_results = surv_eval(shuffled_surv_train, shuffled_surv_test, mirna_prop_train, mirna_prop_test)
  shuffled_sv_eval_results[['Shuffled miRNA']][[i]] = shuffled_mirna_sv_results
  shuffled_concordance_values[i, 'Shuffled miRNA'] = shuffled_mirna_sv_results$concordance
  # /miRNA
  
  # combined
  combined_prop_train = cbind(mrna_prop_train, mirna_prop_train)
  combined_prop_test = cbind(mrna_prop_test, mirna_prop_test)
  combined_full_sv_results = surv_eval(surv_train, surv_test, combined_prop_train, combined_prop_test)
  sv_eval_results[['Combined']][[i]] = combined_sv_results
  concordance_values[i, 'Combined'] = combined_sv_results$concordance
  shuffled_combined_full_sv_results = surv_eval(shuffled_surv_train, shuffled_surv_test, combined_prop_train, combined_prop_test)
  shuffled_sv_eval_results[['Shuffled combined']][[i]] = shuffled_combined_sv_results
  shuffled_concordance_values[i, 'Shuffled combined'] = shuffled_combined_sv_results$concordance
  # /combined
}

save.fig(output.dir(append.extension('concordance_boxplot')))
par(mar=c(9, 4, 4, 2) + 0.1)
boxplot(
  concordance_values,
  col=matplotlib_colors,
  main=paste(c('Survival concordance: ', collapse='')),
  ylab='Survival concordance',
  las=3
)
dev.off()

save.fig(output.dir(append.extension('shuffled_concordance_boxplot')))
par(mar=c(9, 4, 4, 2) + 0.1)
boxplot(
  shuffled_concordance_values,
  col=matplotlib_colors,
  main=paste(c('Shuffled survival concordance: ', collapse='')),
  ylab='Survival concordance',
  las=3
)
dev.off()

save.fig(output.dir(append.extension('combined_concordance_boxplot')))
par(mar=c(9, 4, 4, 2) + 0.1)
boxplot(
  combined_concordance,
  col=matplotlib_colors,
  main=paste(c('Survival concordance: ', collapse='')),
  ylab='Survival concordance',
  las=3
)
dev.off()