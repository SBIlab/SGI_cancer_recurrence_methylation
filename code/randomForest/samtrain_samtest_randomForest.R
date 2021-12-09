library(randomForest)
library(e1071)
library(caret)
library(prediction)
library(performance)
library(ROCR)
library(dplyr)
library(stringr)
library(ranger)
library(binaryLogic)
source("../function/get_column_combi.R")

#data_setup
toget_immune_with_clinical = read.table("../../data/SGI_toget_immune.txt", header = TRUE, sep = "\t", row.names = 1)
TIL_immune_with_clinical = read.table("../../data/SGI_TIL_immune.txt", header = TRUE, sep = "\t", row.names = 1)
PBMC_immune_with_clinical = read.table("../../data/SGI_PBMC_immune.txt", header = TRUE, sep = "\t", row.names = 1)


#column_combination
toget_col = colnames(select(toget_immune_with_clinical, -Recurrence, -Recurrence_int))
toget_colcombi = toget_immune_column_combi()

TIL_col = colnames(select(TIL_immune_with_clinical, -Recurrence, -Recurrence_int))
TIL_colcombi = TIL_immune_column_combi()

PBMC_col = colnames(select(PBMC_immune_with_clinical, -Recurrence, -Recurrence_int))
PBMC_colcombi = PBMC_immune_column_combi()


#variable declearation
toget_perf_rf_list = list()
TIL_perf_rf_list = list()
PBMC_perf_rf_list = list()

#model_function_definition
rf = function(i, x_R, y_R, x_NR, y_NR){
  set.seed(2*i)
  ind_NR = sample(2, nrow(x_NR), replace=TRUE, prob = c(0.7, 0.3))
  ind_R = sample(2, nrow(x_R), replace=TRUE, prob = c(0.7, 0.3))
  
  train_x = bind_rows(filter(x_NR, ind_NR==1), filter(x_R, ind_R==1))
  train_y = c(y_NR[ind_NR==1], y_R[ind_R==1])
  test_x = bind_rows(filter(x_NR, ind_NR==2), filter(x_R, ind_R==2)) 
  test_y = c(y_NR[ind_NR==2], y_R[ind_R==2])
  
  rf_model = randomForest(as.factor(train_y)~., data=train_x, proximity = TRUE, importance = TRUE) 
  
  pred_prob = predict(rf_model, newdata = test_x, type="prob")[,2]
  pred_test = ROCR::prediction(pred_prob, test_y)
  
  roc_value = ROCR::performance(pred_test, measure = "tpr", x.measure = "fpr")
  roc_df = data.frame(roc_value@alpha.values[[1]], roc_value@y.values[[1]], roc_value@x.values[[1]])
  colnames(roc_df) = c(roc_value@alpha.name, roc_value@y.name, roc_value@x.name)
  
  auc = ROCR::performance(pred_test, measure = "auc")@y.values[[1]]
  
  perf_output = data.frame("auc" = auc)
  rownames(perf_output) = i
  output_list = list()
  output_list[[1]] = perf_output
  
  return(output_list)
  
}

#monte-carlo CV
for (n in 1:length(toget_colcombi)){

  gc()
  toget_perf_rf = data.frame()
  TIL_perf_rf = data.frame()
  PBMC_perf_rf = data.frame()
  
  #input_column_parsing
  toget_input_col = toget_col[ strsplit(toget_colcombi[n], "")[[1]] == "1" ]
  if(n<=length(TIL_colcombi)){
    TIL_input_col = TIL_col[ strsplit(TIL_colcombi[n], "")[[1]] == "1" ]
    PBMC_input_col = PBMC_col[ strsplit(PBMC_colcombi[n], "")[[1]] == "1" ]
  }

  toget_x_R = select(filter(toget_immune_with_clinical, Recurrence_int==1), toget_input_col)
  toget_y_R = select(filter(toget_immune_with_clinical, Recurrence_int==1), Recurrence_int)
  toget_x_NR = select(filter(toget_immune_with_clinical, Recurrence_int==0), toget_input_col)
  toget_y_NR = select(filter(toget_immune_with_clinical, Recurrence_int==0), Recurrence_int)
  
  if(n<=length(TIL_colcombi)){
    TIL_x_R = select(filter(TIL_immune_with_clinical, Recurrence_int==1), TIL_input_col)
    TIL_y_R = select(filter(TIL_immune_with_clinical, Recurrence_int==1), Recurrence_int)
    TIL_x_NR = select(filter(TIL_immune_with_clinical, Recurrence_int==0), TIL_input_col)
    TIL_y_NR = select(filter(TIL_immune_with_clinical, Recurrence_int==0), Recurrence_int)

    PBMC_x_R = select(filter(PBMC_immune_with_clinical, Recurrence_int==1), PBMC_input_col)
    PBMC_y_R = select(filter(PBMC_immune_with_clinical, Recurrence_int==1), Recurrence_int)
    PBMC_x_NR = select(filter(PBMC_immune_with_clinical, Recurrence_int==0), PBMC_input_col)
    PBMC_y_NR = select(filter(PBMC_immune_with_clinical, Recurrence_int==0), Recurrence_int)
  }
  
  
  #modeling
  for (i in 1:100){
    toget_output = rf(i, toget_x_R, toget_y_R$Recurrence_int, toget_x_NR, toget_y_NR$Recurrence_int)
    toget_perf_rf = bind_rows(toget_perf_rf, toget_output[[1]])
    
    if(n <= length(TIL_colcombi)){
      
      TIL_output = rf(i, TIL_x_R, TIL_y_R$Recurrence_int, TIL_x_NR, TIL_y_NR$Recurrence_int)
      TIL_perf_rf = bind_rows(TIL_perf_rf, TIL_output[[1]])
      
      PBMC_output = rf(i, PBMC_x_R, PBMC_y_R$Recurrence_int, PBMC_x_NR, PBMC_y_NR$Recurrence_int)
      PBMC_perf_rf = bind_rows(PBMC_perf_rf, PBMC_output[[1]])
      
    }
    
  }
  
  # performance record
  toget_record_df = data.frame(strsplit(toget_colcombi[n], "")[[1]]) %>% t
  colnames(toget_record_df) = toget_col
  toget_perf_rf_list[[n]] = data.frame(toget_record_df, 
                                            auc_mean = mean(toget_perf_rf$auc), auc_sd = sd(toget_perf_rf$auc) )
  if(n <= length(TIL_colcombi)){
    TIL_record_df = data.frame(strsplit(TIL_colcombi[n], "")[[1]]) %>% t
    colnames(TIL_record_df) = TIL_col
    TIL_perf_rf_list[[n]] = data.frame(TIL_record_df,
                                           auc_mean = mean(TIL_perf_rf$auc), auc_sd = sd(TIL_perf_rf$auc) )
    
    PBMC_record_df = data.frame(strsplit(PBMC_colcombi[n], "")[[1]]) %>% t
    colnames(PBMC_record_df) = PBMC_col
    PBMC_perf_rf_list[[n]] = data.frame(PBMC_record_df,
                                             auc_mean = mean(PBMC_perf_rf$auc), auc_sd = sd(PBMC_perf_rf$auc) )
    
  }
  
  
}


output_toget_perf = do.call(rbind, toget_perf_rf_list) %>% as.data.frame %>% arrange(desc(auc_mean))
output_TIL_perf = do.call(rbind, TIL_perf_rf_list) %>% as.data.frame %>% arrange(desc(auc_mean))
output_PBMC_perf = do.call(rbind, PBMC_perf_rf_list) %>% as.data.frame %>% arrange(desc(auc_mean))

dir.create("../../result")
write.table(output_toget_perf, "../../result/Toget_rf_samtrain_samtest.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(output_TIL_perf, "../../result/TIL_rf_samtrain_samtest.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(output_PBMC_perf, "../../result/PBMC_rf_samtrain_samtest.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
