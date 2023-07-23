library(gbm)
load("./model.RData")
cut_point <- 0.892
colname = c('res_d', 'ato_d', 'res_a', 'ato_a', 'chrg_d', 'chrg_a', 
            'pos_d', 'pos_a', 'sec_strc_d', 'sec_strc_a', 
            'res_d_n3', 'res_d_n2', 'res_d_n1', 'res_d_1', 'res_d_2', 'res_d_3', 
            'res_a_n3', 'res_a_n2', 'res_a_n1', 'res_a_1', 'res_a_2', 'res_a_3', 
            'distance', 'atom_type_d', 'atom_type_a', 'res_seq_num_d', 'res_seq_num_a')
test_data <- read.delim("./temp", header=F, sep="", col.names = colname)
original_test <- test_data
test_data <- test_data[ , 1:22]
for(i in c(1:4, 7:22)){
  test_data[, i] <- as.factor(test_data[, i])
}
for(i in 5:6){
  test_data[, i] <- as.integer(test_data[, i])
}
test_data <- test_data[, -7]
y_pred <- matrix(0, dim(test_data)[1], length(models))
for(piece_index in 1:length(models)){
  y_pred[ ,piece_index] <- predict(models[[piece_index]], test_data, n.trees = 5000, type="response")
}
average_prob <- apply((y_pred[,1:10]), 1, mean)
average_prob <- formatC(average_prob, digits = 3, format = "f")
original_test$probability <- average_prob
original_test$class <- "NHB"
original_test$class[original_test$probability > cut_point] <- "SHB"
output <- original_test[, c(1, 3, 23, 28, 26, 24, 27, 25, 29)]
for(i in 1:dim(test_data)[1]){
  output$Donor[i] = paste(paste(as.character(output[i, 1]), as.character(output[i, 5]), sep="_"), as.character(output[i, 6]), sep="@")
  output$Acceptor[i] = paste(paste(as.character(output[i, 2]), as.character(output[i, 7]), sep="_"), as.character(output[i, 8]), sep="@")
}
output <- output[, c(10, 11, 3, 4, 9)]
names(output) <- c("donor residue/atom", "acceptor residue/atom", "R from structure (A)", "Probability of forming a SHB", "Recommended hydrogen bond class (Probability threshold = 0.892)")
output <- output[order(output$Probability, decreasing=TRUE), ]
write.csv(output, file="predict_results.csv",row.names=F)
