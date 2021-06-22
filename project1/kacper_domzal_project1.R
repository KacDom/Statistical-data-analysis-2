library(deal)
library(yeastCC)
library(tidyr)
library(dplyr)
library(ggplot2)
library(dplyr)

#1. Data Preprocessing
  expr <- as.data.frame(t(exprs(yeastCC)[orf800,]))
  cat("Observations:", nrow(expr), "\n")
  cat("Genes:", ncol(expr), "\n")
  
  # Replace missing data with gene expression median
  medianfill <- function(exprCol) {
    exprCol[which(is.na(exprCol))] <- median(exprCol, na.rm=T) 
    return(exprCol)}

  expr <- as.data.frame(apply(expr, 2, medianfill))
  
  # apply() alternative
  gene.iqr <- function(exprCol) {
      quantile(exprCol, c(0.75)) - quantile(exprCol, c(0.25))}
  
  iqr <- apply(expr, 2, gene.iqr)
  
  # Keep only genes with iqr > 1.6
  expr = expr[, iqr > 1.6]   # keep only genes with large variation
  
  # deal package seems not to work if all nodes are continuous.
  # Thus, we need to add a discrete "dummy" node to work around this issue. 
  # Adding this node will not have any influence on the final results.
  expr$dummy <- factor(rep("1", nrow(expr)))

#2. Create prior structure 
  G0  = network(expr, specifygraph=TRUE, inspectprob=TRUE)

  # We don't want any arrows starting from the "dummy" node, thus we construct a list of banned dependencies:
  banlist(G0) <- matrix(c(11,11,11,11,11,11,11,11,11,11,1,2,3,4,5,6,7,8,9,10),ncol=2)
  plot(G0)
  
  
#3. Show local probability distribution
  localprob(G0)
  
  
#4. Compute joint prior distribution
  prior1 <- jointprior(G0, 5)  # equivalent to imaginary sample size = 5

  
#5. Learn the initial network
  G1 <- getnetwork(learn(G0, expr, prior1))
  print(G1$score)


#6. Search for optimal network
  nwSearch <- autosearch(G1, expr, prior1, removecycles=FALSE, trace=FALSE)
  BN <- getnetwork(nwSearch)
  plot(BN)
  BN$score


#7. Calculating the variance for each of 10 genes
  expr_var = apply(X = expr, MARGIN = 2, FUN = var)
  expr_var

 
#8. Perturbating experimental values of each gene 30 times - generating 30 perturbated datasets
  perturb_dataset = function(){
  results = data.frame(matrix(NA, nrow = nrow(expr), ncol = ncol(expr)))
  colnames(results) = colnames(expr)
  perturbated_datasets = list()
  for (data_set_nr in 1:30){
      for (col_index in 1:10){
          noise = rnorm(n=nrow(expr[col_index]), mean = 0, sd = sqrt(expr_var[col_index]/10))
          perturbated = expr[col_index] + noise 
          results[col_index] = perturbated}
      results$dummy = rep(1, nrow(expr))
      rownames(results) = row.names(expr)
      perturbated_datasets[[data_set_nr]] = results}
  return (perturbated_datasets)}
  
  perturbated_datasets = perturb_dataset()

#9. Boxplot of gene YHR143W in each out of 30 perturbated datasets  
  boxplot_of_YHR143W = function(){
    boxplot_result <- data.frame("YHR143W" = rep(1, nrow(expr)))
    for (data_set in 1:30){
        df <- as.data.frame(perturbated_datasets[data_set])
        boxplot_result[data_set] <- df$YHR143W}
        
    boxplot(t(boxplot_result),data=t(boxplot_result), main="Empirical distribution",
       xlab="Experiment no.", ylab="Data distribution")}
  boxplot_of_YHR143W()
  

#10. Repeat steps 2-6 with each of 30 perturbated datasets and plot networtk PBN5

  #functtion for finding optimal networks
  build.optimal.network <- function(exprData, N){
    N0 = getnetwork(learn(G0, exprData, prior1))
    nwasarch = autosearch(N0, exprData, prior1, removecycles=FALSE, trace=FALSE)
    getnetwork(nwasarch)}

  PBNs = list()
  for (dataset in 1:length(perturbated_datasets)){
    PBN <- build.optimal.network(perturbated_datasets[[dataset]], 5)
    PBNs[[dataset]] = PBN}

  plot(PBNs[[5]])


#11. For each edge contained in BN*, calculate the relative frequency of appearance among the PBNi
  # Making a list of edges in BN*
  give_edges_BN_star = function(){
      BN_edges_names = list()
      BN_edges_nr = list()
      edge_counter = 0
      for (i in BN$continuous){
        for (parent_node in BN$nodes[[i]]$parent){
          edge_counter = edge_counter + 1
          node_pairs = c(parent_node, i)
          BN_edges_nr[[length(BN_edges_nr) + 1]] = node_pairs 
          BN_edges_names[[edge_counter]] = paste(BN$nodes[[node_pairs[[1]]]]$name,BN$nodes[[node_pairs[[2]]]]$name,sep="->")}}
     return(list(BN_edges_names, BN_edges_nr))}
  
  # Get edges from a certain network
  read_edges = function(learned_network){
      edges = list()
      for (i in learned_network$continuous){
        for (parent_node in learned_network$nodes[[i]]$parent) {
          edges[[length(edges) + 1]] = c(parent_node, i)}}
      return(edges)}

  calculate_edge_frequency = function(){
      BN_edges_names_nr = give_edges_BN_star()
      # reading edges for each network in perturbed networks
      edges_in_perturbed = list()
      for (i in 1:length(perturbated_datasets)){
        edges_in_perturbed[[i]] = read_edges(PBNs[[i]])}


      # calculating frequency of edges from BN* in each of PBN_i
      edges_frequency = list()
      for (i in 1:length(BN_edges_names_nr[[2]])){
        curr_edge = BN_edges_names_nr[[2]][[i]]
        curr_count = 0
          for (list_of_edges in edges_in_perturbed){
            if (list(curr_edge) %in% list_of_edges){curr_count = curr_count + c(1)}}
        edges_frequency[[i]] = curr_count/length(PBNs)}
      return(list(edges_frequency, BN_edges_names_nr[[1]]))}
    

  edges_frequency_result = calculate_edge_frequency()
  edges_in = as.data.frame(unlist(edges_frequency_result[1]), row.names = unlist(edges_frequency_result[2]))
  colnames(edges_in) = c('frequency')
  
  
  ggplot(edges_in, aes(x = reorder(row.names(edges_in), frequency), y = frequency )) +
    geom_bar(stat = 'identity') + ggtitle('Edge that are in BN* occurrence frequency') + coord_flip(ylim = c(0,1)) +
   xlab('Edge name') + ylab('Edge occurrence frequency') + theme_minimal()


# 12. For each edge not contained in BN*, establish its relative frequency among the PBNi
  edges_not_in_BN_star_frequency = function(){
      edges_nr_in_BN_star = give_edges_BN_star()[[2]]
      edges_frequency_not_in_BN = list()
      result = data.frame()

          for (PBNi in PBNs){
              edges = read_edges(PBNi)
              for (edge in edges){
                  if (!(list(edge) %in% edges_nr_in_BN_star )){
                     current_edge =  paste(BN$nodes[[edge[[1]]]]$name,BN$nodes[[edge[[2]]]]$name,sep="->")
                      if (current_edge %in% row.names(result)){result[current_edge,] = result[current_edge,]+ 1
                      }else{working_df = data.frame('occurrences'= 1, row.names = current_edge)
                        result = rbind(result, working_df)
                            }}}}
      result$frequency = c(0)
      for (index in 1:length(result$frequency)){
          result$frequency[index] = result$occurrences[index] / length(perturbated_datasets)}

      result}

  edges_not_in = edges_not_in_BN_star_frequency() 
  
  ggplot(edges_not_in, aes(x = reorder(row.names(edges_not_in), frequency), y = frequency )) +
    geom_bar(stat = 'identity') + ggtitle('Edge that are NOT in BN* occurrence frequency') + coord_flip(ylim = c(0,1)) +
    xlab('Edge name') + ylab('Edge occurrence frequency') + theme_minimal()
  


 
