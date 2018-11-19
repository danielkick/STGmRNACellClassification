

#' Run Supervised Models
#'
#' This is wrapper funciton for to make it easy to train multiple supervised models at once.
#' It returns a list of all the models run after testing a variety of free parameter values.
#'
#' @title Run Supervised Models
#' @aliases
#'
#' @author Daniel Kick (\email{daniel.r.kick@@protonmail.com})
#'
#' @keywords Clustering
#'
#' @export
#'
#' @return clustering assignments
#'
#' @examples
#'


run_supervised_models <- function(use.seed = 8743436,
                                  input.df = use.input.df, #This df must contain a Cell col
                                  use.cv = 5){

  #use.seed = 8743436
  #input.df = mrna_target

  set.seed(use.seed)
  glmnetm <- train(
    Cell ~ ., input.df,
    metric = "Accuracy",
    method = "glmnet",
    tuneGrid = expand.grid(
      alpha = seq(0, 1, length = 5),
      lambda = seq(0.0001, 1, length = 100)
    ),
    trControl = trainControl(method = "cv", number = use.cv, verboseIter = TRUE)
  )
  # Plot the results
  #plot(glmnetm)


  #kNN
  set.seed(use.seed)
  knnm <- train(
    Cell ~ . ,
    tuneGrid = expand.grid(k = seq(from =1, to =20, by = 1)),
    #tuneLength = 20,
    data = input.df,
    method = "knn",
    trControl = trainControl(method = "cv", number = use.cv, verboseIter = TRUE)
  )

  #plot(knnm)
  #print(max(knnm$results$Accuracy, na.rm = T))

  #LDA
  set.seed(use.seed)
  ldam <- train(
    Cell ~ . ,
    data = input.df,
    method = "lda",
    trControl = trainControl(method = "cv", number = use.cv, verboseIter = TRUE)
  )

  #print(ldam)

  #neural networks
  set.seed(use.seed)
  nnmm <- train(
    Cell ~ . ,
    data = input.df,
    tuneGrid = expand.grid(decay = seq(from = 0.1, to = 1, by = 0.05)),
    method = "multinom",
    trControl = trainControl(method = "cv", number = use.cv, verboseIter = TRUE)
  )

  #plot(nnmm)
  #print(max(nnmm$results$Accuracy, na.rm = T)) #81% accuracy

  set.seed(use.seed)
  nnm <- train(
    Cell ~ . ,
    data = input.df,
    tuneGrid = expand.grid(size = seq(from = 1, to = 12, by = 1),
                           decay = seq(from = 0.1, to = 1.0, by = 0.2)),
    method = "nnet",
    trControl = trainControl(method = "cv", number = use.cv, verboseIter = TRUE)
  )

  #plot(nnm)
  #print(max(nnm$results$Accuracy, na.rm = T))

  #random forest
  max.mtry = 30
  set.seed(use.seed)
  rfm <- train(
    Cell ~ . ,
    #tuneLength = 3,
    tuneGrid = data.frame(mtry = rep(seq(1, max.mtry, by = 1), times = 2), #mtry can be any number from 2 to the number of columns
                          splitrule = rep(c("extratrees", "gini"), each = max.mtry), #the docs make it look like these are the two to use for classification
                          min.node.size = rep(2, each = (max.mtry*2))),
    data = input.df,
    method = "ranger",
    trControl = trainControl(method = "cv", number = use.cv, verboseIter = TRUE)
  )

  #print(rfm)
  #plot(rfm)
  #print(max(rfm$results$Accuracy, na.rm = T))

  #svm
  set.seed(use.seed)
  svm.rad <- train(
    Cell ~ .,
    tuneGrid = expand.grid(sigma = seq(from = 0.001, to = 0.5, by = 0.01),
                           C = seq(from = 0.5, to = 5.5, by = 1)),
    data = input.df,
    method = "svmRadial",
    trControl = trainControl(method = "cv", number = use.cv, verboseIter = TRUE)
  )

  #svm.rad
  #plot(svm.rad)
  #print(max(svm.rad$results$Accuracy, na.rm = T))

  set.seed(use.seed)
  svm.lin <- train(
    Cell ~ .,
    tuneGrid = expand.grid(cost = seq(from = 0.001, to = 0.5, by = 0.01)),
    data = input.df,
    method = "svmLinear2",
    trControl = trainControl(method = "cv", number = use.cv, verboseIter = TRUE)
  )

  #return(list(knnm,ldam,nnmm,nnm,rfm,svm.rad,svm.lin))
  return(list(GLMNet=glmnetm,
              kNN=knnm,
              LDA=ldam,
              NNet=nnm,
              NNet.Multinom=nnmm,
              Rand.Forest=rfm,
              SVM.Radial=svm.rad,
              SVM.Linear=svm.lin))
}
