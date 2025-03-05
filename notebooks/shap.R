tc <- caret::trainControl(method = "repeatedcv",
                   number = 10, repeats = 100)

m <- caret::train(mpg ~ .,
           data = mtcars,
           method = "lm",
           trControl = tc)

m

