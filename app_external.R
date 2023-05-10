# if you don't have these libraries already installed,
# you may need to run install.packages("shiny"), install.packages("ggplot2"),
# and install.packages("dplyr") first.
library(shiny)
library(ggplot2)
library(dplyr)
library(cowplot)
library(magrittr)

#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

# If you are using your own personal training dataset, overwrite the file exampleTrainingDat.csv
# with your own training dataset. Make sure that it has the same column names, all values are separated
# by commas (not spaces or tabs), and all values are numeric (except for the column headings).
# No further changes necessary.
#setwd('/users/rsun3/desktop/kras_calculator-main')
#personalTrainingDat <- read.csv("exampledat.csv", header=T)
personalTrainingDat <- data.frame(bKRAS = rep(0, 4), bNRAS = 0.6, tKRAS = FALSE, tNRAS = c(TRUE, FALSE, FALSE, FALSE),
                                  maxAT = c(0.6, 76.9, 10.6, 14.9))
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#


MLEfun <- function(x, samp) {
  s1 <- x[1]
  s2 <- x[2]
  n <- length(samp)
  ret1 <- sum(log(samp)) - n * (-digamma(s1 + s2) + digamma(s1))
  ret2 <- sum(log(1 - samp)) - n * (-digamma(s1 + s2) + digamma(s2))
  return(c(ret1, ret2))
}


findMLE <- function(dat) {
  solverOut <- nleqslv::nleqslv(c(1, 1), fn=MLEfun, samp=dat)
  neg1 <- ifelse(solverOut$x[1] < 0, TRUE, FALSE)
  neg2 <- ifelse(solverOut$x[2] < 0, TRUE, FALSE)
  
  # negative parameters?
  if (neg1 & neg2) {
    # can't do anything
    return(-1)
  } else if (neg1) {
    solverOut$x[1] <- solverOut$x[2] * mean(dat) / (1 - mean(dat))
  } else if (neg2) {
    solverOut$x[2] <- solverOut$x[1] * (1 - mean(dat)) / mean(dat)
  }
  
  globalMax <- TRUE
  secondD <- c(-(trigamma(solverOut$x[1]) - trigamma(solverOut$x[1] + solverOut$x[2])),
               -(trigamma(solverOut$x[2]) - trigamma(solverOut$x[1] + solverOut$x[2])))
  if (secondD[1] > 0 | secondD[2] > 0 | neg1 | neg2) {
    globalMax <- FALSE
  }
  
  return(list(solverOut = solverOut, params = solverOut$x, secondD = secondD, globalMax = globalMax,
              neg1 = neg1, neg2 = neg2))
}


forceMono <- function(dat) {
  dat <- dat %>% arrange((maxAT))
  
  prevPost <- 10^10
  looking <- FALSE
  t0 <- NULL
  for (row_it in 1:nrow(dat)) {
    currentPost <- dat$posterior[row_it]
    
    if (currentPost <= prevPost & looking == FALSE) {
      # everything is ok
      prevPost <- currentPost
    } else if (currentPost <= prevPost & looking) {
      # found t2
      prevPost <- currentPost
      numToInterpolate <- row_it - t0
      # interpolate
      dat$posterior[t0:row_it] <-  dat$posterior[t0] + 
        (currentPost - dat$posterior[t0]) * 0:numToInterpolate / numToInterpolate
      # solved
      looking <- FALSE
    } else if (currentPost > prevPost & looking == FALSE) {
      # start of issue
      t0 <- row_it - 1
      looking <- TRUE
    } else if (currentPost > prevPost & looking) {
      # in the middle of looking
      next
    }
  }
  
  # if at the end and still looking
  if (looking) {
    numToInterpolate <- row_it - t0
    # interpolate
    dat$posterior[t0:row_it] <-  dat$posterior[t0] + 
      (0 - dat$posterior[t0]) * 0:numToInterpolate / numToInterpolate
  }
  return(dat)
}

# trainDat should have a column tKRAS 0 or 1
# also should have maxAT column
trainTestFunction <- function(trainDat, inputFreq, mu1, params0, params1) {

  # only for personal dataset
  if (is.null(params0)) {
    # train for parameters
    train0out <- findMLE(dat = trainDat %>% filter(tKRAS == FALSE & maxAT > 0) %>% select(maxAT) %>% unlist(.) / 100)
    train1out <- findMLE(dat = trainDat %>% filter(tKRAS == TRUE & maxAT > 0) %>% select(maxAT) %>% unlist(.) / 100)
    if (class(train0out) == "numeric" | class(train1out) == "numeric") {return(-1)}
    
    params0 <- train0out$params
    params1 <- train1out$params
  }
  
  # mu0
  mu0 <- 1 - mu1
  
  # for all possible posteriors
  allPost <- data.frame(maxAT = 1:99) %>%
    mutate(joint0 = mu0 * dbeta(maxAT / 100, shape1 = params0[1],
                                shape2 = params0[2])) %>%
    mutate(joint1 = mu1 * dbeta(maxAT / 100, shape1 = params1[1],
                                shape2 =params1[2])) %>%
    mutate(posterior = joint1 / (joint0 + joint1))
  
  # now make it monotone
  allPost$posterior[which(is.na(allPost$posterior))] <- 0
  allPostMono <- forceMono(allPost) %>%
    mutate(diff = abs(maxAT - inputFreq)) %>%
    arrange(diff)

  return(allPostMono$posterior[1])
}




ui <- fluidPage(
  titlePanel("KRAS False Negative Calculator"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("mafInput", "Your MAF", min = 1, max = 100, value = 1, step = 1),
      radioButtons("trainDataInput", "Training Dataset",
                   choices = c("MDACC", "MGH", "Combined", "Personal"),
                   selected = "MDACC"),
      radioButtons("geneInput", "Gene",
                  choices = c("KRAS", "KRAS and NRAS"),
                  selected = "KRAS")
    ),
    mainPanel(
      plotOutput("posteriorPlot"),
      br(), br(),
      textOutput("calculatedPosterior")
    )
  )
)

server <- function(input, output) {

  toTrain <- reactive({
    if (input$geneInput == "KRAS") {
      if (input$trainDataInput == "MGH") {
        params0 <- c(0.4295405, 1.9639500)
        params1 <- c( 0.3271157, 1.4410169)
        mu1 <- 11 / 194
        chosenDataset <- NULL
      } else if (input$trainDataInput == "MDACC") {
        params0 <- c(0.4127479, 2.2654186)
        params1 <- c(0.6966429, 42.0260599)
        mu1 <- 29 / 237
        chosenDataset <- NULL
      } else if (input$trainDataInput == "Combined") {
        params0 <- c(0.4193147, 2.1048697)
        params1 <- c(0.3179729, 3.4899846)
        mu1 <- 40 / 431
        chosenDataset <- NULL
      } else if (input$trainDataInput == "Personal") {
        chosenDataset <- personalTrainingDat %>% select(tKRAS, tNRAS, bKRAS, bNRAS, maxAT)
        mu1 <- (length(which(chosenDataset$tKRAS == 1))) / nrow(chosenDataset)
        params0 <- NULL; params1 <- NULL
      }
    } else if (input$geneInput == "KRAS and NRAS") {
      if (input$trainDataInput == "MGH") {
        params0 <- c(0.4275095, 1.9290964)
        params1 <- c(0.6285794, 17.9724951)
        mu1 <-  7 / 183
        chosenDataset <- NULL
      } else if (input$trainDataInput == "MDACC") {
        params0 <- c(0.4600886, 1.5114152)
        params1 <- c(0.3253474, 4.3855413)
        mu1 <- 33 / 223
        chosenDataset <- NULL
      } else if (input$trainDataInput == "Combined") {
        params0 <- c(0.4335512, 1.7617013)
        params1 <- c(0.3407908, 4.8926177)
        mu1 <- 40 / 406
        chosenDataset <- NULL
      } else if (input$trainDataInput == "Personal") {
        chosenDataset <- personalTrainingDat %>% filter(bNRAS == 0) %>%
          mutate(tissueKRAS = tKRAS) %>%
          mutate(tKRAS = ifelse(tKRAS == TRUE | tNRAS == TRUE, TRUE, FALSE))
        mu1 <- (length(which(chosenDataset$tKRAS == 1))) / nrow(chosenDataset)
        params0 <- NULL; params1 <- NULL
      }
    }
    return(list(chosenDataset=chosenDataset, mu1=mu1, params0=params0, params1=params1))
  })

  onePosterior <-  reactive({
    # run the training and testing
    trainTestFunction(trainDat=toTrain()$chosenDataset, mu1=toTrain()$mu1, inputFreq=input$mafInput,
                      params0=toTrain()$params0, params1=toTrain()$params1)
  })

  output$posteriorPlot <- renderPlot({
    if (is.null(toTrain())) {
      return()
    }
    pieChartDF <- data.frame(group=c("In Tissue", "Not In Tissue"),
                             prop=c(as.numeric(onePosterior()) * 100,
                                    (1 - as.numeric(onePosterior())) * 100)) %>%
      mutate(ypos = cumsum(prop) - 0.5 * prop) %>%
      mutate(textLabel = paste0(round(prop, 0), "%"))
    ggplot(data = pieChartDF, aes(x="", y=prop, fill=group)) +
      geom_bar(stat="identity", width=1, color="white") +
      scale_fill_manual(values=c("orange", "purple"), name="KRAS") +
      coord_polar("y", start=0) +
      theme_void() +
      geom_text(aes(label = textLabel), size=6, fontface = "bold", color="white", position = position_stack(vjust = 0.5))

  })

  output$calculatedPosterior <- renderText({
    paste0("Your posterior probability is ",
           round(as.numeric(onePosterior()), 4) * 100, "%.")
  })

}

shinyApp(ui = ui, server = server)
