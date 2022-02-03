# if you don't have these libraries already installed,
# you may need to run install.packages("shiny"), install.packages("ggplot2"),
# and install.packages("dplyr") first.
library(shiny)
library(ggplot2)
library(dplyr)
library(cowplot)
library(magrittr)

# to run, put all the files in a folder and change this to be the path to the folder
setwd('/users/rsun3/desktop/krasCalculator')

# load internal data
clia_clean <- read.table("datInternal.txt", header=T, sep=",") %>% 
  filter(!is.na(tKRAS)) %>% 
  mutate(maxAT = pmax(bAPC, bTP53)) %>%
  mutate(Sub = 1:nrow(.))

# load external data
MGHdat <- read.table("datExternal.txt", header=T, sep=',') 

# clean, just the genes we care about
MGHblood <- MGHdat %>%
  filter(!(Alteration == "Not Available" & CDNA == "Not Available" & Gene == "Not Available" & 
             Percentage == "Not Available" & Mutation.Type == "Not Available")) %>%
  filter(Alteration != "AMP") %>%
  mutate(Percentage = as.numeric(as.character(Percentage))) %>%
  tidyr::pivot_wider(., id_cols = Accession.ID, names_from = Gene, values_from = c(Percentage),
                     values_fn = list(Percentage = max), values_fill = list(Percentage = 0)) %>%
  select(Accession.ID, KRAS, BRAF, NRAS, APC, TP53, EGFR, NRAS) 

# get the max information
MGHmax <- MGHdat %>%
  filter(!(Alteration == "Not Available" & CDNA == "Not Available" & Gene == "Not Available" & 
             Percentage == "Not Available" & Mutation.Type == "Not Available")) %>%
  filter(Alteration != "AMP") %>%
  mutate(Percentage = as.numeric(as.character(Percentage))) %>%
  group_by(Accession.ID) %>%
  summarize(All = max(Percentage))

# get the tissue mutation information
MGHtissue <- MGHdat %>% select(Accession.ID, Tissue.Specific.KRAS.NRAS.mutations) %>%
  mutate(tKRAS = ifelse(substr(Tissue.Specific.KRAS.NRAS.mutations, 1, 4) == "KRAS" | 
                          substr(Tissue.Specific.KRAS.NRAS.mutations, 1, 5) == " KRAS", 1, 0)) %>%
  mutate(tNRAS = ifelse(substr(Tissue.Specific.KRAS.NRAS.mutations, 1, 4) == "NRAS" | 
                          substr(Tissue.Specific.KRAS.NRAS.mutations, 1, 5) == " NRAS", 1, 0)) %>%
  filter(tKRAS == 1 | tNRAS == 1)

# put tissue and blood together
# 295 after cleaning
MGHclean <- merge(MGHtissue, MGHblood, by="Accession.ID", all.y=TRUE) %>%
  mutate(tKRAS = ifelse(is.na(tKRAS), 0, 1)) %>%
  mutate(tNRAS = ifelse(is.na(tNRAS), 0, tNRAS)) %>%
  mutate(Sub = 1:nrow(.)) %>%
  merge(., MGHmax, by="Accession.ID") %>%
  select(-Tissue.Specific.KRAS.NRAS.mutations) %>%
  mutate(maxAT = pmax(APC, TP53)) 

# just the blood KRAS = 0
clia_blood0 <- clia_clean %>% filter(bKRAS == 0) 
clia_kn <- clia_clean %>% filter(bKRAS == 0 & bNRAS == 0)
mgh_blood0 <- MGHclean %>% filter(KRAS == 0) 
mgh_kn <- MGHclean %>% filter(KRAS == 0 & NRAS == 0) 


# x is a value of max_AT, xpts are the max ends of the density bins, yvals are the 
# value of the density in that bin.
# expects that bins are closed on the right.
# xpts will have one more element than yvals, the -1
findDens <- function(x, xpts, yvals) {
  densIdx <- max(which(xpts < x))
  return(yvals[densIdx])
}

# give this function a vector of bins and the data, will spit out density.
# the first element in the binVec is the bottom of the first bin, second 
# element is top of first bin/bottom of second, etc.
# the bin is closed at the top, open at the bottom.
# then we also do some adjustment so there are no 0 densities.
calculateDens <- function(binVec, dat) {
  
  # first pass, just standard histogram estimation
  densVec <- rep(NA, length(binVec) - 1)
  n <- length(dat)
  for (dens_it in 1:length(densVec)) {
    densVec[dens_it] <- length(which(dat > binVec[dens_it] & dat <= binVec[dens_it + 1])) / n
  }
  
  # adjustment - if there are zeros, look to the left and right for closest non-zero term, take the
  # min of those.
  # if the first bin is zero, look to right for first non-zero
  if (densVec[1] == 0) {
    tempIdx <- min(which(densVec > 0))
    densVec[1] <- densVec[tempIdx]
  }
  # rest
  for (dens_it in 2:length(densVec)) {
    if (densVec[dens_it] == 0) {
      prevDens <- densVec[dens_it - 1]
      # keep going to the right until you find a non-zero
      tempFill <- c(dens_it)
      temp_it <- dens_it
      while (temp_it < length(densVec) & densVec[temp_it + 1] == 0) {
        temp_it <- temp_it + 1
        tempFill <- c(tempFill, temp_it)
      }
      
      # get the right side value
      # if we hit the right side, just give it fill value of 0
      if (temp_it == length(densVec)) {
        nextDens <- 0
      } else {nextDens <- densVec[temp_it + 1]}
      
      # take the min
      fillValue <- min(prevDens, nextDens)
      
      # fill
      densVec[tempFill] <- fillValue
    } # end if densVec[dens_it] == 0
  }  # end looping through densVec
  
  # now everything has higher density, so divide by total to get sum of 1
  newSum <- sum(densVec)
  returnVec <- densVec / newSum
  
  return(list(returnVec = returnVec, newSum = newSum))
}

# trainDat should have a column tKRAS 0 or 1
# also should have maxAT column
trainTestFunction <- function(binSep = 2.5, trainDat, testDat, mu1=NULL, legTitle=NULL) {
  
  # make bins
  myBins <- seq(from=0, to=100, by=binSep)
  myBins[1] <- -1
  
  # train
  # KRAS not in tissue, should be more density at higher values
  trainKRASneg <- trainDat %>% 
    filter(tKRAS == 0) 
  negDens <- calculateDens(binVec = myBins, dat=trainKRASneg$maxAT)$returnVec
  
  # KRAS in tissue, should be no density at higher values
  trainKRASpos <- trainDat %>% 
    filter(tKRAS == 1) 
  posDens <-  calculateDens(binVec = myBins, dat=trainKRASpos$maxAT)$returnVec
  
  # mu0 and mu1
  if (is.null(mu1)) {
    mu1 <- (length(which(trainDat$tKRAS == 1))) / nrow(trainDat)
  }
  mu0 <- 1 - mu1
  
  # testing 
  testDat <- testDat %>% 
    mutate(numerator = rep(mu1, nrow(.)) * sapply(X=as.list(maxAT), FUN = findDens, 
                                                  xpts=myBins, yvals=posDens)) %>%
    mutate(denominator = numerator + rep(mu0, nrow(.)) * sapply(X=as.list(maxAT), FUN = findDens,
                                                                xpts=myBins, 
                                                                yvals=negDens)) %>%
    mutate(posterior = numerator / denominator) %>%
    # sometimes there is 0/0 if you're way out in the right tail
    mutate(posterior = ifelse(is.na(posterior), 0, posterior)) %>%
    mutate(Sub = 1:nrow(.))
  
  # plot posterior for each subject
  legTitle = ifelse(is.null(legTitle), "KRAS", legTitle)
  returnPlot <- ggplot(data=testDat, aes(x=Sub, y=posterior)) + 
    geom_point(aes(color=as.factor(tKRAS))) + 
    labs(y="Posterior Prob. False Negative", x="Subject ID", color="KRAS") +
    scale_color_manual(labels=c("Not in Tissue", "In Tissue"), values=c("purple", "orange"), 
                       name=legTitle) +
    #ggtitle(paste0("Bin Size ", binSep, "%"))
    ggtitle("") + 
    theme_cowplot()
  
  # get a table summarizing performance
  maxProb <- ceiling(max(testDat$posterior) / 0.05) * 0.05
  if (maxProb == 0) {
    summaryTab <- data.frame(windowLower = 0, windowUpper = 0,
                             nWindow = NA, falseWindow = NA, empiricalProb = NA)
  } else {
    summaryTab <- data.frame(windowLower = seq(from=0, to=maxProb - 0.05, by=0.05),
                             windowUpper = seq(from=0.05, to=maxProb, by=0.05),
                             nWindow = NA, falseWindow = NA, empiricalProb = NA)
  }
  
  for (row_it in 1:nrow(summaryTab)) {
    tempLow <- summaryTab$windowLower[row_it]
    if (tempLow == 0) {tempLow <- -0.1}
    tempHigh <- summaryTab$windowUpper[row_it]
    tempDat <- testDat %>% filter(posterior > tempLow & posterior <= tempHigh)
    summaryTab$nWindow[row_it] <- nrow(tempDat)
    summaryTab$falseWindow[row_it] <- sum(tempDat$tKRAS)
    summaryTab$empiricalProb[row_it] <- ifelse(nrow(tempDat) == 0, 0, sum(tempDat$tKRAS) / nrow(tempDat))
  }
  
  return(list(testDat = testDat, returnPlot = returnPlot, summaryTab = summaryTab))
}


ui <- fluidPage(
  titlePanel("KRAS False Negative Calculator"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("mafInput", "Your MAF", min = 0, max = 100, value = 1, step = 0.1),
      radioButtons("trainDataInput", "Training Dataset",
                   choices = c("CLIA", "MGH", "Combined"),
                   selected = "CLIA"),
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
    if (input$trainDataInput == "MGH") {
      chosenDataset <- mgh_blood0 %>% select(tKRAS, tNRAS, KRAS, NRAS, maxAT) %>%
               set_colnames(c("tKRAS", "tNRAS", "bKRAS", "bNRAS", "maxAT"))
    } else if (input$trainDataInput == "CLIA") {
      chosenDataset <- clia_blood0 %>% select(tKRAS, tNRAS, bKRAS, bNRAS, maxAT)
    } else if (input$trainDataInput == "Combined") {
      chosenDataset <- rbind(mgh_blood0 %>% select(tKRAS, tNRAS, KRAS, NRAS, maxAT) %>%
                     set_colnames(c("tKRAS", "tNRAS", "bKRAS", "bNRAS", "maxAT")),
             clia_blood0 %>% select(tKRAS, tNRAS, bKRAS, bNRAS, maxAT))
    }
    
    if (input$geneInput == "KRAS") {
      return(chosenDataset)
    } else {
      chosenDataset <- chosenDataset %>% filter(bNRAS == 0) %>% 
        mutate(tissueKRAS = tKRAS) %>%
        mutate(tKRAS = ifelse(tKRAS == TRUE | tNRAS == TRUE, TRUE, FALSE))
      return(chosenDataset)
    }
  })
  
  onePosterior <-  reactive({
    newDat <- data.frame(maxAT = rep(as.numeric(input$mafInput), 5), tKRAS = FALSE)
    trainTestFunction(binSep = 1, trainDat=toTrain(), testDat=newDat)$testDat
  })
  
  output$posteriorPlot <- renderPlot({
    if (is.null(toTrain())) {
      return()
    }
    pieChartDF <- data.frame(group=c("In Tissue", "Not In Tissue"), 
                             prop=c(as.numeric(onePosterior()$posterior[1]) * 100, 
                                    (1 - as.numeric(onePosterior()$posterior[1])) * 100)) %>%
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
           round(as.numeric(onePosterior()$posterior[1]), 4) * 100, "%.")
  })
  
}

shinyApp(ui = ui, server = server)