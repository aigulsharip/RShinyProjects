# Install required packages (run once if not already installed)
install.packages("RCurl")
install.packages("randomForest")

# Import libraries
library(shiny)
library(shinythemes)
library(RCurl)
library(randomForest)

# Read data
weather <- read.csv(text = getURL("https://raw.githubusercontent.com/dataprofessor/data/master/weather-weka.csv"))

# Convert the target variable 'play' to a factor for classification
weather$play <- as.factor(weather$play)
weather$outlook <- as.factor(weather$outlook)
weather$windy <- as.factor(weather$windy)
# Build the Random Forest model
model <- randomForest(play ~ ., data = weather, ntree = 500, mtry = 4, importance = TRUE)

####################################
# User Interface (UI)              #
####################################
ui <- fluidPage(theme = shinytheme("united"),
                
                # Page header
                headerPanel('Play Golf?'),
                
                # Input values
                sidebarPanel(
                  HTML("<h3>Input Parameters</h3>"),
                  
                  selectInput("outlook", label = "Outlook:", 
                              choices = list("Sunny" = "sunny", "Overcast" = "overcast", "Rainy" = "rainy"), 
                              selected = "Rainy"),
                  sliderInput("temperature", "Temperature:",
                              min = 64, max = 86,
                              value = 70),
                  sliderInput("humidity", "Humidity:",
                              min = 65, max = 96,
                              value = 90),
                  selectInput("windy", label = "Windy:", 
                              choices = list("Yes" = "TRUE", "No" = "FALSE"), 
                              selected = "TRUE"),
                  
                  actionButton("submitbutton", "Submit", class = "btn btn-primary")
                ),
                
                # Output panel
                mainPanel(
                  tags$label(h3('Status/Output')), # Status/Output Text Box
                  verbatimTextOutput('contents'),
                  tableOutput('tabledata') # Prediction results table
                )
)

####################################
# Server Logic                     #
####################################
server <- function(input, output, session) {
  
  # Reactive input data for prediction
  datasetInput <- reactive({  
    
    # Create a data frame for the input values
    test <- data.frame(
      outlook = factor(input$outlook, levels = levels(weather$outlook)),
      temperature = as.numeric(input$temperature),
      humidity = as.numeric(input$humidity),
      windy = factor(input$windy, levels = levels(weather$windy))
    )
    
    # Generate predictions and probabilities
    Output <- data.frame(
      Prediction = predict(model, test),
      predict(model, test, type = "prob")
    )
    Output
  })
  
  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitbutton > 0) { 
      isolate("Calculation complete.") 
    } else {
      return("Server is ready for calculation.")
    }
  })
  
  # Prediction results table
  output$tabledata <- renderTable({
    if (input$submitbutton > 0) { 
      isolate(datasetInput()) 
    } 
  })
}

####################################
# Create the Shiny App             #
####################################
shinyApp(ui = ui, server = server)
