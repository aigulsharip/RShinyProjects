############################################
# A BMI Calculator is a tool that computes the Body Mass Index (BMI) of an 
# individual. BMI is a measure derived from a person’s weight and height and 
# is used as a quick and simple way to assess whether a person is underweight, 
# normal weight, overweight, or obese.

############################################
install.packages("markdown")
library(shiny)
library(shinythemes)
library(markdown)


####################################
# Server                           #
####################################


shinyServer(function(input, output, session) {
  
  # Input Data
  datasetInput <- reactive({  
    
    bmi <- input$weight/( (input$height/100) * (input$height/100) )
    bmi <- data.frame(bmi)
    names(bmi) <- "BMI"
    print(bmi)
    
  })
  
  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitbutton>0) { 
      isolate("Calculation complete.") 
    } else {
      return("Server is ready for calculation.")
    }
  })
  
  # Prediction results table
  output$tabledata <- renderTable({
    if (input$submitbutton>0) { 
      isolate(datasetInput()) 
    } 
  })
  
}
)
