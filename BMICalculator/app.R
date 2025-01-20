############################################
# A BMI Calculator is a tool that computes the Body Mass Index (BMI) of an 
# individual. BMI is a measure derived from a person’s weight and height and 
# is used as a quick and simple way to assess whether a person is underweight, 
# normal weight, overweight, or obese.

############################################

library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(markdown)

# Function to categorize BMI
BMI_INTERVALS <- list(
  `Severe Thinness` = c(-Inf, 16),
  `Moderate Thinness` = c(16, 17),
  `Mild Thinness` = c(17, 18.5),
  `Normal` = c(18.5, 25),
  `Overweight` = c(25, 30),
  `Obesity I` = c(30, 35),
  `Obesity II` = c(35, 40),
  `Obesity III` = c(40, Inf)
)

categorize_bmi <- function(bmi) {
  for (name in names(BMI_INTERVALS)) {
    range <- BMI_INTERVALS[[name]]
    if (bmi >= range[1] && bmi < range[2]) {
      return(name)
    }
  }
}

# Function to calculate weight change
calculate_weight_change <- function(bmi, height, weight) {
  healthy_min_weight <- 18.5 * ((height / 100) ^ 2)
  healthy_max_weight <- 25 * ((height / 100) ^ 2)
  
  if (bmi < 18.5) {
    return(paste0("Gain at least ", round(healthy_min_weight - weight, 1), " kg"))
  } else if (bmi > 25) {
    return(paste0("Lose at least ", round(weight - healthy_max_weight, 1), " kg"))
  } else {
    return("Your weight is healthy")
  }
}

# Function to create BMI gauge
create_bmi_gauge <- function(bmi) {
  plot_ly(
    type = "indicator",
    mode = "gauge+number+delta",
    value = bmi,
    title = list(text = "BMI Gauge"),
    gauge = list(
      axis = list(range = c(10, 40)),
      steps = list(
        list(range = c(10, 18.5), color = "#FFDDC1"),
        list(range = c(18.5, 25), color = "#C1FFC1"),
        list(range = c(25, 30), color = "#FFFAC1"),
        list(range = c(30, 40), color = "#FFC1C1")
      ),
      threshold = list(
        line = list(color = "red", width = 4),
        value = bmi
      )
    )
  )
}

get_category_color <- function(category) {
  if (category == "Normal") {
    return("green")
  } else if (category %in% c("Overweight", "Mild Thinness")) {
    return("orange")
  } else {
    return("red")
  }
}

# User Interface
ui <- fluidPage(
  theme = shinytheme("spacelab"),
  navbarPage("BMI Calculator",
             
             tabPanel("Home",
                      sidebarPanel(
                        HTML("<h3>Input parameters</h3>"),
                        numericInput("age", 
                                     label = "Age (years)", 
                                     value = 25, 
                                     min = 2, 
                                     max = 120),
                        selectInput("gender", 
                                    label = "Gender",
                                    choices = c("Male", "Female"),
                                    selected = "Male"),
                        sliderInput("height", 
                                    label = "Height (cm)", 
                                    value = 175, 
                                    min = 40, 
                                    max = 250),
                        sliderInput("weight", 
                                    label = "Weight (kg)", 
                                    value = 70, 
                                    min = 20, 
                                    max = 200),
                        actionButton("submitbutton", 
                                     "Calculate", 
                                     class = "btn btn-primary")
                      ),
                      
                      mainPanel(
                        tags$h3("Result"),
                        uiOutput("bmi_output"),
                        plotlyOutput("bmi_gauge"),
                        tags$br(),
                        tags$div(
                          uiOutput("weight_change_output")
                        )
                      )
             ),
             
             tabPanel("About", 
                      titlePanel("About"), 
                      div(includeMarkdown("about.md"), align = "justify")
             )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  datasetInput <- reactive({  
    bmi <- input$weight / ((input$height / 100) ^ 2)
    category <- categorize_bmi(bmi)
    weight_change <- calculate_weight_change(bmi, input$height, input$weight)
    
    data.frame(
      Height = paste0(input$height, " cm"),
      Weight = paste0(input$weight, " kg"),
      Age = input$age,
      Gender = input$gender,
      BMI = round(bmi, 1),
      Category = category,
      WeightChange = weight_change
    )
  })
  
  # Display Text Output
  output$bmi_output <- renderUI({
    if (input$submitbutton > 0) {
      isolate({
        data <- datasetInput()
        category <- data$Category
        color <- get_category_color(category)
        
        HTML(paste0(
          "<h4><b>BMI = ", data$BMI, " kg/m<sup>2</sup></b> (", 
          "<span style='color:", color, "'>", data$Category, "</span>)</h4>"
        ))
      })
    } else {
      HTML("<p>Please enter your height, weight, age, and gender, and click 'Calculate'.</p>")
    }
  })
  
  # Weight Change Output
  output$weight_change_output <- renderUI({
    if (input$submitbutton > 0) {
      isolate({
        data <- datasetInput()
        HTML(paste0(
          "<ul>",
          "<li><b>Healthy BMI range:</b> 18.5 - 25 kg/m<sup>2</sup></li>",
          "<li><b>Healthy weight for your height:</b> ", round(18.5 * ((input$height / 100) ^ 2), 1), " kg - ", 
          round(25 * ((input$height / 100) ^ 2), 1), " kg</li>",
          "<li><b>", data$WeightChange, "</b></li>",
          "</ul>"
        ))
      })
    }
  })
  
  # BMI Gauge Visualization
  output$bmi_gauge <- renderPlotly({
    if (input$submitbutton > 0) {
      isolate({
        bmi <- datasetInput()$BMI
        create_bmi_gauge(bmi)
      })
    }
  })
}

# Create Shiny App
shinyApp(ui = ui, server = server)
