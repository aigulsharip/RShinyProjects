# BMI Calculator

This is a Shiny web application that calculates the Body Mass Index (BMI) of an individual. BMI is a measure derived from a person’s weight and height and is used as a quick and simple way to assess whether a person is underweight, normal weight, overweight, or obese.

## Features

- **Input Parameters**:
  - Age (years)
  - Gender (Male/Female)
  - Height (cm)
  - Weight (kg)

- **Dynamic Outputs**:
  - BMI value with categorized status (e.g., Normal, Overweight, Obese)
  - A BMI gauge visualization.
  - Recommendations on how much weight to gain or lose to reach a healthy BMI range.

- **About Tab**:
  - Provides additional information about BMI and its calculation.

## Screenshot

![BMI Calculator Screenshot](Screenshot%202025-01-20%20at%2010.10.04.png)

## How to Access the App

The app is live and available at: [https://asharip.shinyapps.io/BMICalculator/](https://asharip.shinyapps.io/BMICalculator/)

## How to Run Locally

1. Clone the repository:
   ```bash
   git clone <repository_url>
   ```

2. Open the R project or script in RStudio.

3. Install required packages (if not already installed):
   ```R
   install.packages(c("shiny", "shinythemes", "ggplot2", "plotly", "markdown"))
   ```

4. Run the application:
   ```R
   shiny::runApp()
   ```

5. Access the app in your browser at: [http://localhost:3838](http://localhost:3838)

## Dependencies

This app requires the following R packages:
- `shiny`
- `shinythemes`
- `ggplot2`
- `plotly`
- `markdown`


