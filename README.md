# R Shiny Projects

This repository contains a collection of R Shiny applications. Each app demonstrates specific features and functionality of Shiny for creating interactive web applications in R.

## Applications

### 1. **FirstApp**
- **Description**: A simple app that takes two inputs (Name and Surname) and returns the full name.
- **Features**:
  - Text input fields for Name and Surname.
  - Displays the concatenated full name as output.
- **Screenshot**:
  
<img width="1723" alt="Screenshot 2025-01-20 at 10 22 38" src="https://github.com/user-attachments/assets/eb3aadb8-3457-45db-a9ab-dda6b51ced54" />

### 2. **Dynamic Histogram**
- **Description**: An app that visualizes the `Ozone` levels from the `airquality` dataset using a histogram.
- **Features**:
  - Allows users to adjust the number of bins dynamically.
  - Interactive histogram updates as the number of bins changes.
- **Screenshot**:
  
  <img width="895" alt="Screenshot 2025-01-20 at 10 29 41" src="https://github.com/user-attachments/assets/80a0ec7b-6910-44f7-a5ca-6524c10b5c88" />


### 3. **Weather Prediction App**
- **Description**: A machine learning-based app that predicts whether to play golf or not based on weather conditions using a random forest algorithm.
- **Features**:
  - Inputs for weather parameters such as temperature, humidity, and wind speed.
  - Displays the prediction (Play Golf / Don't Play Golf).
  - Uses a trained random forest model for prediction.
- **Screenshot**:

  <img width="881" alt="Screenshot 2025-01-20 at 10 32 19" src="https://github.com/user-attachments/assets/9e675a18-6b75-4744-9ca0-13b4f6c37509" />


### 4. **BMI Calculator**
- **Description**: An interactive app that calculates the Body Mass Index (BMI) of an individual and provides recommendations.
- **Features**:
  - Inputs for age, gender, height, and weight.
  - Displays BMI with categorized status (e.g., Normal, Overweight, Obese).
  - A dynamic gauge visualization for BMI.
  - Recommendations for weight adjustments.
- **Screenshot**:
  
<img width="1015" alt="Screenshot 2025-01-20 at 10 10 04" src="https://github.com/user-attachments/assets/048493d8-744c-4f10-beb5-0dc2de55ffe5" />

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/aigulsharip/RShinyProjects.git
   ```

2. Open the desired app's folder in RStudio.

3. Install required packages (if not already installed):
   ```R
   install.packages(c("shiny", "shinythemes", "ggplot2", "plotly", "markdown", "randomForest"))
   ```

4. Run the application:
   ```R
   shiny::runApp()
   ```

5. Access the app in your browser at: [http://localhost:127.0.0.1:3399/](http://localhost:(http://127.0.0.1:3399))

