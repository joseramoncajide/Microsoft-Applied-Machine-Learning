
# Applications of Data Science

# Lab 1

# Exploring Time Series Data

## Overview

Time series models are used in a wide range of applications, particularly for forecasting. 

This notebook contains lab material to familarize yourself with key aspects of time series analysis. You will perform analyses on a time series of California dairy data. Specifically, you will explore the structure of the time series and forecast the monthly production of fresh milk in the state of California. 

This exploration is performed in two steps:

- Explore the characteristics of the time series data.
- Decompose the time series of monthly milk production into trend, seasonal components, and remainder components. 
- Apply time series models to the remainder component of the time series.
- Forecast the production of monthly milk produciton for a 12 month period. 

## What you will need

To complete this lab, you will need the following:

- A web browser and Internet connection 
- An Azure Machine Learning workspace
- The lab files for this lab

**Note** To set up the required environment for the lab, follow the instructions in the Setup document from the Setup folder. 


## Load and Examine the Data

As a first step, ensure that you have uploaded the **cadairydata.csv** file as a new dataset in your Azure Machine Learning workspace, and then execute the code below to load the data set into a dataframe.


    options(repr.plot.width=8, repr.plot.height=6)
    options(jupyter.plot_mimetypes = 'image/png')


    to.POSIXct <- function(year, monthNumber){
      ## Function to create a POSIXct time series 
      ## object from a year.month format
      
      ## Create a character vector from the numeric input
      dateStr <- paste(as.character(year), "-",
                       as.character(monthNumber), "-",
                       "01", sep = "")
      ## Return the POSIXct time series object
      as.POSIXct( strptime(dateStr, "%Y-%m-%d"))
    }
    
    order.month <- function(x){
      ## Function to make Month column an ordered factor.
      x <- substr(x, 1, 3) ## Use just the first three letters
      factor(x, 
             levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
             ordered = TRUE)
    }
    
    read.dairy <- function(dataSet = 'cadairydata.csv'){
      library("AzureML")
      ws <- workspace()
      df <- download.datasets(ws, dataSet)
      df$dateTime <- to.POSIXct(df$Year, df$Month.Number)
      dropCols <- c("X","Year.Month", "Month.Number")
      df$Month <- order.month(df$Month)
      df[ , !(names(df) %in% dropCols)]
    }
    
    dairy <- read.dairy()

The code you just executed performs the following steps:

- The data is read from a dataset in your Azure Machine Learning subscription.
- A new column, of type POSIXct, is created. POSIXct is a flexible R data-time class. The **strtime** function formats a text string for conversion to the date-time class. 
- The **Month** column is converted to an ordered R factor class. 
- Some unnessisary columns are removed from the data frame.

Next, execute the code below and examine the head of the data frame.


    head(dairy)


<table>
<thead><tr><th></th><th scope=col>Year</th><th scope=col>Month</th><th scope=col>Cotagecheese.Prod</th><th scope=col>Icecream.Prod</th><th scope=col>Milk.Prod</th><th scope=col>N.CA.Fat.Price</th><th scope=col>dateTime</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>1995</td><td>Jan</td><td>4.37</td><td>51.595</td><td>2.112</td><td>0.9803</td><td>1995-01-01</td></tr>
	<tr><th scope=row>2</th><td>1995</td><td>Feb</td><td>3.695</td><td>56.086</td><td>1.932</td><td>0.8924</td><td>1995-02-01</td></tr>
	<tr><th scope=row>3</th><td>1995</td><td>Mar</td><td>4.538</td><td>68.453</td><td>2.162</td><td>0.8924</td><td>1995-03-01</td></tr>
	<tr><th scope=row>4</th><td>1995</td><td>Apr</td><td>4.28</td><td>65.722</td><td>2.13</td><td>0.8967</td><td>1995-04-01</td></tr>
	<tr><th scope=row>5</th><td>1995</td><td>May</td><td>4.47</td><td>73.73</td><td>2.227</td><td>0.8967</td><td>1995-05-01</td></tr>
	<tr><th scope=row>6</th><td>1995</td><td>Jun</td><td>4.238</td><td>77.994</td><td>2.124</td><td>0.916</td><td>1995-06-01</td></tr>
</tbody>
</table>



The data frame inlcudes **Year**, **Month**, and the **dateTime** column with calendar information. Other columns show the production of dairy products and a benchmark price for milk fat.

Next, examine the class of the **dateTime** column. 


    class(dairy$dateTime)


<ol class=list-inline>
	<li>'POSIXct'</li>
	<li>'POSIXt'</li>
</ol>



As expected, this column is of class **POSIXct**, which inherits from the R **POSIXt** class. 

Now, examine the first 20 values of the **dateTime** column.


    dairy$dateTime[1:20]


     [1] "1995-01-01 UTC" "1995-02-01 UTC" "1995-03-01 UTC" "1995-04-01 UTC"
     [5] "1995-05-01 UTC" "1995-06-01 UTC" "1995-07-01 UTC" "1995-08-01 UTC"
     [9] "1995-09-01 UTC" "1995-10-01 UTC" "1995-11-01 UTC" "1995-12-01 UTC"
    [13] "1996-01-01 UTC" "1996-02-01 UTC" "1996-03-01 UTC" "1996-04-01 UTC"
    [17] "1996-05-01 UTC" "1996-06-01 UTC" "1996-07-01 UTC" "1996-08-01 UTC"


Note that these values inlcude the day and a time zone.

You can perform both arithmatic and logical operations on POSIXct objects. As an example, select the last 12 months of dairy data by executing the next code cell. 


    dairy[dairy$dateTime > '2012-12-01',]


<table>
<thead><tr><th></th><th scope=col>Year</th><th scope=col>Month</th><th scope=col>Cotagecheese.Prod</th><th scope=col>Icecream.Prod</th><th scope=col>Milk.Prod</th><th scope=col>N.CA.Fat.Price</th><th scope=col>dateTime</th></tr></thead>
<tbody>
	<tr><th scope=row>217</th><td>2013</td><td>Jan</td><td>2.282</td><td>62.612</td><td>3.462</td><td>1.7744</td><td>2013-01-01</td></tr>
	<tr><th scope=row>218</th><td>2013</td><td>Feb</td><td>1.941</td><td>66.274</td><td>3.231</td><td>1.6224</td><td>2013-02-01</td></tr>
	<tr><th scope=row>219</th><td>2013</td><td>Mar</td><td>2.15</td><td>75.77</td><td>3.676</td><td>1.6968</td><td>2013-03-01</td></tr>
	<tr><th scope=row>220</th><td>2013</td><td>Apr</td><td>2.217</td><td>76.396</td><td>3.622</td><td>1.7444</td><td>2013-04-01</td></tr>
	<tr><th scope=row>221</th><td>2013</td><td>May</td><td>2.341</td><td>80.899</td><td>3.72</td><td>1.8552</td><td>2013-05-01</td></tr>
	<tr><th scope=row>222</th><td>2013</td><td>Jun</td><td>2.694</td><td>92.091</td><td>3.489</td><td>1.8314</td><td>2013-06-01</td></tr>
	<tr><th scope=row>223</th><td>2013</td><td>Jul</td><td>2.433</td><td>84.026</td><td>3.373</td><td>1.6923</td><td>2013-07-01</td></tr>
	<tr><th scope=row>224</th><td>2013</td><td>Aug</td><td>2.169</td><td>78.382</td><td>3.4</td><td>1.6273</td><td>2013-08-01</td></tr>
	<tr><th scope=row>225</th><td>2013</td><td>Sep</td><td>2.193</td><td>74.871</td><td>3.197</td><td>1.5519</td><td>2013-09-01</td></tr>
	<tr><th scope=row>226</th><td>2013</td><td>Oct</td><td>2.284</td><td>70.011</td><td>3.338</td><td>1.5566</td><td>2013-10-01</td></tr>
	<tr><th scope=row>227</th><td>2013</td><td>Nov</td><td>1.865</td><td>61.394</td><td>3.275</td><td>1.7537</td><td>2013-11-01</td></tr>
	<tr><th scope=row>228</th><td>2013</td><td>Dec</td><td>2.345</td><td>53.653</td><td>3.473</td><td>1.6747</td><td>2013-12-01</td></tr>
</tbody>
</table>



Please try selecting some other subsets of these data using the **POSIXct** column. 

## Create a Time Series Plot

Now that you have examined the data frame, you will create a time series plot of milk production. When creating a time series plot, the **POSIXct** column is used to create the time axis. 

The plot is created using the ggplot2 package. Execute the code in the cell below to create the time series plot of milk production. 


    dairy.plot <- function(df, col = 'Milk.Prod'){
      require(ggplot2)
      ggplot(df, aes_string('dateTime', col)) +
        geom_line() +
        ggtitle(paste('Time series of', col)) +
        xlab('Time in years')
    }
    dairy.plot(dairy)


![png](output_13_0.png)


The production of milk is shown on the vertical scale and the date is on the horizontal scale. 

For most of the time period shown the production of milk increased year over year. However, there is a decline in milk poduction starting in 2009 as a reults of the recession. Also, notice that this time seies exhibits a strong seasonal component, with an annual cycle. 


## Statistical Properties of the Time Series

Having examined the time series of milk production, you will now explore some statistical properties of the time series. 

Autocorrelation is a fundamental property of time series. The **Autocorrelation Function** or ACF provides information on the dependency of the time seires values of previous values. Later in this lab, you will use the results of a ACF analysis to estimate the order of moving average processes.The **Partial autocorrelation Function** or PACF, measures the correlation of the time series with its own lag values. Later in this lab you will use an  PACF to estimate the order of an autoregressive process. 

Execute the code in the following cell and examine the ACF of the milk produciton time series. 


    dairy.acf <- function(df, col = 'remainder', is.df =TRUE){
      if(is.df) temp <- df[, col]
      else temp <- df
      temp = ts(temp, start = 1995, frequency = 12)
      par(mfrow = c(2,1))
      acf(temp, main = paste('ACF of', col))
      pacf(temp, main = paste('PACF of', col))
      par(mfrow = c(1,1))
    }
    dairy.acf(dairy, col = 'Milk.Prod')


![png](output_15_0.png)


Note that the values of the ACF at the various lags decays only slowly. This indicates there is considerable serial correlation betwen the time series values at the various lags, mostly likely from the trend. 

Plotting a histogram provides information on the distribution of values of the time series. Execute the code in the cell below and examine the histogram.


    hist.ts = function(df, col = 'Milk.Prod', bins = 40){
        temp = df[,col]
        breaks = seq(min(temp), max(temp), length.out = (bins + 1))
        hist(temp, breaks = breaks, main = paste('Distribution of ', col), xlab = col)
    }
    hist.ts(dairy)


![png](output_17_0.png)


The histogram of the full milk production time series shows considerable dispursion. Again this behavior is likely the result of the trend. 


### Exercise: Plot Ice Cream Production

In this exercise you will use the **Icecream.Prod** column of the **diary** data frame. Create plots to answer the following questions. 

- Does icecream production have a noticable seasonal component? Can you characterize the trend of ice cream production greater over time as strong or weak? 
- Is the seasonal variation of icecream production noticable in the plot of ACF? Does the ACF plot indicate a strong trend component? 

In the code cell below, use the **dairy.plot** and **dairy.acf** functions to create new plots showing the ice cream production and answer the questions above. 



    dairy.plot(dairy, col="Icecream.Prod")


![png](output_19_0.png)



    dairy.acf(dairy, col = 'Icecream.Prod')


![png](output_20_0.png)



    hist.ts(dairy, col = 'Icecream.Prod')


![png](output_21_0.png)


## Simple Moving Average Decomposition of the Time Series

Time series are typically decomposed into three components: trend, seasonal, and the remainder, or residual. Trend can be modeled by several methods. You will start by decomposing the time series using a simple moving average model. 

The code in the cell below uses moving window method to compute the average of the time series over specified span, or order of the operator. As the moving window operator moves over the data, the average of the values in the windows. Execute the cell to load the function.


    dairy.ma <- function(df, col = 'Milk.Prod', order = 12){
      temp = df[, col]
      end = length(temp) - 1
      out = rep(0, length(temp))
      out[1] = temp[1]
      for(i in 1:end){
        if(i - order <= 1) j = 1 
        else j = j + 1
        out[i + 1] = sum(temp[j:i])/(i - j + 1)
      }
      out
    }

Once the trend has been removed, the seasonal component must be modeled and removed. The function in the cell below computes the seasonal component as a function of the month of the year using a linear model. The **0** in the model formula supresses the intercept term. Since 12 monthly factors are used to model seasonal variation, the model would be over-determined if an intercept was included. 

Execute the code in the cell to load the function.


    dairy.seasons <- function(df, col = 'Milk.Prod'){
      df$y = df[, col]
      fit = lm(y ~ 0 + Month, data = df)
      predict(fit, newdata = df)
    }

Using these functions you will now decompose the time series into its components. The function in the code cell below uses multiplicative decomposition of the time series. The model is transformed to a multiplicative model by taking the log of the time series values. The moving average is computed over a twelve month moving window. 

Execute the code in the cell below and examine the head of the resulting data frame. 


    decomp.dairy <- function(df,  col = 'Milk.Prod', multiplicative = TRUE, order = 12){
      if(multiplicative) {
        temp = log(df[, col])
        df[, col] = temp
      } else { 
        temp = df[, col] 
      }
      trend = dairy.ma(df, col = col, order = order)
      temp = temp - trend
      df[, col] = temp
      seasonal = dairy.seasons(df, col = col)
      remainder = temp - seasonal
      data.frame(trend = trend, seasonal = seasonal, remainder = remainder)
    }
    decomp <- decomp.dairy(dairy, order = 12)
    head(decomp)


<table>
<thead><tr><th></th><th scope=col>trend</th><th scope=col>seasonal</th><th scope=col>remainder</th></tr></thead>
<tbody>
	<tr><th scope=row>1</th><td>0.7476354</td><td>0.02471208</td><td>-0.02471208</td></tr>
	<tr><th scope=row>2</th><td>0.7476354</td><td>-0.04996016</td><td>-0.03911947</td></tr>
	<tr><th scope=row>3</th><td>0.7030956</td><td>0.06783954</td><td>9.862991e-05</td></tr>
	<tr><th scope=row>4</th><td>0.7257416</td><td>0.03872409</td><td>-0.008343715</td></tr>
	<tr><th scope=row>5</th><td>0.7333367</td><td>0.06733682</td><td>-1.813267e-05</td></tr>
	<tr><th scope=row>6</th><td>0.7468004</td><td>0.01530435</td><td>-0.008803681</td></tr>
</tbody>
</table>



The resulting data frame has three components for trend, seasonal and remainder.

The code in the cell below creates visualizations of each of the three components of the decomposition of the time series using ggplot2. run the code in the cell to visualize these components. 

(**Note**: If an error occurs stating that gridExtra can't be found, save the notebook, halt and close it, reopen it, run all the cells above this one, and then re-run the cell below - occasionally newly installed packages fail to load in a timely fashion) 


    decomp.plot <- function(df){
      require(ggplot2)
      install.packages("gridExtra")
      require(gridExtra)
      df$x = 1:nrow(df)
      ycols = c('trend', 'seasonal', 'remainder')
      p <- lapply(ycols, function(y){
                   ggplot(df, aes_string('x', y)) + 
                             geom_line() +
                             ylab(y)
                })
      grid.arrange(p[[1]], p[[2]], p[[3]], nrow = 3)
    }
    decomp.plot(decomp)

    Installing package into '/home/nbcommon/R'
    (as 'lib' is unspecified)


    
    The downloaded source packages are in
    	'/tmp/RtmpEHD0o8/downloaded_packages'



![png](output_29_2.png)


You can see the trend and seasonal components clearly separated in the above plots. The remainder plot looks fairly random, as expected. But, is the remainder actually stationary?  To test for stationarity of the remainder, plot the ACF by executing the code in the cell below. 


    dairy.acf(decomp)


![png](output_31_0.png)


The ACF has 7 significant lag values, indicating the remainder is not, in fact, stationary. 

## Exploring the Multiplicative Model with lowess

Having tried a simple moving average decomposition, you will now use a lowess model to determine the trend. Lowess is a sophisticated non-linear regression. The lowess trend model is combined with a moving window seasonal component model into the R **stl** function.   

The code in the cell below uses **stl** to decompose a time series. Plots are created of the components of the time series. Finally, the columns of the time series decomposition are added as new column to input data frame. Execute the code in the cell below to compute and view a decomposition of the milk production time series. 


    dairy.decomp <- function(df, col = 'Milk.diff', span = 0.5, Mult = TRUE){
      if(Mult) {temp <- ts(log(df[, col]), frequency=12, start=1)
      } else {temp <- ts(df[, col], frequency=24, start=1)}
      span = span * length(temp)  
      dairyFit <- stl(temp, s.window = "periodic", t.window = span)
      plot(dairyFit, main = 'Decompositon of dairy produciton')
      cbind(df, as.data.frame(dairyFit$time.series))
    }
    dairyMult = dairy.decomp(dairy, col = 'Milk.Prod', span = 0.2)


![png](output_33_0.png)


The time series charts show the original time series along with the components of the decomposition. The trend is a bit smoother than was obtained with the simple moving average decomposition. 

The question remains, is the remainder from this decomposition stationary? To find out, execute the code in the cell below.  


    dairy.acf(dairyMult)


![png](output_35_0.png)


The first 4 lag values of the ACF have significant values, indicating that the remainder series is not stationary. Compared to the behavior of the ACF for the simple moving average decomposition, the behavior of the remainder is improved.  

Now, plot the histogram of the remainder by executing the code in the cell below. 


    hist.ts(dairyMult, col = 'remainder')


![png](output_37_0.png)


The distribution of the remainder values is much closer to a Normal distribution than for the original time series you created earlier. This result, combined with the ACF plot, indicates that the stl decomposition is effective.  

You will further, investigate the remainder (non-seasonal residual) component by making a box plot by month of the year. The code in the cell below plots a box plot by month of the remainder component using ggplot2. Execute the code in the cell below to display the box plot. 


    dairy.box <- function(df, col = 'remainder'){
      require(ggplot2)
      p <- ggplot(df, aes_string('Month', col)) +
        geom_boxplot() +
        ggtitle('Variation of remainder component of dairy production by month')
      print(p)
    }
    dairy.box(dairyMult)


![png](output_39_0.png)


The remainder component shows only limited vairation from month to month. The differences are within the interquartile range, indicating that the seasonal model is a reasonably good fit.  

****

### Exercise: Decomposition of Ice Cream Production Time Series

In this exercise you will decompose the **Icecream.Prod** column of the **diary** data frame. Create plots to answer the following questions. 

- Does icecream production have a noticable seasonal compoent or are the values all close to the average over time? Is there a strong seasonal component for icecream production. 
- Does the acf plot indicate that the remainder series is stationary?
- Do the values in the histogram have an approximately normal distribution?
- Does the interquartile range for each month overlap, indicating that the decomposition has produced a reasonably good model of the seasonal variation. 

In the code cell below, use the **dairy.decomp**, **dairy.acf**, **dairy.hist**, and **dairy.box** functions to create new plots of the ice cream production time series and answer the questions above. 

****



    dairyMult = dairy.decomp(dairy, col = 'Icecream.Prod', span = 0.2)


![png](output_41_0.png)



    dairy.acf(dairyMult)


![png](output_42_0.png)



    hist.ts(dairyMult, col = 'remainder')


![png](output_43_0.png)



    dairy.box(dairyMult)


![png](output_44_0.png)


## Moving Average Models

Now that you have explored the decomposition of the time series you, will now construct and test Autoregressive Moving Average (ARMA) models for the remainder of the time series. You will create and test these models in three steps, creating a moving average (MA) model, creating an autoregressive (AR) model and creating an autoregressive moving average (ARMA) model. 

The function in the cell below computes an Autoregressive Integrative Moving Average (ARIMA) model. The summary statistics for the model are printed and the model object returned. By assigning values to the order of each operator different time series models can be specified, as order of MA model, order of Integrative model, and order of AR model. Since the de-trended remainder is being modeled, the **include.mean** argument is set to FALSE in the **arima** function. 

The ACF of the remainder from the **stl** decomposition of the milk production time series had 4 significant lag values. As an inital model, you will now create an MA model of order 4. Execute the code in the cell below to compute the MA(4) model and examine the model summary. 


    model.dairy = function(df, col = 'remainder', order = c(0,0,1)){
      ts = ts(df[, col], frequency = 12, start = 1995)
      dairy.mod = arima(ts, order = order, include.mean = FALSE)
      print(dairy.mod)
      dairy.mod
    }
    ma1 = model.dairy(dairyMult, order = c(0,0,4))

    
    Call:
    arima(x = ts, order = order, include.mean = FALSE)
    
    Coefficients:
             ma1      ma2     ma3      ma4
          0.1379  -0.0150  0.2598  -0.0905
    s.e.  0.0706   0.0732  0.1018   0.0873
    
    sigma^2 estimated as 0.001397:  log likelihood = 425.73,  aic = -841.46


Examine the values of the model coefficients and their standard errors (SE). Notice that the SE of the ma4 coefficient is actually greater than the value of the coefficient itself. This indicates that the value of this coefficient is poorly deterimind and should likely be set to zero. 

The foregoing result indicates that the order of the MA model should be reduced. Generally, the order of an MA model is reduced in unit steps until all the coefficients appear to be significant. The code in the cell below computes an MA(3) model. Run this code and examine the model summary. 


    ma1 = model.dairy(dairyMult, order = c(0,0,3))

    
    Call:
    arima(x = ts, order = order, include.mean = FALSE)
    
    Coefficients:
             ma1     ma2     ma3
          0.1702  0.0124  0.3159
    s.e.  0.0624  0.0659  0.0687
    
    sigma^2 estimated as 0.001403:  log likelihood = 425.15,  aic = -842.3


The small standard error compared to the magnitude of the coefficients indicates that the order of the model is reasonable. 

To test how well this model fits the data, and results in a stationary result, you will plot the ACF of the residuals of the MA(3) model. Run the code in the cell below to plot the ACF of the model result. 


    dairy.acf(ma1$resid[-1], is.df = FALSE)


![png](output_50_0.png)


Note that only the 0 lag of the ACF is significant and that there are no significant lags for the PACS. These observations indicate that the MA(3) model is a good fit. 

## Autoregressive Models

The MA(3) model has been shown to be effective. Next, you will test an autoregressive model. The PACF of the reminder indicates that an AR model might not be the best choice. None the less, a low order AR model might fit these data. To compute an AR(2) model exectue the code in the cell below. 


    ar1 = model.dairy(dairyMult, order = c(2,0,0))

    
    Call:
    arima(x = ts, order = order, include.mean = FALSE)
    
    Coefficients:
             ar1      ar2
          0.1028  -0.0059
    s.e.  0.0663   0.0666
    
    sigma^2 estimated as 0.001515:  log likelihood = 416.61,  aic = -827.23


Examine the values of the coefficients and their stadard errors. Note that the standard error of the second coefficient is of the same magnitude as the coefficient. Clearly, the AR(2) model is over parameterized. 

Next, you will try an AR(1) model, by executing the code in the cell below.


    ar1 = model.dairy(dairyMult, order = c(1,0,0))

    
    Call:
    arima(x = ts, order = order, include.mean = FALSE)
    
    Coefficients:
             ar1
          0.1022
    s.e.  0.0659
    
    sigma^2 estimated as 0.001515:  log likelihood = 416.61,  aic = -829.22


The standard error of the AR(1) model is an order of magnitude less than the value of the coefficient, which is promising.

Next, exectue the code in the cell below to plot the ACF and PACF of the AR(1) model


    dairy.acf(ar1$resid[-1], is.df = FALSE)


![png](output_57_0.png)


Note that only the 0 lag of the ACF is significant and that there are no significant lags for the PACF. These observations indicate that the AR(1) model is a good fit. Compare these results to those of the MA(3) model, noting that they are nearly identical. Evidently, either the MA(3) or AR(1) model is a good choice for these data. 

## Autoregressive Moving Average Models

You have found that both MA(3) and AR(1) models are good fits to the remainder series. You will now investigate the use of autoregressive moving average (ARMA) models on the remainder series. 

As a starting point you will try an ARMA(1,3) model by executing the code in the cell below. 


    arma1 = model.dairy(dairyMult, order = c(1,0,3))

    
    Call:
    arima(x = ts, order = order, include.mean = FALSE)
    
    Coefficients:
              ar1     ma1     ma2     ma3
          -0.1843  0.3359  0.0325  0.2875
    s.e.   0.2222  0.2133  0.0696  0.0844
    
    sigma^2 estimated as 0.001399:  log likelihood = 425.52,  aic = -841.04


In each case, the standard error is of the same order of magnitude as the value of the coefficient, indicating this model is a poor fit to the data.


****

### Exercise: Test Another ARMA Model

In this exercise you will test an ARMA(1,1) model and evaluate its coefficients to answer these questions. 

- How do the standard errors compare to the values of the coefficients?
- Do you think this model fits the data well or is it over parameterize?

Add code to compute the ARMA(1,1) model and plot iots ACFG and PACF to the cell below, and answer the questions above.

****


    arma2 = model.dairy(dairyMult, order = c(1,0,1))

    
    Call:
    arima(x = ts, order = order, include.mean = FALSE)
    
    Coefficients:
             ar1      ma1
          0.3409  -0.2423
    s.e.  0.4001   0.4103
    
    sigma^2 estimated as 0.001514:  log likelihood = 416.65,  aic = -827.29



    dairy.acf(arma2$resid[-1], is.df = FALSE)


![png](output_63_0.png)



## Exploring the Difference Series

Using a difference series is a method to remove trend from a time series. The difference can be computed for any number of lag values, depending on the order of the trend.  In this case you will use a first order difference series to model the trend in the milk production time series. 

Using the code in the cell below, you will compute the first order difference series. Notice that the difference series is  necessarily of length one less than the original series. Execute the code in the cell to compute the difference series. 


    dairy.diff <- function(df, col = 'Milk.Prod', out = 'Milk.Diff'){
      ln <- nrow(df)
      temp <- ts(df[, col], frequency = 12, start = 1995)
      df[2:ln, out] <- diff(temp)
      df <- df[2:ln, ]
      df
    }
    dairyDiff <- dairy.diff(dairy)

Next, compute the stl decomposition of the difference series. Since we are working with a difference series, which has positive and negative values, we use an additive model. No logarithm is taken. Execute the code in the cell below to compute and display the decomposition of the difference series.  


    dairyDiff = dairy.decomp(dairyDiff, col = 'Milk.Diff', span = 0.5, Mult = FALSE)


![png](output_67_0.png)


Examine the results shown above. The difference series is shown in the upper most plot. Notice the small magnitude of the remaing trend indicating that the first order difference model removed most of the trend. However, the seasonal series exhibits a pattern with a 24 month cycle which is a bit odd. 

****

### Exercise: Analysis of the Difference Series

In this exercise you will analyize the remainder of the diference series following these steps. 

- Does the acf plot indicate that the remainder series is stationary?  
- Is the dispersion of these values greater or less than the dispersion obtained directly from the decomposition of the time series (without diferencing). 
- Is the interquartile range of these values greater or less than the interquartile range obtained directly from the decomposition of the time series (without differencing). 

In the code cell below, use the **dairy.acf**, **dairy.hist**, and **dairy.box** functions to create new plots of the difference series and answer the questions above. 
****



    dairy.acf(dairyDiff)


![png](output_69_0.png)



    hist.ts(dairyDiff, col = 'remainder')


![png](output_70_0.png)



    dairy.box(dairyDiff)


![png](output_71_0.png)


## Autoregressive Integrative Moving Average Model

It is clear from the exploration of the ARMA model, that the remainder of the decomposition of the dairy production time series is not stationary. You will now model the remainder series with an autoregressive integrative moving average (ARIMA) model.

Execute the code in the cell below to compute an ARIMA(1,1,1) model.


    arima1 = model.dairy(dairyMult, order = c(1,1,1))

    
    Call:
    arima(x = ts, order = order, include.mean = FALSE)
    
    Coefficients:
             ar1      ma1
          0.1066  -1.0000
    s.e.  0.0663   0.0114
    
    sigma^2 estimated as 0.001521:  log likelihood = 411.71,  aic = -817.42


The standard error of the ar1 coeficient is only about half its value. This model seems to be a reasonable fit. 

Next, plot the ACF and PACF of the model by executing the code in the cell below. 


    dev.off()
    dairy.acf(arima1$resid[-1], is.df = FALSE)


<strong>png:</strong> 2



    Error in plot.xy(xy, type, ...): invalid graphics state
    Traceback:


    1. dairy.acf(arima1$resid[-1], is.df = FALSE)

    2. acf(temp, main = paste("ACF of", col))   # at line 6 of file <text>

    3. plot.acf(acf.out, ...)

    4. plot(x$lag[, i, j], x$acf[, i, j], type = type, xlab = xlab, 
     .     ylab = if (j == 1) ylab else "", ylim = ylim, ...)

    5. plot.default(x$lag[, i, j], x$acf[, i, j], type = type, xlab = xlab, 
     .     ylab = if (j == 1) ylab else "", ylim = ylim, ...)

    6. plot.xy(xy, type, ...)



![png](output_75_2.png)


Note that only the 0 lag of the ACF is significant and that there are no significant lags for the PACF. These observations indicate that the ARIMA(1,1,1) model is a good fit. Compare these results to those of the MA(3) and AR(1) models, noting that they are nearly identical. The ARIMA(1,1,1) model is a good choice for these data as well. 


## Modeling and Forecasting

Now that you have explored the properties of the decomposed time series, you will now compute forecasts of dairy product production. In this exercise you will compute a time series model using the R **forecast** package, and use this model to forecast the next 12 months of dairy product production. 

The R **forecast** package contains the **auto.arima** function which automatically steps through the ARIMA model parameters to find the best fit to the data. The ARIMA model used in the **forecast** package also includes modeling of seasonal differences.  The **auto.arima** function has multiple arguments, specifying the range of parameter values to search. The first argument is a time series object of class **ts**. 

The code in the cell below does the following:

- Creates a time series of class **ts**.
- Automatically finds and computes an ARIMA model.
- Prints a summary of the ARIMA model.

Execute the code in the cell below to compute and print the summary of the ARIMA model. 


    temp <- ts(log(dairy[, 'Milk.Prod']), frequency=12, start=1995)
    install.packages("forecast")
    library(forecast)
    fit.dairy = auto.arima(temp, max.p=3, max.q=3,
          max.P=2, max.Q=2, max.order=5, max.d=2, max.D=1,
          start.p=0, start.q=0, start.P=0, start.Q=0)
    summary(fit.dairy)

    Installing package into '/home/nbcommon/R'
    (as 'lib' is unspecified)
    also installing the dependencies 'timeDate', 'tseries', 'fracdiff', 'RcppArmadillo'
    


    
    The downloaded source packages are in
    	'/tmp/RtmpEHD0o8/downloaded_packages'


    Loading required package: zoo
    
    Attaching package: 'zoo'
    
    The following objects are masked from 'package:base':
    
        as.Date, as.Date.numeric
    
    Loading required package: timeDate
    This is forecast 6.1 
    


    Series: temp 
    ARIMA(0,1,1)(0,1,2)[12]                    
    
    Coefficients:
              ma1     sma1    sma2
          -0.1506  -0.9076  0.1129
    s.e.   0.0743   0.0794  0.0838
    
    sigma^2 estimated as 0.0002547:  log likelihood=577.8
    AIC=-1147.6   AICc=-1147.41   BIC=-1134.12
    
    Training set error measures:
                            ME       RMSE        MAE         MPE    MAPE      MASE
    Training set -0.0003536657 0.01549906 0.01109068 -0.01955938 1.05342 0.2902694
                        ACF1
    Training set 0.005145456


Examine the summary of the model, noticing the following points.

- The model uses an MA(2) model for the seasonal difference. The coefficients of this model, sma1 and sma2, along with their standard errors can be seen in the summary.
- The model of the remainder is and MA(1) model. The coeficient and its standard error can be seen in the summary above. 
- Error metrics, including RMSE, are provided in the summary. Notice that the RMSE is much smaller than the values of the milk production time series indicating good model performance.

The **forecast** function is used to compute the forecast of the next 12 months using the model created using **auto.arima**. Execute the code in the cell below to compute a forecast and print its summary. 


    dairy.fit = forecast(fit.dairy, h=12)
    summary(dairy.fit)

    
    Forecast method: ARIMA(0,1,1)(0,1,2)[12]                   
    
    Model Information:
    Series: temp 
    ARIMA(0,1,1)(0,1,2)[12]                    
    
    Coefficients:
              ma1     sma1    sma2
          -0.1506  -0.9076  0.1129
    s.e.   0.0743   0.0794  0.0838
    
    sigma^2 estimated as 0.0002547:  log likelihood=577.8
    AIC=-1147.6   AICc=-1147.41   BIC=-1134.12
    
    Error measures:
                            ME       RMSE        MAE         MPE    MAPE      MASE
    Training set -0.0003536657 0.01549906 0.01109068 -0.01955938 1.05342 0.2902694
                        ACF1
    Training set 0.005145456
    
    Forecasts:
             Point Forecast    Lo 80    Hi 80    Lo 95    Hi 95
    Jan 2014       1.259101 1.238648 1.279555 1.227820 1.290382
    Feb 2014       1.193990 1.167154 1.220825 1.152948 1.235031
    Mar 2014       1.303803 1.271835 1.335772 1.254912 1.352695
    Apr 2014       1.278470 1.242086 1.314854 1.222825 1.334115
    May 2014       1.307799 1.267480 1.348118 1.246136 1.369462
    Jun 2014       1.256228 1.212325 1.300130 1.189084 1.323371
    Jul 2014       1.253456 1.206240 1.300671 1.181246 1.325665
    Aug 2014       1.243199 1.192889 1.293509 1.166256 1.320141
    Sep 2014       1.198112 1.144887 1.251337 1.116711 1.279512
    Oct 2014       1.232666 1.176677 1.288655 1.147039 1.318293
    Nov 2014       1.209167 1.150544 1.267789 1.119512 1.298821
    Dec 2014       1.252931 1.191789 1.314074 1.159422 1.346440


Much of the summary is the same as before. A 12 month forecast is printed below the model summary. There is a point forecast (the expected value) along with 80 and 95 percent confidence intervals. Notice, that the confidence intervals generally get wider for forecasts further out in time. It is hardly suprising that the forecast has more uncertaintly as time increases from the present. 

Finally, you can plot the forecast by executing the code in the cell below. 


    plot(dairy.fit)


![png](output_81_0.png)


The orignial time series of milk production is show in black in the plot above. The forecast is show in Blue. The 80 and 95 percent confidence intervals are shown in lighter shades of blue-gray.  

****

### Exercise: Forecast Ice Cream Production

In this exercise you will forecast the production of icecream for a 12 month period following these steps. 

- Create a time series object with **frequency = 12** and **start = 1995** from the **Icecream.Prod** column of the **dairy** data frame.
- Fit a model with the **auto.arima** function, following the hint given below.
- Print a summary of the model. What is the order of the MA and AR components of the seasonal and remainder models? Note there will be a drift term, which accounts for linear trend in the time series.
- Compute the forecast of icecream produciton for the next 12 months. 
- Plot the forecast and note the behavior. 

In the code cell below, use the **ts**, **auto.arima**, **summary**, **forecast**, and **plot** functions to create the new model, print the summaries and make the plots. **Hint** use the time series of icecream production as the first argument to the **auto.arima** function. You can copy and paste the other arguments from the milk production model, but set **max.p = 1** to prevent having an over-parameterized model. 

****



    temp2 <- ts(log(dairy[, 'Icecream.Prod']), frequency=12, start=1995)
    fit.dairy2 = auto.arima(temp2, max.p=3, max.q=3,
          max.P=2, max.Q=2, max.order=5, max.d=2, max.D=1,
          start.p=0, start.q=0, start.P=0, start.Q=0)
    summary(fit.dairy2)

    Series: temp2 
    ARIMA(3,0,1)(0,1,2)[12] with drift         
    
    Coefficients:
              ar1     ar2     ar3     ma1     sma1     sma2  drift
          -0.1996  0.1658  0.3746  0.4086  -0.5043  -0.2039  6e-04
    s.e.   0.1818  0.0777  0.0654  0.1981   0.0720   0.0692  2e-04
    
    sigma^2 estimated as 0.001621:  log likelihood=368.05
    AIC=-720.1   AICc=-719.4   BIC=-693.1
    
    Training set error measures:
                          ME       RMSE        MAE        MPE      MAPE      MASE
    Training set 0.001360524 0.03809705 0.02972694 0.02906281 0.7001261 0.7583233
                         ACF1
    Training set -0.004410836



    dairy.fit2 = forecast(fit.dairy2, h=12)
    summary(dairy.fit2)

    
    Forecast method: ARIMA(3,0,1)(0,1,2)[12] with drift        
    
    Model Information:
    Series: temp2 
    ARIMA(3,0,1)(0,1,2)[12] with drift         
    
    Coefficients:
              ar1     ar2     ar3     ma1     sma1     sma2  drift
          -0.1996  0.1658  0.3746  0.4086  -0.5043  -0.2039  6e-04
    s.e.   0.1818  0.0777  0.0654  0.1981   0.0720   0.0692  2e-04
    
    sigma^2 estimated as 0.001621:  log likelihood=368.05
    AIC=-720.1   AICc=-719.4   BIC=-693.1
    
    Error measures:
                          ME       RMSE        MAE        MPE      MAPE      MASE
    Training set 0.001360524 0.03809705 0.02972694 0.02906281 0.7001261 0.7583233
                         ACF1
    Training set -0.004410836
    
    Forecasts:
             Point Forecast    Lo 80    Hi 80    Lo 95    Hi 95
    Jan 2014       4.159970 4.108371 4.211570 4.081055 4.238885
    Feb 2014       4.233898 4.181183 4.286613 4.153278 4.314519
    Mar 2014       4.354554 4.301452 4.407656 4.273341 4.435766
    Apr 2014       4.378331 4.321645 4.435018 4.291636 4.465026
    May 2014       4.417831 4.361133 4.474529 4.331118 4.504543
    Jun 2014       4.517751 4.460791 4.574712 4.430638 4.604865
    Jul 2014       4.462807 4.405474 4.520141 4.375123 4.550491
    Aug 2014       4.409070 4.351736 4.466403 4.321386 4.496754
    Sep 2014       4.319325 4.261906 4.376743 4.231511 4.407139
    Oct 2014       4.259274 4.201827 4.316722 4.171416 4.347133
    Nov 2014       4.116758 4.059310 4.174206 4.028899 4.204617
    Dec 2014       4.010205 3.952740 4.067671 3.922319 4.098092



    plot(dairy.fit2)


![png](output_85_0.png)


## Summary

In this lab you have learned to work with and analyze time series data. Specifically, you have done the following:

- Examined the properties of time series objects.
- Plotted time series data.
- Decomposed time series data into its trend, seasonal, and remainder components.
- Modeled the remainder components as AR, MA, ARMA and ARIMA models. 
- Created and evaluated difference series methods.
- Constructed and evaluated a forecasting model. 
