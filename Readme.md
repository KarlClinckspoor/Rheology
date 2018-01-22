# RheoFC.py

# Description
This is a script that calculates the zero-shear viscosity of flow curves obtained from a rheometer automatically. It has the option to use several models for fitting, which are:

* Linear fitting of the plateau region. Uses a method that selects the best number of points to use by minimizing either the total error of the fit, or the error per point (total error / number of points).
* Nonlinear models: Cross, Carreau, Carreau-Yasuda. At the moment, the script does not vary the region where it fits (it fits the whole curve). However, the script can be easily modified to minimize the error of the viscosity or the error of all parameters by changing which points it considers.

The script can use matplotlib and the uncertainty package to plot the curves together with the data, and then save the plots for a quick way to check if the models are fitting the data well.

# Installation

1. Install python 3.6. Be sure to include python into the $PATH$ variable.
2. Use pip to install the requirements. `pip install -r requirements.txt`
3. Run the script using `python RheoFC.py`

This script uses

* `scipy` for model fitting
* `matplotlib` for plotting data, if needed
* `uncertainties` for error propagation
* `numpy` for arrays and numerical calculations

The following packages are included, but aren't necessarily used by the main script loop.
* `pandas` for data reading and recording.
* `sympy` for symbolically derivating the equations. The `uncertainties` takes care of error propagation, so it is unnecessary to write the propagation explicitly.
 


# Workflow

## Settings

The script checks for a 'settings.dat' file at the same folder as the script itself. If this file is not present, it loads its default settings, coded in the script, and writes a new settings file, correctly formatted. You can then edit this file manually, or edit the settings through the script. Be very careful with editing these. Mind the underlines and capitalization, as several settings are case-sensitive. After changing the settings, they are stored into the settings.dat file.

## Data loading

This script loads plaintext data, separated by commas. It requires a line starting with ';' showing all the parameter names, separated by ';', so that it can find which columns are GP (shear rate) and Eta (viscosity). These .txt files can be converted from .rwd (Haake RheoWin) through the DataManager software, or by configuring the acquisition to save the .txt file together with the .rwd file automatically. This is a sample of how a file should look like: 

    C:\...\file.rwd
    Company / Operator: - / -
    Date / Time / Version: --.--.---- / --:--:-- / HAAKE RheoWin 
    Substance / Sample no: --- / 

    ;GP in 1/s;Tau in Pa;Eta in Pas;T in ºC;t in s;t_seg in s;
    1|1; ;0,06978; ; ;17,86;7,239;
    [...]
    3|1;0,0009986;0,04055;40,61; ;245,4;33,01;   
 
 This specific job had 3 steps, the last of which is a flow curve. This can be seen by the empty values of GP and Eta on the first section (1|...). On the third section (3|...), GP and Eta are present, and the script extracts these. Be careful that sometimes the program will not correctly export all the values, but the script checks if GP and Eta both have the same number of points, and tells you if they don't. This can be fixed by re-exporting the data or adding the value manually in the .txt file.
 
 If your files have a different formatting, consider changing yourself the ExtractData function to suit your needs. As long as it returns two lists containing GP and Eta, of the same length, all the other functions should work without a problem.
 
## Data fitting

There are 4 models that can be chosen to fit your data. Normally, one chooses the linear model together with a nonlinear model. This can be changed by editing the settings file. Fitting using all 4 available models in one go is not implemented at this moment. 

The program can choose automatically which points to be included in the fitting process, but you can also choose them manually, one file at a time.

The automatic linear method uses the following algorithm. Keep in mind that the lists start from 0, so a list with 10 items goes from 0 to 9 (inclusive):
1. Chooses several combinations of two points, fits the data between these two points, stores the parameters and then chooses the one with the smallest error or the smallest error per point.
2. The first point is a number between 0 (first data point) and the total number of data points divided by 3, rounded down. For example, if there are 10 data points, 10//3 = 3, so it will use the first 3 data points, 0, 1 and 2.
3. The last point is at least 3 larger than the first data point, and goes halfway. For example, if there are 10 points, the it will vary, for the first point from 0+3 to 10//2 = 5th point = 4. Should the distance between the first and last data points become smaller than 3, the script ignores it.

So, for example, with 10 data points, the following combinations will be fitted:

    First | Last
        0 | 3
        0 | 4
        1 | 4

Should the first point be 2, the last point would either have to be 5, which is more than half (5 is the 6th point in a list) or, if the second point was kept at 4, the distance between 2 and 4 would be less than 3, which is not allowed.

The nonlinear methods fits one of 3 models into your data. It can minimize the error by varying the first point chosen for the fit and fixing the last point (similar to the linear method), but that is not currently implemented. Implementation is very easy, just change a variable in the `nonlinear_model_auto_fitting` function. In the future, an option to fit using all models will be implemented. 

If there are errors during the fitting, like Overflow errors (the algorithm can't assign a value to a specific parameter, for example), the program will write a log file and also show a brief explanation of what happened, where it happened.

After fitting, it writes only the zero-shear viscosity, and its error, in a .dat file with the name of your model. For example, the file 'test.txt' was fitted automatically with the linear and Carreau models, and had a viscosity of 50 +- 1. The linear.dat file will contain the entry:

    test.txt;50;1;linear_automatic;FP=0;LP=4
and the Carreau.dat file will contain:

    test.txt;50;1;nonlinear_auto_Carreau;

## Plotting

After fitting, one can choose to display the data with the model + error bars. These error bars were obtained through error propagation of the fitting parameters. Each figure will contain the file name, the parameter names and their respective errors. It is not intented to be pretty, just for a quick check.

If you have an interpreter that has inline plotting of graphs, like Spyder, it will show each graph as soon as it is produced. However, if you are running for the normal console, it will show each graph for 0.5s, then close it. Enough to get a look at the data. Choose the options accordingly.

These graphs can also be saved as .png figures if wanted.

If there was an Overflow error, or if the fitting algorithm wasn't able to determine one of the parameters, plotting won't be done and an error will be returned and logged.

## Automatization

There are two `main` functions in the script. The `main()` function uses the settings stored separately and asks you a few questions. The `main_simple()` function just treats all the data present in the current folder without asking anything, using the settings stored in the script. You can change the `if __name__ == '__main__'` section of the script to better suit your needs.