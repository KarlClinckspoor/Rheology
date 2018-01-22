# encoding: utf-8
import matplotlib.pyplot as plt
import glob
import numpy as np
from scipy.optimize import curve_fit
import logging
import traceback

# antes de começar, perguntar se quer analisar tudo usando as configurações anteriores.
# todo: criar uma opção para gerar um arquivo xlsx separando em nome, horário, data e valores.
# todo: sempre checar se o comprimento de GP e Eta são iguais e se só tem valores numéricos em ambos.
# todo: usar um logging bom


# default_settings
settings = {
            'DEBUG': False,
            'INLINE_GRAPHS': False,
            'SORTING_METHOD_LIN': 'by_error_length',  # by_error, by_error_length
            'NL_FITTING_METHOD': 'Carreau',  # Carreau, Cross, Carreau-Yasuda
            'SORTING_METHOD_NL': 'overall',  # eta_0, overall
            'PLOT_GRAPHS': False,
            'SAVE_GRAPHS': False,
            'AUTO_LIN': True,
            'AUTO_NL': True,
            'DO_LIN': True,
            'DO_NL': True,
            'TREAT_ALL':False,
            'EXT':'txt'
            }

valid_options_lin_sorting = ['by_error', 'by_error_length']
valid_options_nl_sorting = ['eta_0', 'overall']
models = ['Carreau', 'Cross', 'Carreau-Yasuda']
Param_names_errs = {'Carreau':'eta_0 eta_inf GP_b n eta_0_err eta_inf_err GP_b_err n_err',
                    'Cross': 'eta_0 eta_inf GP_b n eta_0_err eta_inf_err GP_b_err n_err',
                    'PowerLaw':'k n k_err n_err',
                    'Carreau-Yasuda':'eta_0 eta_inf lambda a n eta_0_err eta_inf_err lambda_err a_err n_err',
                    'Linear':'a a_err'}
# Param_names = {'Carreau':'eta_0 eta_inf GP_b n',
#                'Cross': 'eta_0 eta_inf GP_b n',
#                'PowerLaw':'k n',
#                'Carreau-Yasuda':'eta_0 eta_inf lambda a n',
#                'Linear':'a'}
# Param_names_errors = {'Carreau':'eta_0_err eta_inf_err GP_b_err n_err',
#                      'Cross':'eta_0_err eta_inf_err GP_b_err n_err',
#                      'PowerLaw':'k_err n_err',
#                      'Carreau-Yasuda':'eta_0_err eta_inf_err lambda_err a_err n_err',
#                      'Linear':'a_err'}

def calculate_derivatives():
    """"With the presence of the library "uncertainties", this has become irrelevant, but it's still an interesting
    exercise. This returns no interesting result, other than printing the partial derivatives of the Carreau model."""
    from sympy import symbols, diff
    # Carreau
    eta_inf, eta_0, GP, GP_b, n = symbols('eta_inf eta_0 GP GP_b n', real=True)
    Carr = eta_inf + (eta_0 - eta_inf) / (1 + (GP / GP_b) ** 2) ** (n / 2)
    print(diff(Carr, eta_inf))
    print(diff(Carr, eta_0))
    print(diff(Carr, GP_b))
    print(diff(Carr, n))
    return


# -------------FITTING MODELS-------------------

def fit_Carreau(GP, eta_0, eta_inf, GP_b, n):
    """Eta = eta_inf + (eta_0 - eta_inf) / (1+(GP/GP_b)**2)**(n/2)
    GP_b is a constant with the dimension of time and n is a dimensionless constant"""
    return eta_inf + (eta_0 - eta_inf) / (1 + (GP / GP_b) ** 2) ** (n / 2)



def fit_Cross(GP, eta_0, eta_inf, GP_b, n):
    return eta_inf + (eta_0 - eta_inf) / (1 + (GP / GP_b) ** n)



def fit_PowerLaw(GP, k, n):
    """Power Law: eta = k * GP ** (n-1)"""
    return k * GP ** (n - 1)



def fit_CarreauYasuda(GP, eta_0, eta_inf, lbda, a, n):
    """Carreau-Yasuda: eta(GP) = eta_inf + (eta_0 - eta_inf)(1+(lambda * GP)**a)**((n-1)/a)"""
    return eta_inf + (eta_0 - eta_inf) * (1 + (lbda * GP) ** a) ** ((n - 1) / a)


def fit_lin(x, a, b):
    """Simple function for a linear fit, with a as the linear coefficient and b the angular coefficient."""
    return a + b * x


# ---------------UNCERTAINTY CALCULATIONS------
def carr_uncertainty(GP, eta0, etainf, GPb, n, eta0_err, etainf_err, GPb_err, n_err):
    """Uses the uncertainty package to calculate the Carreau model values. GP
    can be a numpy array, which returns two lists of values and errors, a float64,
    float or int and returns a tuple (val, err)"""
    from uncertainties import ufloat
    f_eta0 = ufloat(eta0, eta0_err)
    f_etainf = ufloat(etainf, etainf_err)
    f_GPb = ufloat(GPb, GPb_err)
    f_n = ufloat(n, n_err)
    Carr = f_etainf + (f_eta0 - f_etainf) / (1 + (GP / f_GPb) ** 2) ** (f_n / 2)

    # Extracts all val +- err pairs if GP is an ndarray
    if type(GP) is np.ndarray:
        Carr_val = [a.nominal_value for a in Carr]
        Carr_err = [a.std_dev for a in Carr]

    # If GP is numeric, separates the two values.
    if (type(GP) is np.float64) or (type(GP) is float) or (type(GP) is int):
        Carr_val = Carr.nominal_value
        Carr_err = Carr.std_dev

    return Carr_val, Carr_err

def cross_uncertainty(GP, eta0, etainf, GPb, n, eta0_err, etainf_err, GPb_err, n_err):
    from uncertainties import ufloat
    f_eta0 = ufloat(eta0, eta0_err)
    f_etainf = ufloat(etainf, etainf_err)
    f_GPb = ufloat(GPb, GPb_err)
    f_n = ufloat(n, n_err)
    Cross = f_etainf + (f_eta0 - f_etainf) / (1 + (GP / f_GPb) ** f_n)

    # Extracts all val +- err pairs if GP is an ndarray
    if type(GP) is np.ndarray:
        Cross_val = [a.nominal_value for a in Cross]
        Cross_err = [a.std_dev for a in Cross]

    # If GP is numeric, separates the two values.
    if (type(GP) is np.float64) or (type(GP) is float) or (type(GP) is int):
        Cross_val = Cross.nominal_value
        Cross_err = Cross.std_dev

    return Cross_val, Cross_err


def CarreauYasuda_uncertainty(GP, eta0, etainf, lbda, a, n, eta0_err, etainf_err, lbda_err, a_err, n_err):
    from uncertainties import ufloat
    f_eta0 = ufloat(eta0, eta0_err)
    f_etainf = ufloat(etainf, etainf_err)
    f_n = ufloat(n, n_err)
    f_lbda = ufloat(lbda, lbda_err)
    f_a = ufloat(a, a_err)
    CY = f_etainf + (f_eta0 - f_etainf) * (1 + (f_lbda * GP) ** f_a) ** ((f_n - 1) / f_a)

    # Extracts all val +- err pairs if GP is an ndarray
    if type(GP) is np.ndarray:
        CY_val = [a.nominal_value for a in CY]
        CY_err = [a.std_dev for a in CY]

    # If GP is numeric, separates the two values.
    if (type(GP) is np.float64) or (type(GP) is float) or (type(GP) is int):
        CY_val = CY.nominal_value
        CY_err = CY.std_dev

    return CY_val, CY_err

# todo: powerlaw model.

# -----------File manipulation----------------
def ExtractData(fname, FC_segment=0):
    """Opens the file fname and extracts the data based on where it finds the word 'Eta' and 'GP', these being
    the Viscosity and the Shear Rate (gamma point). If the file has multiple segments, for example, when multiple
    experiments were done in succession, FC_segment indicates which of those experiments was a Flow Curve."""
    fhand = open(fname, 'r')
    GP = []
    Eta = []
    column_eta = 0
    column_gp = 0
    #FC_segment = '3'
    
    #while FC_segment == 0:
    #    FC_segment = input("What is the segment that has the flow curves? (eg. [1], 2, 3) If you do not know, don't write anything. ")
    #    if FC_segment == '':
    #        print(fhand.read())
    #    elif FC_segment.isnumeric():
    #        break
    #    else:
    #        print('Not a valid number')
    
    for line in fhand:
        if line.startswith(';'):
            column_names = line.rstrip().split(';')
            if settings['DEBUG']:
                print('Debug: column names', column_names)
            for i, column in enumerate(column_names):
                if 'Eta' in column and 'Eta*' not in column:
                    column_eta = i
                    if settings['DEBUG']:
                        print('Debug: Found Eta at', column_eta)
                if 'GP' in column:
                    column_gp = i
                    if settings['DEBUG']:
                        print('Debug: Found GP at', column_gp)
        try:
            GP.append(float(line.replace(',','.').split(';')[column_gp]))
            Eta.append(float(line.replace(',','.').split(';')[column_eta]))
        except:
            pass
        
        #if line.startswith(FC_segment + '|'):
        #    line = line.rstrip()
        #    num, gp, tau, eta, *rest = line.replace(',','.').split(';')
        #    GP.append(float(gp))
        #    Eta.append(float(eta))
        #    #print(line)
    
    fhand.close()
    if len(GP) == 0:
        #print('!!!!No Flow Curve data was found! Re-export the data on file', fname)
        raise ValueError
    #return pd.Series(GP), pd.Series(Eta)
    if settings['DEBUG']:
        print('Debug: Extracted Data: GP:', GP, 'Eta:', Eta)
    return GP, Eta


def record(name, eta0, eta0_err, silent = False, extra = '', fdest_name = 'results.csv'):
    """Writes to a .csv file the fitting results. Default name is results.csv"""
    if not silent:
        print(name + ':', 'Intercept', eta0, '+-', eta0_err, extra)

    with open(fdest_name, 'a', encoding = 'utf-8') as fdest:
        fdest.write(name + ';' + str(eta0) + ';' + str(eta0_err) + ';' + extra +'\n')


def save_settings():
    global settings
    with open('settings.dat', 'w') as settings_file:
        settings_file.write('#This is an automatically created settings file.\n')

        settings_file.write('\n# Inline graphs is meant for those with IDEs like Spyder\n')
        settings_file.write('INLINE_GRAPHS=' + str(settings['INLINE_GRAPHS']))

        settings_file.write('\n# Plot graphs after fitting, with error propagation? Slows down the process greatly.\n')
        settings_file.write('PLOT_GRAPHS=' + str(settings['PLOT_GRAPHS']))

        settings_file.write('\n# Save graphs after fitting? Useless if PLOT_GRAPHS is set to False\n')
        settings_file.write('SAVE_GRAPHS=' + str(settings['SAVE_GRAPHS']))

        settings_file.write('\n# Treat all files in folder?\n')
        settings_file.write('TREAT_ALL=' + str(settings['TREAT_ALL']))

        settings_file.write('\n# Extension of files to look for\n')
        settings_file.write('EXT=' + str(settings['EXT']))

        settings_file.write('\n\n#### Linear Fitting ####')

        settings_file.write('\n# Perform linear fit?\n')
        settings_file.write('DO_LIN=' + str(settings['DO_LIN']))

        settings_file.write("\n# Sorting method of the automatic linear method.\n# Can be 'by_error', minimizing the " +
                            "error, or 'by_error_length', which minimizes the error divided by the total number of " +
                            "points used.\n")
        settings_file.write('SORTING_METHOD_LIN=' + str(settings['SORTING_METHOD_LIN']))

        settings_file.write('\n# Set to True if you want the linear fitting to be done automatically. ' +
                            'False, if to be done manually.\n')
        settings_file.write('AUTO_LIN=' + str(settings['AUTO_LIN']))

        settings_file.write('\n\n#### Non-Linear Fitting ####')

        settings_file.write('\n# Perform non-linear fitting?\n')
        settings_file.write('DO_NL=' + str(settings['DO_NL']))

        settings_file.write("\n# Fitting method. 'Carreau', 'Cross', 'Carreau-Yasuda'\n")
        settings_file.write('NL_FITTING_METHOD=' + str(settings['NL_FITTING_METHOD']))

        settings_file.write("\n# Can be 'overall', minimizing the error of all parameters, or 'eta_0', " +
                            "minimizing the error of only this parameter.\n")
        settings_file.write('SORTING_METHOD_NL=' + str(settings['SORTING_METHOD_NL']))

        settings_file.write('\n# Set to True if you want the non linear fitting to be done automatically\n')
        settings_file.write('AUTO_NL=' + str(settings['AUTO_NL']))

        settings_file.write('\n\n##### Debug #####\n')
        settings_file.write('# Show debug messages\n')
        settings_file.write('DEBUG=' + str(settings['DEBUG']))


def load_settings():
    global settings
    try:
        fhand = open('settings.dat', 'r')
    except FileNotFoundError:
        print('Settings file not found. Loading defaults.')
        print('Creating a new settings file with the defaults')
        save_settings()
        return

    if settings['DEBUG']:
        print('===Builtin settings===')
        for key, value in settings.items():
            print(key, '=', value)
        print('======================')

    for line in fhand:
        if line.startswith('#'):
            continue
        if len(line) < 4:
            continue
        line = line.rstrip()
        line = line.replace(' ', '')
        #print(line, len(line))
        var, value = line.split('=')
        if var not in settings:
            if settings['DEBUG']:
                print("Settings file has an unrecognized parameter", var, "which wasn't loaded")
            continue
        if settings['DEBUG']:  # checks if the default settings in the source code has the Debug enabled
            print(line)
            print('Var =', var, 'val =', value)
        if value.lower() == 'true':
            settings[var] = True
        elif value.lower() == 'false':
            settings[var] = False
        else:
            settings[var] = value
    fhand.close()

    if settings['DEBUG']:
        print('===Current settings===')
        for key, value in settings.items():
            print(key, '=', value)
        print('======================')
    return 


def edit_settings():
    global settings
    counter = 0
    valid_numbers = []
    number_setting_corr = {}
    for key, value in settings.items():
        print(counter, ')', key, '=', value, end='')
        if type(value) == bool:
            print(': Options= True | False')
        if key == 'SORTING_METHOD_LIN':
            print(': Options= ', end='')
            print(*valid_options_lin_sorting, sep=' | ')
        if key == 'NL_FITTING_METHOD':
            print(': Options= ', end='')
            print(*models, sep=' | ')
        if key == 'SORTING_METHOD_NL':
            print(': Options= ', end='')
            print(*valid_options_nl_sorting, sep=' | ')
        if key == 'EXT':
            print(': Options= txt | dat | csv')
        valid_numbers.append(str(counter))
        number_setting_corr[str(counter)] = key
        counter += 1
    print('\n===================\n')
    if settings['DEBUG']:
        print('Valid numbers', valid_numbers)
        print('Correlation', number_setting_corr)
    while True:
        print('Which setting you want to change? Enter "number, new value" to modify, or "quit" to exit.')
        print('Observe the possible values for each setting! They are case sensitive. '
              'Inputting wrong values might break the program. \n')
        choice = input('Input:')
        if choice == 'quit':
            break
        if ',' not in choice:
            print('Invalid input. Place the number, followed by a comma, followed by its value')
            continue
        if len(choice.split(',')) != 2:
            print('Invalid input, must have only one comma')
            continue
        var, val = choice.split(',')
        if var not in valid_numbers:
            print('Invalid number.')
            continue
        real_var = number_setting_corr[var]
        if val.lower() == 'true':
            settings[real_var] = True
            continue
        elif val.lower() == 'false':
            settings[real_var] = False
            continue
        else:
            settings[real_var] = val
        # todo: check for all possible values to avoid inputting wrong settings and messing everything up.
        # if val not in valid_options_nl_sorting:
        #     print('Invalid nonlinear sorting option. Case sensitive! Be very precise.')
        #     continue
        # if val not in valid_options_lin_sorting:
        #     print('Invalid linear sorting option. Case sensitive! Be very precise.')
        #     continue
        # if val not in models:
        #     print('Invalid nonlinear fitting model. Case sensitive! Be very precise.')
        #     continue

    print('===Final settings===')
    print(*[key + '=' + str(value) + '\n' for key, value in settings.items()], sep='')
    print('====================')
    save_settings()
    return


def print_help():
    help_file = open('help', 'r', encoding='utf-8').read()
    print(help_file)


def select_files():
    files = []
    extension = input('What is the file extension? txt, dat, etc:\n')
    allfiles = glob.glob('*.' + extension)
    print(*[str(num) + ')' + file + '\n' for num, file in enumerate(allfiles)], sep='')
    while True:
        file_to_add = input('Which file to add? Number, nothing to continue or "quit" to exit: ')
        if file_to_add == 'quit':
            return []
        elif file_to_add == '':
            break
        else:
            try:
                files.append(allfiles[int(file_to_add)])
            except IndexError:
                print('Invalid value')
            except ValueError:
                print('Enter a number, not text')
    if len(files) == 0:
        print('No file was selected! The program will now quit.')
        return files
    print('====Selected files:====')
    print(*[file + '\n' for file in files], sep='', end='')
    print('=======================')
    return files

# -----------Plotting-------------------------
def PlotData(GP, Eta):
    """Plots the data with labels accompanying each point, for easier attribution.
    Depending on where the program is being run, it is better to leave """
    #labels = [i in range(0, len(GP), 1)]
    plt.xscale('log')
    plt.yscale('log')


    plt.scatter(GP, Eta, marker='*')
    for i in range(0, len(GP)):
        plt.annotate(str(i), (GP[i], Eta[i]))
    
    if not settings['INLINE_GRAPHS']:
        plt.draw()
        plt.pause(0.001)
    
    if settings['INLINE_GRAPHS']:
        plt.show()
    return None


def plot_error_graphs(file, GP, Eta, params=[], model = '', first_point=0, last_point=-1, param_names=''):
    """"Plot fitted model with errorbars, overlaying the data.
    file is the file name only, not the file handle.
    GP, Eta are the extracted data files
    model is the model used for fitting
    params are the fitting parameters, arranged in the following way:
    [par1, par2, par3, par1_err, par2_err, par3_err]
    first and last points are the points used for fitting, must be integers
    param_names contains the parameter names, so that they can be included in the graph"""
    TEXT_FILENAME_X = 0.1
    TEXT_PARAMS_X = 0.3
    TEXT_Y = 0.98
    x = np.logspace(np.log10(GP[0]), np.log10(GP[-1]))
    if settings['DEBUG']:
        print('Debug: x', x)
        print('Debug: params', params)
        print('Debug: GP', GP, 'Eta', Eta)
    
    if model == 'Carreau':
        y, yerr = carr_uncertainty(x, *params)
    if model == 'Cross':
        y, yerr = cross_uncertainty(x, *params)
    if model == 'Carreau-Yasuda':
        y, yerr = CarreauYasuda_uncertainty(x, *params)
    if model == 'Linear':
        y, yerr = np.ones(len(x)) * params[0], np.ones(len(x)) * params[1]  # todo: check if this is right
    # todo: add powerlaw model.

    if not settings['INLINE_GRAPHS']:
        plt.clf()

    plt.xscale('log')
    plt.yscale('log')
    plt.scatter(GP, Eta)
    plt.errorbar(x, y, yerr=yerr)

    plt.annotate(str(first_point+1), (GP[first_point], Eta[first_point]), color = 'red')
    if last_point == -1:
        plt.annotate(str(len(GP)), (GP[last_point], Eta[last_point]), color = 'red')
    else:
        plt.annotate(str(last_point), (GP[last_point], Eta[last_point]), color = 'red')

    # Creates strings to be included in the graphs
    model_param_names = 'Model: ' + model + ' Params: ' + param_names
        
    # This is kinda obscure on purpose, just for fun. What it does is that it transforms
    # each parameter value into 2 significant digits and then joins these as strings
    param_values = ' '.join(list(map(str, [round(x, 2) for x in params[:len(params) // 2]])))  # First half has param values
    param_errors = ' '.join(list(map(str, [round(x, 2) for x in params[len(params) // 2:]])))
    
    #param_values = ' ' * len(model) + ' '.join(list(map(str, params[:len(params) // 2])))  # First half has param values
    #param_errors = ' ' * len(model) + ' '.join(list(map(str, params[len(params) // 2:])))  # Second half has param errors
    
    # Change font to a monospace font.
    # Places strings in the graphs
    plt.figtext(TEXT_FILENAME_X, TEXT_Y, file, size='small')
    plt.figtext(TEXT_FILENAME_X, TEXT_Y - 0.030, model_param_names, size='small')
    plt.figtext(TEXT_FILENAME_X, TEXT_Y - 0.060, param_values, size='small')
    plt.figtext(TEXT_FILENAME_X, TEXT_Y - 0.090, param_errors, size='small')

    if settings['SAVE_GRAPHS']:
        plt.savefig(file[:-4] + '.png')
        print('Figure saved.')
    if not settings['INLINE_GRAPHS']:
        plt.draw()
        plt.pause(0.5)
        plt.clf()
    if settings['INLINE_GRAPHS']:
        plt.show()
    return


# todo: make this generic so that it can accept all fitting methods.
def plot_fitted_data(file, GP, Eta, a, aerr, p1, p2, method=''):
    """Plots the fitted data so that one can evaluate how good the fit is."""
    TEXT_X = 0.6
    TEXT_Y = 1.0
    if settings['DEBUG']:
        print('Debug: GP', GP, 'Eta', Eta)
    x = np.logspace(np.log10(GP[0]), np.log10(GP[-1]))
    if settings['DEBUG']:
        print('Debug: x', x)
    y = np.ones(len(x)) * a
    if settings['DEBUG']:
        print('Debug: y', y)
    yerr = np.ones(len(x)) * aerr

    if not settings['INLINE_GRAPHS']:
        plt.clf()

    plt.xscale('log')
    plt.yscale('log')
    plt.scatter(GP, Eta)
    plt.errorbar(x, y, yerr=yerr)

    plt.annotate(str(p1), (GP[p1], Eta[p1]), color='red')
    plt.annotate(str(p2), (GP[p2], Eta[p2]), color='red')

    plt.figtext(TEXT_X - 0.3, TEXT_Y - 0.05, file)
    plt.figtext(TEXT_X - 0.3, TEXT_Y - 0.1, method)
    plt.figtext(TEXT_X, TEXT_Y - 0.05, 'a = ' + str(a), color='red')
    plt.figtext(TEXT_X, TEXT_Y - 0.10, 'aerr = ' + str(aerr), color='red')

    if settings['SAVE_GRAPHS']:
        plt.savefig(file[:-4] + '.png')
        print('Figure saved.')
    if not settings['INLINE_GRAPHS']:
        plt.draw()
        plt.pause(0.5)
        plt.clf()
    if settings['INLINE_GRAPHS']:
        plt.show()
    return


# -----------Fitting algorithms---------------
def Average_mean(Eta):
    """Function that takes a pandas series and returns the mean and the error."""
    import pandas as pd
    aver = pd.Series(Eta).mean()
    aver_err = pd.Series(Eta).std()
    return aver, aver_err

        
def ManualDataFit(file, model='Linear'):
    """Treats an individual file using a model.
    Available methods:
    Linear
    Carreau
    Cross
    Power Law
    Carreau-Yasuda
    """
    GP, Eta = ExtractData(file)

    PlotData(GP, Eta)

    while True:
        try:
            initialp = int(input('What is the first point you want to use? '))
            finalp = int(input('What is the final point you want to use? '))
            break
        except ValueError:
            print('Not a valid number')
            continue

    GP_arr = np.array(GP[initialp:finalp + 1])
    Eta_arr = np.array(Eta[initialp:finalp + 1])

    if model == 'Mean':
        print('Averaging the points')
        aver, aver_err = Average_mean(Eta[initialp:finalp+1])
        record(file, aver, aver_err, extra=str(initialp)+str(finalp+1))
    
    if model == 'Linear':
        popt, pcov = curve_fit(fit_lin, GP_arr, Eta_arr, p0=(30, 0),
                               bounds=(0, [1000., 0.0001]))
        perr = np.sqrt(np.diag(pcov))

        a = popt[0]
        aerr = perr[0]
        record(file, a, aerr, extra=';first' + str(initialp) + 'final' + 
               str(finalp+1) + ';manual selection')
        #plot_fitted_data(file, GP, Eta, popt[0], perr[0], initialp, finalp, method='manual ' + model)
    
    if model == 'Carreau':
        popt, pcov = curve_fit(fit_Carreau, GP_arr, Eta_arr)
        perr = np.sqrt(np.diag(pcov))
        
        eta_0, eta_inf, GP_b, n = popt[0], popt[1], popt[2], popt[3]
        if settings['DEBUG']:
            print('Carreau fitting, eta_0, eta_inf, GP_b, n', eta_0, eta_inf, GP_b, n)

    if model == 'Carreau-Yasuda':
        popt, pcov = curve_fit(fit_CarreauYasuda, GP_arr, Eta_arr)
        perr = np.sqrt(np.diag(pcov))

        eta_0, eta_inf, lbda, a, n = popt[0], popt[1], popt[2], popt[3], popt[4]
        if settings['DEBUG']:
            print('Carreau-Yasuda fitting, eta_0, eta_inf, lbda, a, n', eta_0, eta_inf, lbda, a, n)

    if model == 'PowerLaw':
        popt, pcov = curve_fit(fit_PowerLaw, GP_arr, Eta_arr)
        perr = np.sqrt(np.diag(pcov))

        k, n = popt[0], popt[1]
        if settings['DEBUG']:
            print('PowerLaw fitting, k, n', k, n)

    if model == 'Cross':
        popt, pcov = curve_fit(fit_Cross, GP_arr, Eta_arr)
        perr = np.sqrt(np.diag(pcov))

        eta_0, eta_inf, GP_b, n = popt[0], popt[1], popt[2], popt[3]
        if settings['DEBUG']:
            print('Cross fitting, eta_0, eta_inf, GP_b, n', eta_0, eta_inf, GP_b, n)

    do_plot = 'Do you want to plot the data? y/[n]'
    #if do_plot == 'y':
    #    plot_fitted_data(file, GP, Eta, a, aerr, method = model + ' manual')
    return (initialp, finalp), popt, perr


def automatic_linear_Fitting(GP, Eta, sorting=settings['SORTING_METHOD_LIN']):
    """Goes through all the files, fits them and selects the best fit according to two algorithms.
    First, it selects two points, a beginning and an end point, the first starting at point 0
    and going to a third of the curve. The second, starting at points to the right,
    going until the middle of the curve.
    Then, it fits the data by fixing the slope at 0 and goes through every possible combination
    of the first and second points.
    It selects the data based on two criteria:
    1. sorting = 'by_error': finds the minimal error. Tends to select less points overall and 
        gives a fitting with a less than ideal representation overall.
    2. sorting = 'by_error_length': divides the error by how many points were used in the fit.
        May result in a higher overall error, but gives a better representation of the curve.
    """
    VISC_LIMIT = 100000. # Upper limit for the viscosity in the fitting.
    
    #segment = input('What is the segment that contains the flow curves for all the data?')
    #GP, Eta = ExtractData(file)  #, FC_segment = segment)
    length = len(GP)
    
    fittings = []
    for first_point in range(0, length // 3, 1):
        for last_point in range(first_point + 3, length // 2, 1):
            GP_arr = np.array(GP[first_point:last_point + 1])
            Eta_arr = np.array(Eta[first_point:last_point + 1])
            popt, pcov = curve_fit(fit_lin, GP_arr, Eta_arr, p0=(30, 0), 
                                   bounds=(0, [VISC_LIMIT, 0.0001]))
            perr = np.sqrt(np.diag(pcov))
            fittings.append((first_point, last_point, popt, perr))
    #                           [0]           [1]      [2]   [3]
    #print('Debug: fittings: ', fittings)
    if sorting == 'by_error':
        fittings.sort(key = lambda x: x[3][0]) # selects the smallest error
    elif sorting == 'by_error_length':
        fittings.sort(key = lambda x: x[3][0] / (x[1]-x[0])) # selects the smallest error divided by
                                                             # the largest length
    else:
        print('Invalid sorting method. Please select change to a valid value.')
    a = fittings[0][2][0]  # popt
    aerr = fittings[0][3][0]  # perr
    if settings['DEBUG']:
        print('Debug: fittings_sorted: ', fittings)
        print('Debug: a: ', a)
        print('Debug: aerr: ', aerr)
    #return fittings[0]
    return (fittings[0][0], fittings[0][1]), a, aerr


def all_models_fitting(file):
    """"Fits all existing models and returns all the possible parameters
    For linear fits: first and last points, intercept, intercept_error
    For the other fits: first point, parameters, errors"""
    GP, Eta = ExtractData(file)  #, FC_segment = segment)
    params = []

    params.append(automatic_linear_Fitting(GP, Eta))
    params.append(nonlinear_model_auto_fitting(GP, Eta, method ='Carreau'))
    params.append(nonlinear_model_auto_fitting(GP, Eta, method ='Cross'))
    params.append(nonlinear_model_auto_fitting(GP, Eta, method='Carreau-Yasuda'))
    """
    lin_first, lin_last, lin_a, lin_a_err = automatic_linear_Fitting(GP, Eta)
    params.append(lin_first, lin_last, lin_a, lin_a_err)
    carr_num, carr_param, carr_err = nonlinear_model_auto_fitting(GP, Eta, method ='Carreau')
    params.append(carr_num, carr_param, carr_err)
    cross_num, cross_param, cross_err = nonlinear_model_auto_fitting(GP, Eta, method ='Cross')
    params.append(cross_num, cross_param, cross_err)
    CY_num, CY_param, CY_err = nonlinear_model_auto_fitting(GP, Eta, method='Carreau-Yasuda')
    params.append(CY_num, CY_param, CY_err)
    """
    # todo: return the parameters separately, for each model.
    return params
    

#todo: add bounds to all fitting functions.
def nonlinear_model_auto_fitting(GP, Eta, method='Carreau', sorting='overall'):
    """Performs three different fittings, starting at the first point (0) and then skipping the first and
    second points. Compares these fittings returns the one with either the smallest 'overall' error, or the
    one with the smallest error in 'eta_0'.
    It returns a tuple with the following components:
    point used (0, 1, 2, etc)
    popt = [eta_0, eta_inf, GP_b, n]: the parameter values
    perr = [eta_0_err, eta_inf_err, GP_b_err, n_err]: the errors of the parameters
    """
    FIRST_POINT_MAX = 0
    fittings = []
    for first_point in range(0, FIRST_POINT_MAX + 1, 1):  # modified so that it does not check for the smallest error.
        GP_arr = np.array(GP[first_point:])
        Eta_arr = np.array(Eta[first_point:])
        if method == 'Carreau':
            #popt, pcov = curve_fit(fit_Carreau, GP_arr, Eta_arr)
            popt, pcov = curve_fit(fit_Carreau, GP_arr, Eta_arr, bounds=(0,np.inf))
            perr = np.sqrt(np.diag(pcov))
            fittings.append((first_point, popt, perr))
            if settings['DEBUG']:
                eta_0, eta_inf, GP_b, n = popt[0], popt[1], popt[2], popt[3]
                print('Carreau fitting, point, eta_0, eta_inf, GP_b, n', first_point, eta_0, eta_inf, GP_b, n)
                print('                                                ', perr[0], perr[1], perr[2], perr[3])
        if method == 'Cross':
            popt, pcov = curve_fit(fit_Cross, GP_arr, Eta_arr)
            perr = np.sqrt(np.diag(pcov))
            fittings.append((first_point, popt, perr))
            if settings['DEBUG']:
                eta_0, eta_inf, GP_b, n = popt[0], popt[1], popt[2], popt[3]
                print('Cross fitting, point, eta_0, eta_inf, GP_b, n', first_point, eta_0, eta_inf, GP_b, n)
                print('                                                ', perr[0], perr[1], perr[2], perr[3])
        if method == 'Carreau-Yasuda':
            popt, pcov = curve_fit(fit_CarreauYasuda, GP_arr, Eta_arr)
            perr = np.sqrt(np.diag(pcov))
            fittings.append((first_point, popt, perr))
            if settings['DEBUG']:
                print('Carreau-Yasuda fitting, point, eta_0, eta_inf, lbda, a, n', first_point, popt[0], popt[1],
                      popt[2], popt[3], popt[4])
                print('                                                ', perr[0], perr[1], perr[2], perr[3], perr[4])
        if method == 'PowerLaw':
            popt, pcov = curve_fit(fit_PowerLaw, GP_arr, Eta_arr)
            perr = np.sqrt(np.diag(pcov))
            fittings.append((first_point, popt, perr))
            if settings['DEBUG']:
                print('PowerLaw fitting, point, k, n', first_point, popt[0], popt[1])
                print('                                                ', perr[0], perr[1])
    if sorting == 'eta_0':
        fittings.sort(key = lambda x: x[2][0])
    if sorting == 'overall':
        fittings.sort(key = lambda x: x[2][0] + x[2][1] + x[2][2] + x[2][2])
    return fittings[0][0], fittings[0][1], fittings[0][2] # first_point, popt, perr



# ------------------ Main Loop ---------------------
# todo: remake the main program loop. Separate it into a function called main, make this call main.
def main_old():
    help = input('Do you want to read the help file? y/[n]')
    if help.lower() == 'y':
        print_help()
    all_or_select = input('Do you want to treat (a)ll files in the folder or a few (s)pecific files?')
    automatic_or_manual = input('Do you want to (s)elect the points or treat them (a)utomatically?')

    if all_or_select.lower() == 'a' or all_or_select.lower() == 'all':  # todo: mudar para all_or_select in valid_choices
        while True:
            ext = input("What is the file extension? Be careful, only select the relevant files?")
            allfilenames = glob.glob('*.' + ext)
            if len(allfilenames) > 0:
                break
            else:
                print('File not found')
        if automatic_or_manual == 's':
            for file in allfilenames:
                ManualDataFit(file)
        elif automatic_or_manual == 'a':
            for file in allfilenames:
                automatic_linear_Fitting(file)

    else:
        while True:
            ext = input('What is the file extension?')
            if ext == 'exit':
                break
            files = glob.glob('*.{}'.format(ext.lower()))
            if len(files) == 0:
                print('No file found with that extension')
                continue
            # pretty printer?
            for index, file in enumerate(files):
                print(index, file)
            while True:
                try:
                    selected_file = files[int(input('Which file number?'))]
                    break
                except ValueError:
                    print('Not a number')
                except IndexError:
                    print('This file number is not present')

            if automatic_or_manual == 's':
                ManualDataFit(selected_file)
            elif automatic_or_manual == 'a':
                automatic_linear_Fitting(selected_file)

            plt.clf()
            if input('Do you want to treat another file? y/[n]') == 'y':
                continue
            else:
                plt.close()
                break


# todo: checar se o método é uma lista ou uma string, e decidir entre fazer passar por todos os modelos.
def main():
    global settings
    load_settings()
    print('====Current settings:====')
    print(*[key + '=' + str(value) + '\n' for key, value in settings.items()], sep='', end='')
    print('=========================')
    do_change = input('Do you want to change the settings? y/[n]')
    if do_change == 'y':
        edit_settings()

    if settings['TREAT_ALL']:
        files = glob.glob('*.' + settings['EXT'])
        if len(files) == 0:
            print('No files with the extension', settings['EXT'], 'found. Please select them manually or change the '
                  'extension accordingly.')
            select_files()

    else:
        files = select_files()

    if len(files) == 0:  # If no files were selected, even manually.
        return

    for file in files:
        try:
            GP, Eta = ExtractData(file)
            GP, Eta = np.array(GP), np.array(Eta)
        except ValueError:
            print('!!!! No Flow Curve data was found! Re-export the data on file', file, 'Skipping.')
            with open('log', 'a') as log:
                log.write('No data found in file ' + file + '\n')
            continue
        if len(GP) != len(Eta):
            print('!!!! GP and Eta have different lengths! Re-export', file, 'or fix the problem manually.')
            print('!!!! Skipping file', file)
            with open('log', 'a') as log:
                log.write('GP and Eta of different lengths on file ' + file + '\n')
            continue

        np.seterr(over='raise') # Makes Overflow errors appear as exceptions, not only as warnings.
        try:
            if settings['DO_LIN']:
                lin_has_error = ''
                try:
                    if settings['AUTO_LIN']:
                        lin_points, a, aerr = automatic_linear_Fitting(GP, Eta, sorting = settings['SORTING_METHOD_LIN'])
                    else:
                        lin_points, a, aerr = ManualDataFit(file, model='Linear')
                except: # todo: debug here and try to find the types of exception that can occur
                    print('!!!!We have encountered an error while processing file', file)
                    print(traceback.format_exc())

                try:
                    if settings['PLOT_GRAPHS']:
                        plot_error_graphs(file[:-4] + '_lin_' + file[-4:], GP, Eta, params=np.array([a,aerr]),
                                          first_point=lin_points[0], last_point=lin_points[1], model = 'Linear',
                                          param_names='eta_0 eta_0_err')
                except:
                    print('Error found while plotting the linear fit')
                    print(traceback.format_exc())
                    lin_has_error = 'error_during_fitting'
                record(file, a, aerr, extra='linear automatic;FP=' + str(lin_points[0]) + 'LP=' + str(lin_points[1]) +
                                            lin_has_error, fdest_name='linear.csv')

            if settings['DO_NL']:
                try:
                    nonlinear_has_error = ''
                    if settings['AUTO_NL']:
                        nl_first, popt, perr = nonlinear_model_auto_fitting(GP, Eta, method=settings['NL_FITTING_METHOD'])
                        #nl_points = (nl_first, -1)
                    else:
                        nl_points, popt, perr = ManualDataFit(file, model=settings['NL_FITTING_METHOD'])
                except FloatingPointError:
                    print('!!!! Overflow detected on one of the parameters. Could not determine all parameters')
                    nonlinear_has_error = ';param_overflow_during_fitting'
                    with open('log', 'a') as log:
                        log.write('Parameter overflow while trying to fit file ' + file + '\n')
                try:
                    if settings['PLOT_GRAPHS']:
                        plot_error_graphs(file[:-4]+'_carr_'+file[-4:], GP, Eta, params=np.concatenate((popt, perr)),
                                          first_point=nl_first, model=settings['NL_FITTING_METHOD'],
                                  param_names=Param_names_errs[settings['NL_FITTING_METHOD']])
                except OverflowError: # todo: write which parameter has overflown
                    print('!!!! Overflow detected on one of the parameters. Could not plot the data')
                    nonlinear_has_error = ';param_overflow_during_fitting'
                try:
                    record(file, popt[0], perr[0], extra= 'nonlinear_auto_' + settings['NL_FITTING_METHOD'] +
                           nonlinear_has_error, fdest_name=settings['NL_FITTING_METHOD']+'.csv')
                except UnboundLocalError:
                    print('Not able to write to file because the subroutine did not return the fitting parameters')
                    record(file, 0, 0, extra='nonlinear_auto_' + settings['NL_FITTING_METHOD'] + ';' +
                                             'unable_to_find_viscosity', fdest_name=settings['NL_FITTING_METHOD']+'.csv')
                    with open('log', 'a') as log:
                        log.write('Unable to find viscosity for file ' + file + '\n')
        except:
            print('!!!!We have encountered an error while processing file', file)
            print(traceback.format_exc())
            with open('log', 'a') as log:
                log.write('Error while processing ' + file + '\n')
                log.write(traceback.format_exc())


def main_simple():
    """"Just gets all the .txt files in the folder, and fits with the linear model and Carreau model"""
    files = glob.glob('*.txt')
    for file in files:
        try:
            lin_has_error = ''
            GP, Eta = ExtractData(file)
            GP, Eta = np.array(GP), np.array(Eta)
        except ValueError:
            print('----------No Flow Curve data was found! Re-export the data on file', file)
            with open('log', 'a') as log:
                log.write('No data found in file ' + file + '\n')
            continue
        if len(GP) != len(Eta):
            print('!!!! GP and Eta have different lengths! Re-export the file or fix the problem manually.')
            print('!!!! Skipping file', file)
            with open('log', 'a') as log:
                log.write('GP and Eta of different lengths on file ' + file + '\n')
            continue
        try:
            # Fitting using the linear model
            lin_points, a, aerr = automatic_linear_Fitting(GP, Eta, sorting=settings['SORTING_METHOD_LIN'])
            try:
                if settings['PLOT_GRAPHS']:
                    plot_error_graphs(file[:-4]+'_lin_'+file[-4:], GP, Eta, params=np.array([a, aerr]),
                                      first_point=lin_points[0], last_point=lin_points[1], model='Linear',
                                      param_names='eta_0 eta_0_err')
            except:
                lin_has_error = 'error_during_fitting'
                print('Error found while plotting the linear fit')
            record(file, a, aerr, extra='linear automatic;FP=' + str(lin_points[0]) + 'LP=' + str(lin_points[1]) +
                   lin_has_error,fdest_name='linear.csv')
            # Fitting using the nonlinear model
            np.seterr(over='raise') # Makes Overflow errors appear as exception, not only as warnings.
            try:
                nonlinear_has_error = ''
                nl_first, popt, perr = nonlinear_model_auto_fitting(GP, Eta, method=settings['NL_FITTING_METHOD'])
                #print('Debug: Carreau ', popt[0], perr[0])
            except FloatingPointError:
                print('------Overflow detected on one of the parameters. Could not determine all parameters')
                nonlinear_has_error = '_param_overflow_during_fitting'
                with open('log', 'a') as log:
                    log.write('Parameter overflow while trying to fit file ' + file + '\n')
            try:  # sometimes a parameter can't be fit and matplotlib gets confused.
                if settings['PLOT_GRAPHS']:
                    plot_error_graphs(file[:-4]+'_carr_'+file[-4:], GP, Eta, params=np.concatenate((popt, perr)), first_point=nl_first,
                              model=settings['NL_FITTING_METHOD'],
                              param_names=Param_names_errs[settings['NL_FITTING_METHOD']])
            except OverflowError:
                print('------Overflow detected on one of the parameters. Could not plot data')
                nonlinear_has_error = '_param_overflow_during_fitting'
            try:
                record(file, popt[0], perr[0], extra='nonlinear_auto_' + settings['NL_FITTING_METHOD'] + nonlinear_has_error,
                   fdest_name='carreau.csv')
            except UnboundLocalError:
                print('Not able to write to file because the subroutine did not return the fitting parameters')
                record(file, 0, 0, extra='nonlinear_auto_' + settings['NL_FITTING_METHOD'] + ';'
                                         + 'unable_to_find_viscosity', fdest_name='carreau.csv')
                with open('log', 'a') as log:
                    log.write('Unable to find viscosity for file ' + file + '\n')
        except ValueError:
            print('The imported file has different lengths of GP and Eta! Re-export the file or fix it manually.')
        except:
            print('!!!!We have encountered an error while processing file', file)
            print(traceback.format_exc())
            with open('log', 'a') as log:
                log.write('Error while processing ' + file + '\n')
                log.write(traceback.format_exc())


if __name__ == '__main__':
    if settings['DEBUG']:
        print(settings)
    #main_simple()
    main()