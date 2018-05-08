import matplotlib.pyplot as plt
import glob
import numpy as np
from scipy.optimize import curve_fit
import traceback
from lmfit import minimize, Parameters
from uncertainties import ufloat
import pandas as pd
import Settings
import sys

# todo: check program with several different settings
# todo: solve the problems with manual fitting

# todo: adjust the plotting function to write the parameters better in the tight_layout
# todo: implement R2' calculations, sorting, and adding them to the plots
#       R'2 = 1 - MQres/MQTcorr; MQTcor = SQTcor / I-1; MQres = SQres / I-p
#       I: degrees of freedom. p: number of parameters in the model
#       SQTcor = sum(yi - mean(y)) ^ 2
#       SQres = chisqr
# todo: check why the R2 on the automatic linear fit is buggy

# todo: Check the results file, if the correct parameters are being recorded
# todo: check if converting GP and Eta to ndarrays at the beginning breaks anything.

# todo: Remove the prev_extracted setting and try to guess this parameter

# todo: remove the debugging setting. Just use the debugging tools.

class Fitter:
    def __init__(self, filename, settings, do_fit=True):
        self.VISC_LIMIT = 10000000
        self.l_first_point = 0
        self.l_last_point = -1
        self.nl_first_point = 0
        self.nl_last_point = -1
        self.manip = FileManip()
        self.filename = filename
        self.settings = settings
        self.model = self.settings.NL_FITTING_METHOD
        self.l_R2 = 0
        self.nl_R2 = 0
        self.wait = float(settings.WAIT)
        self.fixed_fp = settings.FIXED_FP_NL
        if not self.fixed_fp:
            self.max_fp = int(settings.MAX_FP_NL)
        else:
            self.max_fp = 0

        self.lin_done = False
        self.nl_done = False

        if self.model == 'Cross':
            #self.param_names = ['eta_0', 'eta_inf', 'GP_b', 'n']
            self.param_names = 'eta_0 eta_inf GP_b n'
        elif self.model == 'Carreau':
            #self.param_names = ['eta_0', 'eta_inf', 'GP_b', 'n']
            self.param_names = 'eta_0 eta_inf GP_b n'
        elif self.model == 'Carreau-Yasuda':
            #self.param_names = ['eta_0', 'eta_inf', 'lambda', 'a', 'n']
            self.param_names = 'eta_0 eta_inf lambda a n'
        else:
            raise NameError(f'Did not understand model {self.model}')

        self.param_names_lin = ['Int', 'Slp']  # todo: check if this is the correct order.

        if self.settings.DO_LIN:
            self.int = 50
            self.int_err = 0.1
            self.slp = 0
            self.slp_err = 0
        elif self.settings.DO_NL:
            if self.model == 'Carreau' or self.model == 'Cross' or self.model == 'Carreau-Yasuda':
                self.params = [0, 0, 0, 0]
                self.param_errs = [0, 0, 0, 0]
            else:
                raise ValueError(f'Unknown model: {model}')

        try:
            if self.settings.PREV_EXTRACTED:
                self.GP, self.Eta = self.manip.ExtractData_pd(filename)
                self.GP = np.array(self.GP)
                self.Eta = np.array(self.Eta)
            else:
                self.GP, self.Eta = self.manip.ExtractData(filename)
                self.GP = np.array(self.GP)
                self.Eta = np.array(self.Eta)
        except ValueError:
            self.manip.logger(filename, 'Failed to open')
            raise ValueError(f'!!!! No Flow Curve data was found! Re-export the data on file{filename}.')
        except KeyError:
            self.manip.logger(filename, 'Failed to open')
            raise ValueError(f'!!!! No Flow Curve data was found! Re-export the data on file{filename}.')

        if len(self.GP) != len(self.Eta):
            self.manip.logger(self.filename, 'Failed to open')
            raise ValueError(f'!!!! GP and Eta have different lengths. '
                             f'Re-export {filename} or fix the problem manually.')

        if do_fit:
            self.fit()

    def _fit(self):  # Uses fit_curve. Does not provide an R2 value.
        if self.settings.DO_LIN:
            if self.settings.AUTO_LIN:
                self.automatic_lin_fitting(True)
            else:  # todo: plot, save and ask for the required points
                self.manual_fit(0, -1, 'Linear')
        if self.settings.DO_NL:
            if self.settings.AUTO_NL:
                self.automatic_nl_fitting(True)
            else:
                self.manual_fit(0, -1, self.settings.NL_FITTING_METHOD, True)

    def fit(self):
        if self.settings.DO_LIN:
            if self.settings.AUTO_LIN:
                self.automatic_lin_fitting_lm(True)
            else:  # todo: plot, save and ask for the required points
                self.manual_fit(0, -1, 'Linear')
        if self.settings.DO_NL:
            if self.settings.AUTO_NL:
                self.automatic_nl_fitting_lm(True)
            else:
                self.manual_fit(0, -1, self.settings.NL_FITTING_METHOD, True)

    @staticmethod
    def fit_Carreau(GP, eta_0, eta_inf, GP_b, n):
        """Eta = eta_inf + (eta_0 - eta_inf) / (1+(GP/GP_b)**2)**(n/2)
        GP_b is a constant with the dimension of time and n is a dimensionless constant"""
        return eta_inf + (eta_0 - eta_inf) / (1 + (GP / GP_b) ** 2) ** (n / 2)

    @staticmethod
    def fit_Cross(GP, eta_0, eta_inf, GP_b, n):
        return eta_inf + (eta_0 - eta_inf) / (1 + (GP / GP_b) ** n)

    @staticmethod
    def fit_PowerLaw(GP, k, n):
        """Power Law: eta = k * GP ** (n-1)"""
        return k * GP ** (n - 1)

    @staticmethod
    def fit_CarreauYasuda(GP, eta_0, eta_inf, lbda, a, n):
        """Carreau-Yasuda: eta(GP) = eta_inf + (eta_0 - eta_inf)(1+(lambda * GP)**a)**((n-1)/a)"""
        return eta_inf + (eta_0 - eta_inf) * (1 + (lbda * GP) ** a) ** ((n - 1) / a)

    @staticmethod
    def fit_lin(x, a, b):
        """Simple function for a linear fit, with a as the linear coefficient and b the angular coefficient."""
        return a + b * x

    @staticmethod
    def carr_uncertainty(GP, eta0, etainf, GPb, n, eta0_err, etainf_err, GPb_err, n_err):
        """Uses the uncertainty package to calculate the Carreau model values. GP
        can be a numpy array, which returns two lists of values and errors, a float64,
        float or int and returns a tuple (val, err)"""
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

    @staticmethod
    def cross_uncertainty(GP, eta0, etainf, GPb, n, eta0_err, etainf_err, GPb_err, n_err):
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

    @staticmethod
    def carryas_uncertainty(GP, eta0, etainf, lbda, a, n, eta0_err, etainf_err, lbda_err, a_err, n_err):
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

    def lm_curvefit(self, GP, Eta, do_lin=False):
        params = Parameters()
        SStot = sum((Eta - np.mean(Eta)) ** 2)
        if do_lin:  # todo: Check why R2 is very weird here.
            params.add('Int', 50, vary=True, min=0)
            params.add('Slp', 0, vary=False)
            fit = minimize(self.residual_lin, params, args=(GP, Eta))
            slp = fit.params['Slp'].value
            int = fit.params['Int'].value
            slp_err = fit.params['Slp'].stderr
            int_err = fit.params['Int'].stderr
            chisqr = fit.chisqr
            R2 = 1 - fit.chisqr / SStot
            return [slp, int], [slp_err, int_err], R2
        elif self.model == 'Carreau':
            params.add('eta_0', 100, vary=True, min=0)
            params.add('eta_inf', 1, vary=True, min=0)
            params.add('GP_b', 5, vary=True, min=0)
            params.add('n', 1, vary=True, min=0)
            fit = minimize(self.residual, params, args=(GP, Eta))
            params = [fit.params[par].value for par in fit.params]
            param_errs = [fit.params[par].stderr for par in fit.params]
            R2 = 1 - fit.chisqr / SStot
            return params, param_errs, R2
        elif self.model == 'Cross':
            params.add('eta_0', 100, vary=True, min=0)
            params.add('eta_inf', 1, vary=True, min=0)
            params.add('GP_b', 5, vary=True, min=0)
            params.add('n', 1, vary=True, min=0)
            fit = minimize(self.residual, params, args=(GP, Eta))
            params = [fit.params[par].value for par in fit.params]
            param_errs = [fit.params[par].stderr for par in fit.params]
            R2 = 1 - fit.chisqr / SStot
            return params, param_errs, R2
        elif self.model == 'Carreau-Yasuda':
            params.add('eta_0', 100, vary=True, min=0)
            params.add('eta_inf', 1, vary=True, min=0)
            params.add('lbda', 5, vary=True, min=0)
            params.add('a', 1, vary=True, min=0)
            params.add('n', 1, vary=True, min=0)
            fit = minimize(self.residual, params, args=(GP, Eta))
            params = [fit.params[par].value for par in fit.params]
            param_errs = [fit.params[par].stderr for par in fit.params]
            SSres = fit.chisqr
            R2 = 1 - SSres / SStot
            return params, param_errs, R2

    def residual(self, params, x, dataset):
        if self.model == 'Carreau':
            mod = self.fit_Carreau(x, params['eta_0'], params['eta_inf'], params['GP_b'], params['n'])
        elif self.model == 'Cross':
            mod = self.fit_Cross(x, params['eta_0'], params['eta_inf'], params['GP_b'], params['n'])
        elif self.model == 'Carreau-Yasuda':
            mod = self.fit_CarreauYasuda(x, params['eta_0'], params['eta_inf'], params['lbda'], params['a'],
                                         params['n'])
        resid = dataset - mod
        return resid

    def residual_lin(self, params, x, dataset):
        if type(x) == list:
            x = np.array(x)
        mod = params['Int'] + params['Slp'] * x
        resid = dataset - mod
        return resid

    def automatic_lin_fitting_lm(self, save=True):
        length = len(self.GP)
        fittings = []

        # Go through several possible ranges to fit, and fit them, then get the best fit
        for first_point in range(0, length//3, 1):
            for last_point in range(first_point + 3, length // 2, 1):
                GP_arr = np.array(self.GP[first_point:last_point + 1])  # todo: check if this conversion is necessary
                Eta_arr = np.array(self.Eta[first_point:last_point + 1])
                try:
                    #popt, pcov = curve_fit(self.fit_lin, GP_arr, Eta_arr, p0=(30, 0),
                    #                   bounds=(0, [self.VISC_LIMIT, 0.0001]))
                    params, param_errs, R2 = self.lm_curvefit(GP_arr, Eta_arr, do_lin=True)
                except:  # todo: test here and find what types of errors can occur
                    print(f'Error while using linear fit for file {self.filename}')
                    print(traceback.format_exc())
                    self.manip.logger(self.filename, 'Generic')

                #perr = np.sqrt(np.diag(pcov))
                fittings.append((first_point, last_point, params, param_errs, R2))

        if self.settings.LIN_SORTING_METHOD == 'by_error':
            fittings.sort(key=lambda x: np.log(x[3][0]))
        elif self.settings.LIN_SORTING_METHOD == 'by_error_length':
            fittings.sort(key=lambda x: np.log(x[2][1]) / (x[1] - x[0]))
        elif self.settings.LIN_SORTING_METHOD == 'by_R2':
            fittings.sort(key=lambda x: x[4])

        self.l_first_point = fittings[0][0]
        # todo: add variable names to first and last points of linear and nl
        self.l_last_point = fittings[0][1]
        self.int = fittings[0][2][1]
        self.int_err = fittings[0][3][1]
        self.l_R2 = fittings[0][4]
        self.lin_done = True

        if self.settings.DEBUG:
            print('Debug: fittings_sorted: ', fittings)
            print('Debug: a: ', self.int)
            print('Debug: aerr: ', self.int_err)

        if save:
            self.manip.record_fit('linear', self.int, self.int_err, silent=False,
                                  extra=f"{fittings[0][0]};{fittings[0][1]};")
            # todo: check how this was done before, for consistency

        return self.int, self.int_err, self.l_R2

    def automatic_lin_fitting(self, save=True):
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

        length = len(self.GP)
        fittings = []

        # Go through several possible ranges to fit, and fit them, then get the best fit
        for first_point in range(0, length//3, 1):
            for last_point in range(first_point + 3, length // 2, 1):
                GP_arr = np.array(self.GP[first_point:last_point + 1])
                Eta_arr = np.array(self.Eta[first_point:last_point + 1])
                try:
                    popt, pcov = curve_fit(self.fit_lin, GP_arr, Eta_arr, p0=(30, 0),
                                       bounds=(0, [self.VISC_LIMIT, 0.0001]))
                except:  # todo: test here and find what types of errors can occur
                    print(f'Error while using linear fit for file {self.filename}')
                    print(traceback.format_exc())
                    self.manip.logger(self.filename, 'Generic')

                perr = np.sqrt(np.diag(pcov))
                fittings.append((first_point, last_point, popt, perr))

        if self.settings.LIN_SORTING_METHOD == 'by_error':
            fittings.sort(key=lambda x: np.log(x[3][0]))  # gets perr of eta_0
        elif self.settings.LIN_SORTING_METHOD == 'by_error_length':
            fittings.sort( key=lambda x: np.log(x[3][0]) / (x[1] - x[0]) ) # divides perr by last-first

        self.int = fittings[0][2][0]
        self.int_err = fittings[0][3][0]
        self.l_first_point = fittings[0][0]  # todo: add variable names to first and last points of linear and nl
        self.l_last_point = fittings[0][1]
        self.lin_done = True

        if self.settings.DEBUG:
            print('Debug: fittings_sorted: ', fittings)
            print('Debug: a: ', self.int)
            print('Debug: aerr: ', self.int_err)

        if save:
            self.manip.record_fit('linear', self.int, self.int_err, silent=False,
                                  extra=f"{fittings[0][0]};{fittings[0][1]};")

        return self.int, self.int_err

    # todo: change from curve_fit to lm_fit.
    # todo: calculate R2 for all fittings and add it in the end to the class
    # todo: add options to sort by R2.

    def automatic_nl_fitting_lm(self, save=True):
        fittings = []
        try:
            max_range = len(self.GP) // self.max_fp
        except ZeroDivisionError:
            max_range = 1

        for first_point in range(0, max_range, 1):
            GP_arr = np.array(self.GP[first_point:])
            Eta_arr = np.array(self.Eta[first_point:])
            nonlinear_has_error = ''
            try:
                params, param_errs, R2 = self.lm_curvefit(GP_arr, Eta_arr, do_lin=False)
            except FloatingPointError:  # todo: check if these exceptions work
                print('!!!! Overflow detected on one of the parameters. Could not determine all parameters')
                nonlinear_has_error = ';param_overflow_during_fitting'
                self.manip.logger(self.filename, 'Overflow')
            except RuntimeError:
                print('!!!! Overflow detected on one of the parameters. Could not determine all parameters')
                nonlinear_has_error = ';param_overflow_during_fitting'
                self.manip.logger(self.filename, 'Overflow')
            except OverflowError:
                print('!!!! Overflow detected on one of the parameters.')
                self.manip.logger(self.filename, 'Overflow')

            fittings.append((first_point, params, param_errs, R2))

        if self.settings.NL_SORTING_METHOD == 'eta_0':
            fittings.sort(key=lambda x: x[2][0])
        elif self.settings.NL_SORTING_METHOD == 'overall':
            # fittings.sort(key=lambda x: x[2][0] + x[2][1] + x[2][2] + x[2][3])
            fittings.sort(key=lambda x: sum(x[2]))  # sums the errors
        elif self.settings.NL_SORTING_METHOD == 'R2':
            fittings.sort(key=lambda x: x[3])
        else:
            raise ValueError(f'Could not understand the sorting method {self.settings.NL_SORTING_METHOD}')

        self.nl_first_point = fittings[0][0]
        self.params = fittings[0][1]
        self.param_errs = fittings[0][2]
        self.nl_R2 = fittings[0][3]

        if save:  # todo: check here to return a good destination file
            try:
                self.manip.record_fit(
                    self.filename, self.params[0],
                    self.param_errs[0], silent=False,
                    extra=f"{fittings[0][0]};{fittings[0][1]};nonlinear_auto_{self.settings.NL_FITTING_METHOD};"
                          f"{nonlinear_has_error}", fdest_name=self.settings.NL_FITTING_METHOD + '.csv'
                )
            except UnboundLocalError:
                print('Unable to write to file because the subroutine did not return the fitting parameters')
                print(traceback.format_exc())
                self.manip.record_fit(self.filename, 0, 0, extra=f'nonlinear_auto_{self.settings.NL_FITTING_METHOD};'
                                                                 f'unable_to_find_viscosity',
                                      fdest_name=self.settings.NL_FITTING_METHOD + '.csv')
                self.manip.logger(self.filename, 'No Viscosity')

        self.nl_done = True
        return self.nl_first_point, self.params, self.param_errs, self.nl_R2

    def automatic_nl_fitting(self, save=True):
        fittings = []
        try:
            max_range = len(self.GP) // self.max_fp
        except ZeroDivisionError:
            max_range = 1

        for first_point in range(0, max_range, 1):
            GP_arr = np.array(self.GP[first_point:])
            Eta_arr = np.array(self.Eta[first_point:])
            nonlinear_has_error = ''
            try:
                if self.settings.NL_FITTING_METHOD == 'Carreau':
                    popt, pcov = curve_fit(self.fit_Carreau, GP_arr, Eta_arr, bounds=(0, np.inf))
                elif self.settings.NL_FITTING_METHOD == 'Cross':
                    popt, pcov = curve_fit(self.fit_Cross, GP_arr, Eta_arr, bounds=(0, np.inf))
                elif self.settings.NL_FITTING_METHOD == 'Carreau-Yasuda':
                    popt, pcov = curve_fit(self.fit_CarreauYasuda, GP_arr, Eta_arr)
                else:
                    raise ValueError(f'Model not present: {self.settings.NL_FITTING_METHOD}')

            except FloatingPointError:
                print('!!!! Overflow detected on one of the parameters. Could not determine all parameters')
                nonlinear_has_error = ';param_overflow_during_fitting'
                self.manip.logger(self.filename, 'Overflow')
            except RuntimeError:
                print('!!!! Overflow detected on one of the parameters. Could not determine all parameters')
                nonlinear_has_error = ';param_overflow_during_fitting'
                self.manip.logger(self.filename, 'Overflow')
            except OverflowError:
                print('!!!! Overflow detected on one of the parameters.')
                self.manip.logger(self.filename, 'Overflow')

            perr = np.sqrt(np.diag(pcov))
            fittings.append((first_point, popt, perr))

            if self.settings.DEBUG:
                fitting_params_str = ' '.join([str(round(i, 2)) + '+/-' +
                                               str(round(j,2)) for i, j in zip(popt, perr) ])
                                                # 'a+/-aerr b+/-berr ...'
                print(f"{self.settings.NL_FITTING_METHOD} fitting: {fitting_params_str}")

        if self.settings.NL_SORTING_METHOD == 'eta_0':
            fittings.sort(key=lambda x: x[2][0])
        elif self.settings.NL_SORTING_METHOD == 'overall':
            #fittings.sort(key=lambda x: x[2][0] + x[2][1] + x[2][2] + x[2][3])
            fittings.sort(key=lambda x: sum(x[2])) # sums the errors
        else:
            raise ValueError(f'Could not understand the sorting method {self.settings.NL_SORTING_METHOD}')

        self.nl_first_point = fittings[0][0]
        self.params = fittings[0][1]
        self.param_errs = fittings[0][2]

        if save:  # todo: check here to return a good destination file
            try:
                self.manip.record_fit(
                    self.filename, self.params[0],
                    self.param_errs[0], silent=False,
                    extra=f"{fittings[0][0]};{fittings[0][1]};nonlinear_auto_{self.settings.NL_FITTING_METHOD};"
                          f"{nonlinear_has_error}", fdest_name=self.settings.NL_FITTING_METHOD + '.csv'
                )
            except UnboundLocalError:
                print('Unable to write to file because the subroutine did not return the fitting parameters')
                print(traceback.format_exc())
                self.manip.record_fit(self.filename, 0, 0, extra=f'nonlinear_auto_{self.settings.NL_FITTING_METHOD};'
                                                                 f'unable_to_find_viscosity',
                                      fdest_name=self.settings.NL_FITTING_METHOD+'.csv')
                self.manip.logger(self.filename, 'No Viscosity')

        self.nl_done = True
        return self.nl_first_point, self.params, self.param_errs

    # TODO: check if the bounds are correct
    # TODO: increment this function to be able to accept multiple fittings
    def manual_fit(self, first, last, fit_types, save=True):
        GP_arr = np.array(self.GP[first:last + 1])
        Eta_arr = np.array(self.Eta[first:last + 1])
        fittings = []

        for type in fit_types:
            if 'Linear' in type:
                popt, pcov = curve_fit(self.fit_lin, GP_arr, Eta_arr, p0=(30, 0),
                                       bounds=(0, [self.VISC_LIMIT, 0.0001]))
            elif 'Carreau' in type:
                popt, pcov = curve_fit(self.fit_Carreau, GP_arr, Eta_arr, p0=(30, 0),
                                       bounds=(0, np.inf))
            elif 'Cross' in type:
                popt, pcov = curve_fit(self.fit_Carreau, GP_arr, Eta_arr, p0=(30, 0),
                                       bounds=(0, np.inf))
            elif 'Carreau-Yasuda' in type:
                popt, pcov = curve_fit(self.fit_CarreauYasuda, GP_arr, Eta_arr, p0=(30, 0),
                                       bounds=(0, np.inf))
            else:
                raise NameError(f'Could not understand the list fit_types {fit_types}')

            perr = np.sqrt(np.diag(pcov))

            self.params = popt  # Will be continuously overwritten. todo: will this be a problem?
            self.param_errs = perr

            if self.settings.DEBUG:
                # 'a+/-aerr b+/-berr ...'
                fitting_params_str = ' '.join([str(round(i, 2)) + '+/-' +
                                               str(round(j, 2)) for i, j in zip(popt, perr)])
                print(f"{fit_types} fitting: {fitting_params_str}")

            if save:  # todo: check here to return a good destination file
                self.manip.record_fit(self.settings.NL_FITTING_METHOD, self.params,
                                      self.param_errs, silent=False)

            fittings.append((type, popt, perr))

        return fittings

    def plot_error_graphs(self): # todo: If it has both plots, make them side by side
        TEXT_FILENAME_X = 0.1
        TEXT_PARAMS_X = 0.3
        TEXT_Y = 0.98
        x = np.logspace(np.log10(self.GP[0]), np.log10(self.GP[-1]))

        if self.settings.DEBUG:
            print('Debug: x', x)
            print('Debug: params', self.params)
            print('Debug: GP', self.GP, 'Eta', self.Eta)

        if self.nl_done:
            if self.model == 'Carreau':
                y, yerr = self.carr_uncertainty(x, *self.params, *self.param_errs)
            elif self.model == 'Cross':
                y, yerr = self.cross_uncertainty(x, *self.params, *self.param_errs)
            elif self.model == 'Carreau-Yasuda':
                y, yerr = self.carryas_uncertainty(x, *self.params, *self.param_errs)
        if self.lin_done:
            y_l, yerr_l = np.ones(len(x)) * self.int, np.ones(len(x)) * self.int_err
            # Creates a horizontal line with n points

        if self.nl_done and self.lin_done:
            fig, [axn, axl] = plt.subplots(ncols=2, nrows=1, figsize=(12, 4))
        elif self.nl_done and not self.lin_done:
            fig, axn = plt.subplots(ncols=1, nrows=1, figsize=(6, 4))
        elif not self.nl_done and not self.lin_done:
            fig, axl = plt.subplot(ncols=1, nrows=1, figsize=(6, 4))

        if self.nl_done:
            axn.set_xscale('log')
            axn.set_yscale('log')
            axn.plot(self.GP, self.Eta, linewidth=0, marker='o', markersize=5)
            axn.errorbar(x, y, yerr=yerr)
            axn.annotate(str(self.nl_first_point + 1), (self.GP[self.nl_first_point], self.Eta[self.nl_first_point]), color='red')
            if self.nl_last_point == -1:
                axn.annotate(str(len(self.GP)), (self.GP[self.nl_last_point], self.Eta[self.nl_last_point]),
                             color='red')  # todo: check this function
            else:
                axn.annotate(str(self.nl_last_point), (self.GP[self.nl_last_point], self.Eta[self.nl_last_point]), color='red')
            model_param_names = 'Model: ' + self.model + ' Params: ' + self.param_names
            param_text = " ".join([str(round(par, 2)) + '+/-' + str(round(err, 2))
                                   for par, err in zip(self.params, self.param_errs)])

            fig.text(TEXT_FILENAME_X, TEXT_Y, self.filename, size='small')
            fig.text(TEXT_FILENAME_X, TEXT_Y - 0.030, model_param_names, size='small')
            fig.text(TEXT_FILENAME_X, TEXT_Y - 0.060, param_text, size='small')
            fig.text(TEXT_FILENAME_X, TEXT_Y - 0.090, f'R2 = {round(self.nl_R2, 2)}', color='red', size='small')
            # fig.text(TEXT_FILENAME_X, TEXT_Y - 0.060, param_values, size='small')
            # fig.text(TEXT_FILENAME_X, TEXT_Y - 0.090, param_errors, size='small')
            #fig.set_size_inches(6, 4)  # todo: check if this is still necessary
            TEXT_FILENAME_X = 0.5  # Changes it to half in case a linear plot is done.

        if self.lin_done:
            axl.set_xscale('log')
            axl.set_yscale('log')
            axl.plot(self.GP, self.Eta, linewidth=0, marker='o', markersize=5)
            axl.errorbar(x, y_l, yerr=yerr_l)
            axl.annotate(str(self.l_first_point + 1), (self.GP[self.l_first_point], self.Eta[self.l_first_point]), color='red')
            if self.l_last_point == -1:
                axl.annotate(str(len(self.GP)), (self.GP[self.l_last_point], self.Eta[self.l_last_point]),
                             color='red')  # todo: check this function
            else:
                axl.annotate(str(self.l_last_point), (self.GP[self.l_last_point], self.Eta[self.l_last_point]), color='red')
            model_param_names = 'Model: Linear. Params: Intercept'
            param_text = f"int = {self.int}+/-{self.int_err}"

            fig.text(TEXT_FILENAME_X, TEXT_Y, self.filename, size='small')
            fig.text(TEXT_FILENAME_X, TEXT_Y - 0.030, model_param_names, size='small')
            fig.text(TEXT_FILENAME_X, TEXT_Y - 0.060, param_text, size='small')
            fig.text(TEXT_FILENAME_X, TEXT_Y - 0.090, f'R2 = {round(self.l_R2, 2)}', color='red', size='small')
            # fig.text(TEXT_FILENAME_X, TEXT_Y - 0.060, param_values, size='small')
            # fig.text(TEXT_FILENAME_X, TEXT_Y - 0.090, param_errors, size='small')
            #fig.set_size_inches(6, 4)  # todo: check if this is still necessary

        #fig.tight_layout()

        if self.settings.SAVE_GRAPHS:
            fig.savefig(self.filename[:-4] + '.png')
            print('Figure saved.')
        if not self.settings.INLINE_GRAPHS and self.settings.PLOT_GRAPHS:
            plt.draw()
            plt.pause(self.wait)
            #plt.clf()
            plt.close(fig)
        elif self.settings.PLOT_GRAPHS:
            plt.show()
        return


class FileManip:
    #def __init__(self, sett):
    #    self.settings = sett

    @staticmethod
    def ExtractData(fname, FC_segment=0):
        """Opens the file fname and extracts the data based on where it finds the word 'Eta' and 'GP', these being
        the Viscosity and the Shear Rate (gamma point). If the file has multiple segments, for example, when multiple
        experiments were done in succession, FC_segment indicates which of those experiments was a Flow Curve."""
        fhand = open(fname, 'r')
        GP = []
        Eta = []
        column_eta = 0
        column_gp = 0
        # FC_segment = '3'

        # while FC_segment == 0:
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
                # if settings['DEBUG']:
                #     print('Debug: column names', column_names)
                for i, column in enumerate(column_names):
                    if 'Eta' in column and 'Eta*' not in column:
                        column_eta = i
                        #if settings['DEBUG']:
                        #    print('Debug: Found Eta at', column_eta)
                    if 'GP' in column:
                        column_gp = i
                        #if settings['DEBUG']:
                        #    print('Debug: Found GP at', column_gp)
            try:
                GP.append(float(line.replace(',', '.').split(';')[column_gp]))
                Eta.append(float(line.replace(',', '.').split(';')[column_eta]))
            except:
                pass

            # if line.startswith(FC_segment + '|'):
            #    line = line.rstrip()
            #    num, gp, tau, eta, *rest = line.replace(',','.').split(';')
            #    GP.append(float(gp))
            #    Eta.append(float(eta))
            #    #print(line)

        fhand.close()
        if len(GP) == 0:
            # print('!!!!No Flow Curve data was found! Re-export the data on file', fname)
            raise ValueError
        # return pd.Series(GP), pd.Series(Eta)
        # if settings['DEBUG']:
        #     print('Debug: Extracted Data: GP:', GP, 'Eta:', Eta)
        return GP, Eta

    @staticmethod
    def ExtractData_pd(fname):
        """Uses pandas do extract the data if it was exported using the data extraction tool"""
        pd_temp = pd.read_csv(fname, delimiter=';', encoding='latin1', decimal=',')
        pd_temp = pd_temp[pd_temp > 0].dropna()

        col_GP = ''
        col_Eta = ''

        for col in pd_temp.columns:
            if 'GP' in col:
                col_GP = col
                # print('achou GP em', col)
            if 'Eta' in col:
                col_Eta = col
                # print('achou Eta em', col)

        GP = pd_temp[col_GP].tolist()
        Eta = pd_temp[col_Eta].tolist()
        return GP, Eta

    @staticmethod
    def record_fit(name, eta0, eta0_err, silent=False, extra='', fdest_name='results.csv'):
        if not silent:
            print(f"{name}: Intercept={eta0} +- {eta0_err}. Extra={extra}")
            #print(name + ':', 'Intercept', eta0, '+-', eta0_err, extra)

        with open(fdest_name, 'a', encoding='utf-8') as fdest:
            #fdest.write(name + ';' + str(eta0) + ';' + str(eta0_err) + ';' + extra + '\n')
            fdest.write(f"{name};{eta0};{eta0_err};{extra}\n")

    @staticmethod
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

    @staticmethod
    def logger(file, type, extra=''):
        with open('log', 'a') as log:
            if type == 'Overflow':
                log.write(f'Parameter overflow while trying to fit file {file}: {extra}\n')
            if type == 'No Viscosity':
                log.write(f'Unable to find viscosity for file {file}\n')
            if type == 'Failed to open':
                log.write(f'Failed to open file {file}. Re-export the data.')
            else:  # type == 'Generic'
                log.write(f'Error while processing {file}: {extra}\n')


def test():
    settings = Settings.Settings()
    settings.NL_FITTING_METHOD = 'Carreau-Yasuda'
    filename = 'CF_Sac50-3--0.csv'
    fit = Fitter(filename, settings, do_fit=False)
    fit.automatic_nl_fitting_lm(save=True)
    print(fit.model, fit.nl_R2, *fit.params)

    return fit


def main():
    settings = Settings.Settings()
    manip = FileManip()
    settings.print_settings()
    do_change = input('Do you want to change the settings? y/[n]')
    if do_change == 'y':
        settings.edit_settings()

    if settings.TREAT_ALL:
        files = glob.glob(f'*.{settings.EXT}')
        if len(files) == 0:
            print(f'No files with the extension {settings.EXT} found.'
                  f' Please select them manually or change EXT accordingly.')
            files = manip.select_files()

    else:
        files = manip.select_files()

    if len(files) == 0:
        print('No files selected. Quitting.')
        sys.exit()

    for file in files:
        try:
            fit = Fitter(file, settings, do_fit=True)
        except ValueError:  # todo: debug and check what would be needed here.
            print(f'Skipping {file}: Value Error')
            continue
        except KeyError:
            print(f'Skipping {file} Key Error')
            continue
            #print(traceback.format_exc())

        if settings.PLOT_GRAPHS or settings.SAVE_GRAPHS:
            try:
                fit.plot_error_graphs()
            except OverflowError:  # todo: write which parameter has overflown
                print('!!!! Overflow detected on one of the parameters. Could not plot the data')
                nonlinear_has_error = ';param_overflow_during_fitting'
                # todo: log this
            except UnboundLocalError:
                print('Not able to write to file because the subroutine did not return the fitting parameters')
                # todo: log this


            #fit.plot_error_graphs(file[:-4] + '_lin_' + file[-4:], fit.params, fit.first_point, fit.last_point,
            #                      model=fit.settings.NL_FITTING_METHOD, param_names=[''])

                # # except:
                # #     print('Error found while plotting the linear fit')
                # #     print(traceback.format_exc())
                # #     lin_has_error = 'error_during_fitting'
                # record(file, a, aerr,
                #        extra='linear automatic;FP=' + str(lin_points[0]) + 'LP=' + str(lin_points[1]) +
                #              lin_has_error, fdest_name='linear.csv')


            #         if settings['PLOT_GRAPHS']:
            #             plot_error_graphs(file[:-4] + '_carr_' + file[-4:], GP, Eta,
            #                               params=np.concatenate((popt, perr)),
            #                               first_point=nl_first, model=settings['NL_FITTING_METHOD'],
            #                               param_names=Param_names_errs[settings['NL_FITTING_METHOD']])
            #     except OverflowError:  # todo: write which parameter has overflown
            #         print('!!!! Overflow detected on one of the parameters. Could not plot the data')
            #         nonlinear_has_error = ';param_overflow_during_fitting'
            #     try:
            #         record(file, popt[0], perr[0], extra='nonlinear_auto_' + settings['NL_FITTING_METHOD'] +
            #                                              nonlinear_has_error,
            #                fdest_name=settings['NL_FITTING_METHOD'] + '.csv')
            #     except UnboundLocalError:
            #         print('Not able to write to file because the subroutine did not return the fitting parameters')
            #         record(file, 0, 0, extra='nonlinear_auto_' + settings['NL_FITTING_METHOD'] + ';' +
            #                                  'unable_to_find_viscosity',
            #                fdest_name=settings['NL_FITTING_METHOD'] + '.csv')
            #         with open('log', 'a') as log:
        #     #             log.write('Unable to find viscosity for file ' + file + '\n')
        # except:
        #     print('!!!!We have encountered an error while processing file', file)
        #     print(traceback.format_exc())
        #     with open('log', 'a') as log:
        #         log.write('Error while processing ' + file + '\n')
        #         log.write(traceback.format_exc())

if __name__ == '__main__':
    #fit = test()
    main()