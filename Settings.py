import os

# todo: remove debugging feature. Just use the built in debugging tools.

class Settings:
    def __init__(self, debug=False):
        self.valid_settings = ['DEBUG', 'INLINE_GRAPHS', 'LIN_SORTING_METHOD',
                               'NL_FITTING_METHOD', 'NL_SORTING_METHOD',
                               'PLOT_GRAPHS', 'SAVE_GRAPHS', 'AUTO_LIN', 'AUTO_NL', 'DO_LIN',
                               'DO_NL', 'TREAT_ALL', 'EXT', 'PREV_EXTRACTED', 'WAIT', 'FIXED_FP_NL',
                               'MAX_FP_NL']
        self.valid_options_lin_sorting = ['by_error', 'by_error_length', 'by_R2']
        self.valid_options_nl_sorting = ['eta_0', 'overall', 'R2']
        self.models = ['Carreau', 'Cross', 'Carreau-Yasuda']
        if debug:
            self.DEBUG = True
        else:
            self.DEBUG = False

        try:
            self.load_settings()
        except FileNotFoundError:
            print('Settings file not found. Loading defaults.')
            self.DEBUG = False
            self.INLINE_GRAPHS = False
            self.LIN_SORTING_METHOD = 'by_error_length'  # by_error, by_error_length, by_R2
            self.NL_FITTING_METHOD = 'Carreau'  # Carreau, Cross, Carreau-Yasuda
            self.NL_SORTING_METHOD = 'overall'  # eta_0, overall, R2
            self.PLOT_GRAPHS = False
            self.SAVE_GRAPHS = False
            self.AUTO_LIN = True
            self.AUTO_NL = True
            self.DO_LIN = True
            self.DO_NL = True
            self.TREAT_ALL = False
            self.EXT = 'txt'
            self.PREV_EXTRACTED = False
            self.WAIT = '0.5'
            self.FIXED_FP_NL = True
            self.MAX_FP_NL = 2
            print('Creating a new settings file with the defaults')
            self.save_settings()

    def save_settings(self):
        with open('settings.dat', 'w') as settings_file:
            settings_file.write('#This is an automatically created settings file.\n')

            settings_file.write('\n# Inline graphs is meant for those with IDEs like Jupyter Notebooks\n')
            settings_file.write('INLINE_GRAPHS=' + str(self.INLINE_GRAPHS))

            settings_file.write(
                '\n# Plot graphs after fitting, with error propagation? Slows down the process greatly.\n')
            settings_file.write('PLOT_GRAPHS=' + str(self.PLOT_GRAPHS))

            settings_file.write('\n# Save graphs after fitting? Useless if PLOT_GRAPHS is set to False\n')
            settings_file.write('SAVE_GRAPHS=' + str(self.SAVE_GRAPHS))

            settings_file.write('\n# When plotting, time it waits until the next plot is shown,'
                                ' in seconds\n')
            settings_file.write('WAIT=' + str(self.WAIT))

            settings_file.write('\n# Treat all files in folder?\n')
            settings_file.write('TREAT_ALL=' + str(self.TREAT_ALL))

            settings_file.write('\n# Extension of files to look for\n')
            settings_file.write('EXT=' + str(self.EXT))

            settings_file.write('\n# If the important data has been extracted\n')
            settings_file.write('PREV_EXTRACTED=' + str(self.PREV_EXTRACTED))

            settings_file.write('\n\n#### Linear Fitting ####')

            settings_file.write('\n# Perform linear fit?\n')
            settings_file.write('DO_LIN=' + str(self.DO_LIN))

            settings_file.write(
                "\n# Sorting method of the automatic linear method.\n# Can be 'by_error', minimizing the " +
                "error, 'by_error_length', which minimizes the error divided by the total number of " +
                "points used, 'by_R2', which minimizes R squared.\n")
            settings_file.write('LIN_SORTING_METHOD=' + str(self.LIN_SORTING_METHOD))

            settings_file.write('\n# Set to True if you want the linear fitting to be done automatically. ' +
                                'False, if to be done manually.\n')
            settings_file.write('AUTO_LIN=' + str(self.AUTO_LIN))

            settings_file.write('\n\n#### Non-Linear Fitting ####')

            settings_file.write('\n# Perform non-linear fitting?\n')
            settings_file.write('DO_NL=' + str(self.DO_NL))

            settings_file.write("\n# Fitting method. 'Carreau', 'Cross', 'Carreau-Yasuda'\n")
            settings_file.write('NL_FITTING_METHOD=' + str(self.NL_FITTING_METHOD))

            settings_file.write("\n# Can be 'overall', minimizing the error of all parameters,  'eta_0', " +
                                "minimizing the error of only this parameter, or 'R2'.\n")
            settings_file.write('NL_SORTING_METHOD=' + str(self.NL_SORTING_METHOD))

            settings_file.write('\n# Set to True if you want the non linear fitting to be done automatically\n')
            settings_file.write('AUTO_NL=' + str(self.AUTO_NL))

            settings_file.write('\n# Set to true if the first point should be fixed during fitting\n')
            settings_file.write('FIXED_FP_NL=' + str(self.FIXED_FP_NL))

            settings_file.write('\n# If FIXED_FP_NL is true, how much can it travel? This must be in terms of the'
                                'inverse of the length, that is, length/MAX_FP_NL. 1 would be the entire curve (not'
                                'recommended. 2 would be up to half the curve, etc.\n')
            settings_file.write('MAX_FP_NL=' + str(self.MAX_FP_NL))

            settings_file.write('\n\n##### Debug #####\n')
            settings_file.write('# Show debug messages\n')
            settings_file.write('DEBUG=' + str(self.DEBUG))

    def load_settings(self):
        fhand = open('settings.dat', 'r')

        for line in fhand:
            if line.startswith('#'):
                continue
            if len(line) < 4:
                continue
            line = line.rstrip()
            line = line.replace(' ', '')
            param, value = line.split('=')
            if self.DEBUG:
                print(line)
                print(f'Var = {param} val= {param}')
            # todo: think about using getattr

            if (not hasattr(self, param)) and (param in self.valid_settings):  # If param wasn't loaded and is valid

                if value.lower() == 'true':  # Assigns booleans first
                    setattr(self, param, True)
                elif value.lower() == 'false':
                    setattr(self, param, False)
                else:
                    setattr(self, param, value)  # todo: check if the parameter is valid, i.e. Carreau

            elif param not in self.valid_settings:
                print(f"Settings file has an unrecognized parameter {param} which wasn't loaded.")
                continue

        fhand.close()

    def print_settings(self):
        counter = 0
        valid_numbers = []
        number_setting_corr = {}  # To assign a number to a parameter. Ex: 1: 'EXT'

        for param in self.valid_settings:
            sett = getattr(self, param)
            print(f"{counter}) {param} = {sett}", end='')
            if type(sett) == bool:
                print(': Options= True | False')
            elif param == 'LIN_SORTING_METHOD':
                print(': Options= ', end='')
                print(*self.valid_options_lin_sorting, sep=' | ')
            elif param == 'NL_FITTING_METHOD':
                print(': Options= ', end='')
                print(*self.models, sep=' | ')
            elif param == 'NL_SORTING_METHOD':
                print(': Options= ', end='')
                print(*self.valid_options_nl_sorting, sep=' | ')
            elif param == 'EXT':
                print(': Options = txt | dat | csv | etc...')
            elif param == 'WAIT':
                print(': Options = time in seconds')
            elif param == 'MAX_FP_NL':
                print(': Options = 1/n. n=1: whole curve. n=2: half')
            else:
                print('\n')
            valid_numbers.append(str(counter))
            number_setting_corr[str(counter)] = param
            counter += 1

        print('\n===================\n')
        if self.DEBUG:
            print('Valid numbers', valid_numbers)
            print('Correlation', number_setting_corr)

        return valid_numbers, number_setting_corr

    def edit_settings(self):
        """Edits the current settings"""
        while True:
            os.system('cls' if os.name == 'nt' else 'clear')
            valid_numbers, number_setting_corr = self.print_settings()
            print('Which setting you want to change? Enter "number, new value" to modify, or "done" to exit.')
            print('Observe the possible values for each setting! They are case sensitive. '
                  'Inputting wrong values might break the program. \n')
            choice = input('Input:')
            if choice == 'done':
                break
            if ',' not in choice:
                print('Invalid input. Place the number, followed by a comma, followed by its value. Eg: 1,TRUE')
                continue
            if len(choice.split(',')) != 2:
                print('Invalid input, must have only one comma')
                continue

            var, val = choice.split(',')
            if var not in valid_numbers:
                print('Invalid number.')
                continue
            real_var = number_setting_corr[var]  # Changes from a number to the actual parameter
            if val.lower() == 'true':
                setattr(self, real_var, True)
                continue
            elif val.lower() == 'false':
                setattr(self, real_var, False)
                continue
            else:
                setattr(self, real_var, val)

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
        _, _ = self.print_settings()
        self.save_settings()
        return





#
# settings = {
#     'DEBUG': False,
#     'INLINE_GRAPHS': False,
#     'SORTING_METHOD_LIN': 'by_error_length',  # by_error, by_error_length
#     'NL_FITTING_METHOD': 'Carreau',  # Carreau, Cross, Carreau-Yasuda
#     'SORTING_METHOD_NL': 'overall',  # eta_0, overall
#     'PLOT_GRAPHS': False,
#     'SAVE_GRAPHS': False,
#     'AUTO_LIN': True,
#     'AUTO_NL': True,
#     'DO_LIN': True,
#     'DO_NL': True,
#     'TREAT_ALL': False,
#     'EXT': 'txt',
#     'PREV_EXTRACTED': False
# }
#
# valid_options_lin_sorting = ['by_error', 'by_error_length']
# valid_options_nl_sorting = ['eta_0', 'overall']
# models = ['Carreau', 'Cross', 'Carreau-Yasuda']