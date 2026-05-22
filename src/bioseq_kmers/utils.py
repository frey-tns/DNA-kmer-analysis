import argparse
import math
import os
# For shell command string manipulation
import shlex

#################################################
#   Function: Defined the format command line   #
#################################################
def format_command_line(argv):
    """
     Format command line by changing the absolute path to a relative path.

     Args:
         argv (list): Command-line arguments.

     Returns:
         str: Reconstructed command line.
     """
    # Retrieves the folder from which the script is executed
    cwd = os.getcwd()
    # Contains the rebuilt command
    list_cleaned_command = []

    # Executable name only
    program_name = os.path.basename(argv[0])
    list_cleaned_command.append(shlex.quote(program_name))

    # The next argument is a file path
    skip_next = False

    # Iterates through each element of the order
    for arg in argv[1:]:

        # If argument associated with -i or -o
        if skip_next:
            # Converts absolute path to relative path
            rel = os.path.relpath(arg, cwd)
            list_cleaned_command.append(shlex.quote(rel))
            # Returns to the initial state
            skip_next = False

        elif arg in ["-i", "--input", "-o", "--output"]:
            # Keep the current flag
            list_cleaned_command.append(arg)
            # The next argument is a path
            skip_next = True

        else:
            # Normal argument
            list_cleaned_command.append(shlex.quote(arg))

    return " ".join(list_cleaned_command)

########################################
#   Function: positive integer value   #
########################################
def min_interger(min_value):
    """
    Create a validator for integers greater than or equal to `min_value`.

    Args:
        value(str): Command line argument.
    Returns:
        int: Validated positive integer.
    Raises:
        argparse.ArgumentTypeError: If `value` is not an integer greater than or equal to min value.
    """
    def validator(value):
        # Intercept a potential error
        try:
            ## Verify that k length is an integer
            # Convert the argument to the entire file
            value= (int(value))
        except ValueError:
            # If invalid value → Intentionally cause an error suitable for argparse
            raise argparse.ArgumentTypeError(f"Value {value!r} is not an integer.")

        if value < min_value:
            raise argparse.ArgumentTypeError(f"Value {value!r} must be >= {min_value}.")

        return value

    return validator

#####################################
#   Function: scientific notation   #
#####################################
def engineer_mode(log10_value, precision=3):
    """
    Convert a log10 probability into scientific notation string

    """
    # handle infinities
    if log10_value == float("-inf"):
        return "-inf"

    if math.isnan(log10_value):
        return "nan"

    exponent = math.floor(log10_value)
    mantissa = 10 ** (log10_value - exponent)

    return f"{mantissa:.{precision}f}e{exponent}"