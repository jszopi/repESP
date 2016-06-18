import argparse


def input_type(input_charge_type):
    if input_charge_type == 'list':
        result = 'qout'
    elif input_charge_type == 'aim':
        result = 'sumviz'
    else:
        result = 'log'
    return result


parser = argparse.ArgumentParser(add_help=False)

charge_choices = ['list', 'mulliken', 'nbo', 'aim', 'mk', 'chelp', 'chelpg',
                  'hly']
charge_choices_text = "'" + "', '".join(charge_choices[1:-1]) + "' and '" + \
                      charge_choices[-1] + "'."
parser.add_argument(
    "input_charge_type",
    help="type of input charges. To extract charges from a text file, use the "
    "'list' option. Other available options are: " + charge_choices_text,
    choices=charge_choices, metavar="INPUT_CHARGE_TYPE")

parser.add_argument(
    "input_charge_file",
    help="""input file to extract charges from. This can be a Gaussian .log,
    AIMAll .sumviz or a text file. The appropriate format for the text file is
    used throughout this suite of scripts and is very simple:

    1. Values should be space-separated
    2. Atom order should follow that in the template cube file
    3. Line breaks are not significant""",
    metavar="INPUT_CHARGE_FILE")
