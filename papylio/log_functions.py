import inspect
import json
from datetime import datetime
import papylio

def function_arguments(function, function_locals):
    signature_values = inspect.signature(function).parameters.values()
    all_argument_values = {
        parameter.name:
            function_locals[parameter.name] for parameter
        in signature_values if parameter.name is not 'self'
    }

    if list(all_argument_values.keys())==['configuration']:
        all_argument_values = all_argument_values['configuration']

    return all_argument_values

def function_arguments_json(function, function_locals):
    return json.dumps(function_arguments(function, function_locals))

def get_current_datetime():
    current_datetime = datetime.now()
    return current_datetime.strftime("%Y-%m-%d %H:%M:%S")

def add_configuration_to_dataarray(dataarray, function=None, function_locals=None, units=None):
    dataarray.attrs['version'] = papylio.__version__
    if function is not None:
        dataarray.attrs['configuration'] = function_arguments_json(function, function_locals)
    dataarray.attrs['datetime'] = get_current_datetime()
    if units is not None:
        dataarray.attrs['units'] = units
    return dataarray