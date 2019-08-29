#### logging ######
import logging
import sys
# Create the Logger
# logging.basicConfig(filename='rpRanker.log',
#                             filemode='a',
#                             format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
#                             datefmt='%H:%M:%S',
#                             level=logging.INFO)
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
# loggers = logging.getLogger(__name__)
# '''
# # Create a Formatter for formatting the log messages
# logger_formatter = logging.Formatter('%(asctime)s,%(msecs)d %(name)s - %(levelname)s - %(message)s')
#
# # Add the Formatter to the Handler
# logger_handler.setFormatter(logger_formatter)
#
# # Add the Handler to the Logger
# loggers.addHandler(logger_handler)
# '''
# loggers.info('Completed configuring logger()!')

from .rpReader import rpReader
from .rpCache import rpCache
from .rpCofactors import rpCofactors
from .rpSBML import rpSBML
from .rpFBA import rpFBA

try:
    from .rpThermo import rpThermo
    logging_thermo = "Thermodynamics is enabled"
except ModuleNotFoundError:
    logging_thermo = "Thermodynamics is NOT enabled"
from .tools import tools
try:
    from . import component_contribution
    logging_component_contribution = "component_contribution is enabled"
except ModuleNotFoundError:
    logging_component_contribution = "component_contribution is NOT enabled"
#from .component_contribution import *
