"""

 Script name: gwLogger.py

 Description:
    This class contains a custom logger

"""

import logging
import time
import sys
import os

# Now, some additional information

__author__ = "Massimiliano Razzano"
__copyright__ = "Copyright 2016-2020, Massimiliano Razzano"
__credits__ = ["Line for credits"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "M. Razzano"
__email__ = "massimiliano.razzano@unipi.it"
__status__ = "Production"


class GWLogger(logging.Logger):
    """

    This class contains a custom logger

    """

    def __init__(self, name, m_loglevel="INFO",m_file_logging=False,m_log_dir=None):
        # init start time
        self.__name = name
        self.__log_filename = None
        self.__log_dir = m_log_dir
        self.__start_time = str(time.strftime("%y%m%d_%Hh%Mm%Ss", time.localtime()))

        logging.Logger.__init__(self, name)

        base_handler = logging.StreamHandler(sys.stdout)
        base_handler.setLevel(logging.DEBUG)

        formatter = logging.Formatter('%(asctime)s %(levelname)-4s: %(message)s', datefmt='%H:%M:%S')

        base_handler.setFormatter(formatter)
        self.addHandler(base_handler)

        if m_file_logging:
            if not os.path.exists(m_log_dir):
                os.mkdir(m_log_dir)

            #add a file logger
            if m_log_dir is None:
                file_handler_name = self.__name+"_"+self.__start_time+".log"
            else:
                file_handler_name = os.path.join(m_log_dir,self.__name+"_"+self.__start_time+".log")

            try:
                file_handler = logging.FileHandler(file_handler_name)
                self.__log_filename = file_handler_name
                self.addHandler(file_handler)
            except OSError:
                self.error("Cannot open the log file")


        self.set_loglevel(m_loglevel=m_loglevel)

        if m_file_logging:
            self.info("Log file saved to "+file_handler_name)

        self.info('Logger started at ' + self.__start_time + " (LEVEL=" + str(self.getEffectiveLevel()) + ")")

        #self.__stdout_logger = logging.getLogger('STDOUT')
        #self.__sl = StreamToLogger(self.__stdout_logger, logging.INFO)
        #sys.stdout =self.__sl

    def get_log_filename(self):
        """

        :return:
        """

        return self.__log_filename


    def set_loglevel(self,m_loglevel):
        """

        :param m_loglevel:
        :return:

        """

        if m_loglevel=="DEBUG":
            self.setLevel(logging.DEBUG)
        if m_loglevel == "INFO":
            self.setLevel(logging.INFO)

        self.info("Logging level set to "+m_loglevel)


if __name__ == '__main__':
    print("\n****\nTesting Logger...\n****\n")

    my_logger = GWLogger("testlogger")
    my_logger.setLevel(logging.DEBUG)
    my_logger.debug("This is a debug message")

    my_logger.setLevel(logging.INFO)
    my_logger.debug("This is an info message")

    my_logger.setLevel(logging.ERROR)
    my_logger.debug("...and this is an ERROR message!")
