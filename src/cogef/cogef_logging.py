""" handling logging in module cogef """

import logging


def logging_init( log_file = "job_cogef.log", log_level = "DEBUG", stream_level = "INFO"):
    level = logging.getLevelName(log_level)
    slevel = logging.getLevelName(stream_level)

    # init the logging machinery and return a logger
    logger = logging.getLogger()
    logger.setLevel(level)

    # configure file handler
    fhandler  = logging.FileHandler(log_file)
    fhandler.setLevel(level)

    # configure console handler
    shandler = logging.StreamHandler()
    shandler.setLevel(slevel)

    # add formatter
    formatter = logging.Formatter('%(asctime)s - [%(name)8s] - [%(levelname)8s] : %(message)s')
    fhandler.setFormatter(formatter)
    shandler.setFormatter(formatter)

    # add handler to logger
    logger.addHandler(fhandler)
    logger.addHandler(shandler)







