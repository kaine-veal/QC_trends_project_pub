import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path    

def create_logger(log_filename='default.log'):
    """
    Creates a logger that writes to a specific log file in the project's logs directory.
    
    Args:
        log_filename (str): Name of the log file (e.g., 'run_QC.log', 'sampleQC.log')
    
    Returns:
        logging.Logger: Configured logger instance
    """
    # Get the project root: Utils folder is at project root, so just go up one level
    # <repo>/package/utils/logger.py -> <repo>/package
    project_root = Path(__file__).resolve().parent.parent
    
    # Create logs directory if it doesn't exist
    logs_dir = project_root / 'logs'
    logs_dir.mkdir(exist_ok=True)
    
    # Path to the specified log file
    log_file = logs_dir / log_filename
    
    # Create a unique logger name based on the log file to avoid conflicts
    logger_name = f'logger_{log_filename.replace(".", "_")}'
    logger = logging.getLogger(logger_name)
    
    # Prevent duplicate handlers if logger already exists
    if logger.hasHandlers():
        return logger
    
    # Setting the root default logger level to DEBUG to capture everything.
    # DEBUG is the lowest built-in level, meaning it will 'see' every message (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    logger.setLevel(logging.DEBUG) 

    # Creates a stream handler
    stream_handler = logging.StreamHandler() # Creates a streamhandler to print to console
    stream_handler.setLevel(logging.INFO) # Only handle messages DEBUG or higher
    stream_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s") # This defines how log messages should look on the console
    stream_handler.setFormatter(stream_formatter) # Apply the above format to the stream handler

    # Creates a rotating file handler - writes to logs/{log_filename}
    file_handler = RotatingFileHandler(str(log_file),
                                       maxBytes=500000,
                                       backupCount=2)
    file_handler.setLevel(logging.DEBUG) # Only write INFO and above to the log file
    file_formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s") # Same formatting as stream handler
    file_handler.setFormatter(file_formatter) # Apply the above format to the file handler

    # Attach both handlers to the logger
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    return logger

