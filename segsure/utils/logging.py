"""Logging utilities for SegSure."""

import logging
import os
from datetime import datetime
from pathlib import Path
from typing import Optional


def setup_logger(
    name: str,
    output_dir: Optional[str] = None,
    level: int = logging.INFO,
    verbose: bool = False,
) -> logging.Logger:
    """Set up a logger with both console and file handlers.
    
    Parameters
    ----------
    name : str
        Name of the logger.
    output_dir : str, optional
        Directory to save log files. If None, only console output.
    level : int, optional
        Logging level. Default is logging.INFO.
    verbose : bool, optional
        If True, set level to DEBUG. Default is False.
    
    Returns
    -------
    logging.Logger
        Configured logger instance.
    """
    if verbose:
        level = logging.DEBUG
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_format = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    console_handler.setFormatter(console_format)
    
    if logger.handlers:
        logger.handlers.clear()
    
    logger.addHandler(console_handler)
    
    # File handler (if output_dir provided)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        log_file = os.path.join(
            output_dir, f"{name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        )
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_format = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        file_handler.setFormatter(file_format)
        logger.addHandler(file_handler)
    
    return logger


def get_logger(name: str) -> logging.Logger:
    """Get a logger by name.
    
    Parameters
    ----------
    name : str
        Name of the logger.
    
    Returns
    -------
    logging.Logger
        Logger instance.
    """
    return logging.getLogger(name)
