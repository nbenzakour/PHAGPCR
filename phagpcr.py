#!/usr/bin/env python3
"""
Backward-compatible wrapper for PHAGPCR.

This script maintains compatibility with existing usage while delegating
to the modular phagpcr_package implementation.
"""

import sys
from phagpcr_package.main import main

if __name__ == '__main__':
    sys.exit(main())
