"""Plot from bundled data by default; pass --source-data for workbook analysis and plots."""
import sys
from oart_analysis.cli import main

if __name__ == '__main__':
    main(['all',*sys.argv[1:]])
