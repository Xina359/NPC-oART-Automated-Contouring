"""Reproduce the data-based panels of Figure S6 from bundled data or an explicitly selected workbook."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from oart_analysis.cli import main

if __name__ == '__main__':
    main(['plot', '--figure', 'S6', *sys.argv[1:]])
