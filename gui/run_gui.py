#!/usr/bin/env python3

import os
import sys

# Add the parent directory to the path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Import the GUI modules
from gui.gui_pioud_process import PioudProcessGUI

def main():
    app = PioudProcessGUI()
    app.mainloop()
    return 0

if __name__ == "__main__":
    sys.exit(main())