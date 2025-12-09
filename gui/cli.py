#!/usr/bin/env python3
import sys
import os
import importlib.util
import traceback

def main():
    try:
        # Add the current directory to the path to ensure imports work correctly
        current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if current_dir not in sys.path:
            sys.path.insert(0, current_dir)
        
        # Import the GUI module
        from gui.gui_pioud_process import PioudProcessGUI
        
        # Launch the GUI
        app = PioudProcessGUI()
        app.mainloop()
        return 0
    except Exception as e:
        print(f"Error launching GUI: {str(e)}")
        print("Traceback:")
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())