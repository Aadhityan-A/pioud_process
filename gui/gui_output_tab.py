import threading
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from pioud_process.main import pioud_process
from ase.io.extxyz import write_extxyz
from ase.io import read
from ase.visualize import view
import sys
import subprocess
import os

class OutputTab(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)

        # Create progress bar
        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(self, variable=self.progress_var, maximum=100)
        self.progress_bar.pack(fill='x', padx=5, pady=5)

        # Create status label
        self.status_label = ttk.Label(self, text="Ready")
        self.status_label.pack(pady=5)

        # Create input folder selection frame (changed from file to folder)
        input_frame = ttk.Frame(self)
        input_frame.pack(fill='x', padx=5, pady=5)

        ttk.Label(input_frame, text="Input Folder:").pack(side='left')
        self.input_entry = ttk.Entry(input_frame)
        self.input_entry.pack(side='left', fill='x', expand=True, padx=5)
        # Set default to current working directory
        self.input_entry.insert(0, os.getcwd())

        browse_btn = ttk.Button(input_frame, text="Browse", command=self.browse_input_folder)
        browse_btn.pack(side='right')

        # Create centroid selection frame with a single checkbox
        centroid_frame = ttk.LabelFrame(self, text="Centroid Option")
        centroid_frame.pack(fill='x', padx=5, pady=5)

        # Create BooleanVar to track centroid selection
        self.centroid_var = tk.BooleanVar(value=True)  # Default to True

        # Create a single checkbox for centroid selection
        ttk.Checkbutton(
            centroid_frame,
            text="Enable Centroid",
            variable=self.centroid_var,
            onvalue=True,
            offvalue=False
        ).pack(anchor='w', padx=10, pady=5)

        # Create generate button
        self.generate_btn = ttk.Button(self, text="Generate Extended XYZ File", command=self.on_generate)
        self.generate_btn.pack(pady=5)

        # ASE file selection & view controls
        ase_frame = ttk.Frame(self)
        ase_frame.pack(fill='x', padx=5, pady=(10,5))
        ttk.Label(ase_frame, text="ASE File:").pack(side='left')
        self.ase_entry = ttk.Entry(ase_frame)
        self.ase_entry.pack(side='left', fill='x', expand=True, padx=5)
        ttk.Button(ase_frame, text="Browse", command=self.browse_ase_file).pack(side='left', padx=5)
        
        # View structures button
        self.view_btn = ttk.Button(self, text="View Structures", command=self.on_view)
        self.view_btn.pack(pady=5)

        # Folder Browser section with navigation
        browser_frame = ttk.LabelFrame(self, text="Folder Browser")
        browser_frame.pack(fill='both', expand=True, padx=5, pady=5)

        # Navigation frame
        nav_frame = ttk.Frame(browser_frame)
        nav_frame.pack(fill='x', pady=(5, 0))
        
        ttk.Button(nav_frame, text="🏠 Home", command=self.go_home).pack(side='left', padx=2)
        ttk.Button(nav_frame, text="⬆️ Up", command=self.go_up).pack(side='left', padx=2)
        ttk.Button(nav_frame, text="🔄 Refresh", command=self.refresh_folder).pack(side='left', padx=2)

        folder_select_frame = ttk.Frame(browser_frame)
        folder_select_frame.pack(fill='x', pady=(5, 0))
        ttk.Label(folder_select_frame, text="Folder:").pack(side='left')
        self.folder_entry = ttk.Entry(folder_select_frame)
        self.folder_entry.pack(side='left', fill='x', expand=True, padx=5)
        # Set default to current working directory
        self.folder_entry.insert(0, os.getcwd())
        ttk.Button(folder_select_frame, text="Browse", command=self.browse_folder).pack(side='left')

        # Create frame for listbox and scrollbar to fix positioning
        listbox_frame = ttk.Frame(browser_frame)
        listbox_frame.pack(fill='both', expand=True, pady=5)

        # Pack scrollbar first to the right side
        file_scrollbar = ttk.Scrollbar(listbox_frame, orient='vertical')
        file_scrollbar.pack(side='right', fill='y')

        # Then pack listbox with scrollbar command already set
        self.file_listbox = tk.Listbox(listbox_frame, yscrollcommand=file_scrollbar.set)
        self.file_listbox.pack(side='left', fill='both', expand=True)

        # Configure scrollbar to control listbox
        file_scrollbar.config(command=self.file_listbox.yview)
        self.file_listbox.bind('<Double-Button-1>', self.on_file_double_click)

        # Initialize with current working directory
        self.current_folder = os.getcwd()
        self.populate_file_list()

    def browse_input_folder(self):
        """Open a folder dialog to select the input folder."""
        folder_path = filedialog.askdirectory(
            title="Select Input Folder",
            initialdir=os.getcwd()  # Start from current directory
        )
        if folder_path:
            self.input_entry.delete(0, tk.END)
            self.input_entry.insert(0, folder_path)

    def go_home(self):
        """Navigate to home directory."""
        home_dir = os.path.expanduser("~")
        self.current_folder = home_dir
        self.folder_entry.delete(0, tk.END)
        self.folder_entry.insert(0, home_dir)
        self.populate_file_list()

    def go_up(self):
        """Navigate to parent directory."""
        if self.current_folder:
            parent_dir = os.path.dirname(self.current_folder)
            if parent_dir != self.current_folder:  # Prevent infinite loop at root
                self.current_folder = parent_dir
                self.folder_entry.delete(0, tk.END)
                self.folder_entry.insert(0, parent_dir)
                self.populate_file_list()

    def refresh_folder(self):
        """Refresh the current folder contents."""
        if self.current_folder:
            self.populate_file_list()

    def on_generate(self):
        """Handle the 'Generate Extended XYZ File' button click."""
        input_folder = self.input_entry.get().strip()
        if not input_folder:
            messagebox.showwarning("No Input Folder",
                                 "Please select an input folder before generating.")
            return
        
        if not os.path.isdir(input_folder):
            messagebox.showwarning("Invalid Folder",
                                 "The specified path is not a valid directory.")
            return
        
        # Reset previous output and status
        self.clear_output()
        # Launch the PIOUD-SOAP calculation in background
        self.run_calculation(input_folder)

    def clear_output(self):
        """Clear the output text widget and reset progress"""
        # Note: output_text widget doesn't exist in current code, removing reference
        self.progress_var.set(0)
        self.status_label.config(text="Ready")

    def update_progress(self, value):
        """Update the progress bar value"""
        self.progress_var.set(value)

    def update_status(self, text):
        """Update the status label text"""
        self.status_label.config(text=text)

    def run_calculation(self, input_folder, output_file=None, parameters=None):
        """Run the PIOUD-SOAP calculation on a directory and emit an .extxyz."""
        def _task():
            try:
                # Signal start
                self.update_status("Running PIOUD-SOAP…")
                # Get centroid value from radio button selection
                centroid_value = self.centroid_var.get()
                # Call the real entry-point with folder path and selected centroid option
                structures_file_name = pioud_process(input_folder, centroid=centroid_value)

                self.update_status("Completed")
                messagebox.showinfo(
                    "Done", f"Extended-xyz generated:\n{structures_file_name}\nCentroid: {centroid_value}"
                )
            except Exception as e:
                self.update_status("Failed")
                messagebox.showerror("Error", f"Calculation failed:\n{e}")

        threading.Thread(target=_task, daemon=True).start()

    def browse_ase_file(self):
        """Open a file dialog to select the ASE-readable file."""
        file_path = filedialog.askopenfilename(
            title="Select ASE File",
            initialdir=os.getcwd(),  # Start from current directory
            filetypes=[
                ("All Files", "*.*"),
                ("XYZ files", "*.xyz"),
                ("Extended XYZ files", "*.extxyz")
            ]
        )
        if file_path:
            self.ase_entry.delete(0, tk.END)
            self.ase_entry.insert(0, file_path)

    def on_view(self):
        """Read the selected ASE file and launch the viewer."""
        ase_file = self.ase_entry.get().strip()
        if not ase_file:
            messagebox.showwarning("No ASE File", "Please select an ASE file to view structures.")
            return
        try:
            configs = read(ase_file, index=":")
            view(configs)
        except Exception as e:
            messagebox.showerror("View Error", f"Failed to view structures:\n{e}")

    def browse_folder(self):
        folder = filedialog.askdirectory(
            title="Select Folder",
            initialdir=self.current_folder or os.getcwd()  # Use current folder or cwd as initial
        )
        if folder:
            self.folder_entry.delete(0, tk.END)
            self.folder_entry.insert(0, folder)
            self.current_folder = folder
            self.populate_file_list()

    def populate_file_list(self):
        """Populate the file listbox with current folder contents."""
        self.file_listbox.delete(0, tk.END)
        try:
            if not self.current_folder or not os.path.exists(self.current_folder):
                return
            
            entries = os.listdir(self.current_folder)
            # Separate directories and files, show directories first
            dirs = [entry for entry in entries if os.path.isdir(os.path.join(self.current_folder, entry))]
            files = [entry for entry in entries if os.path.isfile(os.path.join(self.current_folder, entry))]
            
            # Add directories first with folder icon
            for dir_name in sorted(dirs):
                self.file_listbox.insert(tk.END, f"📁 {dir_name}")
            
            # Add files with file icon
            for file_name in sorted(files):
                self.file_listbox.insert(tk.END, f"📄 {file_name}")
                
        except Exception as e:
            messagebox.showerror("Error", f"Cannot list directory:\n{e}")

    def on_file_double_click(self, event):
        """Handle double-click on file/folder in the listbox."""
        sel = self.file_listbox.curselection()
        if not sel:
            return
        
        display_name = self.file_listbox.get(sel[0])
        # Remove the icon prefix to get actual name
        name = display_name[2:] if display_name.startswith(('📁 ', '📄 ')) else display_name
        path = os.path.join(self.current_folder, name)
        
        if os.path.isdir(path):
            # Navigate into directory
            self.current_folder = path
            self.folder_entry.delete(0, tk.END)
            self.folder_entry.insert(0, path)
            self.populate_file_list()
        else:
            # Open file with system default application
            self.open_path(path)

    def open_path(self, path):
        """Open file or folder with system default application."""
        try:
            if sys.platform.startswith('darwin'):
                subprocess.call(['open', path])
            elif os.name == 'nt':
                os.startfile(path)
            else:
                subprocess.call(['xdg-open', path])
        except Exception as e:
            messagebox.showerror("Open File Error", f"Cannot open \"{path}\":\n{e}")