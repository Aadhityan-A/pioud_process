# file: gui_pioud_process.py

import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

from pioud_process import pioud_process
from gui.gui_output_tab import OutputTab
from gui.gui_input_tab import InputTab

class PioudProcessGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("pioud_process GUI")
        self.resizable(True, True)  # Changed to allow resizing for better UX

        # Set minimum window size
        self.minsize(800, 600)

        # Create Notebook for tabs
        notebook = ttk.Notebook(self)
        notebook.grid(row=0, column=0, padx=10, pady=10, sticky='nsew')

        # Input Tab (new)
        input_tab = InputTab(notebook)
        notebook.add(input_tab, text="Input Generation")

        # Process AI Tab
        process_frame = ttk.Frame(notebook)
        notebook.add(process_frame, text="Process AI")
        self._build_widgets(process_frame)

        # Output Tab
        output_tab = OutputTab(notebook)
        notebook.add(output_tab, text="Output")

        # Configure grid weights for proper resizing
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

    def _build_widgets(self, parent):
        # === Path Selection ===
        path_frame = ttk.LabelFrame(parent, text="Working Directory")
        path_frame.grid(row=0, column=0, padx=10, pady=5, sticky="ew")

        self.path_var = tk.StringVar()
        ttk.Entry(path_frame, textvariable=self.path_var, width=50).grid(
            row=0, column=0, padx=(5, 2), pady=5
        )
        ttk.Button(path_frame, text="Browse…", command=self._browse_path).grid(
            row=0, column=1, padx=(2, 5), pady=5
        )

        # === descriptor & centroid & list_of_folder ===
        opts_frame = ttk.LabelFrame(parent, text="General Options")
        opts_frame.grid(row=1, column=0, padx=10, pady=5, sticky="ew")

        # descriptor
        ttk.Label(opts_frame, text="Descriptor:").grid(row=0, column=0, padx=5, pady=3, sticky="e")
        self.desc_var = tk.StringVar(value="soap")
        desc_combo = ttk.Combobox(opts_frame, textvariable=self.desc_var, values=["soap", "mbtr"], state="readonly", width=12)
        desc_combo.grid(row=0, column=1, padx=5, pady=3)
        desc_combo.bind("<<ComboboxSelected>>", self._on_descriptor_change)

        # centroid
        self.centroid_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(opts_frame, text="centroid", variable=self.centroid_var).grid(
            row=0, column=2, padx=5, pady=3
        )

        # list_of_folder
        self.list_of_folder_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(opts_frame, text="list_of_folder", variable=self.list_of_folder_var).grid(
            row=0, column=3, padx=5, pady=3
        )

        # === Parameters Frame ===
        self.param_frame = ttk.LabelFrame(parent, text="Parameters")
        self.param_frame.grid(row=2, column=0, padx=10, pady=5, sticky="ew")
        self.param_vars = {}
        self._build_param_fields()

        # === Calculate Button ===
        btn = ttk.Button(parent, text="Calculate", command=self._on_calculate)
        btn.grid(row=3, column=0, pady=10)

        # === Log Output ===
        self.log = scrolledtext.ScrolledText(parent, width=70, height=10, state="disabled")
        self.log.grid(row=4, column=0, padx=10, pady=(0,10))

    def _browse_path(self):
        """Open a directory chooser and set the path."""
        directory = filedialog.askdirectory()
        if directory:
            self.path_var.set(directory)

    def _on_calculate(self):
        """Gather inputs and start the calculation in a background thread."""
        # disable button to avoid double-click
        for child in self.winfo_children():
            if isinstance(child, ttk.Button) and child["text"] == "Calculate":
                child.state(["disabled"])

        # clear and enable log
        self.log.configure(state="normal")
        self.log.delete(1.0, tk.END)

        # start thread
        thread = threading.Thread(target=self._run_task, daemon=True)
        thread.start()

    def _run_task(self):
        """Invoke pioud_process with user inputs, logging progress."""
        try:
            path = self.path_var.get().strip()
            if not path:
                raise ValueError("Please specify a working directory.")

            # build parameter dict
            param_dict = {}
            descriptor = self.desc_var.get().strip()

            for k, var in self.param_vars.items():
                raw_val = var.get() if not isinstance(var, tk.BooleanVar) else var.get()
                # Convert types for SOAP parameters
                if descriptor == "soap":
                    if k in {"r_cut"}:
                        param_dict[k] = float(raw_val)
                    elif k in {"n_max", "l_max", "no_of_jobs", "n_to_select", "n_components"}:
                        param_dict[k] = int(raw_val)
                    elif k in {"periodic"}:
                        param_dict[k] = bool(raw_val)
                    else:
                        param_dict[k] = raw_val
                else:
                    # MBTR type handling
                    param_dict[k] = raw_val

            # fetch other flags
            centroid = self.centroid_var.get()
            list_of_folder = self.list_of_folder_var.get()

            # log parameters
            self._log(f"Calling pioud_process with:\n path={path}\n descriptor={descriptor}"
                      f"\n centroid={centroid}\n list_of_folder={list_of_folder}"
                      f"\n parameters={param_dict}\n")

            # call the function
            output_file = pioud_process(
                path,
                centroid=centroid,
                descriptor=descriptor,
                des_parameters=param_dict,
                list_of_folder=list_of_folder,
            )

            self._log(f"\nSuccess! Output file: {output_file}\n")
            messagebox.showinfo("Done", f"SOAP analysis completed.\nOutput: {output_file}")

        except Exception as e:
            self._log(f"\nError: {e}", err=True)
            messagebox.showerror("Error", str(e))
        finally:
            # re-enable calculate button
            for child in self.winfo_children():
                if isinstance(child, ttk.Button) and child["text"] == "Calculate":
                    child.state(["!disabled"])

    def _on_descriptor_change(self, event=None):
        """Handle descriptor selection changes"""
        # Clear and rebuild parameter fields
        for widget in self.param_frame.winfo_children():
            widget.destroy()
        self._build_param_fields()

    def _build_param_fields(self):
        """Build parameter input fields based on selected descriptor"""
        descriptor = self.desc_var.get()
        self.param_vars.clear()

        if descriptor == "soap":
            # SOAP parameters (all from soap_dict)
            entries = [
                ("r_cut", "5.0", float),
                ("n_max", "5", int),
                ("l_max", "5", int),
                ("n_components", "10", int),
                ("no_of_jobs", "4", int),
                ("n_to_select", "4", int),
            ]
            for i, (name, default, typ) in enumerate(entries):
                ttk.Label(self.param_frame, text=f"{name}:").grid(row=i, column=0, sticky="e", padx=5, pady=2)
                var = tk.StringVar(value=default)
                self.param_vars[name] = var
                ttk.Entry(self.param_frame, textvariable=var, width=10).grid(row=i, column=1, padx=5, pady=2)

            # average option (Combobox)
            ttk.Label(self.param_frame, text="average:").grid(row=0, column=2, sticky="e", padx=5, pady=2)
            self.param_vars["average"] = tk.StringVar(value="inner")
            ttk.Combobox(
                self.param_frame,
                textvariable=self.param_vars["average"],
                values=["inner", "outer"],
                state="readonly",
                width=8,
            ).grid(row=0, column=3, padx=5, pady=2)

            # periodic (Checkbutton)
            ttk.Label(self.param_frame, text="periodic:").grid(row=1, column=2, sticky="e", padx=5, pady=2)
            self.param_vars["periodic"] = tk.BooleanVar(value=True)
            ttk.Checkbutton(self.param_frame, variable=self.param_vars["periodic"]).grid(row=1, column=3, padx=5, pady=2)

            # method (Combobox)
            ttk.Label(self.param_frame, text="method:").grid(row=2, column=2, sticky="e", padx=5, pady=2)
            self.param_vars["method"] = tk.StringVar(value="CUR")
            ttk.Combobox(
                self.param_frame,
                textvariable=self.param_vars["method"],
                values=["FPS", "CUR", "CURFPS"],
                state="readonly",
                width=8,
            ).grid(row=2, column=3, padx=5, pady=2)

            # input_file_name and output_file_name 
            ttk.Label(self.param_frame, text="input_file_name:").grid(row=6, column=2, sticky="e", padx=5, pady=2)
            self.param_vars["input_file_name"] = tk.StringVar(value="ase_configs.xyz")
            ttk.Entry(self.param_frame, textvariable=self.param_vars["input_file_name"], width=20).grid(
                row=6, column=3, padx=5, pady=2
            )

            ttk.Label(self.param_frame, text="output_file_name:").grid(row=7, column=2, sticky="e", padx=5, pady=2)
            self.param_vars["output_file_name"] = tk.StringVar(value="ase_config_for_mace.xyz")
            ttk.Entry(self.param_frame, textvariable=self.param_vars["output_file_name"], width=20).grid(
                row=7, column=3, padx=5, pady=2
            )

        else:
            # MBTR parameters
            # periodic
            ttk.Label(self.param_frame, text="periodic:").grid(row=0, column=0, sticky="e", padx=5, pady=2)
            self.param_vars["periodic"] = tk.BooleanVar(value=True)
            ttk.Checkbutton(self.param_frame, variable=self.param_vars["periodic"]).grid(row=0, column=1, padx=5, pady=2)

            # n_to_select, no_of_jobs, n_components
            entries = [
                ("n_to_select", "2"),
                ("no_of_jobs", "3"),
                ("n_components", "2"),
            ]
            for offset, (name, default) in enumerate(entries, start=1):
                ttk.Label(self.param_frame, text=f"{name}:").grid(row=offset, column=0, sticky="e", padx=5, pady=2)
                var = tk.StringVar(value=default)
                self.param_vars[name] = var
                ttk.Entry(self.param_frame, textvariable=var, width=10).grid(row=offset, column=1, padx=5, pady=2)

            # method
            ttk.Label(self.param_frame, text="method:").grid(row=4, column=0, sticky="e", padx=5, pady=2)
            self.param_vars["method"] = tk.StringVar(value="CURFPS")
            ttk.Combobox(
                self.param_frame,
                textvariable=self.param_vars["method"],
                values=["FPS", "CUR", "CURFPS"],
                state="readonly",
                width=8,
            ).grid(row=4, column=1, padx=5, pady=2)

            # MBTR parameter keys
            mbtr_entries = [
                ("geometry_function", "inverse_distance"),
                ("grid_min", "0"),
                ("grid_max", "1"),
                ("grid_n", "100"),
                ("grid_sigma", "0.1"),
                ("weighting_function", "exp"),
                ("weighting_scale", "0.5"),
                ("weighting_threshold", "1e-3"),
                ("normalization", "none"),
            ]
            for i, (name, default) in enumerate(mbtr_entries, start=5):
                ttk.Label(self.param_frame, text=f"{name}:").grid(row=i, column=0, sticky="e", padx=5, pady=2)
                var = tk.StringVar(value=default)
                self.param_vars[name] = var
                ttk.Entry(self.param_frame, textvariable=var, width=15).grid(row=i, column=1, padx=5, pady=2)

            # input_file_name and output_file_name
            ttk.Label(self.param_frame, text="input_file_name:").grid(row=15, column=2, sticky="e", padx=5, pady=2)
            self.param_vars["input_file_name"] = tk.StringVar(value="ase_configs.xyz")
            ttk.Entry(self.param_frame, textvariable=self.param_vars["input_file_name"], width=20).grid(
                row=15, column=3, padx=5, pady=2
            )

            ttk.Label(self.param_frame, text="output_file_name:").grid(row=16, column=2, sticky="e", padx=5, pady=2)
            self.param_vars["output_file_name"] = tk.StringVar(value="ase_config_for_mace.xyz")
            ttk.Entry(self.param_frame, textvariable=self.param_vars["output_file_name"], width=20).grid(
                row=16, column=3, padx=5, pady=2
            )

    def _log(self, msg: str, err: bool = False):
        """Thread-safe appending to the log widget."""
        # schedule on main thread
        def append():
            self.log.configure(state="normal")
            tag = "ERROR" if err else None
            self.log.insert(tk.END, msg, tag)
            self.log.see(tk.END)
            self.log.configure(state="disabled")

        self.after(0, append)


if __name__ == "__main__":
    app = PioudProcessGUI()
    app.mainloop()