
import os
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
from ase.io import read
from ase.data import atomic_masses, atomic_numbers


class InputTab(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        
        # Initialize variables
        self.cif_file_path = None
        self.current_structure = None
        
        # Create notebook for sub-tabs
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True, padx=5, pady=5)
        
        # Create tabs
        self._create_file_input_tab()
        self._create_pimd_tab()
        self._create_qe_tab()
        self._create_hpcc_tab()
        
    def _create_file_input_tab(self):
        """Create the file input and basic parameters tab"""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="File Input")
        
        # Title
        title_label = ttk.Label(frame, text="Upload CIF File and Set Parameters", 
                               font=('TkDefaultFont', 12, 'bold'))
        title_label.grid(row=0, column=0, columnspan=2, pady=10)
        
        # CIF file selection
        cif_frame = ttk.LabelFrame(frame, text="CIF File Selection")
        cif_frame.grid(row=1, column=0, columnspan=2, sticky='ew', padx=10, pady=5)
        
        self.cif_path_var = tk.StringVar()
        ttk.Entry(cif_frame, textvariable=self.cif_path_var, width=50, state='readonly').grid(
            row=0, column=0, padx=5, pady=5)
        ttk.Button(cif_frame, text="Browse CIF", command=self._browse_cif).grid(
            row=0, column=1, padx=5, pady=5)
        
        # Prefix setting
        prefix_frame = ttk.LabelFrame(frame, text="Calculation Settings")
        prefix_frame.grid(row=2, column=0, columnspan=2, sticky='ew', padx=10, pady=5)
        
        ttk.Label(prefix_frame, text="Prefix:").grid(row=0, column=0, sticky='e', padx=5, pady=5)
        self.prefix_var = tk.StringVar(value='diamond')
        ttk.Entry(prefix_frame, textvariable=self.prefix_var, width=20).grid(
            row=0, column=1, padx=5, pady=5)
        
        # Generate buttons
        button_frame = ttk.Frame(frame)
        button_frame.grid(row=3, column=0, columnspan=2, pady=20)
        
        self.generate_input_btn = ttk.Button(button_frame, text="Generate Input", 
                                           command=self._generate_input, style='Accent.TButton')
        self.generate_input_btn.grid(row=0, column=0, padx=10)
        
        # Instructions
        instructions = ttk.Label(frame, text="""Instructions:
1. Upload a CIF file using the Browse button
2. Set a prefix for your calculation
3. Adjust parameters in other tabs as needed
4. Click 'Generate Input' to create the input file""", 
                               justify='left')
        instructions.grid(row=4, column=0, columnspan=2, padx=10, pady=10)
        
        # Configure grid weights
        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        
    def _create_pimd_tab(self):
        """Create the PIMD parameters tab"""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="PIMD Parameters")
        
        # Title
        title_label = ttk.Label(frame, text="Path Integral Molecular Dynamics Parameters", 
                               font=('TkDefaultFont', 12, 'bold'))
        title_label.grid(row=0, column=0, columnspan=4, pady=10)
        
        # PIMD parameters
        self.pimd_vars = {}
        
        # Row 1: nbeadMD, nblocks
        ttk.Label(frame, text="nbeadMD:").grid(row=1, column=0, sticky='e', padx=5, pady=3)
        self.pimd_vars['nbeadMD'] = tk.StringVar(value='2.0')
        ttk.Entry(frame, textvariable=self.pimd_vars['nbeadMD'], width=10).grid(
            row=1, column=1, padx=5, pady=3)
        
        ttk.Label(frame, text="nblocks:").grid(row=1, column=2, sticky='e', padx=5, pady=3)
        self.pimd_vars['nblocks'] = tk.StringVar(value='2.0')
        ttk.Entry(frame, textvariable=self.pimd_vars['nblocks'], width=10).grid(
            row=1, column=3, padx=5, pady=3)
        
        # Row 2: nstep_block, nunitcells
        ttk.Label(frame, text="nstep_block:").grid(row=2, column=0, sticky='e', padx=5, pady=3)
        self.pimd_vars['nstep_block'] = tk.StringVar(value='2.0')
        ttk.Entry(frame, textvariable=self.pimd_vars['nstep_block'], width=10).grid(
            row=2, column=1, padx=5, pady=3)
        
        ttk.Label(frame, text="nunitcells:").grid(row=2, column=2, sticky='e', padx=5, pady=3)
        self.pimd_vars['nunitcells'] = tk.StringVar(value='1.0')
        ttk.Entry(frame, textvariable=self.pimd_vars['nunitcells'], width=10).grid(
            row=2, column=3, padx=5, pady=3)
        
        # Row 3: iprint, run
        ttk.Label(frame, text="iprint:").grid(row=3, column=0, sticky='e', padx=5, pady=3)
        self.pimd_vars['iprint'] = tk.StringVar(value='1.0')
        ttk.Entry(frame, textvariable=self.pimd_vars['iprint'], width=10).grid(
            row=3, column=1, padx=5, pady=3)
        
        ttk.Label(frame, text="run:").grid(row=3, column=2, sticky='e', padx=5, pady=3)
        self.pimd_vars['run'] = tk.StringVar(value='pioud')
        ttk.Entry(frame, textvariable=self.pimd_vars['run'], width=10).grid(
            row=3, column=3, padx=5, pady=3)
        
        # Row 4: delt, tempMD
        ttk.Label(frame, text="delt:").grid(row=4, column=0, sticky='e', padx=5, pady=3)
        self.pimd_vars['delt'] = tk.StringVar(value='2.0e-3')
        ttk.Entry(frame, textvariable=self.pimd_vars['delt'], width=10).grid(
            row=4, column=1, padx=5, pady=3)
        
        ttk.Label(frame, text="tempMD (K):").grid(row=4, column=2, sticky='e', padx=5, pady=3)
        self.pimd_vars['tempMD'] = tk.StringVar(value='300.0')
        ttk.Entry(frame, textvariable=self.pimd_vars['tempMD'], width=10).grid(
            row=4, column=3, padx=5, pady=3)
        
        # Row 5: gammaMD, delta_force
        ttk.Label(frame, text="gammaMD:").grid(row=5, column=0, sticky='e', padx=5, pady=3)
        self.pimd_vars['gammaMD'] = tk.StringVar(value='0.21')
        ttk.Entry(frame, textvariable=self.pimd_vars['gammaMD'], width=10).grid(
            row=5, column=1, padx=5, pady=3)
        
        ttk.Label(frame, text="delta_force:").grid(row=5, column=2, sticky='e', padx=5, pady=3)
        self.pimd_vars['delta_force'] = tk.StringVar(value='1.0e-4')
        ttk.Entry(frame, textvariable=self.pimd_vars['delta_force'], width=10).grid(
            row=5, column=3, padx=5, pady=3)
        
        # Row 6: delta_harm, restart_pimd
        ttk.Label(frame, text="delta_harm:").grid(row=6, column=0, sticky='e', padx=5, pady=3)
        self.pimd_vars['delta_harm'] = tk.StringVar(value='5.0e-3')
        ttk.Entry(frame, textvariable=self.pimd_vars['delta_harm'], width=10).grid(
            row=6, column=1, padx=5, pady=3)
        
        self.pimd_vars['restart_pimd'] = tk.BooleanVar(value=False)
        ttk.Checkbutton(frame, text="restart_pimd", variable=self.pimd_vars['restart_pimd']).grid(
            row=6, column=2, columnspan=2, sticky='w', padx=5, pady=3)
        
    def _create_qe_tab(self):
        """Create the Quantum ESPRESSO parameters tab"""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="QE Parameters")
        
        # Create scrollable frame
        canvas = tk.Canvas(frame)
        scrollbar = ttk.Scrollbar(frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Title
        title_label = ttk.Label(scrollable_frame, text="Compute SCF Parameters", 
                               font=('TkDefaultFont', 12, 'bold'))
        title_label.grid(row=0, column=0, columnspan=4, pady=10)
        
        self.qe_vars = {}
        
        # Energy Cutoffs & Bands section
        section_label = ttk.Label(scrollable_frame, text="Energy Cutoffs & Bands", 
                                 font=('TkDefaultFont', 10, 'bold'))
        section_label.grid(row=1, column=0, columnspan=4, sticky='w', pady=(10, 5))
        
        # ecutwfc, ecutrho, nbnd
        ttk.Label(scrollable_frame, text="ecutwfc (Ry):").grid(row=2, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['ecutwfc'] = tk.StringVar(value='25.0')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['ecutwfc'], width=10).grid(
            row=2, column=1, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="ecutrho (Ry):").grid(row=2, column=2, sticky='e', padx=5, pady=3)
        self.qe_vars['ecutrho'] = tk.StringVar(value='120.0')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['ecutrho'], width=10).grid(
            row=2, column=3, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="nbnd:").grid(row=3, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['nbnd'] = tk.StringVar(value='20')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['nbnd'], width=10).grid(
            row=3, column=1, padx=5, pady=3)
        
        # Electronic Structure section
        section_label = ttk.Label(scrollable_frame, text="Electronic Structure", 
                                 font=('TkDefaultFont', 10, 'bold'))
        section_label.grid(row=4, column=0, columnspan=4, sticky='w', pady=(10, 5))
        
        # nspin, occupations, smearing
        ttk.Label(scrollable_frame, text="nspin:").grid(row=5, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['nspin'] = tk.StringVar(value='1')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['nspin'], width=10).grid(
            row=5, column=1, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="occupations:").grid(row=5, column=2, sticky='e', padx=5, pady=3)
        self.qe_vars['occupations'] = tk.StringVar(value='smearing')
        occupations_combo = ttk.Combobox(scrollable_frame, textvariable=self.qe_vars['occupations'], 
                                       values=['smearing', 'fixed', 'tetrahedra'], state='readonly', width=8)
        occupations_combo.grid(row=5, column=3, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="smearing:").grid(row=6, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['smearing'] = tk.StringVar(value='fd')
        smearing_combo = ttk.Combobox(scrollable_frame, textvariable=self.qe_vars['smearing'], 
                                    values=['fd', 'mv', 'mp', 'gauss'], state='readonly', width=8)
        smearing_combo.grid(row=6, column=1, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="degauss (Ry):").grid(row=6, column=2, sticky='e', padx=5, pady=3)
        self.qe_vars['degauss'] = tk.StringVar(value='1.0e-3')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['degauss'], width=10).grid(
            row=6, column=3, padx=5, pady=3)
        
        # K-points section
        section_label = ttk.Label(scrollable_frame, text="K-points (automatic)", 
                                 font=('TkDefaultFont', 10, 'bold'))
        section_label.grid(row=7, column=0, columnspan=4, sticky='w', pady=(10, 5))
        
        # K-points grid
        ttk.Label(scrollable_frame, text="K-points X:").grid(row=8, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['kx'] = tk.StringVar(value='4')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['kx'], width=10).grid(
            row=8, column=1, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="Y:").grid(row=8, column=2, sticky='e', padx=5, pady=3)
        self.qe_vars['ky'] = tk.StringVar(value='4')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['ky'], width=10).grid(
            row=8, column=3, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="Z:").grid(row=9, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['kz'] = tk.StringVar(value='4')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['kz'], width=10).grid(
            row=9, column=1, padx=5, pady=3)
        
        # K-point shifts
        ttk.Label(scrollable_frame, text="Shift X:").grid(row=10, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['k_shift_x'] = tk.StringVar(value='0.0')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['k_shift_x'], width=10).grid(
            row=10, column=1, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="Y:").grid(row=10, column=2, sticky='e', padx=5, pady=3)
        self.qe_vars['k_shift_y'] = tk.StringVar(value='0.0')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['k_shift_y'], width=10).grid(
            row=10, column=3, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="Z:").grid(row=11, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['k_shift_z'] = tk.StringVar(value='0.0')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['k_shift_z'], width=10).grid(
            row=11, column=1, padx=5, pady=3)
        
        # Convergence & Mixing section
        section_label = ttk.Label(scrollable_frame, text="Convergence & Mixing", 
                                 font=('TkDefaultFont', 10, 'bold'))
        section_label.grid(row=12, column=0, columnspan=4, sticky='w', pady=(10, 5))
        
        ttk.Label(scrollable_frame, text="conv_thr:").grid(row=13, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['conv_thr'] = tk.StringVar(value='1.0e-8')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['conv_thr'], width=10).grid(
            row=13, column=1, padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="mixing_beta:").grid(row=13, column=2, sticky='e', padx=5, pady=3)
        self.qe_vars['mixing_beta'] = tk.StringVar(value='0.7')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['mixing_beta'], width=10).grid(
            row=13, column=3, padx=5, pady=3)
        
        # Symmetry & Other Options section
        section_label = ttk.Label(scrollable_frame, text="Symmetry & Other Options", 
                                 font=('TkDefaultFont', 10, 'bold'))
        section_label.grid(row=14, column=0, columnspan=4, sticky='w', pady=(10, 5))
        
        self.qe_vars['nosym'] = tk.BooleanVar(value=True)
        ttk.Checkbutton(scrollable_frame, text="nosym", variable=self.qe_vars['nosym']).grid(
            row=15, column=0, sticky='w', padx=5, pady=3)
        
        self.qe_vars['nosym_evc'] = tk.BooleanVar(value=True)
        ttk.Checkbutton(scrollable_frame, text="nosym_evc", variable=self.qe_vars['nosym_evc']).grid(
            row=15, column=1, sticky='w', padx=5, pady=3)
        
        self.qe_vars['tstress'] = tk.BooleanVar(value=True)
        ttk.Checkbutton(scrollable_frame, text="tstress", variable=self.qe_vars['tstress']).grid(
            row=15, column=2, sticky='w', padx=5, pady=3)
        
        ttk.Label(scrollable_frame, text="pseudo_dir:").grid(row=16, column=0, sticky='e', padx=5, pady=3)
        self.qe_vars['pseudo_dir'] = tk.StringVar(value='./')
        ttk.Entry(scrollable_frame, textvariable=self.qe_vars['pseudo_dir'], width=20).grid(
            row=16, column=1, columnspan=2, sticky='w', padx=5, pady=3)
        
        # Pack canvas and scrollbar
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
    def _create_hpcc_tab(self):
        """Create the HPCC submission script tab"""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="HPCC Submit")
        
        # Title
        title_label = ttk.Label(frame, text="HPC Job Submission Script", 
                               font=('TkDefaultFont', 12, 'bold'))
        title_label.grid(row=0, column=0, columnspan=4, pady=10)
        
        self.hpcc_vars = {}
        
        # Scheduler type
        ttk.Label(frame, text="Scheduler:").grid(row=1, column=0, sticky='e', padx=5, pady=3)
        self.hpcc_vars['scheduler'] = tk.StringVar(value='PBS')
        scheduler_combo = ttk.Combobox(frame, textvariable=self.hpcc_vars['scheduler'], 
                                     values=['PBS', 'SLURM'], state='readonly', width=10)
        scheduler_combo.grid(row=1, column=1, padx=5, pady=3)
        
        # Job parameters
        ttk.Label(frame, text="Job name:").grid(row=2, column=0, sticky='e', padx=5, pady=3)
        self.hpcc_vars['job_name'] = tk.StringVar(value='pioud_job')
        ttk.Entry(frame, textvariable=self.hpcc_vars['job_name'], width=15).grid(
            row=2, column=1, padx=5, pady=3)
        
        ttk.Label(frame, text="Wall time:").grid(row=2, column=2, sticky='e', padx=5, pady=3)
        self.hpcc_vars['walltime'] = tk.StringVar(value='24:00:00')
        ttk.Entry(frame, textvariable=self.hpcc_vars['walltime'], width=15).grid(
            row=2, column=3, padx=5, pady=3)
        
        # Resource parameters
        ttk.Label(frame, text="Nodes:").grid(row=3, column=0, sticky='e', padx=5, pady=3)
        self.hpcc_vars['nodes'] = tk.StringVar(value='1')
        ttk.Entry(frame, textvariable=self.hpcc_vars['nodes'], width=15).grid(
            row=3, column=1, padx=5, pady=3)
        
        ttk.Label(frame, text="Cores/node:").grid(row=3, column=2, sticky='e', padx=5, pady=3)
        self.hpcc_vars['cores_per_node'] = tk.StringVar(value='36')
        ttk.Entry(frame, textvariable=self.hpcc_vars['cores_per_node'], width=15).grid(
            row=3, column=3, padx=5, pady=3)
        
        ttk.Label(frame, text="Queue/Partition:").grid(row=4, column=0, sticky='e', padx=5, pady=3)
        self.hpcc_vars['queue'] = tk.StringVar(value='standard')
        ttk.Entry(frame, textvariable=self.hpcc_vars['queue'], width=15).grid(
            row=4, column=1, padx=5, pady=3)
        
        ttk.Label(frame, text="Memory:").grid(row=4, column=2, sticky='e', padx=5, pady=3)
        self.hpcc_vars['memory'] = tk.StringVar(value='120GB')
        ttk.Entry(frame, textvariable=self.hpcc_vars['memory'], width=15).grid(
            row=4, column=3, padx=5, pady=3)
        
        ttk.Label(frame, text="MPI command:").grid(row=5, column=0, sticky='e', padx=5, pady=3)
        self.hpcc_vars['mpi_command'] = tk.StringVar(value='mpirun -np')
        ttk.Entry(frame, textvariable=self.hpcc_vars['mpi_command'], width=15).grid(
            row=5, column=1, padx=5, pady=3)
        
        # Email notifications
        self.hpcc_vars['email_notify'] = tk.BooleanVar(value=False)
        email_check = ttk.Checkbutton(frame, text="Email notifications", 
                                    variable=self.hpcc_vars['email_notify'],
                                    command=self._toggle_email_field)
        email_check.grid(row=5, column=2, sticky='w', padx=5, pady=3)
        
        ttk.Label(frame, text="Email:").grid(row=6, column=0, sticky='e', padx=5, pady=3)
        self.hpcc_vars['email_address'] = tk.StringVar()
        self.email_entry = ttk.Entry(frame, textvariable=self.hpcc_vars['email_address'], 
                                   width=25, state='disabled')
        self.email_entry.grid(row=6, column=1, columnspan=2, sticky='w', padx=5, pady=3)
        
        # Modules
        ttk.Label(frame, text="Modules:").grid(row=7, column=0, sticky='nw', padx=5, pady=3)
        self.modules_text = tk.Text(frame, width=40, height=3)
        self.modules_text.insert('1.0', 'module load quantum-espresso\nmodule load i-pi')
        self.modules_text.grid(row=7, column=1, columnspan=3, sticky='ew', padx=5, pady=3)
        
        # Additional commands
        ttk.Label(frame, text="Setup:").grid(row=8, column=0, sticky='nw', padx=5, pady=3)
        self.setup_text = tk.Text(frame, width=40, height=3)
        self.setup_text.insert('1.0', '# Additional setup commands\n')
        self.setup_text.grid(row=8, column=1, columnspan=3, sticky='ew', padx=5, pady=3)
        
        # Generate script button
        self.generate_script_btn = ttk.Button(frame, text="Generate Submit Script", 
                                            command=self._generate_script, style='Accent.TButton')
        self.generate_script_btn.grid(row=9, column=0, columnspan=4, pady=20)
        
        # Configure grid weights
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(3, weight=1)
        
    def _browse_cif(self):
        """Open a file dialog to select a CIF file"""
        file_path = filedialog.askopenfilename(
            title="Select CIF File",
            filetypes=[("CIF files", "*.cif"), ("All files", "*.*")]
        )
        if file_path:
            self.cif_file_path = file_path
            self.cif_path_var.set(file_path)
            try:
                # Validate the file by attempting to read it
                self.current_structure = read(file_path, format='cif')
                messagebox.showinfo("Success", "CIF file loaded successfully!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to read CIF file:\n{e}")
                self.cif_file_path = None
                self.current_structure = None
                self.cif_path_var.set("")
                
    def _toggle_email_field(self):
        """Enable/disable email field based on checkbox"""
        if self.hpcc_vars['email_notify'].get():
            self.email_entry.config(state='normal')
        else:
            self.email_entry.config(state='disabled')
            
    def _generate_input(self):
        """Generate PIOUD input file"""
        if not self.current_structure:
            messagebox.showwarning("No Structure", "Please load a CIF file first.")
            return
            
        try:
            prefix = self.prefix_var.get().strip()
            if not prefix:
                messagebox.showwarning("No Prefix", "Please specify a prefix.")
                return
                
            # Create directory if it doesn't exist
            if not os.path.exists(prefix):
                os.makedirs(prefix)
                
            # Extract structure information
            cell = self.current_structure.get_cell()
            symbols = self.current_structure.get_chemical_symbols()
            frac_positions = self.current_structure.get_scaled_positions()
            
            # Collect unique species
            unique_species = []
            for s in symbols:
                if s not in unique_species:
                    unique_species.append(s)
                    
            # Generate atomic species lines
            atomic_species_lines = []
            for s in unique_species:
                Z = atomic_numbers[s]
                mass = atomic_masses[Z]
                pp = f"{s}.UPF"
                atomic_species_lines.append(f"{s:<2} {mass:8.6f} {pp}")
                
            # Generate atomic positions
            atomic_positions_lines = []
            for s, pos in zip(symbols, frac_positions):
                atomic_positions_lines.append(
                    f"{s:<2} {pos[0]:8.6f} {pos[1]:8.6f} {pos[2]:8.6f}"
                )
                
            # Build the input file
            nat = len(symbols)
            ntyp = len(unique_species)
            
            lines = [
                "BEGIN",
                "BEGIN_PIMD_INPUT",
                "&dynamics",
                f"  nbeadMD={self.pimd_vars['nbeadMD'].get()}d0",
                f"  nblocks={self.pimd_vars['nblocks'].get()}d0",
                f"  nstep_block={self.pimd_vars['nstep_block'].get()}d0",
                f"  restart_pimd={'.true.' if self.pimd_vars['restart_pimd'].get() else '.false.'}",
                f"  nunitcells={self.pimd_vars['nunitcells'].get()}d0",
                f"  iprint={self.pimd_vars['iprint'].get()}d0",
                f"  run='{self.pimd_vars['run'].get()}'",
                f"  delt={float(self.pimd_vars['delt'].get()):.3e}d0",
                f"  tempMD={float(self.pimd_vars['tempMD'].get()):.1f}d0",
                f"  gammaMD={float(self.pimd_vars['gammaMD'].get()):.2f}d0",
                f"  delta_force={float(self.pimd_vars['delta_force'].get()):.3e}d0",
                f"  delta_harm={float(self.pimd_vars['delta_harm'].get()):.3e}d0",
                "/",
                "END_PIMD_INPUT",
                "BEGIN_ENGINE_INPUT",
                "&CONTROL",
                f'  prefix = "{prefix}"',
                f'  pseudo_dir = "{self.qe_vars["pseudo_dir"].get()}"',
                f"  tstress = {'.true.' if self.qe_vars['tstress'].get() else '.false.'}",
                "/",
                "&SYSTEM",
                "  ibrav= 0,",
                f"  nat= {nat},",
                f"  ntyp= {ntyp},",
                f"  ecutwfc= {float(self.qe_vars['ecutwfc'].get())},",
                f"  nosym= {'.true.' if self.qe_vars['nosym'].get() else '.false.'}",
                f"  nosym_evc= {'.true.' if self.qe_vars['nosym_evc'].get() else '.false.'}",
                f"  nspin= {int(self.qe_vars['nspin'].get())},",
                f"  ecutrho= {float(self.qe_vars['ecutrho'].get())},",
                f"  nbnd= {int(self.qe_vars['nbnd'].get())},",
                f"  occupations='{self.qe_vars['occupations'].get()}'",
                f"  smearing='{self.qe_vars['smearing'].get()}'",
                f"  degauss= {float(self.qe_vars['degauss'].get()):.3e},",
                "/",
                "&ELECTRONS",
                f"  conv_thr = {float(self.qe_vars['conv_thr'].get()):.3e},",
                f"  mixing_beta = {float(self.qe_vars['mixing_beta'].get())},",
                "/",
                "&IONS",
                "/",
                "ATOMIC_SPECIES",
            ]
            lines += atomic_species_lines
            lines += [
                "BEGIN_POSITIONS",
                "ATOMIC_POSITIONS crystal",
            ]
            lines += atomic_positions_lines
            lines += [
                "END_POSITIONS",
                "CELL_PARAMETERS angstrom",
            ]
            for vec in cell:
                lines.append(f" {vec[0]:.14f} {vec[1]:.14f} {vec[2]:.14f}")
            lines += [
                "K_POINTS automatic",
                f"{int(self.qe_vars['kx'].get())} {int(self.qe_vars['ky'].get())} {int(self.qe_vars['kz'].get())} "
                f"{float(self.qe_vars['k_shift_x'].get())} {float(self.qe_vars['k_shift_y'].get())} {float(self.qe_vars['k_shift_z'].get())}",
                "END_ENGINE_INPUT",
                "END",
            ]
            
            contents = "\n".join(lines)
            out_fname = f"{prefix}_input.in"
            filepath = os.path.join(prefix, out_fname)
            
            with open(filepath, 'w') as f:
                f.write(contents)
                
            # Show success message with run commands
            success_msg = f"""✅ Input file generated: {filepath}

To run locally (serial):
pioud.x -inp {filepath} > {prefix}_output.out

To run locally (parallel):
mpirun -n 4 pioud.x -inp {filepath} > {prefix}_output.out"""
            
            messagebox.showinfo("Success", success_msg)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate input file:\n{e}")
            
    def _generate_script(self):
        """Generate HPC submission script"""
        prefix = self.prefix_var.get().strip()
        if not prefix:
            messagebox.showwarning("No Prefix", "Please set a prefix first.")
            return
            
        input_file = f"{prefix}/{prefix}_input.in"
        if not os.path.exists(input_file):
            messagebox.showwarning("No Input File", 
                                 f"Input file {input_file} does not exist. Generate it first.")
            return
            
        try:
            # Calculate total cores
            total_cores = int(self.hpcc_vars['nodes'].get()) * int(self.hpcc_vars['cores_per_node'].get())
            output_file = f"{prefix}_output.out"
            
            # Create script content based on scheduler type
            if self.hpcc_vars['scheduler'].get() == 'PBS':
                script_content = [
                    "#!/bin/bash",
                    f"#PBS -N {self.hpcc_vars['job_name'].get()}",
                    f"#PBS -l select={self.hpcc_vars['nodes'].get()}:ncpus={self.hpcc_vars['cores_per_node'].get()}:mem={self.hpcc_vars['memory'].get()}",
                    f"#PBS -l walltime={self.hpcc_vars['walltime'].get()}",
                    f"#PBS -q {self.hpcc_vars['queue'].get()}"
                ]
                
                if self.hpcc_vars['email_notify'].get() and self.hpcc_vars['email_address'].get():
                    script_content.extend([
                        "#PBS -m abe",
                        f"#PBS -M {self.hpcc_vars['email_address'].get()}"
                    ])
                    
                script_content.extend([
                    "",
                    "# Change to submission directory",
                    "cd $PBS_O_WORKDIR",
                    ""
                ])
                
            else:  # SLURM
                script_content = [
                    "#!/bin/bash",
                    f"#SBATCH --job-name={self.hpcc_vars['job_name'].get()}",
                    f"#SBATCH --nodes={self.hpcc_vars['nodes'].get()}",
                    f"#SBATCH --ntasks-per-node={self.hpcc_vars['cores_per_node'].get()}",
                    f"#SBATCH --time={self.hpcc_vars['walltime'].get()}",
                    f"#SBATCH --partition={self.hpcc_vars['queue'].get()}",
                    f"#SBATCH --mem={self.hpcc_vars['memory'].get()}"
                ]
                
                if self.hpcc_vars['email_notify'].get() and self.hpcc_vars['email_address'].get():
                    script_content.extend([
                        "#SBATCH --mail-type=ALL",
                        f"#SBATCH --mail-user={self.hpcc_vars['email_address'].get()}"
                    ])
                    
                script_content.extend([
                    "",
                    "# Change to submission directory",
                    "cd $SLURM_SUBMIT_DIR",
                    ""
                ])
                
            # Add module loading commands
            script_content.append("# Load required modules")
            script_content.extend(self.modules_text.get('1.0', 'end-1c').splitlines())
            script_content.append("")
            
            # Add additional setup commands
            script_content.append("# Setup environment")
            script_content.extend(self.setup_text.get('1.0', 'end-1c').splitlines())
            script_content.append("")
            
            # Add the command to run the job
            script_content.append("# Run the job")
            script_content.append(f"{self.hpcc_vars['mpi_command'].get()} {total_cores} pioud.x -inp {prefix}_input.in > {output_file}")
            script_content.append("")
            
            # Write the script to a file
            script_filename = f"{prefix}/{prefix}_submit.sh"
            with open(script_filename, 'w') as f:
                f.write('\n'.join(script_content))
                
            # Show success message
            submit_cmd = 'qsub' if self.hpcc_vars['scheduler'].get() == 'PBS' else 'sbatch'
            success_msg = f"""✅ Submission script generated: {script_filename}

To submit: {submit_cmd} {script_filename}"""
            
            messagebox.showinfo("Success", success_msg)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate submission script:\n{e}")