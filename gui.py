#!/usr/bin/env python3
"""
EZ-ConvPhase GUI Wrapper
A streamlined GUI for the ez-convphase pipeline with ConvPhase parameter control
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import os
import subprocess
import threading
import json
import signal
from pathlib import Path
import sys

class EZConvPhaseGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("EZ-ConvPhase Pipeline GUI")
        self.root.geometry("800x650")
        
        # Variables for file paths and parameters
        self.input_dir = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.script_dir = tk.StringVar()
        
        # ConvPhase parameters
        self.iterations = tk.IntVar(value=10000)
        self.thinning = tk.IntVar(value=10)
        self.burnin = tk.IntVar(value=1000)
        self.phase_threshold = tk.DoubleVar(value=0.9)
        self.allele_threshold = tk.DoubleVar(value=0.9)
        self.jobs = tk.IntVar(value=4)
        
        # File handling options
        self.skip_existing = tk.BooleanVar(value=True)
        self.keep_intermediate = tk.BooleanVar(value=True)
        
        # Process control
        self.process = None
        self.is_running = False
        
        self.create_widgets()
        self.detect_script_location()
        
    def create_widgets(self):
        # Create notebook for tabs
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill="both", expand=True, padx=10, pady=5)
        
        # Main tab
        main_frame = ttk.Frame(notebook)
        notebook.add(main_frame, text="Main Pipeline")
        
        # Parameters tab
        params_frame = ttk.Frame(notebook)
        notebook.add(params_frame, text="Parameters")
        
        # Parameter guide tab
        guide_frame = ttk.Frame(notebook)
        notebook.add(guide_frame, text="Parameter Guide")
        
        # Logs tab
        logs_frame = ttk.Frame(notebook)
        notebook.add(logs_frame, text="Logs & Output")
        
        self.create_main_tab(main_frame)
        self.create_params_tab(params_frame)
        self.create_guide_tab(guide_frame)
        self.create_logs_tab(logs_frame)
        
    def create_main_tab(self, parent):
        # Title
        title_label = ttk.Label(parent, text="EZ-ConvPhase Pipeline", 
                               font=("Arial", 16, "bold"))
        title_label.pack(pady=10)
        
        # Input/Output section
        io_frame = ttk.LabelFrame(parent, text="Input/Output Configuration", padding=10)
        io_frame.pack(fill="x", padx=10, pady=5)
        
        # Input directory
        ttk.Label(io_frame, text="Input Directory (FASTA files):").grid(row=0, column=0, sticky="w", pady=2)
        input_frame = ttk.Frame(io_frame)
        input_frame.grid(row=1, column=0, columnspan=2, sticky="ew", pady=2)
        ttk.Entry(input_frame, textvariable=self.input_dir, width=60).pack(side="left", fill="x", expand=True)
        ttk.Button(input_frame, text="Browse", command=self.browse_input_dir).pack(side="right", padx=(5,0))
        
        # Output directory
        ttk.Label(io_frame, text="Output Directory:").grid(row=2, column=0, sticky="w", pady=2)
        output_frame = ttk.Frame(io_frame)
        output_frame.grid(row=3, column=0, columnspan=2, sticky="ew", pady=2)
        ttk.Entry(output_frame, textvariable=self.output_dir, width=60).pack(side="left", fill="x", expand=True)
        ttk.Button(output_frame, text="Browse", command=self.browse_output_dir).pack(side="right", padx=(5,0))
        
        # Script directory
        ttk.Label(io_frame, text="EZ-ConvPhase Scripts Directory:").grid(row=4, column=0, sticky="w", pady=2)
        script_frame = ttk.Frame(io_frame)
        script_frame.grid(row=5, column=0, columnspan=2, sticky="ew", pady=2)
        ttk.Entry(script_frame, textvariable=self.script_dir, width=60).pack(side="left", fill="x", expand=True)
        ttk.Button(script_frame, text="Browse", command=self.browse_script_dir).pack(side="right", padx=(5,0))
        
        io_frame.columnconfigure(0, weight=1)
        
        # Processing options
        proc_frame = ttk.LabelFrame(parent, text="Processing Options", padding=10)
        proc_frame.pack(fill="x", padx=10, pady=5)
        
        # Parallel jobs
        ttk.Label(proc_frame, text="Parallel jobs:").grid(row=0, column=0, sticky="w", padx=5)
        ttk.Spinbox(proc_frame, from_=1, to=32, textvariable=self.jobs, width=10).grid(row=0, column=1, sticky="w", padx=5)
        ttk.Label(proc_frame, text="(0 = use all available cores)").grid(row=0, column=2, sticky="w", padx=5)
        
        # File handling options
        file_frame = ttk.LabelFrame(parent, text="File Handling", padding=10)
        file_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Checkbutton(file_frame, text="Skip files with existing output", 
                       variable=self.skip_existing).pack(anchor="w")
        
        ttk.Checkbutton(file_frame, text="Keep intermediate files (clean, mask, phased)", 
                       variable=self.keep_intermediate).pack(anchor="w")
        
        # File type info
        info_frame = ttk.LabelFrame(parent, text="Supported File Types", padding=10)
        info_frame.pack(fill="x", padx=10, pady=5)
        
        info_text = ("Supported formats: .fa, .fasta, .fas\n"
                    "All sequences must be aligned (same length)\n" 
                    "Headers can be: sample|Species or sample.Species format")
        ttk.Label(info_frame, text=info_text, justify="left").pack()
        
        # Control buttons
        control_frame = ttk.Frame(parent)
        control_frame.pack(fill="x", padx=10, pady=10)
        
        self.run_button = ttk.Button(control_frame, text="Run Pipeline", 
                                   command=self.run_pipeline, style="Accent.TButton")
        self.run_button.pack(side="left", padx=5)
        
        self.stop_button = ttk.Button(control_frame, text="Stop", 
                                    command=self.stop_pipeline, state="disabled")
        self.stop_button.pack(side="left", padx=5)
        
        ttk.Button(control_frame, text="Validate Setup", 
                  command=self.validate_setup).pack(side="left", padx=5)
        
        ttk.Button(control_frame, text="Import Config", 
                  command=self.import_config).pack(side="right", padx=5)
        
        ttk.Button(control_frame, text="Export Config", 
                  command=self.export_config).pack(side="right", padx=5)
        
        ttk.Button(control_frame, text="Open Output Folder", 
                  command=self.open_output_folder).pack(side="right", padx=5)
        
        # Progress bar
        self.progress = ttk.Progressbar(parent, mode='indeterminate')
        self.progress.pack(fill="x", padx=10, pady=5)
        
    def create_params_tab(self, parent):
        # ConvPhase MCMC Parameters
        mcmc_frame = ttk.LabelFrame(parent, text="MCMC Parameters", padding=10)
        mcmc_frame.pack(fill="x", padx=10, pady=5)
        
        # Iterations
        ttk.Label(mcmc_frame, text="Number of iterations:").grid(row=0, column=0, sticky="w", pady=2)
        ttk.Entry(mcmc_frame, textvariable=self.iterations, width=15).grid(row=0, column=1, sticky="w", padx=5)
        ttk.Label(mcmc_frame, text="(Number of MCMC iterations)").grid(row=0, column=2, sticky="w", padx=5)
        
        # Thinning
        ttk.Label(mcmc_frame, text="Thinning interval:").grid(row=1, column=0, sticky="w", pady=2)
        ttk.Entry(mcmc_frame, textvariable=self.thinning, width=15).grid(row=1, column=1, sticky="w", padx=5)
        ttk.Label(mcmc_frame, text="(Sample every nth iteration)").grid(row=1, column=2, sticky="w", padx=5)
        
        # Burn-in
        ttk.Label(mcmc_frame, text="Burn-in:").grid(row=2, column=0, sticky="w", pady=2)
        ttk.Entry(mcmc_frame, textvariable=self.burnin, width=15).grid(row=2, column=1, sticky="w", padx=5)
        ttk.Label(mcmc_frame, text="(Iterations to discard at start)").grid(row=2, column=2, sticky="w", padx=5)
        
        # Threshold Parameters
        thresh_frame = ttk.LabelFrame(parent, text="Certainty Thresholds", padding=10)
        thresh_frame.pack(fill="x", padx=10, pady=5)
        
        # Phase threshold (-p)
        ttk.Label(thresh_frame, text="Phase threshold (-p):").grid(row=0, column=0, sticky="w", pady=2)
        ttk.Entry(thresh_frame, textvariable=self.phase_threshold, width=15).grid(row=0, column=1, sticky="w", padx=5)
        ttk.Label(thresh_frame, text="(Phase certainty from 0 to 1)").grid(row=0, column=2, sticky="w", padx=5)
        
        # Allele threshold (-q)
        ttk.Label(thresh_frame, text="Allele threshold (-q):").grid(row=1, column=0, sticky="w", pady=2)
        ttk.Entry(thresh_frame, textvariable=self.allele_threshold, width=15).grid(row=1, column=1, sticky="w", padx=5)
        ttk.Label(thresh_frame, text="(Genotype certainty from 0 to 1)").grid(row=1, column=2, sticky="w", padx=5)
        
        # Parameter presets
        preset_frame = ttk.LabelFrame(parent, text="Parameter Presets", padding=10)
        preset_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Button(preset_frame, text="Fast (1000 iter)", 
                  command=lambda: self.set_preset("fast")).pack(side="left", padx=5)
        ttk.Button(preset_frame, text="High Quality (10000+ iter)", 
                  command=lambda: self.set_preset("high")).pack(side="left", padx=5)
        
    def create_guide_tab(self, parent):
        # Parameter explanation
        help_frame = ttk.LabelFrame(parent, text="Parameter Guide", padding=10)
        help_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        help_text = """ConvPhase Parameter Guide:

• Iterations: More iterations = more accurate results but longer runtime
  - Fast: 1,000 iterations (quick testing)
  - High quality: 10,000+ iterations (publication grade)

• Thinning: Sample every nth iteration to reduce autocorrelation
  - Default: 10 (keeps every 10th sample)

• Burn-in: Initial iterations discarded to allow convergence
  - Usually 10-20% of total iterations

• Phase threshold (-p): Minimum certainty for phase assignment
  - 0.9 = 90% certainty required (recommended)

• Allele threshold (-q): Minimum certainty for genotype calls
  - 0.9 = 90% certainty required (recommended)

Processing Options:

• Parallel jobs: Number of loci to process simultaneously
  - More jobs = faster processing but higher memory usage
  - Set to 0 to use all available CPU cores

• Skip existing: Avoids reprocessing completed files
  - Useful for resuming interrupted runs

• Keep intermediate: Retains cleaned and phased files
  - Helpful for debugging pipeline issues
  - Disable to save disk space

Pipeline Stages:

1. Sanitize: Remove problematic characters, create mask
2. Phase: Run ConvPhase on cleaned sequences  
3. Rebuild: Restore original alignment length using mask

The pipeline automatically handles header formats like:
• sample|Species (MolD format)
• sample.Species (HapView format)
• sample_a/sample_b (haplotype suffixes)
"""
        
        help_label = ttk.Label(help_frame, text=help_text, justify="left")
        help_label.pack(anchor="w")
        
    def create_logs_tab(self, parent):
        # Log display
        log_frame = ttk.LabelFrame(parent, text="Pipeline Output", padding=5)
        log_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.log_text = scrolledtext.ScrolledText(log_frame, height=20, wrap=tk.WORD)
        self.log_text.pack(fill="both", expand=True)
        
        # Results summary
        results_frame = ttk.LabelFrame(parent, text="Results Summary", padding=5)
        results_frame.pack(fill="x", padx=10, pady=5)
        
        self.results_text = tk.Text(results_frame, height=4, wrap=tk.WORD)
        self.results_text.pack(fill="x")
        
        # Log controls
        log_control_frame = ttk.Frame(parent)
        log_control_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Button(log_control_frame, text="Clear Logs", 
                  command=self.clear_logs).pack(side="left", padx=5)
        ttk.Button(log_control_frame, text="Save Logs", 
                  command=self.save_logs).pack(side="left", padx=5)
        ttk.Button(log_control_frame, text="View Results", 
                  command=self.analyze_results).pack(side="left", padx=5)
        
        # Status
        self.status_var = tk.StringVar(value="Ready")
        self.status_label = ttk.Label(log_control_frame, textvariable=self.status_var)
        self.status_label.pack(side="right", padx=5)
        
    def detect_script_location(self):
        """Try to automatically detect the script directory"""
        possible_locations = [
            "./scripts",
            "../scripts", 
            "./ez-convphase/scripts",
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts")
        ]
        
        for loc in possible_locations:
            if os.path.exists(os.path.join(loc, "run_pipeline_parallel.sh")):
                self.script_dir.set(os.path.abspath(loc))
                break
                
    def browse_input_dir(self):
        directory = filedialog.askdirectory(title="Select Input Directory with FASTA files")
        if directory:
            self.input_dir.set(directory)
            # Auto-suggest output directory
            if not self.output_dir.get():
                parent_dir = os.path.dirname(directory)
                suggested_output = os.path.join(parent_dir, "ezconvphase_out")
                self.output_dir.set(suggested_output)
                
    def browse_output_dir(self):
        directory = filedialog.askdirectory(title="Select Output Directory")
        if directory:
            self.output_dir.set(directory)
            
    def browse_script_dir(self):
        directory = filedialog.askdirectory(title="Select EZ-ConvPhase Scripts Directory")
        if directory:
            self.script_dir.set(directory)
            
    def set_preset(self, preset_type):
        """Set parameter presets"""
        if preset_type == "fast":
            self.iterations.set(1000)
            self.thinning.set(5)
            self.burnin.set(100)
        elif preset_type == "high":
            self.iterations.set(10000)
            self.thinning.set(10)
            self.burnin.set(1000)
            
    def validate_setup(self):
        """Validate that all required components are available"""
        issues = []
        
        # Check input directory
        if not self.input_dir.get():
            issues.append("Input directory not selected")
        elif not os.path.exists(self.input_dir.get()):
            issues.append("Input directory does not exist")
        else:
            # Check for FASTA files
            fasta_files = []
            for ext in ["*.fa", "*.fasta", "*.fas"]:
                fasta_files.extend(Path(self.input_dir.get()).glob(ext))
            if not fasta_files:
                issues.append("No FASTA files found in input directory")
                
        # Check script directory
        if not self.script_dir.get():
            issues.append("Scripts directory not selected")
        elif not os.path.exists(self.script_dir.get()):
            issues.append("Scripts directory does not exist")
        else:
            required_scripts = ["run_pipeline_parallel.sh", "sanitize_for_convphase.py", "rebuild_from_mask.py"]
            for script in required_scripts:
                if not os.path.exists(os.path.join(self.script_dir.get(), script)):
                    issues.append(f"Missing required script: {script}")
                    
        # Check ConvPhase availability
        try:
            result = subprocess.run(["convphase", "-h"], capture_output=True, text=True)
            if result.returncode != 0:
                issues.append("ConvPhase not found in PATH")
        except FileNotFoundError:
            issues.append("ConvPhase not found in PATH")
            
        # Check parallel availability
        try:
            result = subprocess.run(["parallel", "--version"], capture_output=True, text=True)
            if result.returncode != 0:
                issues.append("GNU parallel not found in PATH")
        except FileNotFoundError:
            issues.append("GNU parallel not found in PATH")
            
        # Report results
        if issues:
            messagebox.showerror("Validation Failed", 
                               "Setup validation failed:\n\n" + "\n".join(f"• {issue}" for issue in issues))
        else:
            messagebox.showinfo("Validation Passed", 
                              "All components are properly configured!")
            
    def log_message(self, message):
        """Add message to log display"""
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.root.update_idletasks()
        
    def run_pipeline(self):
        """Run the EZ-ConvPhase pipeline"""
        if self.is_running:
            return
            
        # Validate required fields
        if not all([self.input_dir.get(), self.output_dir.get(), self.script_dir.get()]):
            messagebox.showerror("Error", "Please specify input directory, output directory, and scripts directory")
            return
            
        self.is_running = True
        self.run_button.config(state="disabled")
        self.stop_button.config(state="normal")
        self.progress.start()
        self.status_var.set("Running...")
        
        # Run in separate thread to prevent GUI freezing
        thread = threading.Thread(target=self._run_pipeline_thread)
        thread.daemon = True
        thread.start()
        
    def _run_pipeline_thread(self):
        """Run pipeline in separate thread"""
        try:
            self.log_message("Starting EZ-ConvPhase pipeline...")
            self.log_message(f"Input: {self.input_dir.get()}")
            self.log_message(f"Output: {self.output_dir.get()}")
            self.log_message(f"Scripts: {self.script_dir.get()}")
            
            # Prepare environment variables
            env = os.environ.copy()
            env.update({
                "JOBS": str(self.jobs.get()),
                "ITER_VAL": str(self.iterations.get()),
                "THIN_VAL": str(self.thinning.get()),
                "BURN_VAL": str(self.burnin.get())
            })
            
            # Build command
            script_path = os.path.join(self.script_dir.get(), "run_pipeline_parallel.sh")
            cmd = ["bash", script_path, self.input_dir.get(), self.output_dir.get()]
            
            self.log_message(f"Command: {' '.join(cmd)}")
            self.log_message("=" * 50)
            
            # Run process
            self.process = subprocess.Popen(
                cmd,
                env=env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True,
                preexec_fn=os.setsid if os.name != 'nt' else None  # Create process group
            )
            
            # Stream output
            for line in iter(self.process.stdout.readline, ''):
                if not self.is_running:  # Check if stopped
                    break
                self.root.after(0, lambda l=line.rstrip(): self.log_message(l))
                
            self.process.wait()
            
            if self.process.returncode == 0:
                self.root.after(0, lambda: self.log_message("Pipeline completed successfully!"))
                self.root.after(0, lambda: self.status_var.set("Completed"))
            else:
                self.root.after(0, lambda: self.log_message(f"Pipeline failed with exit code {self.process.returncode}"))
                self.root.after(0, lambda: self.status_var.set("Failed"))
                
        except Exception as e:
            self.root.after(0, lambda: self.log_message(f"Error: {str(e)}"))
            self.root.after(0, lambda: self.status_var.set("Error"))
            
        finally:
            self.root.after(0, self._pipeline_finished)
            
    def _pipeline_finished(self):
        """Clean up after pipeline finishes"""
        self.is_running = False
        self.run_button.config(state="normal")
        self.stop_button.config(state="disabled")
        self.progress.stop()
        self.process = None
        
    def stop_pipeline(self):
        """Stop the running pipeline and all parallel processes"""
        if self.process and self.is_running:
            self.is_running = False
            self.log_message("Stopping pipeline...")
            
            try:
                # Kill the process group (includes all parallel jobs)
                if os.name != 'nt':  # Unix-like systems
                    os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
                    # Also kill any remaining parallel processes
                    try:
                        subprocess.run(["pkill", "-f", "parallel"], check=False)
                        subprocess.run(["pkill", "-f", "convphase"], check=False)
                    except:
                        pass
                else:  # Windows
                    self.process.terminate()
                    
                self.log_message("Pipeline stopped by user")
                self.status_var.set("Stopped")
            except Exception as e:
                self.log_message(f"Error stopping pipeline: {str(e)}")
                
    def open_output_folder(self):
        """Open the output folder in file manager"""
        if self.output_dir.get() and os.path.exists(self.output_dir.get()):
            if sys.platform == "win32":
                os.startfile(self.output_dir.get())
            elif sys.platform == "darwin":
                subprocess.run(["open", self.output_dir.get()])
            else:
                subprocess.run(["xdg-open", self.output_dir.get()])
        else:
            messagebox.showwarning("Warning", "Output directory does not exist yet")
            
    def clear_logs(self):
        """Clear the log display"""
        self.log_text.delete(1.0, tk.END)
        
    def save_logs(self):
        """Save logs to file"""
        filename = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
        )
        if filename:
            with open(filename, 'w') as f:
                f.write(self.log_text.get(1.0, tk.END))
            messagebox.showinfo("Success", f"Logs saved to {filename}")
            
    def analyze_results(self):
        """Analyze and summarize pipeline results"""
        if not self.output_dir.get() or not os.path.exists(self.output_dir.get()):
            messagebox.showwarning("Warning", "Output directory does not exist")
            return
            
        try:
            results = self._scan_results()
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, results)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to analyze results: {str(e)}")
            
    def _scan_results(self):
        """Scan output directory for results summary"""
        output_path = Path(self.output_dir.get())
        
        # Count files in each stage
        clean_files = len(list((output_path / "clean").glob("*.fasta"))) if (output_path / "clean").exists() else 0
        phased_files = len(list((output_path / "phased").glob("*.fasta"))) if (output_path / "phased").exists() else 0
        rebuilt_files = len(list((output_path / "rebuilt").glob("*.fasta"))) if (output_path / "rebuilt").exists() else 0
        
        # Check for errors in logs
        error_count = 0
        if (output_path / "logs").exists():
            for log_file in (output_path / "logs").glob("*.log"):
                try:
                    with open(log_file, 'r') as f:
                        content = f.read()
                        if "ERROR" in content or "FAILED" in content:
                            error_count += 1
                except:
                    pass
                    
        summary = f"""Pipeline Results Summary:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Processed Files:
• Cleaned: {clean_files} files
• Phased: {phased_files} files  
• Rebuilt: {rebuilt_files} files

Status:
• Errors detected: {error_count} files
• Success rate: {(rebuilt_files / max(clean_files, 1)) * 100:.1f}%

Output Structure:
• Clean alignments: {output_path}/clean/
• Mask files: {output_path}/mask/
• Phased sequences: {output_path}/phased/
• Final results: {output_path}/rebuilt/
• Processing logs: {output_path}/logs/
"""
        return summary

    def export_config(self):
        """Export current configuration to JSON"""
        config = {
            "input_dir": self.input_dir.get(),
            "output_dir": self.output_dir.get(),
            "script_dir": self.script_dir.get(),
            "parameters": {
                "iterations": self.iterations.get(),
                "thinning": self.thinning.get(),
                "burnin": self.burnin.get(),
                "phase_threshold": self.phase_threshold.get(),
                "allele_threshold": self.allele_threshold.get(),
                "jobs": self.jobs.get(),
                "skip_existing": self.skip_existing.get(),
                "keep_intermediate": self.keep_intermediate.get()
            }
        }
        
        filename = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if filename:
            with open(filename, 'w') as f:
                json.dump(config, f, indent=2)
            messagebox.showinfo("Success", f"Configuration saved to {filename}")
            
    def import_config(self):
        """Import configuration from JSON"""
        filename = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if filename:
            try:
                with open(filename, 'r') as f:
                    config = json.load(f)
                    
                # Set paths
                if "input_dir" in config:
                    self.input_dir.set(config["input_dir"])
                if "output_dir" in config:
                    self.output_dir.set(config["output_dir"])
                if "script_dir" in config:
                    self.script_dir.set(config["script_dir"])
                    
                # Set parameters
                if "parameters" in config:
                    params = config["parameters"]
                    for key, var in {
                        "iterations": self.iterations,
                        "thinning": self.thinning,
                        "burnin": self.burnin,
                        "phase_threshold": self.phase_threshold,
                        "allele_threshold": self.allele_threshold,
                        "jobs": self.jobs,
                        "skip_existing": self.skip_existing,
                        "keep_intermediate": self.keep_intermediate
                    }.items():
                        if key in params:
                            var.set(params[key])
                            
                messagebox.showinfo("Success", f"Configuration loaded from {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load configuration: {str(e)}")

def main():
    root = tk.Tk()
    
    # Set up styling
    style = ttk.Style()
    if "clam" in style.theme_names():
        style.theme_use("clam")
    
    app = EZConvPhaseGUI(root)
    
    # Center window
    root.update_idletasks()
    x = (root.winfo_screenwidth() - root.winfo_width()) // 2
    y = (root.winfo_screenheight() - root.winfo_height()) // 2
    root.geometry(f"+{x}+{y}")
    
    root.mainloop()

if __name__ == "__main__":
    main()
                