import os, shlex, subprocess
import tkinter as tk
from tkinter import filedialog, messagebox

def main():
    root = tk.Tk(); root.withdraw()

    in_dir = filedialog.askdirectory(title="Choose input FASTA folder")
    if not in_dir:
        return

    out_dir = filedialog.askdirectory(title="Choose output folder (Cancel for default)")
    if not out_dir:
        out_dir = os.path.abspath(os.path.join(in_dir, os.pardir, "ezconvphase_out"))

    # repo root = one level up from this file
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    script = os.path.join(repo, "scripts", "run_pipeline_parallel.sh")

    cmd = f'/bin/bash {shlex.quote(script)} {shlex.quote(in_dir)} {shlex.quote(out_dir)}'
    try:
        subprocess.check_call(cmd, shell=True)
        messagebox.showinfo("ez-convphase", f"Done!\nOutput: {out_dir}")
    except subprocess.CalledProcessError as e:
        messagebox.showerror("ez-convphase", f"Pipeline failed (code {e.returncode}).")

if __name__ == "__main__":
    main()
