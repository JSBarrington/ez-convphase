import os, sys, subprocess, shlex

def main():
    if len(sys.argv) not in (2,3):
        print("Usage: ezconvphase <input_folder> [output_folder]", file=sys.stderr)
        sys.exit(2)

    in_dir = os.path.abspath(sys.argv[1])
    out_dir = os.path.abspath(sys.argv[2]) if len(sys.argv) == 3 else None

    if not os.path.isdir(in_dir):
        print(f"[ERROR] Not a directory: {in_dir}", file=sys.stderr)
        sys.exit(2)

    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
    script = os.path.join(repo, "scripts", "run_pipeline_parallel.sh")
    if not os.path.exists(script):
        print(f"[ERROR] Missing pipeline script at: {script}", file=sys.stderr)
        sys.exit(2)

    cmd = f'/bin/bash {shlex.quote(script)} {shlex.quote(in_dir)}'
    if out_dir:
        cmd += f' {shlex.quote(out_dir)}'

    print(f"[ez-convphase] {cmd}")
    sys.exit(subprocess.call(cmd, shell=True))
