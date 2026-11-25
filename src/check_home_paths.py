"""Release check: fail if any file holds a local path. Same expression as the sister project's release script (grigorovia-mito-insertions/code/release/export_prereg.py), written in pieces so that this file does not match itself.
usage: check_home_paths.py [FOLDER]   default: every tracked or untracked-but-not-ignored file of this project; with FOLDER: every file under it (.git excluded). Exit 1 and a list if anything is found."""
import os,re,subprocess,sys
PAT=re.compile(("/"+r"Users/[^/\s]+|"+"/"+r"home/[^/\s]+|"+"/"+r"private/tmp|"+"/"+r"opt/homebrew").encode())
def files():
    if len(sys.argv)>1:
        for d,sub,fs in os.walk(sys.argv[1]):
            sub[:]=[x for x in sub if x!=".git"]
            for f in fs: yield os.path.join(d,f)
    else:
        root=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..")
        for f in subprocess.run(["git","ls-files","-co","--exclude-standard","-z"],cwd=root,capture_output=True,text=True,check=True).stdout.split("\0"):
            if f and os.path.isfile(os.path.join(root,f)): yield os.path.join(root,f)
bad=0;n=0
for f in files():
    n+=1; m=PAT.findall(open(f,"rb").read())
    if m: bad+=1; print(os.path.relpath(f),len(m),sorted({x.decode("utf8","replace") for x in m})[:3])
print(f"{n} files checked, {bad} with local paths"); sys.exit(1 if bad else 0)
