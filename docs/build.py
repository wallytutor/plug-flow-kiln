# -*- coding: utf-8 -*-
import subprocess

filename = "theory.tex"

latex_cmd = [
    "xelatex.exe",
    "-synctex=1",
    "-shell-escape",
    "-interaction=nonstopmode",
    filename
]

bibtex_cmd = [
    "bibtex.exe",
    filename
]

subprocess.run(latex_cmd)
subprocess.run(latex_cmd)
subprocess.run(bibtex_cmd)
subprocess.run(latex_cmd)
