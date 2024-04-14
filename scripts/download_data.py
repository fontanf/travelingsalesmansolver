import argparse
import gdown
import os
import shutil
import pathlib


gdown.download(id="1fnpRMRA2-AWhYf9HMPxs7AblZViZb_cr", output="data.7z")
os.system("7z x data.7z")
pathlib.Path("data.7z").unlink()
shutil.copytree("traveling_salesman", "data", dirs_exist_ok=True)
