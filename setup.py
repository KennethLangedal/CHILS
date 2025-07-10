import os
import subprocess
import platform
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py

class CustomBuild(build_py):
    def run(self):
        # Determine the library name based on the OS
        if platform.system() == "Windows":
            subprocess.check_call(["make", "CC=gcc", "libCHILS.dll"])
            lib_name = "libCHILS.dll"
        elif platform.system() == "Darwin":            
            # Find the correct gcc compiler from brew
            try:
                brew_prefix = subprocess.check_output(["brew", "--prefix", "gcc"], text=True).strip()
                gcc_bin_dir = os.path.join(brew_prefix, "bin")
                compilers = sorted([f for f in os.listdir(gcc_bin_dir) if f.startswith("gcc-")])
                if not compilers:
                    raise FileNotFoundError("No gcc compiler found in brew directory")
                cc_path = os.path.join(gcc_bin_dir, compilers[-1])
            except (subprocess.CalledProcessError, FileNotFoundError):
                # Fallback to system gcc if brew command fails
                cc_path = "gcc"

            # Build the library with the correct compiler
            subprocess.check_call(["make", f"CC={cc_path}", "libCHILS.so"])
            lib_name = "libCHILS.so"
        else:
            subprocess.check_call(["make", "CC=gcc", "libCHILS.so"])
            lib_name = "libCHILS.so"

        # The target directory for the library is inside the build folder
        target_dir = os.path.join(self.build_lib, "chils")
        self.mkpath(target_dir)
        
        # The source library is in the root of the sdist temporary directory
        self.copy_file(lib_name, os.path.join(target_dir, lib_name))

        super().run()


setup(
    name='chils',
    version='1.0.0',
    author='Kenneth Langedal',
    packages=find_packages(where="python"),
    package_dir={"" : "python"},
    cmdclass={
        'build_py': CustomBuild,
    },
    package_data={
        'chils': ['*.so', '*.dll'],
    },
    include_package_data=True,
    has_ext_modules=lambda : True
)
