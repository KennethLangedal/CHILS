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
            subprocess.check_call(["make", "CC=$(find $(brew --prefix gcc)/bin -name \"gcc-*\" | head -n 1)", "libCHILS.so"])
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
